// Route a SimpleRouteJson with tscircuit's autorouter (@tscircuit/capacity-autorouter) and save its traces.
//
//   node route.mjs SRJ.json OUT.json [--seconds S] [--pipeline NAME]
//
// The package is resolved from $TSCIRCUIT_AUTOROUTER (a node_modules/@tscircuit/capacity-autorouter directory)
// or the usual node resolution. Network access is disabled: some pipelines can ask a remote solution cache,
// and board geometry must not leave the machine.
import fs from "node:fs"
import path from "node:path"
import { pathToFileURL } from "node:url"

const args = process.argv.slice(2)
const opt = (name, dflt) => {
  const i = args.indexOf(name)
  return i >= 0 ? args.splice(i, 2)[1] : dflt
}
const seconds = Number(opt("--seconds", "600"))
const pipeline = opt("--pipeline", "AutoroutingPipelineSolver")
const liftFlag = args.indexOf("--lift-pair-cap")
const liftPairCap = liftFlag >= 0
if (liftPairCap) args.splice(liftFlag, 1)
let liftedPairCap = null
const [srjPath, outPath] = args

globalThis.fetch = async () => {
  throw new Error("network disabled for benchmark")
}

const pkgDir = process.env.TSCIRCUIT_AUTOROUTER
const mod = pkgDir
  ? await import(pathToFileURL(path.join(pkgDir, "dist/index.js")).href)
  : await import("@tscircuit/capacity-autorouter")
const version = JSON.parse(
  fs.readFileSync(path.join(pkgDir ?? "node_modules/@tscircuit/capacity-autorouter", "package.json"), "utf8"),
).version

const srj = JSON.parse(fs.readFileSync(srjPath, "utf8"))
const Solver = mod[pipeline]
const t0 = performance.now()
const solver = new Solver(srj)
let steps = 0
let lastLog = t0
let timedOut = false
let error = null
try {
  while (!solver.solved && !solver.failed) {
    solver.step()
    steps++
    // NetToPointPairsSolver takes one iteration per connection but keeps BaseSolver's cap of 1000, so a board
    // with more connections fails before routing; --lift-pair-cap raises that one cap (recorded in the result)
    const pairs = solver.netToPointPairsSolver
    if (liftPairCap && pairs && !pairs.solved && pairs.MAX_ITERATIONS < srj.connections.length + 10) {
      pairs.MAX_ITERATIONS = srj.connections.length + 10
      liftedPairCap = pairs.MAX_ITERATIONS
    }
    const now = performance.now()
    if (now - t0 > seconds * 1000) {
      timedOut = true
      break
    }
    if (now - lastLog > 10000) {
      lastLog = now
      const stage = solver.getCurrentPhase?.() ?? solver.activeSubSolver?.constructor?.name ?? "?"
      console.error(`${((now - t0) / 1000).toFixed(0)} s, ${steps} steps, ${stage}`)
    }
  }
} catch (e) {
  error = String(e?.stack ?? e)
}
const routeSeconds = (performance.now() - t0) / 1000
let traces = []
let outputError = null
try {
  traces = solver.getOutputSimpleRouteJson().traces ?? []
} catch (e) {
  outputError = String(e?.message ?? e)
}
const result = {
  autorouter: "@tscircuit/capacity-autorouter",
  version,
  pipeline,
  solver: solver.getSolverName?.() ?? solver.constructor.name,
  seconds_limit: seconds,
  lifted_pair_cap: liftedPairCap,
  route_seconds: routeSeconds,
  steps,
  solved: Boolean(solver.solved),
  failed: Boolean(solver.failed),
  timed_out: timedOut,
  solver_error: solver.error ?? null,
  exception: error,
  output_error: outputError,
  phase: solver.getCurrentPhase?.() ?? null,
  traces,
}
fs.writeFileSync(outPath, JSON.stringify(result))
console.log(
  JSON.stringify({ ...result, traces: traces.length, exception: error?.split("\n")[0] ?? null }),
)
process.exit(0)
