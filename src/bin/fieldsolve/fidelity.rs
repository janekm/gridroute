//! Fidelity policies for geometry iteration: `fieldsolve run SPEC --fast` (screening) or
//! `--fidelity screening|refined|validated`. A plain `fieldsolve run` is level `spec`: the spec as given.
//! `validated` keeps the spec's mesh and outputs and only replaces the energy/port stop by a stricter
//! objective rule (0.1 dB over 8 periods); `refined` and `screening` coarsen as below.
//!
//! A policy rewrites the spec's numerical effort and nothing of its physics: the whole board stays in the
//! domain (radiation needs the ground plane), materials, loss, ports, PML and air are untouched. It
//! - coarsens the fine mesh (`cell`, and `cell_z` with it) by `cell_scale`, while every port and the spec's
//!   `"fidelity": {"protect": [...]}` features (feed, gaps, shorting stubs) keep the spec's own cell size
//!   through `"refine"` boxes;
//! - drops volume fields, current frames and the scene export, and keeps far-field frequencies only in the
//!   objective's band (plus `objective.harmonics`);
//! - stops each excitation when the objective's components have stayed within a spread of `tol_db` over
//!   `window_periods` periods of the band's lowest frequency, after the source and a minimum settling time
//!   of `min_settle_periods` periods (the energy/port decay rules still apply if they come first; a run that
//!   reaches `max_steps` is reported as incomplete).
//!
//! Every override is recorded (`"fidelity"` in result.json) together with the actual mesh, time step,
//! excitation count, stop reasons and objective stability. Policy numbers can be overridden in the spec's
//! `"fidelity"` object (`cell_scale`, `protect_margin`, `tol_db`, `window_periods`, `min_settle_periods`), and
//! `"fidelity": {"set": {"screening": {"key": value}, ...}}` sets any spec key for one level: the hook for
//! model corrections (edge singularity, conformal metal) that matter most on coarse cells.
//! The measured accuracy of the default policies is in docs/fieldsolve-fast-mode.md.

use serde_json::{json, Value};

#[derive(Clone, Debug, PartialEq)]
pub struct Policy {
    pub level: String,
    pub cell_scale: f64,
    pub protect_margin: f64,
    pub tol_db: f64,
    pub window_periods: f64,
    pub min_settle_periods: f64,
}

impl Policy {
    pub fn defaults(level: &str) -> Result<Policy, String> {
        let p = |cell_scale, tol_db, window_periods, min_settle_periods| Policy { level: level.to_string(), cell_scale, protect_margin: 0.25, tol_db, window_periods, min_settle_periods };
        match level {
            // stop rules from replays on the MTH PIFA's full-mesh history (docs/fieldsolve-fast-mode.md)
            "screening" => Ok(p(2.0, 0.2, 4.0, 8.0)),
            "refined" => Ok(p(1.5, 0.1, 6.0, 8.0)),
            "validated" => Ok(p(1.0, 0.1, 8.0, 8.0)),
            "spec" => Ok(p(1.0, 0.0, 0.0, 0.0)),
            _ => Err(format!("unknown fidelity {level:?}: use screening, refined or validated")),
        }
    }
}

/// Resolve `level` against `spec`: the rewritten spec (with `"fidelity_record"` for result.json).
pub fn resolve(spec: &Value, level: &str) -> Result<Value, String> {
    let mut pol = Policy::defaults(level)?;
    let user = spec.get("fidelity").cloned().unwrap_or(json!({}));
    let mut overrides: Vec<Value> = Vec::new();
    let mut user_set = Vec::new();
    if level != "spec" {
        for (k, slot) in [("cell_scale", &mut pol.cell_scale), ("protect_margin", &mut pol.protect_margin), ("tol_db", &mut pol.tol_db), ("window_periods", &mut pol.window_periods), ("min_settle_periods", &mut pol.min_settle_periods)] {
            if let Some(v) = user.get(k).and_then(|v| v.as_f64()) {
                *slot = v;
                user_set.push(k);
            }
        }
    }
    let mut s = spec.clone();
    // a key set to `to` (Null removes it), recorded when it changes
    let set = |s: &mut Value, ov: &mut Vec<Value>, key: &str, to: Value, why: &str| {
        let from = s.get(key).cloned().unwrap_or(Value::Null);
        if from != to {
            ov.push(json!({"key": key, "from": from, "to": to, "why": why}));
            if to.is_null() {
                s.as_object_mut().unwrap().remove(key);
            } else {
                s[key] = to;
            }
        }
    };
    // the objective: the spec's, or the far-field band at the first port
    if s.get("objective").and_then(|o| o.get("band")).is_none() {
        if level == "spec" {
            // nothing to track: record the level only
        } else {
            let ff = s.pointer("/farfield/freqs").cloned().ok_or("--fast needs an \"objective\": {\"port\", \"band\"} (or far-field frequencies to use as the band)")?;
            set(&mut s, &mut overrides, "objective", json!({"band": ff}), "no objective in the spec: the far-field frequencies as the band, first port");
        }
    }
    if level == "validated" {
        // the full model and outputs, with strict objective convergence
        set(&mut s, &mut overrides, "stop", json!("objective"), "stop on strict objective convergence");
        let mut o = s["objective"].clone();
        o["tol_db"] = json!(pol.tol_db);
        o["window_periods"] = json!(pol.window_periods);
        o["min_settle_periods"] = json!(pol.min_settle_periods);
        set(&mut s, &mut overrides, "objective", o, "stop rule parameters of the policy");
    } else if level != "spec" {
        let cell = s.get("cell").and_then(|v| v.as_f64()).unwrap_or(0.1);
        set(&mut s, &mut overrides, "cell", json!(cell * pol.cell_scale), "coarser fine mesh");
        if let Some(cz) = s.get("cell_z").and_then(|v| v.as_f64()) {
            set(&mut s, &mut overrides, "cell_z", json!(cz * pol.cell_scale), "coarser fine mesh through the board");
        }
        // protected features keep the spec's cell size: every port, plus the spec's list
        let mut refine = s.get("refine").and_then(|v| v.as_array()).cloned().unwrap_or_default();
        refine.push(json!({"ports": true, "cell": cell, "margin": pol.protect_margin}));
        for p in user.get("protect").and_then(|v| v.as_array()).cloned().unwrap_or_default() {
            let mut e = if p.is_string() { json!({"pads": [p]}) } else if p.as_array().is_some_and(|a| a.len() == 4) { json!({"box": p}) } else { p.clone() };
            if e.get("cell").is_none() {
                e["cell"] = json!(cell);
            }
            if e.get("margin").is_none() {
                e["margin"] = json!(pol.protect_margin);
            }
            refine.push(e);
        }
        set(&mut s, &mut overrides, "refine", json!(refine), "protected features keep the spec's cell size");
        for k in ["fields", "frames"] {
            set(&mut s, &mut overrides, k, Value::Null, "not needed to score the candidate");
        }
        set(&mut s, &mut overrides, "scene", json!(false), "no scene export");
        let band: Vec<f64> = s.pointer("/objective/band").and_then(|v| v.as_array()).map(|a| a.iter().filter_map(|x| x.as_f64()).collect()).unwrap_or_default();
        let harm: Vec<f64> = s.pointer("/objective/harmonics").and_then(|v| v.as_array()).map(|a| a.iter().filter_map(|x| x.as_f64()).collect()).unwrap_or_default();
        if let Some(ff) = s.get("farfield").cloned() {
            let keep: Vec<f64> = band.iter().chain(&harm).copied().collect();
            let mut f2 = ff.clone();
            f2["freqs"] = json!(keep);
            set(&mut s, &mut overrides, "farfield", f2, "far field only at the objective's frequencies");
        }
        set(&mut s, &mut overrides, "stop", json!("objective"), "stop on objective stability");
        let mut o = s["objective"].clone();
        o["tol_db"] = json!(pol.tol_db);
        o["window_periods"] = json!(pol.window_periods);
        o["min_settle_periods"] = json!(pol.min_settle_periods);
        set(&mut s, &mut overrides, "objective", o, "stop rule parameters of the policy");
    }
    // per-level spec switches, e.g. {"set": {"screening": {"conformal": true}}}: model corrections that help
    // most on coarse cells (edge singularity, conformal metal) are switched on per level from the spec,
    // without this module knowing their names; recorded like every other override
    if let Some(m) = user.pointer(&format!("/set/{level}")).and_then(|v| v.as_object()) {
        for (k, v) in m {
            set(&mut s, &mut overrides, k, v.clone(), "per-level setting from the spec's fidelity.set");
        }
    }
    s["fidelity_record"] = json!({
        "level": level,
        "policy": {"cell_scale": pol.cell_scale, "protect_margin_mm": pol.protect_margin, "tol_db": pol.tol_db, "window_periods": pol.window_periods, "min_settle_periods": pol.min_settle_periods},
        "user_overrides": user_set,
        "overrides": overrides,
        "physics_kept": "whole board, materials, loss model, ports, PML and air as in the spec",
    });
    Ok(s)
}

/// result.json's `"fidelity"`: the record plus what the run actually did, and the status it earns:
/// the level when every excitation stopped on a convergence rule, else "incomplete".
pub fn finish_record(record: &Value, res: &Value) -> Value {
    let mut r = record.clone();
    if r.is_null() {
        return r;
    }
    let stops: Vec<String> = res["runs"].as_array().map(|a| a.iter().map(|x| x["stop"].as_str().unwrap_or("").to_string()).collect()).unwrap_or_default();
    let complete = !stops.is_empty() && stops.iter().all(|s| s != "max_steps");
    r["status"] = json!(if complete { r["level"].as_str().unwrap_or("").to_string() } else { "incomplete".to_string() });
    r["actual"] = json!({
        "mesh": res["mesh"],
        "dt": res["dt"],
        "excitations": stops.len(),
        "stops": stops,
        "steps": res["runs"].as_array().map(|a| a.iter().map(|x| x["steps"].clone()).collect::<Vec<_>>()),
        "seconds": res["runs"].as_array().map(|a| a.iter().map(|x| x["seconds"].as_f64().unwrap_or(0.0)).sum::<f64>()),
        "objective_worst_db": res.pointer("/objective/worst_db"),
        "objective_stability_db": res.pointer("/objective/stability_db"),
    });
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec() -> Value {
        json!({"cell": 0.05, "fields": {"freqs": [2.44e9]}, "frames": {"layers": ["F.Cu"]},
               "farfield": {"freqs": [2.4e9, 2.44e9, 2.48e9, 4.88e9, 7.32e9], "step_deg": 5},
               "objective": {"port": "C20", "band": [2.4e9, 2.44e9, 2.48e9], "harmonics": [4.88e9]},
               "fidelity": {"protect": ["C20.1", [153.9, 75.5, 156.2, 77.0]]},
               "refine": [{"box": [0, 0, 1, 1], "cell": 0.02}]})
    }

    /// Every override is recorded with its old and new value; the physics keys are left alone.
    #[test]
    fn screening_records_every_override() {
        let s = resolve(&spec(), "screening").unwrap();
        let rec = &s["fidelity_record"];
        assert_eq!(rec["level"], "screening");
        let keys: Vec<&str> = rec["overrides"].as_array().unwrap().iter().map(|o| o["key"].as_str().unwrap()).collect();
        for k in ["cell", "refine", "fields", "frames", "scene", "farfield", "stop", "objective"] {
            assert!(keys.contains(&k), "{k} not recorded: {keys:?}");
        }
        for o in rec["overrides"].as_array().unwrap() {
            let k = o["key"].as_str().unwrap();
            assert_eq!(o["to"], s.get(k).cloned().unwrap_or(Value::Null), "{k}: recorded value is not what the spec got");
        }
        assert_eq!(s["cell"], json!(0.1));
        assert!(s.get("fields").is_none() && s.get("frames").is_none());
        assert_eq!(s["farfield"]["freqs"], json!([2.4e9, 2.44e9, 2.48e9, 4.88e9]));
        assert_eq!(s["stop"], "objective");
        assert_eq!(s["objective"]["tol_db"], json!(0.2));
        // ports and protected features keep the spec's 0.05 mm; the user's own refine box is kept
        let r = s["refine"].as_array().unwrap();
        assert_eq!(r[0], json!({"box": [0, 0, 1, 1], "cell": 0.02}));
        assert!(r.iter().any(|e| e["ports"] == json!(true) && e["cell"] == json!(0.05)));
        assert!(r.iter().any(|e| e["pads"] == json!(["C20.1"]) && e["cell"] == json!(0.05)));
        assert!(r.iter().any(|e| e["box"] == json!([153.9, 75.5, 156.2, 77.0])));
        for k in ["air", "pml", "loss", "dielectrics", "extent", "ports", "elements", "max_steps"] {
            assert!(!keys.contains(&k), "{k} is physics and must not be changed");
        }
    }

    #[test]
    fn spec_changes_nothing_validated_only_the_stop_and_user_values_win() {
        let s = resolve(&spec(), "spec").unwrap();
        let mut plain = s.clone();
        plain.as_object_mut().unwrap().remove("fidelity_record");
        assert_eq!(plain, spec());
        assert!(s["fidelity_record"]["overrides"].as_array().unwrap().is_empty());
        let v = resolve(&spec(), "validated").unwrap();
        let keys: Vec<&str> = v["fidelity_record"]["overrides"].as_array().unwrap().iter().map(|o| o["key"].as_str().unwrap()).collect();
        assert_eq!(keys, ["stop", "objective"]);
        assert_eq!(v["objective"]["tol_db"], json!(0.1));
        assert_eq!((v["cell"].clone(), v["fields"].clone(), v["farfield"].clone()), (spec()["cell"].clone(), spec()["fields"].clone(), spec()["farfield"].clone()));
        let mut sp = spec();
        sp["fidelity"]["cell_scale"] = json!(1.5);
        sp["fidelity"]["tol_db"] = json!(0.2);
        let s = resolve(&sp, "screening").unwrap();
        assert_eq!(s["cell"], json!(0.05 * 1.5));
        assert_eq!(s["objective"]["tol_db"], json!(0.2));
        assert_eq!(s["fidelity_record"]["user_overrides"], json!(["cell_scale", "tol_db"]));
        assert!(resolve(&json!({"cell": 0.1}), "screening").is_err(), "no objective and no far field");
        // per-level pass-through, recorded
        let mut sp = spec();
        sp["fidelity"]["set"] = json!({"screening": {"conformal": true}, "validated": {"conformal": false}});
        let s = resolve(&sp, "screening").unwrap();
        assert_eq!(s["conformal"], json!(true));
        assert!(s["fidelity_record"]["overrides"].as_array().unwrap().iter().any(|o| o["key"] == "conformal" && o["to"] == json!(true)));
        assert_eq!(resolve(&sp, "validated").unwrap()["conformal"], json!(false));
        assert!(resolve(&sp, "refined").unwrap().get("conformal").is_none());
        assert!(resolve(&spec(), "quick").is_err());
    }

    #[test]
    fn status_is_incomplete_after_max_steps() {
        let rec = json!({"level": "screening"});
        let ok = json!({"runs": [{"stop": "objective", "steps": 10, "seconds": 1.0}], "mesh": {}, "dt": 1e-13});
        assert_eq!(finish_record(&rec, &ok)["status"], "screening");
        let bad = json!({"runs": [{"stop": "objective"}, {"stop": "max_steps"}]});
        assert_eq!(finish_record(&rec, &bad)["status"], "incomplete");
        assert_eq!(finish_record(&rec, &bad)["actual"]["excitations"], 2);
    }
}
