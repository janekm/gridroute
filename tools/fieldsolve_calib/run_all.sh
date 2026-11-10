#!/bin/zsh
# Fast-mode calibration batch, one GPU job at a time: full runs (the spec as given), screening on every
# candidate, refined on six. FS is the fieldsolve binary (default target/release/fieldsolve of this checkout).
cd "$(dirname "$0")"
FS=${FS:-${0:a:h}/../../target/release/fieldsolve}
C=(base trim1 trim2 ext1 ext15 feedL05 feedL10 feedR03 gnd+05 gnd+10 gnd-05 trim1_feedL05)
for c in $C; do
  [ -f runs/$c/result.json ] || { echo "$(date +%T) full $c"; $FS run $c.json --out runs/$c > runs/$c.log 2>&1; }
done
for c in $C; do
  [ -f runs/${c}_screen/result.json ] || { echo "$(date +%T) screening $c"; $FS run $c.json --fast --out runs/${c}_screen > runs/${c}_screen.log 2>&1; }
done
for c in base trim2 ext15 feedL10 gnd+10 gnd-05; do
  [ -f runs/${c}_refined/result.json ] || { echo "$(date +%T) refined $c"; $FS run $c.json --fidelity refined --out runs/${c}_refined > runs/${c}_refined.log 2>&1; }
done
echo "$(date +%T) all done"
