#!/usr/bin/env bash
# End-to-end test against the installed SolidCraft (solidcraft-cli, headless, no window).
# Runs every library material through the real card validator, applies materials, reads the
# report, checks undo, and reopens the saved design with no plugin involved.
# Usage: tests/cli-test.sh   (after ./build.sh)
set -euo pipefail
cd "$(dirname "$0")/.."
pkg=$(pwd)/target/package/gen-design-engine
[ -f "$pkg/plugin.wasm" ] || { echo "run ./build.sh first" >&2; exit 1; }
work=$(mktemp -d /tmp/gde-cli-test.XXXXXX)
trap '[ "${KEEP:-}" ] || rm -rf "$work"' EXIT
cp -r "$pkg" "$work/gen-design-engine"
cd "$work"
fail() { echo "FAIL: $*" >&2; exit 1; }

# 1. Inspect (no code runs) and take the approval hash.
echo '[{"command":"plugin.inspect","params":{"path":"gen-design-engine"}}]' > inspect.json
solidcraft-cli run inspect.json --results --json-out inspect.out --quiet
approval=$(jq -r '.[0].approval' inspect.out)
[ "$(jq -r '.[0].within_host_limits' inspect.out)" = true ] || fail "package is over the host limits: $(jq -c '.[0].limits_error' inspect.out)"
[ ${#approval} -ge 32 ] || fail "no approval hash: $(cat inspect.out)"
echo "inspect ok, approval ${approval:0:16}..."

run_step() { # entry params-json
  jq -nc --arg a "$approval" --arg e "$1" --argjson p "$2" '{command:"plugin.run",params:{path:"gen-design-engine",approval:$a,entry:$e,params:$p}}'
}

# 2. Every library material must pass SolidCraft's own card check (material.create).
ids=$(jq -r '.runtime.entries[] | select(.name=="apply-material") | .input_schema.properties.material.oneOf[].const' gen-design-engine/plugin.json)
{
  echo '['
  i=0
  for id in $ids; do
    echo "{\"command\":\"solid.box\",\"params\":{\"length\":10,\"width\":10,\"height\":10,\"body_name\":\"B_$i\"}},"
    run_step apply-material "{\"material\":\"$id\",\"bodies\":[\"B_$i\"]}"; echo ','
    i=$((i+1))
  done
  run_step material-report '{"safety_factor":2}'; echo ','
  echo '{"command":"inspect.measure","params":{}}'
  echo ']'
} > all.json
solidcraft-cli run all.json --results --json-out all.out --save all.sc --quiet || fail "run failed: $(head -c 2000 all.out 2>/dev/null)"
n=$(echo "$ids" | wc -w)
created=$(jq '[.[] | (.result? // .) | select(type=="object" and .created==true)] | length' all.out)
[ "$created" = "$n" ] || fail "created $created of $n materials"
report=$(jq -c '[.[] | (.result? // .) | select(type=="object" and has("stress_analysis"))][0]' all.out)
[ "$(echo "$report" | jq '.bodies | map(select(.allowable_mpa != null)) | length')" = "$n" ] || fail "report misses strengths: $report"
echo "all $n library materials created and assigned; report has allowable stress for each"
echo "$report" | jq -r '.bodies[] | "  \(.body)  \(.material)  \(.mass_g) g  allowable \(.allowable_mpa) MPa (\(.allowable_basis))"'

# 3. Reapplying reuses the material; a one-step undo removes only the assignment.
{
  echo '['
  echo '{"command":"solid.box","params":{"length":10,"width":20,"height":5,"body_name":"Plate"}},'
  run_step apply-material '{"material":"steel-s235jr"}'; echo ','
  run_step apply-material '{"material":"steel-s235jr","bodies":["Plate"]}'; echo ','
  echo '{"command":"inspect.measure","params":{"bodies":["Plate"]}},'
  echo '{"command":"edit.undo","params":{}},'
  echo '{"command":"inspect.measure","params":{"bodies":["Plate"]}}'
  echo ']'
} > undo.json
solidcraft-cli run undo.json --results --json-out undo.out --quiet || fail "undo run failed: $(head -c 2000 undo.out 2>/dev/null)"
[ "$(jq -r '.[1].result.created' undo.out)" = true ] && [ "$(jq -r '.[2].result.created' undo.out)" = false ] || fail "second apply should reuse"
mass=$(jq -r '.[3].bodies[0].mass_g' undo.out)
python3 -I -c "import sys; sys.exit(0 if abs($mass - 7.85) < 1e-6 else 1)" || fail "steel plate mass $mass, expected 7.85 g"
[ "$(jq -r '.[5].bodies[0].material' undo.out)" = "Steel S235JR" ] || fail "undo of the reapply changed the material"
echo "reuse + undo ok (plate 7.85 g)"

# 4. Hostile input is refused and changes nothing.
{
  echo '['
  echo '{"command":"solid.box","params":{"length":1,"width":1,"height":1,"body_name":"X"}},'
  run_step apply-material '{"material":"unobtainium"}'
  echo ']'
} > bad.json
if solidcraft-cli run bad.json --results --json-out bad.out --quiet 2>bad.err; then fail "unknown material was accepted"; fi
grep -qs "not in the library" bad.err bad.out || fail "refusal does not say why: $(cat bad.err)"
echo "bad input refused: $(grep -hos 'not in the library[^;]*' bad.err bad.out | head -n 1)"

# 5. The saved design reopens and keeps its materials with no plugin present.
rm -rf gen-design-engine
echo '[{"command":"inspect.measure","params":{}}]' > reopen.json
solidcraft-cli run reopen.json --in all.sc --results --json-out reopen.out --quiet
[ "$(jq '[.[0].bodies[] | select(.material | test("Aluminium|Steel|Stainless|Titanium|printed"))] | length' reopen.out)" = "$n" ] \
  || fail "reopened design lost materials"
echo "saved design reopens without the plugin; all $n materials kept"
echo "PASS"
