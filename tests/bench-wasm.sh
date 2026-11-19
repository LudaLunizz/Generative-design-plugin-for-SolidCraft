#!/usr/bin/env bash
# Times the solver-benchmark entry inside the installed SolidCraft (WebAssembly, headless CLI)
# and the same runs natively, and finds the fuel one run needs by bisection.
# Usage: tests/bench-wasm.sh   (after ./build.sh); CASES="nx ny its solver;..." CG_TOL=1e-4 PARTS=1
set -euo pipefail
cd "$(dirname "$0")/.."
root=$(pwd)
pkg=$root/target/package/gen-design-engine
[ -f "$pkg/plugin.wasm" ] || { echo "run ./build.sh first" >&2; exit 1; }
cargo build -q --release --example bench
work=$(mktemp -d /tmp/gde-bench.XXXXXX)
trap '[ "${KEEP:-}" ] || rm -rf "$work"' EXIT
cd "$work"

# Copy the package with a given fuel limit and print its approval hash.
package() { # fuel
  rm -rf gen-design-engine && cp -r "$pkg" gen-design-engine
  jq -c --argjson f "$1" '.runtime.limits.fuel = $f' "$pkg/plugin.json" > gen-design-engine/plugin.json
  echo '[{"command":"plugin.inspect","params":{"path":"gen-design-engine"}}]' > i.json
  solidcraft-cli run i.json --results --json-out i.out --quiet
  jq -r '.[0].approval' i.out
}

# Run one benchmark; print "ms<TAB>status<TAB>result-json".
run() { # approval nx ny iterations [solver]
  jq -nc --arg a "$1" --argjson nx "$2" --argjson ny "$3" --argjson it "$4" --arg s "${5:-multigrid}" --argjson tol "${CG_TOL:-1e-6}" \
    '[{command:"plugin.run",params:{path:"gen-design-engine",approval:$a,entry:"solver-benchmark",params:{nx:$nx,ny:$ny,iterations:$it,solver:$s,cg_tol:$tol,picture:false}}}]' > r.json
  if solidcraft-cli run r.json --results --timings --json-out r.out --quiet 2> r.err; then
    printf '%s\t%s\t%s\n' "$(jq '.timings[0].ms' r.out)" ok "$(jq -c '.results[0].result // .results[0]' r.out)"
  else
    printf '%s\t%s\t%s\n' "-" "failed" "$(tr '\n' ' ' < r.err | cut -c1-200)"
  fi
}

echo "== Time: WebAssembly in SolidCraft vs native (fuel 1e11, deadline 30 s)"
approval=$(package 100000000000)
printf '%-10s %5s %-9s %10s %10s %7s %8s %s\n' grid its solver wasm_ms native_ms ratio cg_total status
# CASES: "nx ny iterations solver" entries separated by ";"
IFS=';' read -ra cases <<< "${CASES:-60 30 20 jacobi;120 60 5 jacobi;60 30 100 multigrid;120 60 100 multigrid;160 80 100 multigrid;240 120 100 multigrid}"
for case in "${cases[@]}"; do
  set -- $case
  IFS=$'\t' read -r ms status result < <(run "$approval" "$1" "$2" "$3" "$4")
  native=$("$root/target/release/examples/bench" "$1" "$2" "$3" cantilever "$4" | head -n 1 | sed -n 's/.*: \([0-9]*\) ms native.*/\1/p')
  cg=$(echo "$result" | jq -r '.cg_iterations_total // "-"' 2>/dev/null || echo -)
  ratio=$( [ "$ms" != "-" ] && python3 -I -c "print(f'{$ms/max($native,1):.2f}')" || echo -)
  printf '%-10s %5s %-9s %10s %10s %7s %8s %s\n' "${1}x${2}" "$3" "$4" "$ms" "$native" "$ratio" "$cg" "$status"
  [ "$status" = ok ] || echo "    $result"
done

FUEL_SOLVER=${FUEL_SOLVER:-multigrid}
echo "== Fuel for 60x30, 5 iterations, $FUEL_SOLVER (bisection)"
lo=1000000; hi=100000000000
products=0
for _ in $(seq 1 14); do
  mid=$(python3 -I -c "import math; print(int(math.sqrt($lo*$hi)))")
  a=$(package "$mid")
  IFS=$'\t' read -r ms status result < <(run "$a" 60 30 5 "$FUEL_SOLVER")
  if [ "$status" = ok ]; then hi=$mid; products=$(echo "$result" | jq '.element_products'); else lo=$mid; fi
done
echo "fuel needed: between $lo and $hi; element products: $products; fuel per element product: $(python3 -I -c "print(round($hi/$products))")"

# PARTS=1: fuel of each part on its own (fuel per call = fuel(repeat=R) / R, measured at two
# repeat counts so the fixed start-up cost drops out).
if [ "${PARTS:-}" ]; then
  fuel_of() { # part repeat
    local lo=100000 hi=100000000000 mid a
    for _ in $(seq 1 16); do
      mid=$(python3 -I -c "import math; print(int(math.sqrt($lo*$hi)))")
      a=$(package "$mid")
      jq -nc --arg a "$a" --arg part "$1" --argjson r "$2" \
        '[{command:"plugin.run",params:{path:"gen-design-engine",approval:$a,entry:"solver-benchmark",params:{nx:120,ny:60,part:$part,repeat:$r}}}]' > p.json
      if solidcraft-cli run p.json --results --json-out p.out --quiet 2>/dev/null; then hi=$mid; else lo=$mid; fi
    done
    echo "$hi"
  }
  echo "== Fuel per part at 120x60 (one call)"
  for part in apply vcycle setup filter; do
    case $part in apply|filter) r1=20; r2=60 ;; vcycle) r1=5; r2=15 ;; setup) r1=2; r2=6 ;; esac
    f1=$(fuel_of $part $r1); f2=$(fuel_of $part $r2)
    python3 -I -c "per=($f2-$f1)/($r2-$r1); print(f'{\"$part\":8} {per:14,.0f} fuel per call  ({per/7200:,.0f} per cell)')"
  done
fi
