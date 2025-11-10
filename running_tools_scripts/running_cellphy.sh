#!/usr/bin/env bash
set -euo pipefail

# --- Tools & opts ---
CELLPHY="${CELLPHY:-$HOME/tools/cellphy-0.9.3/cellphy.sh}"
MODEL="${MODEL:-GT10}"
THREADS=1

# --- Fixed base path to the input of simulation data ---
ROOT="$HOME/Documents/promotion/projects/cancer_models/benchmark_gt10/simulation_data_CellPhy"

SCEN="${1%/}"
MODEL="${2:-$MODEL}"
THREADS=1 #"{3:-$THREADS}"

INPUT_DIR="$ROOT/$SCEN/true_haplotypes_dir"
OUTPUT_DIR="$ROOT/$SCEN/cellphy0.9.3"

# ---------------------------------
# Validatee iqtree3 and INPUT_DIR

[[ -x "$CELLPHY" ]] || { echo "Error: cellphy not found at $CELLPHY" >&2; exit 1; }

[[ -d "$INPUT_DIR" ]] || { echo "Error: INPUT_DIR not found: $INPUT_DIR" >&2; exit 1; }

# Ensure output dir exists
mkdir -p "$OUTPUT_DIR"

shopt -s nullglob
for f in "$INPUT_DIR"/true_hap.*; do
    [[ -f "$f" ]] || continue
    
    base="$(basename "$f")" # e.g. true_hap.0042
    id="${base#true_hap.}"
    
    #32-bit random seed
    if seed="$(od -An -N4 -tu4 /dev/urandom 2>/dev/null | tr -d ' ')"; then
        : #ok
    else
        seed="$(( (RANDOM<<16) ^ RANDOM ))" #fallback
    fi

    prefix="${OUTPUT_DIR}/${base}.cellphy93.gt10" #e.g., true_hap.0042.iqtree3.gt10FO
    
    # Skip if run already finished
    if [[ -e "${prefix}.iqtree" ]]; then
        echo "Skipping $base (found ${prefix}.raxml.log)"
        continue
    fi
    
    echo "Running: $base | seed=$seed | prefix=$prefix"
    "$CELLPHY" RAXML -msa "$f" --model GT10 --seed "$seed" --prefix "$prefix" --threads "$THREADS"
done

echo "Done!"
