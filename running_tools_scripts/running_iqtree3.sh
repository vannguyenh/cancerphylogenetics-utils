#!/usr/bin/env bash
set -euo pipefail

# --- Tools & opts ---
IQTREE="${IQTREE:-$HOME/tools/build-iqtree3/iqtree3}"
MODEL="${MODEL:-GT10+FO}"
THREADS=1

# --- Fixed base path to the input of simulation data ---
ROOT="$HOME/Documents/promotion/projects/cancer_models/benchmark_gt10/simulation_data_CellPhy"

SCEN="${1%/}"
MODEL="${2:-$MODEL}"
THREADS=1 #"{3:-$THREADS}"

INPUT_DIR="$ROOT/$SCEN/true_haplotypes_dir"
OUTPUT_DIR="$ROOT/$SCEN/iqtree3"

# ---------------------------------
# Validatee iqtree3 and INPUT_DIR

[[ -x "$IQTREE" ]] || { echo "Error: iqtree3 not found at $IQTREE" >&2; exit 1; }

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

    prefix="${OUTPUT_DIR}/${base}.iqtree3.gt10FO" #e.g., true_hap.0042.iqtree3.gt10FO
    
    # Skip if run already finished
    if [[ -e "${prefix}.iqtree" ]]; then
        echo "Skipping $base (found ${prefix}.iqtree)"
        continue
    fi
    
    echo "Running: $base | seed=$seed | prefix=$prefix"
    "$IQTREE" -s "$f" -m "$MODEL" -seed "$seed" --prefix "$prefix" -T "$THREADS"
done

echo "Done!"
