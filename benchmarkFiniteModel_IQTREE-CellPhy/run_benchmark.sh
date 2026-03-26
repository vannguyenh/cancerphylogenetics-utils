#!/bin/bash
################################################################################
# Full benchmark pipeline: IQ-TREE vs CellPhy genotype models
#
# Steps:
#   1. Generate CellCoal parameter files (9 scenarios)
#   2. Run CellCoal simulations
#   3. Convert output data for GT16 (IQ-TREE and CellPhy encodings)
#   4. Run tree inference (IQ-TREE GT16, IQ-TREE GT10, CellPhy GT16, CellPhy GT10)
#   5. Evaluate inferred trees against true trees
#
# Usage:
#   bash run_benchmark.sh          # Run all steps
#   bash run_benchmark.sh step3    # Run from step 3 onwards
################################################################################

set -euo pipefail

# ============================================================
# TOOL PATHS — adjust these to your system
# ============================================================
CELLCOAL="/Users/u7875558/tools/cellcoal/bin/cellcoal-1.2.0"
IQTREE="/Users/u7875558/tools/build-iqtree3/iqtree3"
CELLPHY="/Users/u7875558/tools/cellphy/cellphy.sh"
GT16_MAP="/Users/u7875558/tools/cellphy/bin/GT16.map"

# Working directory (where this script lives)
BENCHMARK_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$BENCHMARK_DIR"

# Number of replicates (must match generate_params.py)
NUM_REPS=20

# Scenarios
SCENARIOS=(
    "S1_ADO0.00_ERR0.00"
    "S2_ADO0.00_ERR0.05"
    "S3_ADO0.00_ERR0.10"
    "S4_ADO0.10_ERR0.00"
    "S5_ADO0.10_ERR0.05"
    "S6_ADO0.10_ERR0.10"
    "S7_ADO0.25_ERR0.00"
    "S8_ADO0.25_ERR0.05"
    "S9_ADO0.25_ERR0.10"
)

# Start from a specific step? (default: step1)
START_STEP="${1:-step1}"

# ============================================================
# STEP 1: Generate CellCoal parameter files
# ============================================================
if [[ "$START_STEP" == "step1" ]]; then
    echo "========================================"
    echo "STEP 1: Generating parameter files"
    echo "========================================"
    python3 generate_params.py
    echo ""
    START_STEP="step2"
fi

# ============================================================
# STEP 2: Run CellCoal simulations
# ============================================================
if [[ "$START_STEP" == "step2" ]]; then
    echo "========================================"
    echo "STEP 2: Running CellCoal simulations"
    echo "========================================"
    for SCENARIO in "${SCENARIOS[@]}"; do
        RESULTS_DIR="results_${SCENARIO}"
        PARAM_FILE="params/params_${SCENARIO}.txt"

        if [[ -d "$RESULTS_DIR" ]]; then
            echo "SKIP: $RESULTS_DIR already exists"
            continue
        fi

        echo "Running CellCoal for $SCENARIO ..."
        $CELLCOAL -F"$PARAM_FILE"
        echo "Done: $RESULTS_DIR"
    done
    echo ""
    START_STEP="step3"
fi

# ============================================================
# STEP 3: Convert CellCoal output for GT16
# ============================================================
if [[ "$START_STEP" == "step3" ]]; then
    echo "========================================"
    echo "STEP 3: Converting data for GT16"
    echo "========================================"
    for SCENARIO in "${SCENARIOS[@]}"; do
        RESULTS_DIR="results_${SCENARIO}"

        for REP in $(seq 1 $NUM_REPS); do
            for SITES in snv full; do
                # IQ-TREE GT16 encoding (M,R,W,S,Y,K for heterozygotes)
                python3 convert_cellcoal_to_phylip.py \
                    --input_dir "$RESULTS_DIR" \
                    --replicate "$REP" \
                    --tool iqtree \
                    --sites "$SITES"

                # CellPhy GT16 encoding (1,2,3,4,5,6 for heterozygotes)
                python3 convert_cellcoal_to_phylip.py \
                    --input_dir "$RESULTS_DIR" \
                    --replicate "$REP" \
                    --tool cellphy \
                    --sites "$SITES"
            done
        done
        echo "Converted $SCENARIO ($NUM_REPS replicates, both encodings, snv+full)"
    done
    echo ""
    START_STEP="step4"
fi

# ============================================================
# STEP 4: Run tree inference
# ============================================================
if [[ "$START_STEP" == "step4" ]]; then
    echo "========================================"
    echo "STEP 4: Running tree inference"
    echo "========================================"
    for SCENARIO in "${SCENARIOS[@]}"; do
        RESULTS_DIR="results_${SCENARIO}"

        for REP in $(seq 1 $NUM_REPS); do
            REP_FMT=$(printf "%04d" "$REP")

            for SITES in snv full; do
                # --- Input files ---
                # IQ-TREE GT16: uses converted phylip with 16-symbol encoding
                GT16_IQTREE="${RESULTS_DIR}/phylip_gt16_iqtree_${SITES}/rep${REP_FMT}_GT16_iqtree_${SITES}.phy"

                # CellPhy GT16: uses converted phylip with CellPhy encoding + GT16.map
                GT16_CELLPHY="${RESULTS_DIR}/phylip_gt16_cellphy_${SITES}/rep${REP_FMT}_GT16_cellphy_${SITES}.phy"

                # GT10 (IQ-TREE only): uses snv_hap/full_hap directly (IUPAC encoded)
                if [[ "$SITES" == "snv" ]]; then
                    # snv_hap files have a site indices line (line 2) that must
                    # be removed for PHYLIP compatibility
                    SNV_HAP_RAW="${RESULTS_DIR}/snv_haplotypes_dir/snv_hap.${REP_FMT}"
                    GT10_CLEAN_DIR="${RESULTS_DIR}/phylip_gt10_${SITES}"
                    mkdir -p "$GT10_CLEAN_DIR"
                    GT10_INPUT="${GT10_CLEAN_DIR}/rep${REP_FMT}_GT10_${SITES}.phy"
                    if [[ ! -f "$GT10_INPUT" ]]; then
                        # Keep header (line 1) and sequences (line 3+), skip site indices (line 2)
                        sed '2d' "$SNV_HAP_RAW" > "$GT10_INPUT"
                    fi
                else
                    GT10_INPUT="${RESULTS_DIR}/full_haplotypes_dir/full_hap.${REP_FMT}"
                fi

                # --- Output directories ---
                IQTREE_GT16_DIR="${RESULTS_DIR}/inference_iqtree_gt16_${SITES}"
                IQTREE_GT10_DIR="${RESULTS_DIR}/inference_iqtree_gt10_${SITES}"
                CELLPHY_GT16_DIR="${RESULTS_DIR}/inference_cellphy_gt16_${SITES}"
                CELLPHY_GT10_DIR="${RESULTS_DIR}/inference_cellphy_gt10_${SITES}"
                mkdir -p "$IQTREE_GT16_DIR" "$IQTREE_GT10_DIR" "$CELLPHY_GT16_DIR" "$CELLPHY_GT10_DIR"

                # --- IQ-TREE GT16 ---
                IQTREE_GT16_PREFIX="${IQTREE_GT16_DIR}/rep${REP_FMT}"
                if [[ ! -f "${IQTREE_GT16_PREFIX}.treefile" ]]; then
                    echo "IQ-TREE GT16 ($SITES): $SCENARIO rep $REP"
                    $IQTREE -s "$GT16_IQTREE" \
                        --seqtype GT \
                        -m GT16+FO+E \
                        --prefix "$IQTREE_GT16_PREFIX" \
                        -nt 1 -quiet
                fi

                # --- IQ-TREE GT10 ---
                IQTREE_GT10_PREFIX="${IQTREE_GT10_DIR}/rep${REP_FMT}"
                if [[ ! -f "${IQTREE_GT10_PREFIX}.treefile" ]]; then
                    echo "IQ-TREE GT10 ($SITES): $SCENARIO rep $REP"
                    $IQTREE -s "$GT10_INPUT" \
                        --seqtype GT \
                        -m GT10+FO+E \
                        --prefix "$IQTREE_GT10_PREFIX" \
                        -nt 1 -quiet
                fi

                # --- CellPhy GT16 ---
                # Use RAXML mode with +M{GT16.map} so raxml-ng knows the
                # CellPhy 28-character GT16 encoding (1-6 for het, !@#$%&^ for rev-het)
                CELLPHY_GT16_PREFIX="${CELLPHY_GT16_DIR}/rep${REP_FMT}"
                if [[ ! -f "${CELLPHY_GT16_PREFIX}.raxml.bestTree" ]]; then
                    echo "CellPhy GT16 ($SITES): $SCENARIO rep $REP"
                    $CELLPHY RAXML --search \
                        --msa "$GT16_CELLPHY" \
                        --model "GT16+FO+E+M{${GT16_MAP}}" \
                        --blopt nr_safe \
                        --prefix "$CELLPHY_GT16_PREFIX" \
                        --threads 1 > /dev/null 2>&1
                fi

                # --- CellPhy GT10 ---
                CELLPHY_GT10_PREFIX="${CELLPHY_GT10_DIR}/rep${REP_FMT}"
                if [[ ! -f "${CELLPHY_GT10_PREFIX}.raxml.bestTree" ]]; then
                    echo "CellPhy GT10 ($SITES): $SCENARIO rep $REP"
                    $CELLPHY RAXML --search \
                        --msa "$GT10_INPUT" \
                        --model GT10+FO+E \
                        --prefix "$CELLPHY_GT10_PREFIX" \
                        --threads 1 > /dev/null 2>&1
                fi
            done

        done
        echo "Finished inference: $SCENARIO"
    done
    echo ""
    START_STEP="step5"
fi

# ============================================================
# STEP 5: Evaluate results
# ============================================================
if [[ "$START_STEP" == "step5" ]]; then
    echo "========================================"
    echo "STEP 5: Evaluating results"
    echo "========================================"
    echo "TODO: Run evaluation script (compute nRF, branch length correlation, etc.)"
    echo "      python3 evaluate_benchmark.py"
    echo ""
    echo "The evaluation script needs to:"
    echo "  - Compare inferred trees against true trees (trees_dir/trees.XXXX)"
    echo "  - Compute normalized Robinson-Foulds distance"
    echo "  - Extract estimated parameters (ADO rate, error rate)"
    echo "  - Measure runtime from log files"
fi

echo ""
echo "========================================"
echo "BENCHMARK PIPELINE COMPLETE"
echo "========================================"