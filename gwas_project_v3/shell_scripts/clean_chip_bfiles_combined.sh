#!/bin/bash

# This script performs chip-specific SNP filtering (--geno)
# followed by chip-specific sample filtering (--mind).
# The final output for each chip will be <chip_name>_cleaned.bed/bim/fam.

# Navigate to the parent directory where all chip_* folders are located
# This script assumes it's being run from /faststorage/project/populationgenomics/students/estherhelga/GWAS_project/
# If not, cd to it first.
# cd /faststorage/project/populationgenomics/students/estherhelga/GWAS_project/

# Define chip basenames and their specific --geno (SNP missingness) thresholds
declare -A chip_geno_thresholds
chip_geno_thresholds["chip_HTS_iSelect_HD"]="0.05"
chip_geno_thresholds["chip_Illumina_GSAs"]="0.05"
chip_geno_thresholds["chip_OmniExpress"]="0.05"
chip_geno_thresholds["chip_OmniExpress_plus"]="0.05"
chip_geno_thresholds["chip_unknown"]="0.20"

# Define chip basenames and their specific --mind (sample missingness) thresholds
declare -A chip_mind_thresholds
chip_mind_thresholds["chip_HTS_iSelect_HD"]="0.05"
chip_mind_thresholds["chip_Illumina_GSAs"]="0.05"
chip_mind_thresholds["chip_OmniExpress"]="0.05"
chip_mind_thresholds["chip_OmniExpress_plus"]="0.05"
chip_mind_thresholds["chip_unknown"]="0.1" # Maybe more lenient for unknown if quality is expected to be lower

# Default thresholds if a chip is not listed above (optional, can be removed if all chips are always listed)
DEFAULT_GENO_THRESHOLD="0.05"
DEFAULT_MIND_THRESHOLD="0.02"

echo "Starting chip-specific SNP (--geno) and Sample (--mind) filtering..."

# Iterate through the chip names (assuming geno_thresholds has all chip names)
for CHIP_BASE in "${!chip_geno_thresholds[@]}"; do
    echo "Processing: ${CHIP_BASE}"

    CURRENT_GENO_THRESHOLD="${chip_geno_thresholds[$CHIP_BASE]:-$DEFAULT_GENO_THRESHOLD}"
    CURRENT_MIND_THRESHOLD="${chip_mind_thresholds[$CHIP_BASE]:-$DEFAULT_MIND_THRESHOLD}"

    # Intermediate filename for after --geno filtering
    GENO_FILTERED_OUT="${CHIP_BASE}_geno_filtered"
    # Final filename after --mind filtering
    FINAL_CLEANED_OUT="${CHIP_BASE}_cleaned"

    # Check if the chip-specific directory exists
    if [ -d "${CHIP_BASE}" ]; then
        # Navigate into the chip directory
        cd "${CHIP_BASE}" || { echo "ERROR: Failed to cd into ${CHIP_BASE}. Exiting."; exit 1; }

        # Check if the initial input bfiles exist
        if [ -f "${CHIP_BASE}.bed" ] && [ -f "${CHIP_BASE}.bim" ] && [ -f "${CHIP_BASE}.fam" ]; then
            echo "  Input files found: ${CHIP_BASE}.bed/bim/fam"
            
            snp_count_initial=$(wc -l < "${CHIP_BASE}.bim")
            sample_count_initial=$(wc -l < "${CHIP_BASE}.fam")
            echo "  Initial SNPs: ${snp_count_initial}, Initial Samples: ${sample_count_initial}"

            # Step 1: Apply --geno (SNP missingness filter)
            echo "  Applying --geno ${CURRENT_GENO_THRESHOLD}..."
            plink --bfile "${CHIP_BASE}" \
                  --geno "${CURRENT_GENO_THRESHOLD}" \
                  --make-bed \
                  --out "${GENO_FILTERED_OUT}"

            if [ ! -f "${GENO_FILTERED_OUT}.bed" ]; then
                echo "  ERROR: Failed to create ${GENO_FILTERED_OUT}.bed. Check PLINK log. Skipping further steps for ${CHIP_BASE}."
                cd ..
                continue # Skip to the next chip
            fi
            
            snp_count_after_geno=$(wc -l < "${GENO_FILTERED_OUT}.bim")
            sample_count_after_geno=$(wc -l < "${GENO_FILTERED_OUT}.fam") # Should be same as initial samples at this stage
            echo "  SNPs after --geno: ${snp_count_after_geno} (Removed: $((snp_count_initial - snp_count_after_geno)))"
            echo "  Samples after --geno: ${sample_count_after_geno}"


            # Step 2: Apply --mind (Sample missingness filter) to the output of Step 1
            echo "  Applying --mind ${CURRENT_MIND_THRESHOLD} to ${GENO_FILTERED_OUT}..."
            plink --bfile "${GENO_FILTERED_OUT}" \
                  --mind "${CURRENT_MIND_THRESHOLD}" \
                  --make-bed \
                  --out "${FINAL_CLEANED_OUT}"

            if [ ! -f "${FINAL_CLEANED_OUT}.bed" ]; then
                echo "  ERROR: Failed to create ${FINAL_CLEANED_OUT}.bed. Check PLINK log. Skipping cleanup for ${CHIP_BASE}."
                # Attempt to clean up intermediate files even on error if the _geno_filtered files exist
                if [ -f "${GENO_FILTERED_OUT}.bed" ]; then
                    echo "  Cleaning up intermediate files: ${GENO_FILTERED_OUT}.*"
                    rm -f "${GENO_FILTERED_OUT}.bed" "${GENO_FILTERED_OUT}.bim" "${GENO_FILTERED_OUT}.fam" "${GENO_FILTERED_OUT}.log"
                fi
                cd ..
                continue # Skip to the next chip
            fi

            snp_count_final=$(wc -l < "${FINAL_CLEANED_OUT}.bim") # Should be same as after geno
            sample_count_final=$(wc -l < "${FINAL_CLEANED_OUT}.fam")
            echo "  SNPs after --mind (in ${FINAL_CLEANED_OUT}.bim): ${snp_count_final}"
            echo "  Samples after --mind (in ${FINAL_CLEANED_OUT}.fam): ${sample_count_final} (Removed: $((sample_count_after_geno - sample_count_final)))"
            echo "  Successfully created ${FINAL_CLEANED_OUT}.bed/bim/fam"

            # Optional: Clean up intermediate files from the --geno step
            echo "  Cleaning up intermediate files: ${GENO_FILTERED_OUT}.*"
            rm -f "${GENO_FILTERED_OUT}.bed" "${GENO_FILTERED_OUT}.bim" "${GENO_FILTERED_OUT}.fam" "${GENO_FILTERED_OUT}.log"

        else
            echo "  WARNING: Initial input files for ${CHIP_BASE} not found in directory ${CHIP_BASE}. Skipping."
        fi
        
        # Navigate back to the parent directory
        cd ..
    else
        echo "WARNING: Directory ${CHIP_BASE} not found. Skipping."
    fi
    echo "------------------------------------"
done

echo "Chip-specific SNP and Sample filtering complete."
echo "Final cleaned files are named <chip_name>_cleaned.bed/bim/fam within their respective chip folders."
echo "These are now ready to be listed in your merge_list.txt for the next merging step."