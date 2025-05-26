############################################################
# 0.  ONE-TIME ENVIRONMENT & PROJECT SET-UP
############################################################
conda create -n gwas \
  -c bioconda -c conda-forge \
  plink plink2 vcftools bcftools samtools gcta admixture \
  python=3.10 numpy pandas matplotlib seaborn scikit-learn jupyter \
  r-base r-ggplot2 r-data.table r-dplyr r-optparse r-ggpubr r-readr

cd /faststorage/project/populationgenomics/students/estherhelga/GWAS_project/

ln -s  ~/populationgenomics/project_data/GWAS
cp ~/populationgenomics/project_data/GWAS/{eye_color.txt,metadata.txt,gwas_data.{fam,bim,bed}} .

############################################################
# 1.  PHENOTYPE EXPLORATION (shell helpers only)
############################################################
awk '{print $2}' eye_color.txt | sort | uniq -c | sort -nr    # category counts

############################################################
# 2.  CHIP-BASED EXTRACTION
############################################################
# (for each chip in {HTS_iSelect_HD, Illumina_GSAs, OmniExpress, OmniExpress_plus, unknown})
plink --bfile gwas_data \
      --keep chip_group_<CHIP>.txt \
      --make-bed \
      --out  chip_<CHIP>

############################################################
# 3.  CHIP-SPECIFIC QC  (inside clean_chip_bfiles.sh)
############################################################
# Known chips
plink --bfile chip_<CHIP> \
      --geno 0.05  --make-bed --out chip_<CHIP>_geno05
plink --bfile chip_<CHIP>_geno05 \
      --mind 0.02 --make-bed --out chip_<CHIP>_cleaned

# Unknown-chip group
plink --bfile chip_unknown \
      --geno 0.20 --make-bed --out chip_unknown_geno20
plink --bfile chip_unknown_geno20 \
      --mind 0.05 --make-bed --out chip_unknown_cleaned

############################################################
# 4.  MERGE CLEANED CHIP DATA  (OmniExpress excluded)
############################################################
echo "chip_HTS_iSelect_HD_cleaned"   >  merge_list.txt
echo "chip_Illumina_GSAs_cleaned"    >> merge_list.txt
echo "chip_OmniExpress_plus_cleaned" >> merge_list.txt
echo "chip_unknown_cleaned"          >> merge_list.txt

plink --merge-list merge_list.txt \
      --make-bed \
      --out  allchips_qc1_merged

############################################################
# 5.  KEEP ONLY ID’s WITH BOTH GENOTYPE & PHENOTYPE
############################################################
awk '{print $2}' gwas_data.fam  | sort -u > gwas_iids.txt
awk '{print $1}' eye_color.txt  | sort -u > pheno_iids.txt
comm -12 gwas_iids.txt pheno_iids.txt          > iids_common_to_both.txt
awk 'NR==FNR{a[$1];next} $2 in a{print $1,$2}' \
        iids_common_to_both.txt gwas_data.fam  > individuals_to_keep_for_gwas.txt

plink --bfile allchips_qc1_merged \
      --keep individuals_to_keep_for_gwas.txt \
      --make-bed \
      --out  merged_data_core_cohort

############################################################
# 6.  HETEROZYGOUS-HAPLOID → MISSING
############################################################
plink --bfile merged_data_core_cohort \
      --set-hh-missing \
      --make-bed \
      --out  merged_data_core_cohort_hhcleaned

############################################################
# 7.  SEX CHECK & FILTER
############################################################
plink --bfile merged_data_core_cohort_hhcleaned \
      --check-sex \
      --out  core_cohort_sex_check_results

# (remove 16 ambiguous samples; list created in R)
plink --bfile merged_data_core_cohort_hhcleaned \
      --remove core_cohort_individuals_to_remove_sex_issues.txt \
      --make-bed \
      --out  merged_data_core_cohort_sexcleaned

############################################################
# 8.  FIRST LD-PRUNE → PCA → OUTLIER REMOVAL
############################################################
plink --bfile merged_data_core_cohort_sexcleaned \
      --exclude range high_ld_regions_exclude.txt \
      --make-bed \
      --out merged_LD_excluded

plink --bfile merged_LD_excluded \
      --indep-pairwise 50 5 0.2 \
      --out pruning_results

plink --bfile merged_LD_excluded \
      --extract pruning_results.prune.in \
      --make-bed \
      --out merged_data_final_pruned_for_pca

plink --bfile merged_data_final_pruned_for_pca \
      --pca 10 \
      --out final_pca_results

# (remove 21 PCA outliers)
plink --bfile merged_data_core_cohort_sexcleaned \
      --remove pca_outliers_to_remove.txt \
      --make-bed \
      --out merged_data_core_cohort_no_pca_outliers

############################################################
# 9.  SECOND LD-PRUNE → SECOND PCA  (post-outlier)
############################################################
plink --bfile merged_data_core_cohort_no_pca_outliers \
      --exclude range high_ld_regions_exclude.txt \
      --make-bed \
      --out no_outliers_LD_excluded

plink --bfile no_outliers_LD_excluded \
      --indep-pairwise 50 5 0.2 \
      --out no_outliers_pruning

plink --bfile no_outliers_LD_excluded \
      --extract no_outliers_pruning.prune.in \
      --make-bed \
      --out merged_data_final_pruned_no_outliers_for_pca

plink --bfile merged_data_final_pruned_no_outliers_for_pca \
      --pca 30 \
      --out final_pca_results_no_outliers

############################################################
# 10.  RELATEDNESS (PI_HAT) & REMOVE CLOSE RELATIVES
############################################################
plink --bfile merged_data_final_pruned_no_outliers_for_pca \
      --genome \
      --out relatedness_check_results_no_outliers

plink --bfile merged_data_core_cohort_no_pca_outliers \
      --remove related_individuals_to_remove.txt \
      --make-bed \
      --out final_analysis_cohort_unrelated

############################################################
# 11.  SNP-LEVEL QC  (MAF + HWE)
############################################################
plink --bfile final_analysis_cohort_unrelated \
      --maf 0.01 \
      --make-bed \
      --out final_analysis_cohort_maf_filtered   # –356 459 SNPs

plink --bfile final_analysis_cohort_maf_filtered \
      --hwe 1e-6 \
      --make-bed \
      --out final_analysis_cohort_hwe_filtered   # –414 SNPs

############################################################
# 12.  GWAS RUNS
############################################################
# Initial (4 PCs + sex)
plink --bfile final_analysis_cohort_hwe_filtered \
      --pheno gwas_pheno_1196.txt \
      --covar gwas_covariates_4pcs_sex_1194.txt \
      --linear \
      --allow-no-sex \
      --out gwas_eyecolor_linear_4pcs_sex

# Re-run with 7 PCs + sex (after re-PCA)
plink --bfile final_analysis_cohort_hwe_filtered \
      --pheno gwas_pheno_1196.txt \
      --covar gwas_covariates_7pcs_sex_1194.txt \
      --linear \
      --allow-no-sex \
      --out gwas_eyecolor_linear_correct_7pcs_sex

# Core-cluster (1124) GWAS with 6 PCs + sex
plink --bfile gwas_core_1124_finalqc \
      --pheno gwas_core_pheno_1124.txt \
      --covar gwas_core_covar_6pcs_sex_1122.txt \
      --linear \
      --allow-no-sex \
      --out gwas_eyecolor_CORE_SUBGROUP_6pcs_sex

# Super-core GWAS (example with 4 PCs + sex)
plink --bfile gwas_super_core_N1048_finalqc \
      --pheno gwas_super_core_pheno_N1048.txt \
      --covar gwas_super_core_covar_4pcs_sex_N1046.txt \
      --linear \
      --allow-no-sex \
      --out gwas_eyecolor_SUPERCORE_4pcs_sex

############################################################
# 13.  CONDITIONAL GWAS  (rs12913832 as covariate)
############################################################
plink --bfile final_analysis_cohort_hwe_filtered \
      --pheno gwas_pheno_1196.txt \
      --covar gwas_covariates_7pcs_sex_1194.txt \
      --condition rs12913832 \
      --linear \
      --allow-no-sex \
      --out gwas_cond_rs12913832

############################################################
# 14.  TOP-SNP GENOTYPE EXTRACTION  (for phenotype plot)
############################################################
echo "rs12913832" > top_snp_for_recode.txt

plink --bfile gwas_core_1124 \
      --snp rs12913832 \
      --recode A \
      --out top_snp_genotypes_for_supercore_analysis





