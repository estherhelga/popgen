
# 1: create a new enviroment on the cluster containing the nessecary libraries, python and ect. 

conda create -n gwas \
  -c bioconda -c conda-forge \
  plink plink2 vcftools bcftools samtools gcta admixture \
  python=3.10 numpy pandas matplotlib seaborn scikit-learn jupyter \
  r-base r-ggplot2 r-data.table r-dplyr r-optparse r-ggpubr r-readr

----------------------------------------------------------------------------------

# 2: Make a project folder under your user. This is where we will keep the GWAS data. You can CD into the folder using the following command: 

cd /faststorage/project/populationgenomics/students/estherhelga/GWAS_project/

# 2.1: If you want, you can create a symlink to the data, or copy the data into your folder:

  # symlink
  
  ln -s ~/populationgenomics/project_data/GWAS

  # Specific files for eye color (copy to your project folder)

  cp /home/estherhelga/populationgenomics/project_data/GWAS/eye_color.txt .
  cp /home/estherhelga/populationgenomics/project_data/GWAS/metadata.txt .
  cp /home/estherhelga/populationgenomics/project_data/GWAS/gwas_data.fam .
  cp /home/estherhelga/populationgenomics/project_data/GWAS/gwas_data.bim .
  cp /home/estherhelga/populationgenomics/project_data/GWAS/gwas_data.bed .

# 3: Explore the eyecolor data file.
  
  # How many unique eye color categories are present? We can use the following command:

  awk '{print $2}' eye_color.txt | sort | uniq -c | sort -nr

  # this gives us the following:  
    491 brown
    374 hazel/brown-green
    328 blue
    188 blue-grey
    168 green
    123 blue-green
     92 dark_brown
     39 blue-green-grey
     19 dark_blue
     16 green-gray
     13 amber-brown
      8 blue-green-gold
      2 green_with_brown_freckles
      2 blue/green/grey_-_changes_with_lighting_and_clothing
      2 black
      1 red/blood
      1 hazel_green
      1 blue-green_heterochromia
    

# data analysis and thoughts
  Next task is to actually perform the GWAS on our cleaned data. Essentially, we are now running GWAS for the eyecolour phenotype. We use the --assoc command from PLINK to test to see if we find any eye color assosiated SNPs. 

  However, now GWAS design comes into play. On the one hans we can do a Binary GWAS on light vs dark eyes, or we could do a categorical GWAS study using the 6 categories of colors. 

  The pros of a binary study are:
    - Increased Statistical Power (Potentially): By collapsing categories, you create two larger groups. If the primary genetic signal distinguishes broadly between "lightness" and "darkness" (which is true for major genes like HERC2/OCA2), this approach can maximize your power to detect these major loci.

    - Simplicity: The analysis is straightforward (logistic regression). Interpretation is simpler: variants associated with being "light-eyed" or "dark-eyed."

    - Reduced Misclassification (Potentially): Self-reported distinction between "light" and "dark" might be more reliable than distinguishing between, say, "green" and "hazel" or "blue" and "gray." This can reduce noise.

    - Good for Identifying Major Effect Genes: This approach is well-suited for finding genes with strong, overarching effects on melanin production leading to a general light/dark phenotype.

  However, there are also some cons when it comes to a binary study. My main worry would be oversimlification of biology, since eye color genetics are more complex than a simple binary trait. Grouping distinct phenotypes likem green and blue together might mask specific genetic signals. 

  Other cons are:
    - Loss of Information/Nuance: You lose the ability to detect genes that might specifically influence shades within the light category (e.g., blue vs. green) or within the dark category (e.g., hazel vs. dark brown). Green eyes, for example, have a distinct genetic basis partially separate from just "lightness."

    - Arbitrary Cutoffs/Grouping: The definition of "light" and "dark" and where to place intermediate colors like hazel can be somewhat arbitrary and influence results. If hazel has distinct genetic determinants, lumping it with brown might obscure those.

    - Reduced Power for Specific Modifiers: If a gene modifies blue to green, but doesn't strongly differentiate "light" from "dark" overall, its effect might be diluted and missed.

  
  Then we are the categorical study. Here the pros and cons are a bit different. 

  Pros:
    - Captures More Phenotypic Variation: This approach has the potential to identify genes associated with more subtle variations in eye color (e.g., genes influencing green vs. blue, or hazel vs. brown).

    - More Biologically Granular Insights: If successful, you could understand the genetic architecture at a finer scale.

    - Allows for Ordered Trait Analysis: You could potentially treat this as an ordered categorical trait (e.g., using ordered logistic regression), assuming a continuum of "lightness" or melanin content. This can be more powerful than unordered multinomial logistic regression if the ordering assumption is valid.

  Cons:
    - Reduced Statistical Power per Category: Each category will have a smaller sample size. If some categories are rare, you'll have very little power to detect associations specific to them.

    - Increased Potential for Misclassification: Self-reporting 6 distinct categories is much harder. The distinction between "green" and "hazel," or "blue" and "gray," or "light brown" and "hazel" can be highly subjective and variable between individuals. This noise can significantly reduce power.

    - More Complex Analysis:
      -  Unordered: You might use multinomial logistic regression, comparing each category to a reference (e.g., blue). This increases the number of tests and can be harder to interpret.
      - Ordered: Ordered logistic regression makes assumptions about the equal "distance" or ordering of categories, which might not perfectly reflect the underlying biology.
      - Quantitative (treating 1-6 as a scale): This is another option, but assumes the categories are equally spaced on some underlying quantitative scale of "color," which is a strong assumption.

    - Interpretation Challenges: Pinpointing which specific comparison (e.g., blue vs. green, or hazel vs. brown) is driving a signal can be more complex.

    - Multiple Testing Burden: If you analyze categories pairwise or look at multiple contrasts, you'll need to correct for more tests.


  Just from these pros and cons alone, we need to do some critical thinking in relation to our data and color groupings. 

  For now, lets imagine we exclude the two self reports of blue_green_heterochromia and red/blood, since they seem to be very different from the rest. 

  Our current color scale would look like: 
  Total N ≈ 1865 (This is a decent overall sample size for a GWAS, but category sizes matter)
    0: Blue (N ≈ 535)
        blue (328), dark_blue (19), blue-grey (188)
        Comment: Good. This is a well-defined and large category.
    1: Blue-Green (N ≈ 171) (123+39+8+1(if blue-green heterochromia kept as blue-green)+2(changes))
        blue-green (123), blue-green-grey (39), blue-green-gold (8), (blue-green_heterochromia (1)), (blue/green/grey_-_changes... (2))
        Comment: This is an intermediate category. Its size is okay. The key question is its distinctness from "Green" and "Blue."
    2: Green (N ≈ 184)
        green (168), green-gray (16)
        Comment: Good size. "Green-gray" fits well here.
    3: Hazel (N ≈ 377)
        hazel/brown-green (374), hazel_green (1), green_with_brown_freckles (2)
        Comment: Good size. Dominated by hazel/brown-green. The inclusion of green_with_brown_freckles (N=2) here is key. If it went to "Green," it wouldn't change Ns much but conceptually blurs lines slightly. This definition seems reasonable for hazel.
    4: Brown (N ≈ 504)
        brown (491), amber-brown (13)
        Comment: Good. "Amber-brown" is a good fit.
    5: Dark Brown (N ≈ 94) (92+2 (if black kept))
        dark_brown (92), black (2)
        Comment: This is your smallest category. For a 6-category analysis, N=94 is on the lower side for having strong independent power, especially if comparing it to another specific category

  Option 1: Binary (Light vs. Dark)
    A) Strict Light vs. Dark (Hazel as Dark):
        Light (N ≈ 890): Blue (535) + Blue-Green (171) + Green (184)
        Dark (N ≈ 975): Hazel (377) + Brown (504) + Dark Brown (94)
        Pros: Simple, powerful for major loci, good balance.
        Cons: Loses nuance. blue-green is light, hazel (often with green flecks) is classified dark.
    B) Alternative Binary (Hazel as Intermediate/Lighter):
        Light (N ≈ 1267): Blue (535) + Blue-Green (171) + Green (184) + Hazel (377)
        Dark (N ≈ 598): Brown (504) + Dark Brown (94)
        Pros: "Dark" becomes more unequivocally dark. Hazel often has significant light components.
        Cons: Less balanced. May be too broad a "Light" category. The traditional view often puts hazel as more distinct from clear light colors.
    Recommendation for Binary: Option 1A is more conventional and likely what you initially intended. It provides a strong contrast.
      
  Option 2: 3-Category Model (Blue, Green/Intermediate, Brown/Hazel)
  This is a very common and robust approach.
    1. Blue (N ≈ 535): blue, dark_blue, blue-grey
    2. Green/Intermediate (N ≈ 355): blue-green, blue-green-grey, blue-green-gold, changes..., green, green-gray. (Essentially your category 1 + category 2)
    3. Brown/Hazel (N ≈ 975): hazel/brown-green, hazel_green, green_with_brown_freckles, brown, amber-brown, dark_brown, black. (Essentially your category 3 + 4 + 5)
  Pros: Good sample sizes in each. Clearer distinctions than 6 categories. Captures the three major poles often discussed.
  Cons: Still loses some of the finer detail (e.g., is "blue-green" distinct from "green"? Is "hazel" distinct from "brown"?)

  Option 3: 4-Category Model (e.g., Blue, Green, Hazel, Brown)
  This tries to preserve more key distinctions.
    1. Blue (N ≈ 535): As above.
    2. Green/Blue-Green (N ≈ 355): Your current Blue-Green (171) + Green (184). This combines the "pure green" and "blue-green" which are often genetically related and can be hard to distinguish perfectly by self-report.
    3. Hazel (N ≈ 377): As your current Hazel.
    4. Brown (N ≈ 598): Your current Brown (504) + Dark Brown (94). (Combines all browns).
  Pros: All categories have good N. Keeps Hazel distinct. Simplifies the "green spectrum." Simplifies browns.
  Cons: Does blue-green truly belong with green more than as its own entity or with blue?

  Option 4: 5-Category Model (Refining your 6-category by merging Dark Brown)
  This is very close to your current model but addresses the smallest category.
    0: Blue (N ≈ 535)
    1: Blue-Green (N ≈ 171)
    2: Green (N ≈ 184)
    3: Hazel (N ≈ 377)
    4: Brown (N ≈ 598): brown, amber-brown, dark_brown, black (Your current category 4 + 5)
  Pros: Maintains most of your granularity. All categories have N > 170.
  Cons: Still relies on the self-reported distinction between "Blue-Green" and "Green" being robust. If these are noisy, their separation might not yield distinct genetic signals and could reduce power compared to merging them.



  From all of this I think we will start by going for a 4-category model: 

    | Category             | Self-Reported Colors                                                                                                        | Approximate Count (N) |
    | -------------------- | --------------------------------------------------------------------------------------------------------------------------- | --------------------- |
    | **Blue**             | blue, dark_blue, blue-grey                                                                                                  | ≈ 535                 |
    | **Green/Blue-Green** | blue-green, blue-green-grey, blue-green-gold, blue/green/grey_changes_with_lighting_and_clothing, green, green-gray         | ≈ 355                 |
    | **Hazel**            | hazel/brown-green, hazel_green, green_with_brown_freckles                                                                   | ≈ 377                 |
    | **Brown**            | brown, amber-brown, dark_brown, black                                                                                       | ≈ 600                 |

  # 6.1: Creating the new phenotype file
    The original eye_color.txt file had the self reported colors, so using R, I made a new file called eye_color_updated.txt that has the "IID", "FID", "phenotype". This R script is called 4_model_scale_creation

# 5.2: Metadata exploration

    After exploring the metadata file, We found that: 
      The metadata file contains the following columns:
      user, 	build, 	chip, 	chip_version, 	inferred_sex, 	source

      We have 4 different chips:
      HTS iSelect HD (587 entries),    Illumina GSAs (248 entires),      OmniExpress (291 entries),	 OmniExpress plus (484 entries).

      We have a lot of different companies, some who use the same chips, some who dont, and some who never report their chip. In total, we have 1610 chip reports, and 399 that did not report a chip.

    From this info, we can now group our data by chip, and do QC on each individual chip type. 

  # Chip grouping

    First, we split our main dataset into 4 chip-based groups:

      HTS iSelect HD

      Illumina GSAs

      OmniExpress

      OmniExpress Plus

    Using the split_by_chip.R on the cluster, and the command Rscript to run it, we create 4 individual chip txt files, that we can then use alongside PLINK to --keep and --remove during the individual QC steps. The .txt files contain the Family ID (FID), and the Individual ID (IID), which in our case are the same. But PLINK needs both to operate the --keep and the --remove.

    We might also consider doing one more group, which includs the 399 entries that did not report a chip, and do a QC on that group also. This might tell us something. 

  # Chip specific genotype data

    Using the following command for each chip group, we extract each subgroups genotyping data:
      plink --bfile gwas_data \
      --keep chip_group_<CHIP>.txt \
      --make-bed \
      --out chip_<CHIP>

    For the group of unknown chips, which was made at the same time as the other groups were made, we first rename it, before using the same plink command on it:

    mv chip_group_.txt chip_group_unknown.txt

    plink --bfile gwas_data \
      --keep chip_group_unknown.txt \
      --make-bed \
      --out chip_unknown

  # Prune the datafiles so that they only include the SNPs from each chip. 

    To do this, we make an .sh script that runs trhough our folders and creates new bfiles that use the --geno 0.05. However for the Unknown one, we do it manually with a geno threshold of 0.2.

    Script is called clean_chip_bfiles.sh, and we make it executable using: 
      chmod +x clean_chip_bfiles.sh
    And we run it using:
      ./clean_chip_bfiles.sh

    Filter SNPs by missingness (--geno) for each chip.
Filter samples by missingness (--mind) on the SNP-filtered data for each chip.
The output of this two-step process per chip will then be ready for merging.

Unknown
SNPs after --geno: 146398 (Removed: 1230255)
Samples after --mind (in chip_unknown_cleaned.fam): 303 (Removed: 96)

Illumina_GSAs
SNPs after --geno: 507262 (Removed: 869391)
Samples after --mind (in chip_Illumina_GSAs_cleaned.fam): 239 (Removed: 9)

OmniExpress
SNPs after --geno: 612572 (Removed: 764081)
Samples after --mind (in chip_OmniExpress_cleaned.fam): 280 (Removed: 11)

HTS_iSelect_HD
SNPs after --geno: 569460 (Removed: 807193)
Samples after --mind (in chip_HTS_iSelect_HD_cleaned.fam): 584 (Removed: 3)

OmniExpress_plus
SNPs after --geno: 910557 (Removed: 466096)
Samples after --mind (in chip_OmniExpress_plus_cleaned.fam): 474 (Removed: 10)


# Step 0: Merge the Datasets
We already did sample and marker quality on the chip specific subsets.
Now we need to merge.

Making merge_list.txt

  echo "chip_HTS_iSelect_HD/chip_HTS_iSelect_HD_cleaned" > merge_list.txt
  echo "chip_Illumina_GSAs/chip_Illumina_GSAs_cleaned" >> merge_list.txt
  echo "chip_OmniExpress/chip_OmniExpress_cleaned" >> merge_list.txt
  echo "chip_OmniExpress_plus/chip_OmniExpress_plus_cleaned" >> merge_list.txt
  echo "chip_unknown/chip_unknown_cleaned" >> merge_list.txt

Using PLINK to merge using the txt file

  plink --merge-list merge_list.txt \
        --make-bed \
        --out merged_data_step0

Number of Individuals and Variants after merge: 1880 people, 1325823 variants.

Command: 
plink --bfile merged_data_step0 --set-hh-missing --make-bed --out merged_data_step0_hhCleaned

Why we are doing this:
  PLINK identified "heterozygous haploid genotypes" in the merged dataset (merged_data_step0.log warning). This means it found genotype calls showing two different alleles (e.g., A/G) on chromosomes where an individual is expected to have only one copy (e.g., non-pseudoautosomal regions of the X chromosome in males, or the Y chromosome). Such calls are biologically implausible and usually represent genotyping errors, issues with pseudo-autosomal region definitions, or mis-sexed samples.
Purpose of the command:
  The --set-hh-missing command specifically targets these erroneous heterozygous haploid genotype calls and converts them to 'missing' (0/0) in the new output dataset (merged_data_step0_hhCleaned). The purpose is to clean these likely genotyping artifacts from the data, preventing them from skewing subsequent quality control calculations (like allele frequencies, HWE tests) and the final association analysis, thereby improving the overall data quality and reliability of the results. It does not remove entire individuals or SNPs, but rather corrects specific erroneous genotype calls.


# Sex Inference (Sex Check)
  Purpose:
    To genetically infer the sex of each individual using X chromosome heterozygosity (and Y chromosome call rates if Y SNPs are present and reliable).
    To identify individuals whose genetically inferred sex is ambiguous or conflicts with any recorded sex (though your project description says sex is missing, your merge log showed some sex codes in the .fam).
    To update the .fam file with inferred sexes or flag problematic samples for removal.

  Plink Command:
    plink --bfile merged_data_step0_hhCleaned \
      --check-sex \
      --out sex_check_results

Objective:
  To genetically infer the sex of individuals in the merged dataset (merged_data_step0_hhCleaned), identify discrepancies or ambiguities, and filter samples accordingly to ensure data quality for downstream GWAS of eye color.
Methodology:
  Initial Dataset: The analysis began with the merged dataset merged_data_step0_hhCleaned.bed/bim/fam, which had undergone per-chip quality control (SNP and sample missingness), merging of 5 distinct chip datasets, and correction of heterozygous haploid genotypes (--set-hh-missing). This dataset contained 1,880 individuals and 1,325,823 SNPs.
Preliminary X/Y Chromosome SNP Check:
  A command-line check (awk) was performed to count X and Y chromosome SNPs in merged_data_step0_hhCleaned.bim. This revealed 32,118 X-chromosome SNPs and a substantial number of Y-chromosome SNPs, indicating sufficient markers were present for an initial sex check.
Genetic Sex Inference (--check-sex):
  PLINK's --check-sex command was run on merged_data_step0_hhCleaned to infer sex based on X chromosome heterozygosity (F-statistic) and Y chromosome calls.
    Output file: sex_check_results_on_pre_global_geno.sexcheck.
Data Integration in R:
  The .sexcheck output was loaded into R.
  A mapping file (iid_to_chip.txt), previously generated to link individual IDs (IIDs) to their original genotyping chip platform, was also loaded.
  These two datasets were merged in R (sexcheck_data_merged) to allow for chip-specific analysis of sex check results.
Analysis of PLINK's "PROBLEM" Status by Chip Type:
  Initial examination of PLINK's STATUS column revealed 301 samples flagged as "PROBLEM".
  Analysis of these "PROBLEM" samples by chip type showed a striking pattern: all 280 samples genotyped on the chip_OmniExpress platform were flagged as "PROBLEM". This accounted for the vast majority of "PROBLEM" samples.
  Other chips (chip_HTS_iSelect_HD, chip_Illumina_GSAs, chip_OmniExpress_plus, chip_unknown) had very few "PROBLEM" samples (4, 1, 2, and 14 respectively).
Investigation of OmniExpress Chip Data:
  X/Y SNP Content on OmniExpress: A crucial diagnostic step involved checking the X and Y chromosome SNP content specifically within the pre-merge, per-chip cleaned file for OmniExpress (chip_OmniExpress/chip_OmniExpress_cleaned.bim). This revealed 0 X-chromosome SNPs and 0 Y-chromosome SNPs on this specific chip dataset. For comparison, other chips like chip_HTS_iSelect_HD contained >17,000 X-SNPs and >2,000 Y-SNPs.
  R Analysis of OmniExpress F and SNPSEX: For all 280 OmniExpress samples in the sexcheck_data_merged R dataframe, the F-statistic was NaN (Not a Number) and the PLINK-inferred sex (SNPSEX) was 0 (undetermined).
Conclusion for OmniExpress Samples:
  The absence of X and Y chromosome SNPs on the chip_OmniExpress platform made it impossible for PLINK to genetically infer sex for these 280 individuals. This was the direct cause of their "PROBLEM" status and F=NaN/SNPSEX=0 results.
Strategy for OmniExpress Samples:
  Given that the target phenotype (eye color) is primarily autosomal and not strongly sex-dependent, a decision was made to retain these 280 OmniExpress samples in the dataset.
  Their sex will be treated as "unknown" (typically coded '0' in the .fam file).
  The PLINK flag --allow-no-sex will be used during the GWAS association testing to ensure their inclusion.
Refined Filtering Criteria for Non-OmniExpress Samples:
  For the remaining individuals (not on chip_OmniExpress), a more nuanced filtering strategy was adopted, focusing on the F-statistic rather than solely on PLINK's "PROBLEM" status (which could be triggered by PEDSEX=0 even if genetic sex was clear).
  F-statistic thresholds were defined (e.g., female_F_max = 0.2, male_F_min = 0.9 – note: actual thresholds used should be based on visual inspection of the F-statistic histogram for non-OmniExpress samples).
  Non-OmniExpress individuals were identified for removal if:
    Their F-statistic was NaN.
    Their F-statistic fell into the ambiguous range (e.g., >0.2 AND <0.9).
  A list of these individuals (FID and IID) was generated and saved to individuals_to_remove_sex_issues.txt.
Key Findings & Outcomes:
  A significant batch effect or chip-specific limitation was identified: the chip_OmniExpress platform data lacked X/Y chromosome SNPs, preventing genetic sex inference for all 280 samples processed on it.
  These 280 OmniExpress samples will be retained with "unknown" sex for the autosomal eye color GWAS, utilizing PLINK's --allow-no-sex option. This preserves sample size while acknowledging the data limitation.
  For non-OmniExpress samples, a targeted removal list (individuals_to_remove_sex_issues.txt) was created to exclude only those with genuinely ambiguous or failed genetic sex determination based on their F-statistics. This avoids unnecessary removal of samples where PLINK's "PROBLEM" status was due to PEDSEX=0 but genetic sex was otherwise clear.
  This iterative QC process, combining PLINK output with R-based diagnostics, allowed for a more informed and tailored approach to handling sex data irregularities in a multi-chip cohort.

  Now, we move on to PLINK where we remove the individuals we marked as unsusable, before moving on to global SNP missingness. 

    Plink command:
      plink --bfile merged_data_step0_hhCleaned \
      --remove individuals_to_remove_sex_issues.txt \
      --make-bed \
      --out merged_data_after_sex_removal



# Global Marker Genotyping Call Rate (SNP missingness)

For this, we used the following command, with the global missingness threshhold of 0.05. However, that only left us with aproxx 62 thousand variants out of the initial 1.3 million, which feels very strict. Thus, we want to try doing some R analysis to see if we can find a more appropreate threshold. 

    plink --bfile merged_data_after_sex_removal \
      --geno 0.05 \
      --make-bed \
      --out merged_data_step1_geno_final 

      Try geno (from R analysis below)
      0.3
      0.5
      0.6

  # R analysis and file generation
    Using PLINK, we generate SNP missingness statistics using:
      plink --bfile merged_data_after_sex_removal \
      --missing \
      --out snp_missing_stats

    We then use this file in R to do some analysis. See global_missingness_analysis.rmd

    After looking at the plots, we are going to use the strict threshold for inital analysis, with the posibility of changing the trheshold at a later time. Some interesting thresholds would be: 
      0.3
      0.5
      0.6


# LD Pruning (Preparation for PCA & Relatedness)
  Objective: 
    To create a subset of SNPs from your current dataset (merged_data_step1_geno_final) that are in approximate linkage equilibrium (LD) with each other. This pruned set is ideal for PCA and relatedness estimation, as these methods assume (or perform better with) SNPs that are not highly correlated due to LD.

  Method: Window-based LD Pruning using PLINK (--indep-pairwise)
    This method works by:
      Calculating LD (typically r²) between pairs of SNPs within a defined window.
      If a pair of SNPs has r² greater than a specified threshold, one SNP from the pair is removed.
      The window then slides along the chromosome.

  Step 1: Excluding Regions of High and Complex LD (Optional but Recommended)
    We are looking to exclude the MHC (found on Chr6), since it is known to shwocase highly complex LD patterns. 
      To do this, we first had to find out the genome build for our data, which using rs positions and the NIH website <https://www.ncbi.nlm.nih.gov/snp/>, we found to be GRCh37. 
      Next, we need to find the specific coordinates for the regions we want to exclude, which to start with is just the MHC on chromosome 6. 
        This time, using NCBI and the Genome Referance Consortium, we see that on our specific build, the MHC coordinates on chromosme 6 are:
          28,477,797 - 33,448,354.
          A commonly used broad range for exclusion is chr6:25000000-35000000, however I am not sure which one we should go with. 

            To check how many SNPs we have in the broad reagion, we can use the following command:
              awk '$1 == "6" && $4 >= 25000000 && $4 <= 35000000 {print $0}' merged_data_step1_geno_final.bim | wc -l

        For now, I think only exlcuding the MHC is fine, however, other areas of note are: 
          The 8p23 inversion on chr8. (chr 8 : 7–13 Mb?)
          chr 11 : 45–57 Mb
          chr 17 : 40–46 Mb

      
    After the analysis of the data, as well as finding the coordinates of the MHC within our specific genome build, we created a textfile called high_ld_regions_exclude.txt, whcih includes the chromosome and range, as well as the MHC label. We need this file for the following PLINK command:

      plink --bfile merged_data_step1_geno_final \
      --exclude range high_ld_regions_exclude.txt \
      --make-bed \
      --out merged_data_step1_geno_final_noMHC
    
    This plink command creates a new set of bfiles, that we then can use for further pruning before checking relatedness and pop structure. 

  Step 2: Perform LD Pruning with --indep-pairwise

    Objective:
      To identify a subset of SNPs from merged_data_step1_geno_final_noMHC that are in approximate linkage equilibrium (i.e., not highly correlated with each other due to physical proximity on the chromosome). This command doesn't create a new dataset directly, but rather generates lists of SNPs to keep or remove.

      plink --bfile merged_data_step1_geno_final_noMHC \
      --indep-pairwise 50 5 0.2 \
      --out ld_pruning_results

        --indep-pairwise 50 5 0.2:
          This is the core flag for window-based LD pruning.
          50: This is the window size in SNPs. PLINK will consider a sliding window of 50 SNPs at a time.
          5: This is the SNP count to shift the window by at each step. After processing a window, it moves 5 SNPs forward to define the next window.
          0.2: This is the r² (squared correlation coefficient) threshold. For any pair of SNPs within the current window, if their r² value is greater than 0.2, one SNP from that pair will be marked for removal. PLINK typically keeps the SNP with the higher MAF from such a pair, or if MAFs are equal, the one that appears first in the .bim file.
          Parameter Choice:
          50 5 0.2 are very common default parameters used in many studies and tutorials. They provide a reasonable balance, aiming to remove strong LD without being overly aggressive.
          You can adjust these:
          Smaller window size (e.g., 20) or larger r² threshold (e.g., 0.5) would result in less pruning (more SNPs kept).
          Larger window size (e.g., 100) or smaller r² threshold (e.g., 0.1) would result in more pruning (fewer SNPs kept).
          For your project, 50 5 0.2 is a perfectly good starting point.

            Principal Component Analysis (PCA) / Admixture / Population Structure: You want SNPs to be as independent as possible. A lower r² (e.g., 0.1, 0.2) is common. 0.2 means you're removing SNPs that share about 20% or more of their variance with another SNP in the window.
            IBD Estimation / Relatedness: Moderate LD can sometimes be tolerated or even informative, but very high LD can skew estimates. 0.2 to 0.5 might be used.

      Looking at the .log file from the pruning results, we can see how many variants were removed from each chromosome. In total over all 22 autosomal chromosome pairs we removed  19087 of 62091 variants.
    

  Step 3: Create a New PLINK Fileset with Only the Pruned-In SNPs 
    
    Objective:
      To create a new PLINK binary fileset (.bed, .bim, .fam) that contains all your current individuals but only the ~43,004 SNPs that were identified as being in approximate linkage equilibrium by the --indep-pairwise command (i.e., those listed in ld_pruning_results.prune.in).

    PLINK command:

      plink --bfile merged_data_step1_geno_final_noMHC \
      --extract ld_pruning_results.prune.in \
      --make-bed \
      --out merged_data_step2_pruned


# Population Stratification (PCA) & Batch Effect Visualization
  We can do these two "together" or at the same time, as they are close together. PCA is the primary tool used to visualize both popuæation stratification and potential batch effects. 
    - Population Stratification: PCA identifies the major axes of genetic variation in your dataset. These axes often correspond to ancestral differences or population substructure.
    - Batch Effect Visualization: By plotting the samples on the principal components (e.g., PC1 vs PC2) and then coloring the points by a suspected batch variable (like 'chip type' in your case), you can visually assess if samples from different batches cluster separately along these PCs. If they do, and this separation isn't easily explained by known population structure, it suggests a batch effect.

    Using PLINK for our calculations: 
      plink --bfile merged_data_step2_pruned \
      --pca 10 \
      --out pca_results_initial


      If we feel like looking at more than the topp 10 PCs, we can change the --pca command. 

      --out pca_results_initial: Basename for the output files. PLINK will create:
        pca_results_initial.eigenvec: This is the most important file. It contains the principal components (eigenvectors) for each individual. It will have columns like: FID, IID, PC1, PC2, ..., PC10.
        pca_results_initial.eigenval: Contains the eigenvalues for each component, which indicate the amount of variance explained by each PC. Useful for creating a scree plot (though not strictly required by your project description for the main plot).
        pca_results_initial.log: Standard PLINK log.


    Next, we load the eigen files into R, where we do the PCA visualiations. Please refer to the following r script for this: 

    During this step, we discovered a crutial data mismatch, whcih we had not found earlier. We had some individuals in the eye_color data that we did not have genotyping data for, as well as individuals we had genotyping data for but no reported phenotype in the eye_color data. These mismatches explained the troubles we had during the R analysis, where a lot of our individuals were coming up as NAs. This was due to the mismatch. 

      50 iids_in_pheno_not_gwas.txt
      191 iids_in_gwas_not_pheno.txt
      1818 iids_common_to_both.txt


  # Updates to the whole pipeline:

    Create the PLINK "keep" file (FID and IID for common individuals)

        awk 'NR==FNR{keep[$1]; next} $2 in keep {print $1, $2}' iids_common_to_both.txt gwas_data.fam > individuals_to_keep_for_gwas.txt

    Filter the Merged Genetic Data to Keep Only Common Individuals

        plink --bfile merged_data_step0_hhCleaned \
      --keep individuals_to_keep_for_gwas.txt \
      --make-bed \
      --out merged_data_core_cohort

    Sex Inference (Sex Check) on the Core Cohort

        plink --bfile merged_data_core_cohort \
      --check-sex \
      --out core_cohort_sex_check_results

      Analysis done in R using the sex_individuals_remove.Rmd, Where we create a new core_cohort_individuals_to_remove_sex_issues.txt file. 

    Remove Individuals with Ambiguous Sex (from Core Cohort)

      plink --bfile merged_data_core_cohort \
      --remove core_cohort_individuals_to_remove_sex_issues.txt \
      --make-bed \
      --out merged_data_core_cohort_sexcleaned

    Global Marker Genotyping Call Rate (SNP missingness) on Sex-Cleaned Core Cohort

      plink --bfile merged_data_core_cohort_sexcleaned \
      --geno 0.05 \
      --make-bed \
      --out merged_data_core_cohort_geno_filtered

    Exclude High LD Regions (e.g., MHC)

      plink --bfile merged_data_core_cohort_geno_filtered \
      --exclude range high_ld_regions_exclude.txt \
      --make-bed \
      --out merged_data_core_cohort_noMHC

    Perform LD Pruning (Identify SNPs to keep)

      plink --bfile merged_data_core_cohort_noMHC \
      --indep-pairwise 50 5 0.2 \
      --out core_cohort_ld_pruning_results

          Pruning complete.  19102 of 62091 variants removed.

    Create the Final LD-Pruned Dataset:

      plink --bfile merged_data_core_cohort_noMHC \
      --extract core_cohort_ld_pruning_results.prune.in \
      --make-bed \
      --out merged_data_final_pruned_for_pca


    Population Stratification (PCA) & Batch Effect Visualization (Run PCA)

      plink --bfile merged_data_final_pruned_for_pca \
      --pca 10 \
      --out final_pca_results

  # PCA analysis
      Based on the PCA plots (especially the one colored by chip type and the one just showing general structure):
      Population Structure: Do the samples form one homogenous cloud, or are there distinct clusters or gradients? If so, what might these correspond to (e.g., European ancestries, other ancestries if present in openSNP data)?
      Batch Effects: Do samples cluster primarily by chip type on any of the top PCs? If chip_OmniExpress (or any other chip) forms a separate cluster from the others, this is strong evidence of a batch effect. You'd mention that the PCA reveals chip-specific clustering.
      Outliers: Are there any individual samples that are extreme outliers on the PCA plot? These might warrant further investigation (though QC should have caught most problematic samples).

        PCA: PC1 vs PC2 (Population Structure - Uncolored Plot):
          Observation: You see a main dense "cloud" of points near the origin (0,0), but also distinct "arms" or gradients extending away, particularly towards the positive PC1 direction, and separating further along PC2 in the upper-right and lower-right quadrants.
          Interpretation (Population Structure): This is a classic picture of significant population substructure.
          PC1 (46.1% variance): This axis captures the largest source of genetic variation among your samples. Given that openSNP data often has a strong European component but also includes individuals of other ancestries, PC1 frequently separates individuals along a major axis of European variation (e.g., North-West Europe vs. South/East Europe) or sometimes separates Europeans from non-Europeans. The large spread along this axis means there are substantial genetic differences captured here.
          PC2 (21.3% variance): This captures the second largest amount of variation. It further differentiates groups. It might separate different European subgroups or distinguish individuals with non-European ancestry (e.g., East Asian, African, Ashkenazi Jewish – whose genetic profiles differ from the main European cluster) from the main European cluster(s).
          The Dense Cloud: Likely represents the most numerous ancestry group in your subset (often Northern/Western European descent in openSNP).
          The "Arms": Represent individuals who are genetically distinct from the main cloud along these primary axes. They could be from different European populations, individuals with non-European ancestry, or potentially admixed individuals.
          Conclusion: Your dataset is not genetically homogeneous. It contains substantial population structure, likely reflecting diverse ancestries present in the openSNP cohort. This structure must be controlled for in your GWAS to avoid spurious associations. Including PC1 and PC2 (and possibly more PCs) as covariates is essential.


        PCA: PC1 vs PC2, Colored by Chip Type:
          Observation: As you noted, the different colors (chip types) are largely intermingled throughout the structure revealed in the first plot. There isn't one color forming its own isolated island separate from the others along PC1 or PC2.
          Interpretation (Batch Effects): This is generally very good news. It indicates that the major population structure patterns (captured by PC1 and PC2) are not primarily driven by which chip the sample was run on. If there were strong, systematic technical differences between chips affecting many of the ~43k common SNPs used for PCA, you would expect to see distinct clusters based purely on color. The mixing suggests that genetic variation (ancestry) is a much stronger signal than chip type for these top PCs after your QC steps.
          Conclusion: While originating from different genotyping platforms, the samples do not show strong clustering by chip type along the primary axes of genetic variation (PC1/PC2) after quality control, suggesting major batch effects are not confounding the observed population structure in this view.

        PCA: PC1 vs PC2, Colored by Mapped Eye Color:
          Observation: You see clear non-random distribution. Blue (0) is concentrated in the main dense cloud. Green/Hazel (1, 2) are also mostly in/near that cloud. Brown (3) is present in the cloud but is much more prevalent in the "arms" extending towards positive PC1. Your observation about darker eyes "stretching out" is accurate. The two NAs are also plotted somewhere based on their genetics.
          Interpretation (Phenotype-Structure Correlation): This strongly indicates that eye color is correlated with the genetic ancestry captured by PCA.
          The main cluster (low PC1/PC2), rich in blue eyes, likely represents populations with high frequencies of alleles associated with blue eyes (e.g., common variants near HERC2/OCA2 prevalent in Europeans, especially Northern Europeans).
          The "arms" (higher PC1), rich in brown eyes, likely represent populations where alleles for brown eyes are much more common (e.g., Southern Europeans, non-Europeans).
          This is biologically expected. Eye color frequencies vary significantly across global populations.
          Conclusion: Eye color phenotype is significantly correlated with the population structure identified by PC1 and PC2. This underscores the critical importance of using PCs as covariates in the GWAS. Failing to do so would lead to highly inflated results and false positives, simply rediscovering that different ancestry groups have different eye colors, rather than finding specific causal variants within those populations.
              
        The PCA clearly reveals substantial population stratification within your openSNP cohort, likely reflecting diverse European and potentially non-European ancestries. The first two principal components capture over 67% of the genetic variance in the LD-pruned dataset. Encouragingly, samples do not cluster strongly by genotyping chip along these main PCs, suggesting major batch effects are not driving the primary structure after QC. However, there is a clear correlation between the observed population structure (PCs) and the eye color phenotype, with blue eyes predominating in the main genetic cluster and brown eyes more prevalent in genetically distinct groups. This highlights the necessity of including principal components (e.g., PC1, PC2, and potentially more) as covariates in the subsequent GWAS to control for confounding due to ancestry.


      Bjarke says that the "elbow shape" we see is usually unwanted. But why? Thats for us to find out. 

      If it is ancestry, it does look pretty similar to these: https://www.nature.com/articles/s41598-021-97129-2/figures/1 

        It’s almost never genuine population structure
          * In a genome‑wide human data set that has been LD‑pruned, the first PC usually explains ≤ 5 % of the variance. In your plot PC 1 already captures ≈ 46 % and PC 2 another ≈ 21 %—far too much. Those huge eigenvalues and the characteristic fan/“elbow” that splits into two arms are classic fingerprints of a single long‑range LD block or a common chromosomal inversion dominating the PCA rather than broad ancestry differences. https://pmc.ncbi.nlm.nih.gov/articles/PMC7750941/

        Why an inversion/long‑range‑LD block gives a ‘fan’?
          SNPs inside the block are in very strong mutual LD, so they behave like one giant marker.

          Individuals that are homozygous for the reference orientation pile up on one side of PC 1, those homozygous for the inverted orientation on the other, and heterozygotes lie in between. When PC 2 also picks up variation from the same region, the three genotype classes stretch into a triangle or V.

          Price et al. 2008 first warned that these regions can swamp genome‑wide PCs and create apparent “structure” where none exists. https://pubmed.ncbi.nlm.nih.gov/18606306/ 

        Why it matters for your GWAS
          If you feed these PCs into your association model you are not “correcting for population stratification”; you are in effect adjusting for (and possibly generating artefactual signals from) one inversion locus.

          Conversely, if you leave the LD block in, variants inside it are very likely to pop up as genome‑wide‑significant hits for eye colour simply because they co‑vary with those PCs, not because they influence the trait.

          So, one by one, we can add new blocks of LD regions or inversions and see if it improves our PCA. 

          So now, when looking at liteature, it seems that a common culprit of a V shape PCA could be due to the 8p23.1 inversion, so looking up the coordinates of that one in our build on ncbi, we find:
              chr8:8093169-11898980 - the usual range used is 8000000 - 12000000. 

          Another block could be the 17q21.31 Inversion Region (MAPT locus) "The locus structure was first described by Stefansson et al., with inversion breakpoints defined at chr17:43,628,944–44,571,603 (GRCh37–hg19 reference assembly) https://pubmed.ncbi.nlm.nih.gov/15654335/ ":
                chr17:43,628,944–44,571,603 - So the range we can use is either the specific one or we rounf down and up - 43000000 - 45000000

          Lactase Persistence Region (LCT gene and flanking regions):
            Chromosome: 2
            Approximate Size: The region of strong selection and LD around LCT can be quite extended (hundreds of kb to over 1Mb).
            Reason: Strong positive selection for lactase persistence in European populations has led to long haplotypes and high LD in this region. While not an inversion, its LD pattern can influence PCA.
            Search: "LCT gene region LD coordinates [your build]" or find the LCT gene and define a generous window around it (e.g., +/- 500kb to 1Mb from the gene's start/end).

            | Chr    | hg19 start–end (Mb) | Why it misbehaves                    | Source                                         |
            | ------ | ------------------- | ------------------------------------ | ---------------------------------------------- |
            | **2**  | **86–100.5**        | *LCT* sweep + segmental duplications | Price et al. 2008 ([reich.hms.harvard.edu][1]) |
            | 2      | 183–190             | Recombination cold spot              | Price et al. 2008 ([PMC][2])                   |
            | **11** | **46–57**           | Extended LD plateau                  | Price et al. 2008 ([PMC][2])                   |
            | **1**  | **48–52**           | Large tandem duplications            | Price et al. 2008 ([PMC][2])                   |


          I am starting by running the following PLINK commands again, after adding the inversion and LD regions from centrer for statistical genetics (https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)). However, This did not change the PCA shape.  
              # Updated exclusion list is in high_ld_regions_exclude.txt

              #: Exclude High LD Regions (MHC + new ones)
              plink --bfile merged_data_core_cohort_geno_filtered \
                    --exclude range high_ld_regions_exclude.txt \
                    --make-bed \
                    --out merged_data_core_cohort_multi_LDexcluded 
                    # New name, e.g., _multi_LDexcluded

              #Perform LD Pruning 
              plink --bfile merged_data_core_cohort_multi_LDexcluded \
                    --indep-pairwise 50 5 0.2 \
                    --out core_cohort_multi_LDexcluded_pruning_results 
                    # New name

              #Create the Final LD-Pruned Dataset
              plink --bfile merged_data_core_cohort_multi_LDexcluded \
                    --extract core_cohort_multi_LDexcluded_pruning_results.prune.in \
                    --make-bed \
                    --out merged_data_final_pruned_for_pca_v2 
                    # New name, e.g., _v2

              #Run PCA
              plink --bfile merged_data_final_pruned_for_pca_v2 \
                    --pca 10 \
                    --out final_pca_results_v2 
                    # New name


        After finding out that the shape was not due to LD or inverted regions, the next steps to try are:

            Project onto 1000 Genomes Reference Panel (Recommended Next Step for Interpretation):
              This is the standard way to interpret what your PCs mean in terms of global ancestry.
              You take your ~1800 individuals and project them onto the PC space defined by the diverse 1000 Genomes Project samples (which have known population labels like CEU, YRI, CHB, etc.).
              How to do it:
              Download 1000 Genomes Phase 3 PLINK files (ensure good overlap of SNPs with your data).
              Merge your data with 1000 Genomes data, keeping only common SNPs.
              Run PCA on this combined dataset.
              Plot the PCA results, coloring 1000 Genomes samples by their known population labels, and overlay your study samples.
              This will show you where your samples fall relative to known reference populations, helping you label the arms/clusters in your own PCA (e.g., "This arm clusters with 1000G Europeans," "This group clusters with 1000G East Asians").
              PLINK can do this, or specialized tools like smartpca from EIGENSOFT, or custom R scripts.
            Consider ADMIXTURE (Alternative/Complementary):
              Run ADMIXTURE on your LD-pruned dataset (merged_data_final_pruned_for_pca_v2). Try different values of K (number of ancestral populations, e.g., K=2, 3, 4, 5).
              Visualize the ancestry proportions for each individual. This can also help delineate population groups.
              You can then color your original PCA plot by the dominant ADMIXTURE component for each individual.






# Sample Relatedness

# Final Sample Genotyping Call Rate (Optional Post-Relatedness Check)
# Minor Allele Frequency (MAF) Filter

# Hardy-Weinberg Equilibrium (HWE) Filter
