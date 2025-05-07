
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




# Population Stratification (PCA) & Batch Effect Visualization

# Sample Relatedness

# Final Sample Genotyping Call Rate (Optional Post-Relatedness Check)

# Minor Allele Frequency (MAF) Filter

# Hardy-Weinberg Equilibrium (HWE) Filter
