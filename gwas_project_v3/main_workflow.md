# General project description and requirements

The last part of this course focus on applying the knowledge acquired during these 10 weeks into real genetic datasets. We encourage you to work on the project in class Monday and Thursday where one (or several of us) will be there to guide you and answer questions.

The project is mandatory and needs to be handed in as a report.

The requirements of the report are:

It should be at most 10 pages (including references and figures)

It must be divided into sections: abstract, introduction, results/discussion, conclusion, and references.

It must be cohesive and coherent.

Source code must be provided (appended or linked/GitHub repository).

It must be in a PDF format.

Chosen project: GWAS of eye color or height

Deadline: likely 14th june

Submission: The project needs to be submitted through Brightspace. The name of your report must state your name and the chosen project:

MYNAME_archaic.pdf, MYNAME_xchromosome.pdf or MYNAME_GWAS.pdf

# Project choice description: GWAS of eye color or height

Listed relevant papers:
    Genome-wide association study in almost 195,000 individuals identifies 50 previously unidentified genetic loci for eye color: https://www.science.org/doi/10.1126/sciadv.abd1239

    A saturated map of common genetic variants associated with human height | Nature: https://www.nature.com/articles/s41586-022-05275-y

In this project, you will be looking at GWAS data from openSNP, which is a website where users of direct-to-customer genetic tests can share their personal data with other users. The phenotypes we will be looking at are self-reported eye color and height. When looking at the data, you should be aware that:

The data comes from different companies that use different chips so there are markers that are missing from some individuals because they were not present on the chip used by their company.
The gender information is missing from the file and by default plink will ignore the phenotype of individuals without gender information. So you have to use “--allow-no-sex” option in plink.

Investigate the following
A. Do QC. Are there any closely related individuals in the sample?

B. Do a PCA plot. What does it tell you about the samples?

C. The files eye_color.txt and height.txt contain the self-reported eye color and height of the individuals in the data. Do a GWAS on one of these traits. There are 12 eye color categories, and you can group some of them together to create a binary phenotype. How many significant loci do you find?

D. Do additional analyses using Plink, GCTA, R, or any other tool you might find relevant. The list below contains some suggestions for possible analyses, but you can also come up with your ideas

Suggestions for further analyses:

Use mixed model for GWAS.
Do imputation (either of the whole genome or the region around the most significant SNP) and see if you can then find variants with lower p-values.
If you use half of the data set to calculate a polygenic score, how well does that score predict height on the other half?
Find a trained height PRS on the internet. How well does it predict the height in this data set?
Test for epistasis.
What are the distribution of phenotypes for each of the genotypes at th most significant SNP? If you want to analyse it in R you can use the "--recode A" together with the "--snp" and "--window" option in plink to get the variants around a specific SNP written to a text file that it is easy to load in R.
How many of the significant variants found in the largest published GWAS study can you replicate those hits in this data set?
Make association tests where you condition on the most significant variant (you can use the --condition option in plink)


Data for the project can be found in this folder on the cluster:
~/populationgenomics/project_data/GWAS

Next will be a workflow / log.

# Create a new enviroment on the cluster containing the nessecary libraries, python and ect. 

conda create -n gwas \
  -c bioconda -c conda-forge \
  plink plink2 vcftools bcftools samtools gcta admixture \
  python=3.10 numpy pandas matplotlib seaborn scikit-learn jupyter \
  r-base r-ggplot2 r-data.table r-dplyr r-optparse r-ggpubr r-readr

# Make a project folder under your user. This is where we will keep the GWAS data. You can CD into the folder using the following command: 

cd /faststorage/project/populationgenomics/students/estherhelga/GWAS_project/

# Create a symlink to the data, and copy the data into your folder:

*Create symlink*
ln -s ~/populationgenomics/project_data/GWAS
  
*Extract from symlink:*
cp /home/estherhelga/populationgenomics/project_data/GWAS/eye_color.txt .
cp /home/estherhelga/populationgenomics/project_data/GWAS/metadata.txt .
cp /home/estherhelga/populationgenomics/project_data/GWAS/gwas_data.fam .
cp /home/estherhelga/populationgenomics/project_data/GWAS/gwas_data.bim .
cp /home/estherhelga/populationgenomics/project_data/GWAS/gwas_data.bed .

# Explore the eyecolor data file
  
*How many unique eye color categories are present? We can use the following command:*

awk '{print $2}' eye_color.txt | sort | uniq -c | sort -nr

*This gives us the following:*
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
    

# Data analysis and thoughts

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

*Creating the new phenotype file*
    The original eye_color.txt file had the self reported colors, so using R, I made a new file called eye_color_updated.txt that has the "IID", "FID", "phenotype". This R script is called '4_model_scale_creation'

# Metadata exploration

*After exploring the metadata file, we found that:*

The metadata file contains the following columns:
user, 	build, 	chip, 	chip_version, 	inferred_sex, 	source

*We have 4 different chips:*
HTS iSelect HD (587 entries),    Illumina GSAs (248 entires),      OmniExpress (291 entries),	 OmniExpress plus (484 entries).

We have a lot of different companies, some who use the same chips, some who dont, and some who never report their chip. In total, we have 1610 chip reports, and 399 that did not report a chip. We mark that as 'unknown'.

From this info, we can now group our data by chip, and do QC on each individual chip type. 

# Chip grouping

First, we split our main dataset into 4 chip-based groups:

    HTS iSelect HD

    Illumina GSAs

    OmniExpress

    OmniExpress Plus

Using the split_by_chip.R on the cluster, and the command Rscript to run it, we create 4 individual chip txt files, that we can then use alongside PLINK to --keep and --remove during the individual QC steps. The .txt files contain the Family ID (FID), and the Individual ID (IID), which in our case are the same. But PLINK needs both to operate the --keep and the --remove.

We might also consider doing one more group, which includs the 399 entries that did not report a chip, and do a QC on that group also. This might tell us something. MEaning we now have these 5 chip groups:

5 chip-based groups:

    HTS iSelect HD

    Illumina GSAs

    OmniExpress

    OmniExpress Plus

    Unknown

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

*Objective*
The initial chip-specific datasets (chip_<CHIPNAME>.bed/bim/fam) were created by subsetting individuals based on their reported genotyping chip. However, the associated .bim files at this stage still listed all SNPs present in the original merged gwas_data fileset. Many of these SNPs were not actually assayed on each specific chip, leading to artificial 100% missingness for those SNPs within a given chip-specific dataset. This step aims to:

    Remove SNPs not genuinely genotyped or poorly genotyped on each specific chip (SNP-level QC using --geno).
    Remove individuals with high missing genotype rates relative to the SNPs assayed on their specific chip (sample-level QC using --mind).

*Method*
This was achieved using a Bash script (clean_chip_bfiles.sh) that automates PLINK commands for each of the five chip-specific datasets (chip_HTS_iSelect_HD, chip_Illumina_GSAs, chip_OmniExpress, chip_OmniExpress_plus, chip_unknown).

To perform this step:

Prepare the script: *Ensure clean_chip_bfiles.sh is present in your project directory*. The script should contain PLINK commands to:

    First, filter SNPs using --geno (e.g., 0.05 for known chips, 0.2 for chip_unknown).
    Second, filter samples using --mind on the SNP-filtered data (e.g., 0.02 for known chips, 0.05 for chip_unknown).

The output for each chip group should be named chip_<CHIPNAME>_cleaned.

Make the script executable:
    *chmod +x clean_chip_bfiles.sh*

Run the script:
    *./clean_chip_bfiles.sh*

The output of this two-step process per chip will then be ready for merging.

*Outcome*

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


# Merge the Datasets

*Purpose*
To combine the five cleaned, chip-specific PLINK filesets (chip_<CHIPNAME>_cleaned) into a single, unified dataset for subsequent cohort-wide quality control and analysis.

*Method*
PLINK's --merge-list function was used.

To perform this step:
Create a merge list file (merge_list.txt):
This file specifies the basenames of the PLINK filesets to be merged (one per line, path relative to the working directory). From your main project directory, run:

echo "chip_HTS_iSelect_HD/chip_HTS_iSelect_HD_cleaned" > merge_list.txt
echo "chip_Illumina_GSAs/chip_Illumina_GSAs_cleaned" >> merge_list.txt
echo "chip_OmniExpress/chip_OmniExpress_cleaned" >> merge_list.txt
echo "chip_OmniExpress_plus/chip_OmniExpress_plus_cleaned" >> merge_list.txt
echo "chip_unknown/chip_unknown_cleaned" >> merge_list.txt


Verify the contents: cat merge_list.txt

Run PLINK to merge:
plink --merge-list merge_list.txt \
      --make-bed \
      --out allchips_qc1_merged

*Outcome*

A new PLINK binary fileset: allchips_qc1_merge.bed/bim/fam.

This dataset contains:
Individuals: 1880 (sum of individuals from all cleaned chip filesets).
Variants: 1,325,823 (union of all unique SNPs from the cleaned chip filesets).
The associated merged_data_step0.log file should be checked for any merge warnings (e.g., regarding SNP strand or position mismatches).

# Making sure every reported eye color has genotype data, and that every genotype data has a reported eye color; if only has one = removed

*Purpose*
To ensure that the individuals included in the downstream quality control and GWAS have both genetic data (from gwas_data.fam) and phenotype data (from the original eye_color.txt). This step identifies and retains only those individuals present in both datasets.

*Method*
Command-line tools (awk, sort, comm) were used to compare lists of Individual IDs (IIDs) from the genetic data and the phenotype data. A PLINK --keep command was then used to filter the merged genetic dataset.

To perform this step:

Extract Unique IIDs from Original Genetic Data (gwas_data.fam):

This file contains all individuals initially present in the genetic dataset.
    *awk '{print $2}' gwas_data.fam | sort -u > gwas_iids.txt*

Expected output: gwas_iids.txt containing a sorted list of unique IIDs (e.g., 2009 IIDs).

Extract Unique IIDs from Original Phenotype Data (eye_color.txt):

This file contains all individuals for whom original phenotype data was provided.

(Adjust the awk command based on your eye_color.txt structure: no header, IID is the first column, space-separated assumed here).
    *awk '{print $1}' eye_color.txt | sort -u > pheno_iids.txt*

Expected output: pheno_iids.txt containing a sorted list of unique IIDs (e.g., 1868 IIDs).

Identify IIDs Common to Both Datasets:

This identifies individuals for whom both types of data are available.
    *comm -12 gwas_iids.txt pheno_iids.txt > iids_common_to_both.txt*

Expected output: iids_common_to_both.txt containing only IIDs present in both lists (e.g., 1818 IIDs).

Create a PLINK "Keep" File (FID and IID format):

PLINK's --keep function requires both Family ID (FID) and IID. We use gwas_data.fam to retrieve the FIDs for the common IIDs.

    awk 'NR==FNR{keep[$1]; next} $2 in keep {print $1, $2}' iids_common_to_both.txt gwas_data.fam > individuals_to_keep_for_gwas.txt

Expected output: individuals_to_keep_for_gwas.txt with two columns (FID and IID) for each of the common individuals (e.g., 1818 lines).

Filter the Merged Genetic Dataset:

This step applies the filter to your main merged genetic data file created in the previous section.

Input file: allchips_qc1_merged.bed/bim/fam.

plink --bfile allchips_qc1_merged \
      --keep individuals_to_keep_for_gwas.txt \
      --make-bed \
      --out merged_data_core_cohort

*1696* people remaining.

*Outcome*

A new PLINK binary fileset: merged_data_core_cohort.bed/bim/fam.
This dataset now contains only individuals for whom both original genetic data and original phenotype data were available (1696 individuals).

The number of SNPs in merged_data_core_cohort.bim will be the same as in allchips_qc1_merged.bim (1,325,823 variants), as only individuals were filtered at this stage.

This merged_data_core_cohort fileset forms the basis for all subsequent QC and analysis steps.

# Heterozygous Haploid Genotype Correction

*Purpose*
To identify and correct erroneous "heterozygous haploid" genotype calls within the merged_data_core_cohort dataset. These occur when a genotype call shows two different alleles (e.g., A/G) on a chromosome where an individual is expected to have only one copy (e.g., non-pseudoautosomal regions of the X chromosome in males, or the Y chromosome). Such calls are biologically implausible and typically represent genotyping errors or issues with pseudo-autosomal region definitions. Setting these calls to 'missing' cleans the data and prevents these artifacts from skewing downstream QC and analysis.

*Method*
The PLINK command --set-hh-missing is used. This command does not remove individuals or SNPs but converts the specific problematic genotype calls to missing (0/0).

To perform this step:
Run PLINK to correct heterozygous haploid genotypes:
Input file: merged_data_core_cohort.bed/bim/fam
plink --bfile merged_data_core_cohort \
      --set-hh-missing \
      --make-bed \
      --out merged_data_core_cohort_hhcleaned

*Outcome*
A new PLINK binary fileset: merged_data_core_cohort_hhcleaned.bed/bim/fam.
This dataset will have the same number of individuals and SNPs as the input (merged_data_core_cohort).

The merged_data_core_cohort_hhcleaned.log file will indicate if any heterozygous haploid genotypes were found and set to missing. (You would have seen warnings about this in the log of the initial merge that created allchips_qc1_merged if they were present across the full set of 1880 individuals. This step now specifically cleans them within your defined core cohort of 1696 individuals).
This merged_data_core_cohort_hhcleaned fileset is now ready for sex inference.

# Sex Inference (Sex Check)
  
*Purpose*
To genetically infer the sex of individuals in the merged_data_core_cohort_hhcleaned dataset using X chromosome heterozygosity (F-statistic).
To identify individuals with ambiguous or failed genetic sex determination.
To make informed decisions about handling these individuals, particularly considering chip-specific data limitations (i.e., the OmniExpress chip).
To filter out individuals (excluding OmniExpress samples) whose genetic sex data suggests broader genotyping quality issues.

*Methodology Overview*
This involved running PLINK's --check-sex, followed by detailed analysis in R to interpret the results in the context of chip-specific data quality, and finally, filtering the dataset based on these findings.

To perform this step:
Run PLINK --check-sex:
Input file: merged_data_core_cohort_hhcleaned.bed/bim/fam

plink --bfile merged_data_core_cohort_hhcleaned \
      --check-sex \
      --out core_cohort_sex_check_results

Output files: core_cohort_sex_check_results.sexcheck (main results), .log, .nosex.
Note: The .fam file loaded at this stage (from merged_data_core_cohort_hhcleaned) contains the individuals common to both genetic and phenotype data. The number of males/females/ambiguous reported by PLINK here will be based on the sex codes present in this specific input .fam file.


Detailed Analysis in R (using sex_analysis.Rmd):

Load core_cohort_sex_check_results.sexcheck.
Load iid_to_chip.txt (mapping IIDs to their original chip platform).
Merge these two dataframes.

Key Findings from R Analysis (to be summarized from your actual R output):
OmniExpress Chip: Confirmed that all individuals genotyped on chip_OmniExpress (e.g., N=~280, use your actual count) had F=NaN and SNPSEX=0 (undetermined). This was previously traced back to the chip_OmniExpress_cleaned.bim file containing no X or Y chromosome SNPs.
Non-OmniExpress Samples:
Visualize the F-statistic distribution for these samples.
Define F-statistic thresholds to categorize individuals as genetically male (e.g., F > 0.8), female (e.g., F < 0.2), or ambiguous (e.g., 0.2 ≤ F ≤ 0.8 or F = NaN). Specify the exact thresholds you decided upon.
Create Removal List: Generate a file named core_cohort_individuals_to_remove_sex_issues.txt containing the FID and IID of only non-OmniExpress individuals who fall into the "ambiguous" or "NaN F-statistic" categories.

*Lists 16 individuals to be removed.*

Remove Identified Problematic Non-OmniExpress Individuals:
Input file: merged_data_core_cohort_hhcleaned.bed/bim/fam

plink --bfile merged_data_core_cohort_hhcleaned \
      --remove core_cohort_individuals_to_remove_sex_issues.txt \
      --make-bed \
      --out merged_data_core_cohort_sexcleaned


*Final Decisions and Rationale for Sex Handling:*

OmniExpress Samples: Due to the confirmed absence of X/Y chromosome SNPs on the chip_OmniExpress platform after per-chip cleaning, genetic sex could not be inferred for these individuals. These individuals (e.g., N=~280) will be retained in the dataset. Their sex will be treated as "unknown" (typically coded '0' in the .fam file) for the GWAS. The PLINK flag --allow-no-sex will be utilized during association testing to include them. This decision prioritizes sample size for the autosomal eye color analysis.

Non-OmniExpress Samples: Individuals from other chip platforms whose F-statistic was NaN or fell into a pre-defined ambiguous range (e.g., between 0.2 and 0.8, state your thresholds) were removed. This step acts as an additional quality control measure, as unresolved sex based on X-chromosome data (for chips that do have X-SNPs) can sometimes indicate broader sample quality issues.

Use of Sex as a Covariate: Given that eye color is not considered a strongly sex-dependent trait, and a significant portion of the cohort (OmniExpress samples) will have unknown genetic sex, sex will not be used as a covariate in the primary GWAS model. This simplifies the analysis and treats all individuals consistently in this regard.

Project Requirement for --allow-no-sex: This aligns with the project instructions, which noted that gender information might be missing and the --allow-no-sex option in PLINK would be necessary.

*Outcome*
A new PLINK binary fileset: merged_data_core_cohort_sexcleaned.bed/bim/fam.
This dataset contains individuals from the core cohort after removing a small number of non-OmniExpress samples with problematic genetic sex inference. All OmniExpress samples are retained. The sex codes in the .fam file reflect the best available information (genetically inferred for some, original PEDSEX, or '0' for unknown/OmniExpress).
This fileset is now ready for global SNP-level QC.

*1680 people remaining.*

# [OMITTED] Global Marker Genotyping Call Rate (SNP missingness)

For this, we used the following command, with the global missingness threshhold of 0.05. However, that only left us with aproxx 62 thousand variants out of the initial 1.3 million, which feels very strict. Thus, we want to try doing some R analysis to see if we can find a more appropreate threshold. 

    plink --bfile merged_data_after_sex_removal \
      --geno 0.05 \
      --make-bed \
      --out merged_data_step1_geno_final 

      Try geno (from R analysis below)
      0.3
      0.5
      0.6

R analysis and file generation
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

*Purpose*
To generate a subset of SNPs that are in approximate linkage equilibrium (LD). This is crucial because methods like Principal Component Analysis (PCA) and Identity-by-Descent (IBD) for relatedness estimation perform best with, or assume, SNPs that are not highly correlated due to physical proximity. This process involves two main parts: excluding known large regions of complex LD, and then applying a window-based pruning algorithm.

*Method*

*Part 1: Exclusion of Known Long-Range High-LD / Inversion Regions*

Rationale: Certain genomic regions, such as the Major Histocompatibility Complex (MHC) and common chromosomal inversions (e.g., on 8p23.1, 17q21.31), and regions under strong selection (e.g., LCT locus), exhibit complex and extensive LD patterns. These can disproportionately influence PCA results, potentially masking true genome-wide ancestry signals or creating artifactual clusters. Excluding SNPs within these pre-defined regions is a common step to mitigate their impact.

To perform this sub-step:
Identify Genome Build: Confirmed as GRCh37 based on rsID position checks.
Define Exclusion Regions: A file named high_ld_regions_exclude.txt was created. This file lists the chromosome, start, and end positions (GRCh37 coordinates) for regions to be excluded. It initially included the MHC and was subsequently updated to include other known regions like the 8p23.1 inversion, 17q21.31 inversion, LCT region, and potentially others identified from literature (e.g., UMich list).
Example content of high_ld_regions_exclude.txt (ensure this matches your final list):
Chr Start      End        Label
6   25000000   35000000   MHC
8   8000000    12000000   Inv8p23_1 
17  43000000   45000000   Inv17q21_31
2   135000000  137000000  LCT_region 
 ... other regions ...

Apply Exclusion with PLINK:
Input file: merged_data_core_cohort_sexcleaned.bed/bim/fam

plink --bfile merged_data_core_cohort_sexcleaned \
      --exclude range high_ld_regions_exclude.txt \
      --make-bed \
      --out merged_data_core_cohort_LDregions_excluded

*Outcome* of Part 1:
A new PLINK fileset: merged_data_core_cohort_LDregions_excluded.bed/bim/fam.
This dataset has SNPs from the specified complex LD regions removed. The merged_data_core_cohort_LDregions_excluded.log file details the number of SNPs excluded.
--exclude range: 43898 variants excluded.
--exclude range: 1281925 variants remaining.

*Part 2: Window-Based LD Pruning*

Rationale: After excluding large known complex regions, this step further refines the SNP set by removing SNPs that are still in LD with nearby markers based on a sliding window approach and an r² threshold.

To perform this sub-step:
Identify SNPs for Pruning (--indep-pairwise):
This command generates a list of SNPs to keep.
Input file: merged_data_core_cohort_LDregions_excluded.bed/bim/fam

plink --bfile merged_data_core_cohort_LDregions_excluded \
      --indep-pairwise 50 5 0.2 \
      --out core_cohort_LDregions_excluded_pruning_results

Parameters: 50 (window size in SNPs), 5 (step size for window), 0.2 (r² threshold for pruning).
Output files: core_cohort_LDregions_excluded_pruning_results.prune.in (list of SNPs to keep) and .prune.out (list of SNPs removed).
From the log (core_cohort_LDregions_excluded_pruning_results.log):
Pruned 75427 variants from chromosome 1, leaving 28547.
Pruned 73693 variants from chromosome 2, leaving 27034.
Pruned 63387 variants from chromosome 3, leaving 24060.
Pruned 56823 variants from chromosome 4, leaving 22950.
Pruned 54453 variants from chromosome 5, leaving 21449.
Pruned 53544 variants from chromosome 6, leaving 20218.
Pruned 52391 variants from chromosome 7, leaving 19551.
Pruned 46982 variants from chromosome 8, leaving 17381.
Pruned 40433 variants from chromosome 9, leaving 14370.
Pruned 49763 variants from chromosome 10, leaving 18959.
Pruned 46936 variants from chromosome 11, leaving 17392.
Pruned 44216 variants from chromosome 12, leaving 16605.
Pruned 34893 variants from chromosome 13, leaving 12870.
Pruned 29732 variants from chromosome 14, leaving 11293.
Pruned 28404 variants from chromosome 15, leaving 10706.
Pruned 29218 variants from chromosome 16, leaving 11818.
Pruned 26992 variants from chromosome 17, leaving 10781.
Pruned 26326 variants from chromosome 18, leaving 10582.
Pruned 18334 variants from chromosome 19, leaving 8274.
Pruned 21886 variants from chromosome 20, leaving 8993.
Pruned 12522 variants from chromosome 21, leaving 5136.
Pruned 13316 variants from chromosome 22, leaving 5562.
Pruned 25181 variants from chromosome 23, leaving 6937.
Pruned 1412 variants from chromosome 24, leaving 220.
Pruned 72 variants from chromosome 25, leaving 60.
Pruned 3218 variants from chromosome 26, leaving 623.
Pruning complete.  929554 of 1281925 variants removed.

Create the Final LD-Pruned Dataset:
This command uses the list of SNPs to keep to create the new, pruned fileset.
Input files: merged_data_core_cohort_LDregions_excluded.bed/bim/fam and core_cohort_LDregions_excluded_pruning_results.prune.in

plink --bfile merged_data_core_cohort_LDregions_excluded \
      --extract core_cohort_LDregions_excluded_pruning_results.prune.in \
      --make-bed \
      --out merged_data_final_pruned_for_pca

Overall Outcome of LD Pruning:
A new PLINK binary fileset: merged_data_final_pruned_for_pca.bed/bim/fam.
This dataset contains all individuals from merged_data_core_cohort_LDregions_excluded but only a subset of SNPs that are in approximate linkage equilibrium and outside the specified major complex LD regions.

--extract: 352371 variants remaining.

*Explanation for large pruning ratio (~72%)*

Does PLINK account for merging different chips with different SNP sets?

For LD Calculation: When PLINK calculates LD (r²) between a pair of SNPs for the --indep-pairwise command, it only uses individuals who are genotyped for both SNPs in that pair.
If many individuals are missing genotypes for one or both SNPs in a pair (which is common in merged multi-chip data where not all SNPs are on all chips), those individuals are excluded from that specific r² calculation.
If too few individuals have data for a pair, the LD calculation might be unreliable or skipped for that pair.

Impact: This means that for SNPs unique to smaller chips (and thus having high missingness across the full cohort), their LD relationships might be calculated based on a smaller effective sample size or might be harder to accurately assess. However, the pruning algorithm itself still functions: if a SNP is found to be in high LD (r² > 0.2) with another kept SNP based on the available data, it will be pruned.

PLINK doesn't inherently "know" which chip a SNP came from during this step. It just sees genotypes (0, 1, 2, or missing) and calculates LD.


Why such a large number of SNPs initially (1.28 million) and why so many pruned?

Nature of LD Pruning: LD pruning aims to find a set of "tagging" SNPs that are relatively independent. In dense regions, many SNPs can be highly correlated, so only a few are needed to represent the common variation. With 1.28 million SNPs, many of which will be in LD with each other (especially if coverage is dense from multiple chips in some regions), a reduction to 352k is not entirely unexpected for typical human genome LD patterns and r²=0.2.


This merged_data_final_pruned_for_pca dataset is the appropriate input for Principal Component Analysis (PCA) and subsequent relatedness checks.


# Population Stratification (PCA) & Batch Effect Visualization
  
*Purpose*
To identify major axes of genetic variation (population structure) within the core analysis cohort (merged_data_final_pruned_for_pca).
To visually assess potential batch effects by observing the distribution of samples from different chip platforms on the principal components (PCs).
To understand the relationship between population structure and the eye color phenotype.
To generate PCs for use as covariates in the downstream GWAS to control for population stratification.

*Methodology - PCA Calculation*
Input Dataset for PCA:
The LD-pruned dataset merged_data_final_pruned_for_pca.bed/bim/fam (containing 1680 individuals and ~352,371 SNPs after LD pruning and exclusion of complex LD regions) was used.

Run PCA using PLINK:

plink --bfile merged_data_final_pruned_for_pca \
      --pca 10 \
      --out final_pca_results

This command calculates the top 10 principal components.
Output files:
final_pca_results.eigenvec: Contains the eigenvectors (PC values) for each individual (FID, IID, PC1, PC2, ..., PC10).
final_pca_results.eigenval: Contains the eigenvalues for each calculated PC.
final_pca_results.log: PLINK log file.

*Methodology - Visualization and Initial Interpretation in R*
(R Markdown script: updated_PCA_analysis.Rmd)
Data Loading and Preparation:
Loaded final_pca_results.eigenvec and final_pca_results.eigenval.
Calculated and confirmed percentage of variance explained by each PC (PC1 ≈ 2.5%, PC2 ≈ 1% relative to N individuals).
Merged PCA data with chip information (iid_to_chip.txt), mapped eye color phenotypes (eye_color_updated.txt), original self-reported eye colors (eye_color.txt), and inferred sex information (core_cohort_sex_check_results.sexcheck) into a comprehensive dataframe (pca_with_metadata).

*Key Visualizations Generated*
Scree plot showing variance explained by each PC.
PC1 vs PC2 (and other combinations like PC1 vs PC3, PC3 vs PC4):
Uncolored (to show overall structure).
Colored by ChipType.
Colored by phenotype_labeled (mapped 0-3 eye color).
Colored by ReportedEyeColor_Factor (original self-reported).
Colored by Sex_Labeled.
Faceted plots (e.g., by ChipType while coloring by phenotype).

*Summary of Initial PCA Interpretation (incorporating your notes)*
Population Structure: The PCA plots (e.g., PC1 vs PC2) reveal substantial population stratification. A main dense cluster is observed, with distinct "arms" or gradients extending, particularly along PC1. This "V-shape" or "elbow" pattern, even after excluding several known long-range LD regions, suggests complex underlying ancestry within the openSNP cohort. The relatively low variance explained by PC1 and PC2 after correction is more typical for broad ancestry components than for LD artifacts. The observed structure appears similar to PCA plots of diverse human populations (e.g., reference to Nature Scientific Reports paper).
Batch Effects (Chip Type): Samples from different chip types appear largely intermingled on the plots of the top PCs (e.g., PC1 vs PC2, PC1 vs PC3). No single chip forms a completely isolated cluster distinct from the ancestry-driven patterns. This suggests that after the extensive QC, major batch effects related to chip platform are not the primary drivers of the observed structure in these top PCs for this SNP set.
Phenotype Correlation: A clear correlation exists between PC scores and eye color. For instance, individuals with brown eyes are more prevalent in the "arms" of the PCA plot that diverge from the main cluster, while blue eyes are concentrated in the main cluster. This highlights that eye color is strongly associated with the underlying genetic ancestry captured by the PCs.
Overall Conclusion from Visualization: The data is not genetically homogeneous and contains significant population structure related to ancestry. This structure correlates with the phenotype of interest (eye color). It is therefore crucial to use an adequate number of PCs as covariates in the GWAS to control for this confounding.
Addressing Your Current Question: Outlier Removal Based on PCA (e.g., +/- 6 SDs)
Yes, after visualizing the PCA, it's a common subsequent QC step to identify and potentially remove individuals who are extreme outliers on the principal components. The rationale is that these individuals might represent:
Mislabeled samples.
Individuals with very distinct ancestry not well-represented by the rest of the cohort (which could overly influence association tests or not be the target population).
Rare genetic variation leading to an extreme position.
Poor sample quality that wasn't caught by earlier filters but manifests as an odd genetic profile.
Is using +/- 6 Standard Deviations (SDs) from the mean of PC1 and PC2 correct?
Yes, this is a standard and well-accepted heuristic for identifying PCA outliers.
Typically, you would apply this to the first few informative PCs (e.g., PC1, PC2, and maybe up to PC4 or PC5, or however many PCs you decide to use as covariates based on the scree plot or further interpretation).
The idea is that individuals falling far outside the main distribution (6 SDs is quite extreme) are genetically very different from the bulk of your samples along that particular axis of variation.

*Next Steps (How to Implement Outlier Removal)*
Decide which PCs to use for outlier detection: Usually PC1 and PC2 are primary. You might extend this to PC3, PC4, etc., if they also show clear structure and you plan to use them as covariates.
Calculate Mean and SD for each selected PC in R:
Using your pca_with_metadata dataframe (or pca_eigenvectors if you prefer working with just the PCs initially):

PC1 and PC2
mean_pc1 <- mean(pca_with_metadata$PC1, na.rm = TRUE)
sd_pc1 <- sd(pca_with_metadata$PC1, na.rm = TRUE)
mean_pc2 <- mean(pca_with_metadata$PC2, na.rm = TRUE)
sd_pc2 <- sd(pca_with_metadata$PC2, na.rm = TRUE)

Define outlier thresholds
pc1_lower_bound <- mean_pc1 - (6 * sd_pc1)
pc1_upper_bound <- mean_pc1 + (6 * sd_pc1)
pc2_lower_bound <- mean_pc2 - (6 * sd_pc2)
pc2_upper_bound <- mean_pc2 + (6 * sd_pc2)

cat("PC1 Mean:", mean_pc1, "SD:", sd_pc1, "Bounds:", pc1_lower_bound, "to", pc1_upper_bound, "\n")
cat("PC2 Mean:", mean_pc2, "SD:", sd_pc2, "Bounds:", pc2_lower_bound, "to", pc2_upper_bound, "\n")

Identify Outliers in R:
outlier_pca_samples <- pca_with_metadata %>%
  filter(PC1 < pc1_lower_bound | PC1 > pc1_upper_bound |
         PC2 < pc2_lower_bound | PC2 > pc2_upper_bound)

cat("\nNumber of potential PCA outliers identified:", nrow(outlier_pca_samples), "\n")
if (nrow(outlier_pca_samples) > 0) {
  cat("Details of PCA outliers:\n")
  print(outlier_pca_samples %>% select(FID, IID, PC1, PC2, ChipType)) # Show some relevant info
}

Create a Removal List for PLINK:
If nrow(outlier_pca_samples) > 0, create a file (e.g., pca_outliers_to_remove.txt) with the FID and IID of these individuals.
if (nrow(outlier_pca_samples) > 0) {
  write.table(outlier_pca_samples[, c("FID", "IID")], "pca_outliers_to_remove.txt",
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  cat("Written PCA outliers to pca_outliers_to_remove.txt\n")
} else {
  cat("No PCA outliers identified with +/- 6SD criteria on PC1/PC2.\n")
}

Remove Outliers using PLINK:
You need to decide which dataset to remove them from. This should be a dataset before LD pruning for PCA, but after your core cohort definition and sex cleaning, because you want these outliers removed from the set that will eventually go into MAF/HWE filtering and GWAS.
A good candidate would be merged_data_core_cohort_geno_filtered (the one that has ~62k SNPs after global geno, or if you skipped global geno, then merged_data_core_cohort_sexcleaned which has ~1.28M SNPs).
Let's assume you apply it to merged_data_core_cohort_geno_filtered (if you kept that global --geno step) or merged_data_core_cohort_sexcleaned (if you omitted the global --geno). For consistency with your latest decision to omit global --geno before LD pruning, let's use merged_data_core_cohort_sexcleaned.

Previous output was merged_data_core_cohort_sexcleaned

Input for this PLINK command:

INPUT_BFILE_FOR_OUTLIER_REMOVAL="merged_data_core_cohort_sexcleaned"

plink --bfile ${INPUT_BFILE_FOR_OUTLIER_REMOVAL} \
      --remove pca_outliers_to_remove.txt \
      --make-bed \
      --out data_after_pca_outlier_removal

IMPORTANT: If you remove PCA outliers, you would then need to re-run the exclusion of high LD regions and LD pruning, and then re-run PCA on this newly filtered data_after_pca_outlier_removal dataset to get your final PCs for use as covariates. This is because removing individuals can change the LD structure and thus the PCs.
  



  
  
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

From updated:

    Population Stratification (PCA) & Batch Effect Visualization (Run PCA)

      plink --bfile merged_data_final_pruned_for_pca \
      --pca 10 \
      --out final_pca_results

      Based on the PCA plots (especially the one colored by chip type and the one just showing general structure):
      Population Structure: Do the samples form one homogenous cloud, or are there distinct clusters or gradients? If so, what might these correspond to (e.g., European ancestries, other ancestries if present in openSNP data)?
      Batch Effects: Do samples cluster primarily by chip type on any of the top PCs? If chip_OmniExpress (or any other chip) forms a separate cluster from the others, this is strong evidence of a batch effect. You'd mention that the PCA reveals chip-specific clustering.
      Outliers: Are there any individual samples that are extreme outliers on the PCA plot? These might warrant further investigation (though QC should have caught most problematic samples).

        PCA: PC1 vs PC2 (Population Structure - Uncolored Plot):
          Observation: You see a main dense "cloud" of points near the origin (0,0), but also distinct "arms" or gradients extending away, particularly towards the positive PC1 direction, and separating further along PC2 in the upper-right and lower-right quadrants.
          Interpretation (Population Structure): This is a classic picture of significant population substructure.
          PC1: This axis captures the largest source of genetic variation among your samples. Given that openSNP data often has a strong European component but also includes individuals of other ancestries, PC1 frequently separates individuals along a major axis of European variation (e.g., North-West Europe vs. South/East Europe) or sometimes separates Europeans from non-Europeans. The large spread along this axis means there are substantial genetic differences captured here.
          PC2: This captures the second largest amount of variation. It further differentiates groups. It might separate different European subgroups or distinguish individuals with non-European ancestry (e.g., East Asian, African, Ashkenazi Jewish – whose genetic profiles differ from the main European cluster) from the main European cluster(s).
          The Dense Cloud: Likely represents the most numerous ancestry group in your subset (often Northern/Western European descent in openSNP).
          The "Arms": Represent individuals who are genetically distinct from the main cloud along these primary axes. They could be from different European populations, individuals with non-European ancestry, or potentially admixed individuals.
          Conclusion: Your dataset is not genetically homogeneous. It contains substantial population structure, likely reflecting diverse ancestries present in the openSNP cohort. This structure must be controlled for in your GWAS to avoid spurious associations. Including PC1 and PC2 (and possibly more PCs) as covariates is essential.


        PCA: PC1 vs PC2, Colored by Chip Type:
          Observation: As you noted, the different colors (chip types) are largely intermingled throughout the structure revealed in the first plot. There isn't one color forming its own isolated island separate from the others along PC1 or PC2.
          Interpretation (Batch Effects): This is generally very good news. It indicates that the major population structure patterns (captured by PC1 and PC2) are not primarily driven by which chip the sample was run on. If there were strong, systematic technical differences between chips affecting many of the common SNPs used for PCA, you would expect to see distinct clusters based purely on color. The mixing suggests that genetic variation (ancestry) is a much stronger signal than chip type for these top PCs after your QC steps.
          Conclusion: While originating from different genotyping platforms, the samples do not show strong clustering by chip type along the primary axes of genetic variation (PC1/PC2) after quality control, suggesting major batch effects are not confounding the observed population structure in this view.

        PCA: PC1 vs PC2, Colored by Mapped Eye Color:
          Observation: You see clear non-random distribution. Blue (0) is concentrated in the main dense cloud. Green/Hazel (1, 2) are also mostly in/near that cloud. Brown (3) is present in the cloud but is much more prevalent in the "arms" extending towards positive PC1. Your observation about darker eyes "stretching out" is accurate. The two NAs are also plotted somewhere based on their genetics.
          Interpretation (Phenotype-Structure Correlation): This strongly indicates that eye color is correlated with the genetic ancestry captured by PCA.
          The main cluster (low PC1/PC2), rich in blue eyes, likely represents populations with high frequencies of alleles associated with blue eyes (e.g., common variants near HERC2/OCA2 prevalent in Europeans, especially Northern Europeans).
          The "arms" (higher PC1), rich in brown eyes, likely represent populations where alleles for brown eyes are much more common (e.g., Southern Europeans, non-Europeans).
          This is biologically expected. Eye color frequencies vary significantly across global populations.
          Conclusion: Eye color phenotype is significantly correlated with the population structure identified by PC1 and PC2. This underscores the critical importance of using PCs as covariates in the GWAS. Failing to do so would lead to highly inflated results and false positives, simply rediscovering that different ancestry groups have different eye colors, rather than finding specific causal variants within those populations.
              
        The PCA clearly reveals substantial population stratification within your openSNP cohort, likely reflecting diverse European and potentially non-European ancestries.
     Encouragingly, samples do not cluster strongly by genotyping chip along these main PCs, suggesting major batch effects are not driving the primary structure after QC. However, there is a clear correlation between the observed population structure (PCs) and the eye color phenotype, with most colors mixing in the main genetic cluster and brown eyes more prevalent in genetically distinct groups. This highlights the necessity of including principal components (e.g., PC1, PC2, and potentially more) as covariates in the subsequent GWAS to control for confounding due to ancestry.


      If it is ancestry, it does look pretty similar to these: https://www.nature.com/articles/s41598-021-97129-2/figures/1 

        
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
            


# Sample Relatedness

# Final Sample Genotyping Call Rate (Optional Post-Relatedness Check)

# Minor Allele Frequency (MAF) Filter

# Hardy-Weinberg Equilibrium (HWE) Filter
