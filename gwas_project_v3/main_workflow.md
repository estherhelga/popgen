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
echo "chip_OmniExpress_plus/chip_OmniExpress_plus_cleaned" >> merge_list.txt
echo "chip_unknown/chip_unknown_cleaned" >> merge_list.txt

[REMOVED] echo "chip_OmniExpress/chip_OmniExpress_cleaned" >> merge_list.txt

Verify the contents: cat merge_list.txt

Run PLINK to merge:
plink --merge-list merge_list.txt \
      --make-bed \
      --out allchips_qc1_merged

*Outcome*

A new PLINK binary fileset: allchips_qc1_merge.bed/bim/fam.

This dataset contains:
1306452 variants and 1600 people pass filters and QC.
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

*Outcome*

*1306452 variants and 1449 people pass filters and QC.*

A new PLINK binary fileset: merged_data_core_cohort.bed/bim/fam.
This dataset now contains only individuals for whom both original genetic data and original phenotype data were available.

The number of SNPs in merged_data_core_cohort.bim will be the same as in allchips_qc1_merged.bim, as only individuals were filtered at this stage.

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

*1306452 variants and 1449 people pass filters and QC.*

A new PLINK binary fileset: merged_data_core_cohort_hhcleaned.bed/bim/fam.
This dataset will have the same number of individuals and SNPs as the input (merged_data_core_cohort).

The merged_data_core_cohort_hhcleaned.log file will indicate if any heterozygous haploid genotypes were found and set to missing. (You would have seen warnings about this in the log of the initial merge that created allchips_qc1_merged if they were present across the full set of 1880 individuals. This step now specifically cleans them within your defined core cohort of 1449 individuals).
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

[OMNIEXPRESS WAS REMOVED] OmniExpress Chip: Confirmed that all individuals genotyped on chip_OmniExpress (e.g., N=~280, use your actual count) had F=NaN and SNPSEX=0 (undetermined). This was previously traced back to the chip_OmniExpress_cleaned.bim file containing no X or Y chromosome SNPs.

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

--remove: 1433 people remaining.

*Final Decisions and Rationale for Sex Handling:*

[OMNIEXPRESS WAS REMOVED] OmniExpress Samples: Due to the confirmed absence of X/Y chromosome SNPs on the chip_OmniExpress platform after per-chip cleaning, genetic sex could not be inferred for these individuals. These individuals (e.g., N=~280) will be retained in the dataset. Their sex will be treated as "unknown" (typically coded '0' in the .fam file) for the GWAS. The PLINK flag --allow-no-sex will be utilized during association testing to include them. This decision prioritizes sample size for the autosomal eye color analysis.

Non-OmniExpress Samples: Individuals from other chip platforms whose F-statistic was NaN or fell into a pre-defined ambiguous range (e.g., between 0.2 and 0.8, state your thresholds) were removed. This step acts as an additional quality control measure, as unresolved sex based on X-chromosome data (for chips that do have X-SNPs) can sometimes indicate broader sample quality issues.

Use of Sex as a Covariate: Given that eye color is not considered a strongly sex-dependent trait, and a significant portion of the cohort (OmniExpress samples) will have unknown genetic sex, sex will not be used as a covariate in the primary GWAS model. This simplifies the analysis and treats all individuals consistently in this regard.

Project Requirement for --allow-no-sex: This aligns with the project instructions, which noted that gender information might be missing and the --allow-no-sex option in PLINK would be necessary.

*Outcome*
A new PLINK binary fileset: merged_data_core_cohort_sexcleaned.bed/bim/fam.
This dataset contains individuals from the core cohort after removing a small number of non-OmniExpress samples with problematic genetic sex inference. All OmniExpress samples are retained. The sex codes in the .fam file reflect the best available information (genetically inferred for some, original PEDSEX, or '0' for unknown/OmniExpress).
This fileset is now ready for global SNP-level QC.

*1433 people remaining.*

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
--exclude range: 43199 variants excluded.
--exclude range: 1263253 variants remaining.

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
Pruned 74349 variants from chromosome 1, leaving 28058.
Pruned 72801 variants from chromosome 2, leaving 26547.
Pruned 62596 variants from chromosome 3, leaving 23694.
Pruned 56187 variants from chromosome 4, leaving 22600.
Pruned 53705 variants from chromosome 5, leaving 21026.
Pruned 52901 variants from chromosome 6, leaving 19883.
Pruned 51643 variants from chromosome 7, leaving 19257.
Pruned 46216 variants from chromosome 8, leaving 17136.
Pruned 39725 variants from chromosome 9, leaving 14112.
Pruned 48994 variants from chromosome 10, leaving 18597.
Pruned 46165 variants from chromosome 11, leaving 17142.
Pruned 43522 variants from chromosome 12, leaving 16309.
Pruned 34373 variants from chromosome 13, leaving 12593.
Pruned 29329 variants from chromosome 14, leaving 11123.
Pruned 28026 variants from chromosome 15, leaving 10538.
Pruned 28816 variants from chromosome 16, leaving 11589.
Pruned 26642 variants from chromosome 17, leaving 10521.
Pruned 25854 variants from chromosome 18, leaving 10424.
Pruned 18030 variants from chromosome 19, leaving 8071.
Pruned 21635 variants from chromosome 20, leaving 8820.
Pruned 12324 variants from chromosome 21, leaving 5081.
Pruned 13092 variants from chromosome 22, leaving 5484.
Pruned 25181 variants from chromosome 23, leaving 6937.
Pruned 1412 variants from chromosome 24, leaving 220.
Pruned 72 variants from chromosome 25, leaving 60.
Pruned 3218 variants from chromosome 26, leaving 623.
Pruning complete.  916808 of 1263253 variants removed.

Create the Final LD-Pruned Dataset:
This command uses the list of SNPs to keep to create the new, pruned fileset.
Input files: merged_data_core_cohort_LDregions_excluded.bed/bim/fam and core_cohort_LDregions_excluded_pruning_results.prune.in

plink --bfile merged_data_core_cohort_LDregions_excluded \
      --extract core_cohort_LDregions_excluded_pruning_results.prune.in \
      --make-bed \
      --out merged_data_final_pruned_for_pca

*--extract: 346445 variants remaining.*

Overall Outcome of LD Pruning:
A new PLINK binary fileset: merged_data_final_pruned_for_pca.bed/bim/fam.
This dataset contains all individuals from merged_data_core_cohort_LDregions_excluded but only a subset of SNPs that are in approximate linkage equilibrium and outside the specified major complex LD regions.

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
The LD-pruned dataset merged_data_final_pruned_for_pca.bed/bim/fam (containing xxxx individuals and xxx.xxx SNPs after LD pruning and exclusion of complex LD regions) was used.

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
A good candidate would be merged_data_core_cohort_sexcleaned (which has ~1.28M SNPs).
Let's assume you apply it to merged_data_core_cohort_geno_filtered (if you kept that global --geno step) or merged_data_core_cohort_sexcleaned (if you omitted the global --geno). For consistency with your latest decision to omit global --geno before LD pruning, let's use merged_data_core_cohort_sexcleaned.

Previous output was merged_data_core_cohort_sexcleaned

Input for this PLINK command:

INPUT_BFILE_FOR_OUTLIER_REMOVAL="merged_data_core_cohort_sexcleaned"

plink --bfile merged_data_core_cohort_sexcleaned \
      --remove pca_outliers_to_remove.txt \
      --make-bed \
      --out merged_data_core_cohort_no_pca_outliers

*1306452 variants and 1412 people pass filters and QC.*

IMPORTANT: If you remove PCA outliers, you would then need to re-run the exclusion of high LD regions and LD pruning, and then re-run PCA on this newly filtered data_after_pca_outlier_removal dataset to get your final PCs for use as covariates. This is because removing individuals can change the LD structure and thus the PCs.

*merged_data_core_cohort_no_pca_outliers*


# Second LD pruning post PCA outlier Removal

*Purpose*
After removing PCA outliers from the `merged_data_core_cohort_sexcleaned` dataset (resulting in `merged_data_core_cohort_no_pca_outliers`), it is necessary to re-derive Principal Components (PCs) that accurately reflect the population structure of this refined cohort. To do this, we need to again prepare an appropriate subset of SNPs that are in approximate linkage equilibrium. This involves re-applying the exclusion of known complex LD regions and then performing window-based LD pruning on the outlier-removed dataset.


*Method*

This process mirrors the initial LD pruning but is applied to the dataset that has had PCA outliers removed.

*Part 1: Exclusion of Known Long-Range High-LD / Inversion Regions (on Outlier-Removed Data)*

Rationale: To mitigate the disproportionate influence of known complex genomic regions (MHC, inversions, LCT, etc.) on the subsequent PCA calculation for the refined cohort.

To perform this sub-step:
Ensure the `high_ld_regions_exclude.txt` file (listing chromosome, start, end positions for regions to exclude) is available.

Apply Exclusion with PLINK:
Input file: `merged_data_core_cohort_no_pca_outliers.bed/bim/fam` (This is your dataset after removing individuals identified as PCA outliers from `merged_data_core_cohort_sexcleaned`).

plink --bfile merged_data_core_cohort_no_pca_outliers \
      --exclude range high_ld_regions_exclude.txt \
      --make-bed \
      --out merged_data_core_cohort_no_pca_outliers_LDregions_excluded

*1263253 variants and 1412 people pass filters and QC.*

*Outcome of Part 1:*
A new PLINK fileset: `merged_data_core_cohort_no_pca_outliers_LDregions_excluded.bed/bim/fam`.
This dataset has SNPs from the specified complex LD regions removed from the outlier-filtered cohort. The `merged_data_core_cohort_no_pca_outliers_LDregions_excluded.log` file will detail the number of SNPs excluded. Check this log to confirm (it should be a similar number of SNPs excluded as in the first LD region exclusion, but applied to the slightly smaller set of individuals).

*Part 2: Window-Based LD Pruning (on Outlier-Removed, Region-Excluded Data)*

Rationale: To further refine the SNP set by removing SNPs that are in LD with nearby markers, ensuring the SNPs used for the new PCA are relatively independent.

To perform this sub-step:
Identify SNPs for Pruning (`--indep-pairwise`):
This command generates a list of SNPs to keep based on the specified LD parameters.
Input file: `merged_data_core_cohort_no_pca_outliers_LDregions_excluded.bed/bim/fam`

plink --bfile merged_data_core_cohort_no_pca_outliers_LDregions_excluded \
      --indep-pairwise 50 5 0.2 \
      --out core_cohort_no_pca_outliers_LDregions_excluded_pruning_results

Pruned 74901 variants from chromosome 1, leaving 27506.
Pruned 73314 variants from chromosome 2, leaving 26034.
Pruned 63077 variants from chromosome 3, leaving 23213.
Pruned 56567 variants from chromosome 4, leaving 22220.
Pruned 54033 variants from chromosome 5, leaving 20698.
Pruned 53236 variants from chromosome 6, leaving 19548.
Pruned 52049 variants from chromosome 7, leaving 18851.
Pruned 46526 variants from chromosome 8, leaving 16826.
Pruned 40010 variants from chromosome 9, leaving 13827.
Pruned 49344 variants from chromosome 10, leaving 18247.
Pruned 46500 variants from chromosome 11, leaving 16807.
Pruned 43865 variants from chromosome 12, leaving 15966.
Pruned 34651 variants from chromosome 13, leaving 12315.
Pruned 29543 variants from chromosome 14, leaving 10909.
Pruned 28229 variants from chromosome 15, leaving 10335.
Pruned 29087 variants from chromosome 16, leaving 11318.
Pruned 26836 variants from chromosome 17, leaving 10327.
Pruned 26076 variants from chromosome 18, leaving 10202.
Pruned 18195 variants from chromosome 19, leaving 7906.
Pruned 21810 variants from chromosome 20, leaving 8645.
Pruned 12454 variants from chromosome 21, leaving 4951.
Pruned 13178 variants from chromosome 22, leaving 5398.
Pruned 25300 variants from chromosome 23, leaving 6818.
Pruned 1420 variants from chromosome 24, leaving 212.
Pruned 74 variants from chromosome 25, leaving 58.
Pruned 3227 variants from chromosome 26, leaving 614.
Pruning complete.  923502 of 1263253 variants removed.

Parameters: `50` (window size in SNPs), `5` (step size for window), `0.2` (r² threshold for pruning).
Output files:
`core_cohort_no_pca_outliers_LDregions_excluded_pruning_results.prune.in` (list of SNPs to keep)
`core_cohort_no_pca_outliers_LDregions_excluded_pruning_results.prune.out` (list of SNPs removed)
`core_cohort_no_pca_outliers_LDregions_excluded_pruning_results.log` (PLINK log file).
*Action: Review this log file to see the number of variants pruned per chromosome and the total number remaining. Compare these numbers to your first LD pruning run; they should be similar, but potentially slightly different due to the removal of outlier individuals.*

Create the Final LD-Pruned Dataset for Re-PCA:
This command uses the list of SNPs to keep (`.prune.in` file) to create the new, pruned fileset.
Input files: `merged_data_core_cohort_no_pca_outliers_LDregions_excluded.bed/bim/fam` and `core_cohort_no_pca_outliers_LDregions_excluded_pruning_results.prune.in`.

plink --bfile merged_data_core_cohort_no_pca_outliers_LDregions_excluded \
      --extract core_cohort_no_pca_outliers_LDregions_excluded_pruning_results.prune.in \
      --make-bed \
      --out merged_data_final_pruned_no_outliers_for_pca

*339751 variants and 1412 people pass filters and QC.*


*Note: Ensure the input filename in `--extract` matches the output from the `--indep-pairwise` step exactly.*

*Overall Outcome of Second LD Pruning:*
A new PLINK binary fileset: `merged_data_final_pruned_no_outliers_for_pca.bed/bim/fam`.
This dataset contains individuals from the outlier-removed cohort (`merged_data_core_cohort_no_pca_outliers`) but only a subset of SNPs (from `merged_data_core_cohort_no_pca_outliers_LDregions_excluded`) that are in approximate linkage equilibrium and outside the specified major complex LD regions.
The `merged_data_final_pruned_no_outliers_for_pca.log` file will state the final number of variants remaining.
This `merged_data_final_pruned_no_outliers_for_pca` dataset is the appropriate input for re-running the Principal Component Analysis to generate PCs for the outlier-removed cohort.
            

# Second PCA (Post Outlier Removal) - Generating Final Covariates

*Purpose*
Having removed PCA outliers and subsequently re-pruned the dataset for linkage disequilibrium (`merged_data_final_pruned_no_outliers_for_pca`), the next step is to perform Principal Component Analysis (PCA) again. This new PCA will generate eigenvectors (principal components, PCs) that accurately reflect the population structure of the refined cohort (i.e., the individuals remaining after outlier removal). These PCs will be used as covariates in the downstream GWAS to control for population stratification.

*Methodology - PCA Calculation*

Input Dataset for PCA:
The LD-pruned dataset created after outlier removal and re-pruning:
`merged_data_final_pruned_no_outliers_for_pca.bed/bim/fam`

This dataset contains:
*   Individuals: The original N=XXXX minus the number of PCA outliers removed.
*   SNPs: The set of SNPs remaining after excluding high-LD regions and performing window-based LD pruning on the outlier-removed cohort = XXXXXX

Run PCA using PLINK:

plink --bfile merged_data_final_pruned_no_outliers_for_pca \
      --pca 10 \
      --out final_pca_results_no_outliers

This command calculates the top 10 principal components for the refined cohort.
*   `--bfile merged_data_final_pruned_no_outliers_for_pca`: Specifies the input dataset.
*   `--pca 10`: Instructs PLINK to calculate the top 10 PCs. You can adjust this number if desired (e.g., 20), but 10 is a common starting point.
*   `--out final_pca_results_no_outliers`: Specifies the prefix for the output files.

Output files:
*   `final_pca_results_no_outliers.eigenvec`: Contains the eigenvectors (PC values) for each individual in the refined cohort (FID, IID, PC1, PC2, ..., PC10). **This is the key file you will use for GWAS covariates.**
*   `final_pca_results_no_outliers.eigenval`: Contains the eigenvalues for each calculated PC. This tells you how much variance each PC explains.
*   `final_pca_results_no_outliers.log`: PLINK log file.

*Methodology - Visualization and Review in R (Recommended)*
It is highly recommended to briefly visualize these new PCs, similar to how you analyzed the first PCA.

(Open R with `PCA_analysis_post_outlier_removal.Rmd`)

Data Loading and Preparation:
*   Load `final_pca_results_no_outliers.eigenvec` and `final_pca_results_no_outliers.eigenval`.
*   Calculate and review the percentage of variance explained by each new PC.
*   Merge this new PCA data with your relevant metadata (chip information, mapped eye color phenotypes, original self-reported eye colors, sex information) for the *remaining individuals*.

Key Visualizations to Generate (for review):
*   Scree plot showing variance explained by each PC (from `.eigenval`).
*   PC1 vs PC2 (and other combinations like PC1 vs PC3):
    *   Uncolored (to show overall structure).
    *   Colored by ChipType.
    *   Colored by `phenotype_labeled` (your 4-category eye color).
    *   Colored by Sex_Labeled.

*Expected Outcome and Interpretation of Second PCA:*
*   The general population structure observed should be similar to your first PCA, but potentially "tighter" or more refined, as the extreme outliers have been removed.
*   The variance explained by the top PCs might shift slightly.
*   Confirm that no new, unexpected gross outliers have appeared. The goal here is not to do another round of outlier removal unless something is clearly wrong, but to ensure the PCs look reasonable for the refined cohort.
*   These `final_pca_results_no_outliers.eigenvec` are now ready to be used as covariates in your association testing to account for the observed population structure in your final analysis cohort.


TIME TO GO BACK AND REMOVE OMINEXPRESS IN THE VERY BEGINNING! :) DONE


# Sample Relatedness Check (IBD Estimation)

*Purpose*
To identify pairs of individuals in the analysis cohort who are closely related (e.g., duplicates, parent-offspring, full siblings). Standard GWAS association tests assume samples are unrelated. Including closely related individuals can inflate test statistics and lead to false positives due to shared genetic material beyond what's expected by chance in a random population sample. This step aims to estimate identity-by-descent (IBD) to quantify relatedness and remove one individual from each related pair to ensure sample independence.

*Methodology*

PLINK's IBD estimation capabilities will be used. This is typically performed on an LD-pruned set of SNPs, as high LD can inflate IBD estimates. The dataset `merged_data_final_pruned_no_outliers_for_pca` (which was used for your final PCA) is suitable for this.

Input Dataset for IBD Estimation:
*   `merged_data_final_pruned_no_outliers_for_pca.bed/bim/fam`
    *   This dataset contains: 1412 individuals (after OmniExpress removal and one round of PCA outlier removal) and 339,751 LD-pruned SNPs.
    *   It's crucial to use an LD-pruned set to avoid LD confounding IBD estimates.

*Step 1: Initial IBD Calculation with PLINK*

Run PLINK to calculate IBD sharing:
This command calculates pairwise IBD estimates (PI_HAT, proportion IBD) for all pairs of individuals.

plink --bfile merged_data_final_pruned_no_outliers_for_pca \
      --genome \
      --out relatedness_check_results_no_outliers

Parameters:
*   `--bfile merged_data_final_pruned_no_outliers_for_pca`: Specifies the input LD-pruned dataset.
*   `--genome`: Instructs PLINK to perform IBD estimation. By default, it uses a method of moments approach. It considers parameters like minimum MAF for SNPs used in IBD estimation (default usually 0.01 or 0.05).
*   `--out relatedness_check_results_no_outliers`: Specifies the prefix for the output files.

Output files:
*   `relatedness_check_results_no_outliers.genome`: This is the primary output file. It lists pairs of individuals and their estimated IBD sharing statistics, including Z0, Z1, Z2 (probabilities of sharing 0, 1, or 2 alleles IBD), and PI_HAT (overall proportion of genome IBD).
*   `relatedness_check_results_no_outliers.log`: PLINK log file.

*Step 2: Identify and Handle Related Individuals in R*

(Run in R Studio `relatedness_analysis.Rmd`)

Define Relatedness Thresholds:
Commonly used PI_HAT thresholds for identifying related pairs:
Duplicates / Monozygotic Twins: PI_HAT > 0.9 (very close to 1)
Parent-Offspring or Full Siblings: PI_HAT > 0.4 (typically around 0.5)
Second-degree relatives (e.g., half-sibs, grandparent-grandchild, avuncular): PI_HAT > 0.1875 (typically around 0.25)
Third-degree relatives (e.g., first cousins): PI_HAT > 0.09375 (typically around 0.125)
For GWAS, it's common to remove one individual from pairs that are second-degree relatives or closer. A PI_HAT threshold of 0.1875 or 0.2 is often used as a cutoff to identify pairs needing resolution. Let's use PI_HAT > 0.1875.

Based on the number of relations per individual plot, the following algorithm runs. There were 20000+ related pairs, that generated 229 unique related individuals across all pairs.

This R script that creates 'related_individuals_to_remove.txt' implements a greedy algorithm: it finds the person involved in the most high-PI_HAT relationships, adds them to a removal list, removes all their relationships, and repeats until no high-PI_HAT relationships remain. This is a common approach.

"Number of pairs with PI_HAT > 0.1875: 21382"
"Written 216 individuals to related_individuals_to_remove.txt"

plink --bfile merged_data_core_cohort_no_pca_outliers \
      --remove related_individuals_to_remove.txt \
      --make-bed \
      --out final_analysis_cohort_unrelated

*1306452 variants and 1196 people pass filters and QC.*

# Final Check of Individual Genotyping Call Rates (F_MISS)

*Purpose*
To assess the distribution of global individual genotyping call rates in the `final_analysis_cohort_unrelated` dataset and to make a final decision on sample inclusion based on this metric.

*Methodology & Findings*
Individual missingness statistics were calculated using PLINK:

plink --bfile final_analysis_cohort_unrelated \
      --missing \
      --out final_cohort_imissing_stats

The resulting *final_cohort_imissing_stats.imiss* file was loaded into R (final_check_individual_genotyping_call_rates.rmd)for examination.

The distribution of F_MISS (proportion of missing SNPs per individual across the ~1.3 million SNPs in the merged dataset) was plotted. The histogram revealed a multi-modal distribution, with distinct clusters of individuals showing different overall missingness rates (e.g., peaks around 30-35% missing, 55-63% missing, and 85-90% missing). Approximately 263 individuals exhibited global missingness rates exceeding 85%.

*Rationale for Retaining All Individuals at This Stage*
Despite some individuals showing high global missingness rates, a decision was made to retain all 1196 individuals from the final_analysis_cohort_unrelated dataset for the following reasons:
- Effective Per-Chip QC: Robust individual call rate filtering (--mind) was performed on a per-chip basis before the datasets were merged. This ensured that individuals had high-quality genotypes for the SNPs actually present on their respective arrays. The high global F_MISS values observed post-merging are therefore primarily interpreted as a consequence of differing SNP content across the various genotyping platforms, rather than an indication of poor individual sample quality for assayed SNPs. Individuals from chips with less SNP overlap with the total merged SNP pool will naturally show higher global F_MISS.
- Preservation of Sample Size and Diversity: Removing a large number of individuals (such as the 263 with >85% F_MISS) solely based on global missingness would significantly reduce sample size and potentially remove individuals from less common chip platforms, thereby losing genetic and phenotypic diversity that could be valuable.
- Information Content for Overlapping SNPs: Even individuals with high global F_MISS contribute valid genotype information for the subset of SNPs (e.g., 10-15% for the highest missingness group) for which they do have data. These overlapping SNPs are often the more common ones.

*Downstream Analysis Capabilities:*
- Standard GWAS association models in PLINK (e.g., logistic/linear regression) perform analysis on a per-SNP basis and can typically handle individuals with missing genotypes for the SNP being tested (those individuals are simply excluded from the test for that specific SNP).
- Should more advanced analyses be pursued later, methods like mixed models (e.g., in GCTA or SAIGE) are robust to some levels of missing data.
- Furthermore, genotype imputation could be considered in future analyses. While individuals with very sparse genome-wide data can be challenging to impute accurately, imputation could potentially recover information for these samples if a suitable reference panel is used. Retaining them keeps this option open.

*Conclusion for this Step*
No individuals were removed based on global F_MISS thresholds at this stage. The final_analysis_cohort_unrelated dataset, containing 1196 individuals, will proceed to SNP-level quality control. The observed F_MISS distribution underscores the heterogeneous nature of the merged multi-chip dataset.



# Minor Allele Frequency (MAF) Filter

*Purpose*
To remove SNPs with very low minor allele frequencies (MAF) from the `final_analysis_cohort_unrelated` dataset. SNPs with extremely low MAF typically have very little statistical power to detect associations in a GWAS of the current sample size. They can also sometimes be more prone to genotyping errors, and tests involving them can be less stable. Filtering by MAF reduces the multiple testing burden by removing these less informative SNPs and focuses the analysis on variants common enough to yield potentially meaningful results with the available sample size.

*Methodology*

The PLINK command `--maf <threshold>` is used to filter out SNPs where the frequency of the minor allele is below the specified threshold. A common threshold is 0.01 (1%), meaning SNPs where the rarer allele is present in less than 1% of the chromosomes in the sample will be removed. Some studies use 0.05 (5%), especially with smaller sample sizes or for initial exploratory analyses. Given your sample size of 1196, a MAF of 0.01 is a reasonable starting point.

We ran the following PLINK command, which calculated allele frequencies without filtering any out:

plink --bfile final_analysis_cohort_unrelated \
      --freq \
      --out current_cohort_allele_freqs

We then used the output in R to visualize the distribution of MAF, so we would know which threshold would serve us best. We see that by far the biggest bin in the histogram is below 0.01, so thats the threshold we went with. You can see the R script under maf_threshold_calcs.rmd, and the histrogram under MAF_exclusion_threshold_histogram.png

Input Dataset for MAF Filtering:
*   `final_analysis_cohort_unrelated.bed/bim/fam`
    *   This dataset contains: 1196 individuals and 1,306,452 SNPs (after OmniExpress removal, PCA outlier removal, and relatedness filtering).

*Step 1: Apply MAF Filter using PLINK*

Run PLINK to filter SNPs based on MAF:

plink --bfile final_analysis_cohort_unrelated \
      --maf 0.01 \
      --make-bed \
      --out final_analysis_cohort_maf_filtered

Parameters:
--bfile final_analysis_cohort_unrelated: Specifies the input dataset.
--maf 0.01: Sets the minimum MAF threshold. SNPs with MAF below 0.01 will be removed. You can adjust this value (e.g., to 0.05) if desired, but 0.01 is a common choice.
--make-bed: Instructs PLINK to output the filtered dataset in binary PED format.
--out final_analysis_cohort_maf_filtered: Specifies the prefix for the output files.

Outcome
A new PLINK binary fileset: final_analysis_cohort_maf_filtered.bed/bim/fam.
This dataset will contain the same 1196 individuals but will have fewer SNPs than the input dataset, as rare variants have been removed.
The final_analysis_cohort_maf_filtered.log file will detail:
The number of SNPs initially present.
The number of SNPs removed due to the MAF filter.
The number of SNPs remaining.

Action: Review the PLINK log file (final_analysis_cohort_maf_filtered.log) to note how many SNPs were removed and how many remain. This gives an indication of the proportion of rare variants in your dataset.

Considerations and Rationale for MAF Threshold Choice (e.g., 0.01):
Statistical Power: Detecting associations for very rare variants requires extremely large sample sizes or very large effect sizes. With N=1196, power for variants with MAF < 0.01 is generally low unless the effect size is substantial.
Genotyping Accuracy: Very rare variants can sometimes be enriched for genotyping errors.
Allele Frequency Estimation: Frequencies for very rare alleles are estimated with less precision.
Multiple Testing: Removing these low-MAF SNPs reduces the total number of tests performed in the GWAS.

Choosing a MAF threshold is a balance. Too strict (e.g., MAF 0.10) might remove potentially interesting, less common variants. Too lenient (e.g., MAF 0.001) might retain many SNPs with little power, increasing computational burden and multiple testing. For a study of this size, MAF 0.01 is a standard and defensible choice. MAF 0.05 is also common and more conservative (removes more SNPs).

*356459 variants removed due to minor allele threshold(s)*
*949993 variants and 1196 people pass filters and QC*

# Hardy-Weinberg Equilibrium (HWE) Filter

*Purpose*
To identify and remove SNPs that show a significant deviation from Hardy-Weinberg Equilibrium (HWE). HWE describes a principle stating that allele and genotype frequencies in a population will remain constant from generation to generation in the absence of other evolutionary influences (like mutation, selection, genetic drift, gene flow, and non-random mating). Significant deviation from HWE for a particular SNP can indicate:
*   Genotyping errors for that SNP.
*   Population stratification (though PCA covariates aim to handle broader stratification).
*   True biological selection at or near the locus.
*   Issues with copy number variations affecting genotyping.

In GWAS QC, the primary concern is often genotyping error, so SNPs deviating significantly from HWE are typically removed.

*Methodology*

PLINK's `--hwe <p-value threshold>` command is used to calculate HWE p-values for each SNP and filter out those below the specified significance threshold. The test is usually performed on a specific group of individuals, ideally **unrelated controls** if conducting a case-control study. If you don't have distinct cases and controls (as with eye color, which is a trait, not a disease status), or if defining controls is complex, the test can be performed on **all founders** or **all individuals** in the current sample set, while being mindful that true population structure can cause HWE deviations. Given you've handled relatedness and will use PCA for structure, applying it to all current individuals is a common approach for trait GWAS.

We ran the following PLINK command to calculate out HWE p values for all our SNPS, which we then took the output file from and visualized the different thresholds in R:

plink --bfile final_analysis_cohort_maf_filtered \
      --hardy \
      --out current_cohort_hwe_stats

From the visualization we chose to stick with 1e-6 as our p value. 

Input Dataset for HWE Filtering:
*   `final_analysis_cohort_maf_filtered.bed/bim/fam`
    *   This dataset contains: 1196 individuals and 949,993 SNPs (after MAF filtering).

*Step 1: Apply HWE Filter using PLINK*

A common p-value threshold for HWE filtering in GWAS is `1e-6`. This is a relatively stringent threshold, aiming to remove only SNPs with strong evidence of HWE deviation.

Run PLINK to filter SNPs based on HWE:

plink --bfile final_analysis_cohort_maf_filtered \
      --hwe 1e-6 \
      --make-bed \
      --out final_analysis_cohort_hwe_filtered


Parameters:
--bfile final_analysis_cohort_maf_filtered: Specifies the input dataset.
--hwe 1e-6: Sets the HWE p-value significance threshold. SNPs with a p-value below 1e-6 (i.e., more significant deviation) will be removed.

By default, PLINK applies the HWE test to founders, or all individuals if founder status isn't clear or relevant (which is likely your case here). You can add options like --hwe-all to force it on all, or specify a case/control context if you had one. For your trait data, PLINK's default behavior of testing across all unrelated individuals is generally appropriate.
--make-bed: Instructs PLINK to output the filtered dataset in binary PED format.
--out final_analysis_cohort_hwe_filtered: Specifies the prefix for the output files.

Outcome
A new PLINK binary fileset: final_analysis_cohort_hwe_filtered.bed/bim/fam.
This dataset will contain the same 1196 individuals but may have fewer SNPs than the input dataset, as SNPs significantly deviating from HWE have been removed.
The final_analysis_cohort_hwe_filtered.log file will detail:
The number of SNPs initially present.
The number of SNPs removed due to the HWE filter.
The number of SNPs remaining.

Action: Review the PLINK log file (final_analysis_cohort_hwe_filtered.log) to note how many SNPs were removed and how many remain. Typically, only a small fraction of SNPs are removed by HWE filtering if prior QC has been good and the MAF filter has been applied.

Considerations for HWE Filtering:
Threshold Choice: 1e-6 is common. More lenient thresholds (e.g., 1e-4 or 1e-3) would remove more SNPs but might also remove SNPs deviating due to minor population structure rather than just error. Stricter thresholds (e.g., 1e-10) remove fewer.
Population Context: HWE assumptions are best met in large, randomly mating populations. Your openSNP cohort is diverse. While PCA covariates will adjust for major population structure in the association analysis, some residual fine-scale structure or admixture could still cause SNPs to deviate from HWE. The stringent 1e-6 threshold helps mitigate removing SNPs for these reasons, focusing more on potential genotyping errors.
Case/Control Studies: In case-control studies, HWE is typically only tested in the control group, as allele frequencies (and thus HWE) in cases might be skewed by the disease association itself for true causal variants. For a quantitative trait like eye color, testing on all (unrelated) individuals is standard.


*949579 variants and 1196 people pass filters and QC.* 


# GWAS design

With our phenotypes and our 4 category eye color model a solid and "simple" start point is for us to treat this as a quantitative trait, which opens the door to Linear Regression. This essentially is done by telling PLINK to treat our phenotype scale, which is coded as 0-3, as a quantitative score, which PLINK can easily use to run linear regression on for each SNP: Phenotype_Score ~ Genotype + Covariates 

This model is very easy to implement, and due to previous liteature suggesting that there "somewhat" a linear relationship between the genetic effect and the ordered phenotype (melanin amount = eye color shade), we can assume that this model would have enough power. However, this assumes that the "distance" between Blue and Green/Blue-Green is the same as the "distance" between Green/Blue-Green and Hazel, which we cant confirm is biologically true. Violations of this assumption can reduce power or lead to misinterpretation.  

This can be our initial GWAS design, and for later additional analysis, one option would be to convert our data so that it is compatible with PLINK2, which would allow us to do Ordinal Logistic Regression, which is more statistically fitting for our model, due to the ordered categorical data. It does not assume equal spacing, but it respects the order (our order being light to dark). We find this appropriate due to previous studies which suggest that the "lightness" of eyecolor is tied to the individuals melanin production, meaning that more melanin production equals more eye pigment, which means darker eyes.  

If we wanted to start at a more natural point for our further analysis, we could skip assuming order of the eye colors. This would allow us to do Multinomial Logistic Regression. This compares each eyecolor category to a reference category. However, we feel that since we have solid evidence suggesting that there is an order to eyecolor, this is not the model we want to go with. 


# GWAS 4 category Quantitative Trait Linear Model

The Model: For each of your ~950,000 SNPs, PLINK is fitting a linear regression model:
Y = β₀ + β₁ * SNP + β₂ * PC1 + β₃ * PC2 + β₄ * PC3 + β₅ * PC4 + ε

Where:
Y is your eye color score (0-3).
SNP is the genotype of the individual at that SNP (coded as 0, 1, or 2 copies of the tested allele).
PC1 to PC4 are the values of the principal components for that individual.
β₀ is the intercept.
β₁ is the effect size of the SNP on eye color – this is what you're most interested in.
β₂ to β₅ are the effects of the PCs.
ε is the error term.

The Hypothesis Test: For each SNP, PLINK tests the null hypothesis that β₁ = 0 (i.e., the SNP has no effect on eye color, after accounting for PCs). The p-value tells you the probability of observing an effect as large as or larger than β₁ if the null hypothesis were true.

*Purpose*
To perform a genome-wide association study to identify genetic variants associated with the 4-category eye color phenotype. A linear regression model will be used, treating the eye color categories as a quantitative score (0-3). The model will include the top 4 Principal Components (PCs) and sex as covariates to control for population stratification and potential sex-specific effects, respectively.

*Input Files*
1.  **Genotypes:** `final_analysis_cohort_hwe_filtered.bed/bim/fam`
    *   Contains: 1196 individuals, 949,579 SNPs. This is the fully QCed dataset.
2.  **Phenotypes (Original Full):** `eye_color_updated.txt`
    *   Contains: FID, IID, and the 4-category eye color phenotype score (0-3) for potentially more than 1196 individuals.
3.  **Principal Components (Eigenvectors):** `final_pca_results_no_outliers.eigenvec`
    *   Contains: FID, IID, PC1, PC2, ..., PC10 for the 1196 individuals.
4.  **FAM file for Sex Information:** `final_analysis_cohort_hwe_filtered.fam`
    *   Contains: FID, IID, and sex codes (1=male, 2=female) for the 1196 individuals.

*Methodology - File Preparation (R)*

To ensure clarity and that only the 1196 individuals in the final QCed dataset are used with corresponding phenotypes and covariates, new subsetted phenotype and combined covariate files will be created using R. the R script XXX creates the new phenotype file: *gwas_pheno_1196.txt* and the new combined file: *gwas_covariates_pcs_sex_1196.txt*

In these files, Sex was coded as 0 for males (FAM code 1) and 1 for females (FAM code 2). Two individuals had an undefined sex code (0) in the final FAM file and were consequently assigned an NA for the sex covariate, despite previous sex checks and filtering. These two individuals were therefore excluded from the GWAS analysis by PLINK due to missing covariate data, resulting in an effective sample size of 1194 for the association tests.

*Methodology - GWAS Execution (PLINK)*
A linear regression model will be run for each SNP using PLINK 1.9.

PLINK Command:

plink --bfile final_analysis_cohort_hwe_filtered \
      --pheno gwas_pheno_1196.txt \
      --covar gwas_covariates_pcs_sex_1196.txt \
      --linear \
      --allow-no-sex \
      --out gwas_eyecolor_linear_4pcs_sex


Breakdown of PLINK command:
--bfile final_analysis_cohort_hwe_filtered: Specifies the final QCed genotype dataset (1196 individuals, 949,579 SNPs).
--pheno gwas_pheno_1196.txt: Specifies the subsetted phenotype file containing the 4-category eye color score for the 1196 individuals.
--covar gwas_covariates_pcs_sex_1196.txt: Specifies the combined covariate file (top 4 PCs and sex) for the 1196 individuals.
--linear: Instructs PLINK to perform linear regression. The model for each SNP is:
EyeColorScore ~ SNP_Genotype + PC1 + PC2 + PC3 + PC4 + Sex + Intercept.
The additive genetic model (ADD) is tested by default for SNPs.
--allow-no-sex: Included for consistency, although sex is now explicitly used as a covariate and should be defined for all individuals in the covariate file. PLINK's primary use of this flag is for phenotype inclusion when sex in the FAM file is '0'; here, the covariate file dictates sex usage in the model.
--out gwas_eyecolor_linear_4pcs_sex: Specifies the prefix for the output files.


# GWAS Results Visualization: Manhattan and QQ Plots

*Purpose*
To visualize the results of the genome-wide association study (GWAS) for eye color. Two main plots will be generated:
1.  **Manhattan Plot:** Displays the strength of association (-log10 P-value) for each SNP across the genome, organized by chromosome and position. It helps to quickly identify chromosomal regions harboring SNPs with significant associations.
2.  **QQ (Quantile-Quantile) Plot:** Compares the observed distribution of GWAS P-values against the expected uniform distribution under the null hypothesis of no association. It helps to assess overall inflation or deflation of test statistics, which can indicate issues like uncorrected population stratification, cryptic relatedness, or other systematic biases.

*Input File*
*   `gwas_eyecolor_linear_4pcs_sex.assoc.linear`: The output file from the PLINK linear regression, containing association statistics for each SNP (CHR, SNP, BP, P, etc.).

*Methodology - Plotting in R*

We started by loading our data into the R script gwas_linear_result_visualizations.rmd, and generated a QQ plot, our lamda and a Manhattan plot.

*Discussion*

When we visualized our resluts, and calculated our initial lamda, we already could see problems. With a lambda of almost 10, an unreadable QQ plot and a Manhattan plot with only one hit, we already know that something is causing problems. 

The high lambda would, under normal circumstances, cause many peaks or hits to manifest on a Manhattan plot. However, when visualizing our data on a manhattan plot, we only had one hit/peak. This is often called the "winner takes it all" effect, and is usually caused by Extreme Stratification for a Highly Differentiated Trait. Eye color is a trait with extremely strong genetic differentiation that aligns with major ancestral populations (e.g., HERC2/OCA2 variants are at very different frequencies in European vs. African vs. East Asian populations).

The hit in the Manhattan plot is actually a pretty big indicator of what went wrong. The hit is located on chromosome 15, which contains the HERC2/OCA2 region. The variants on this locus are major determinants of blue vs. brown eyes and have vastly different frequencies across global populations. Thus, this is an indicate that we have not correctly controlled for population stratification. 

So, our next step is essentually to try and better control for population stratification.



# Addressing Population Stratification in GWAS

*Purpose*
The initial GWAS run (using 4 PCs + sex as covariates) revealed severe genomic inflation (Lambda ≈ 9.7), as evidenced by the QQ plot. This indicates that population stratification is not adequately controlled. The following steps aim to improve stratification control by including more principal components (PCs) in the regression model.

*Underlying Theory*
Population stratification occurs when allele frequencies and phenotype distributions differ systematically across underlying ancestral subgroups within the study sample. If not accounted for, SNPs that are merely markers of ancestry (rather than being causally related to the trait) can show spurious associations. Principal Components (PCs) derived from genome-wide SNP data capture major axes of genetic variation, which often correspond to ancestral differences. Including these PCs as covariates in the GWAS regression model statistically adjusts for these broad ancestral effects, reducing false positives due to stratification. The goal is to include enough PCs to capture the relevant population structure that correlates with the trait, thereby reducing the genomic inflation factor (lambda) to an acceptable level (ideally close to 1.0).

*Methodology - Iterative GWAS with Increased PCs*

We will re-run the GWAS, systematically increasing the number of PCs used as covariates, and monitor the lambda value after each run. The PCs are taken from the `final_pca_results_no_outliers.eigenvec` file (which contains up to 10 PCs for the 1196 individuals who passed all QC prior to the GWAS run).

**Key Files (already prepared or used previously):**
*   **Genotypes:** `final_analysis_cohort_hwe_filtered.bed/bim/fam` (1196 individuals, 949,579 SNPs)
*   **Phenotypes (Subsetted):** `gwas_pheno_1196.txt` (1196 individuals, created in previous R step)
*   **Full PCA Eigenvectors:** `final_pca_results_no_outliers.eigenvec` (1196 individuals, 10 PCs)
*   **FAM file for Sex:** `final_analysis_cohort_hwe_filtered.fam` (1196 individuals)

**Step 1: Prepare New Covariate Files with More PCs (Using R)**

We will create covariate files for, e.g., 10 PCs + Sex. We might also consider trying 6, 8, 15, or 20 PCs if needed, based on your scree plot and the results of these iterations. We did this in R, with some simple commands and data filtering. The new file being used as covariates is *gwas_covariates_10pcs_sex_1194.txt*


**Step 2: Re-run GWAS with 10 PCs + Sex**

PLINK Command:

plink --bfile final_analysis_cohort_hwe_filtered \
      --pheno gwas_pheno_1196.txt \
      --covar gwas_covariates_10pcs_sex_1194.txt \
      --linear \
      --allow-no-sex \
      --out gwas_eyecolor_linear_10pcs_sex

*Action:* Note the new output prefix gwas_eyecolor_linear_10pcs_sex.

**Step 3: Evaluate Lambda and QQ Plot for the New Results**

Use the R script previously provided for visualization, but point it to the new results file:
gwas_results_file <- "gwas_eyecolor_linear_10pcs_sex.assoc.linear"

Using 10 PCs instead of 4, we see a drastic inprovement in our lambda (now 2.811), however, our QQ plot essentially looks unchanged. We will try using more PCs to see if that yields a positiver result, if not, we have to move on to a different approach. 

However, turns out that we are using the wrong PC data as our covariates. The PCA data was calculated before the final filtering of related individuals, which means that they *might* be pulling/effecting the way the PCs are supposed to correct our data when they act as covariates. So before anything, we should try and recalculate the PCAs, this time 30 of them, then do the GWAS again and see if that has an effect. 

# Addressing Population Stratification in GWAS (Revised Approach)

*Purpose*
The initial GWAS runs (using 4 and subsequently 10 PCs + sex as covariates) revealed severe and persistent genomic inflation (Lambda remaining high, e.g., 2.811 with 10 PCs), with a QQ plot shape indicative of uncontrolled population stratification. A critical review revealed that the Principal Components (PCs) used were derived from a cohort of 1412 individuals (post-PCA outlier removal but pre-relatedness filtering), while the GWAS was performed on a subsequent cohort of 1196 individuals (after relatedness filtering). Using PCs derived from a different sample set than the analysis set is suboptimal for stratification control.

The following steps aim to:
1.  Generate new PCs based on the correct final set of 1196 unrelated individuals.
2.  Use an appropriate number of these new PCs to control for population stratification.
3.  Re-run the GWAS and evaluate inflation.

*Underlying Theory*
Principal Components (PCs) must be calculated on the specific set of (preferably unrelated) individuals included in the association analysis to accurately capture and correct for population structure within that cohort. Using PCs from a different or larger ancestral pool can lead to incomplete or inappropriate adjustment.

*Methodology - Generating Correct PCs and Iterative GWAS*

**Key Files:**
*   **Current Genotype Dataset for GWAS (but PCs need re-derivation):** `final_analysis_cohort_hwe_filtered.bed/bim/fam` (1196 individuals, 949,579 SNPs). This is the dataset *after* relatedness removal and all SNP QC (MAF, HWE).
*   **Phenotypes (Subsetted for 1196):** `gwas_pheno_1196.txt`
*   **FAM file for Sex (for 1196):** `final_analysis_cohort_hwe_filtered.fam`

**Step 1: Prepare an LD-Pruned Dataset for the Correct 1196 Individuals**

The dataset `final_analysis_cohort_hwe_filtered` contains all SNPs. For PCA, we need an LD-pruned version of *these 1196 individuals*.

Input: final_analysis_cohort_hwe_filtered (1196 individuals, ~950k SNPs)

1a. Exclude known high-LD regions

plink --bfile final_analysis_cohort_hwe_filtered \
      --exclude range high_ld_regions_exclude.txt \
      --make-bed \
      --out gwas_1196_ldregions_excluded

1b. Perform window-based LD pruning

plink --bfile gwas_1196_ldregions_excluded \
      --indep-pairwise 50 5 0.2 \
      --out gwas_1196_pruning_results

plink --bfile gwas_1196_ldregions_excluded \
      --extract gwas_1196_pruning_results.prune.in \
      --make-bed \
      --out gwas_1196_final_ld_pruned_for_pca

Output: gwas_1196_final_ld_pruned_for_pca.bed/bim/fam. This file now contains the 1196 individuals and an LD-pruned set of SNPs derived specifically from them. Note the number of SNPs remaining.

**Step 2: Re-calculate PCA on the Correct LD-Pruned 1196-Individual Dataset (Generate e.g., 30 PCs)**

plink --bfile gwas_1196_final_ld_pruned_for_pca \
      --pca 30 \
      --out gwas_1196_pca_30pcs

Outputs: gwas_1196_pca_30pcs.eigenvec (PCs for 1196 individuals) and gwas_1196_pca_30pcs.eigenval.

**Step 3: Create a New Scree Plot (up to 30 PCs for the 1196 Individuals)**

Use R with gwas_1196_pca_30pcs.eigenval to visualize variance explained. Examine this scree plot carefully to decide how many PCs (N_new_pcs) to test as covariates (e.g., start with 10, then try 15, 20 based on this plot). We chose 7.

**Step 4: Prepare New Covariate File (using N_NEW_PCS_TO_TRY from gwas_1196_pca_30pcs.eigenvec + Sex)**

Useing our covariate_creation_for_popstrat R script, we made a new file using that would be used as covariates, from the 7 PCs chosen plus sex. The output file is gwas_covariates_7pcs_sex_1194.txt

**Step 5: Re-run GWAS with the Corrected and Potentially Increased Number of PCs + Sex**

plink --bfile final_analysis_cohort_hwe_filtered \
      --pheno gwas_pheno_1196.txt \
      --covar gwas_covariates_7pcs_sex_1194.txt \
      --linear \
      --allow-no-sex \
      --out gwas_eyecolor_linear_correct_7pcs_sex

**Step 6: Evaluate Lambda and QQ Plot (Iteratively)**

Using our gwas_linear_result_visualizations.rmd we will visualize the QQ plot again and calculate our lambda


# Restricting Analysis to a More Homogeneously Ancestral Subgroup

*Purpose*
Previous GWAS attempts on the full cohort of 1196 individuals, even with up to 20 Principal Components (PCs) as covariates, resulted in unacceptably high genomic inflation (Lambda ≈ 1.973) and a QQ plot shape indicating persistent, uncontrolled population stratification. This is likely due to the significant ancestral diversity within the openSNP cohort and the strong correlation between ancestry and eye color.

To achieve better control of population stratification, the analysis will now be restricted to a more ancestrally homogeneous subgroup identified from the principal components. This approach aims to reduce the extent of confounding by focusing on a set of individuals who are genetically more similar.

*Input Files for this Section*
1.  **PCA Eigenvectors (for 1196 individuals, up to 30 PCs):** `pca_results_1196ind_30pcs.eigenvec`
2.  **Phenotype Data (for 1196 individuals):** `gwas_pheno_1196.txt` (or your original `eye_color_updated.txt` which PLINK can subset)
3.  **Fully QCed Genotype Data (for 1196 individuals, pre-GWAS):** `final_analysis_cohort_hwe_filtered.bed/bim/fam`

*Methodology*

**Step 1: Visualize Current PCA Results to Define the Core Subgroup (R)**

The first step is to plot the principal components from `pca_results_1196ind_30pcs.eigenvec` to visually identify the main, dense cluster of individuals. We will also color by phenotype to confirm your earlier observation about the "arms."

*Decision on Core Subgroup Definition*
Based on visualizing the PCA plots (`pca_results_1196ind_30pcs.eigenvec`) and iteratively applying standard deviation (SD) based thresholds, a core subgroup was defined by selecting individuals within +/- 2.0 SD of the mean for both PC1 and PC2. This resulted in the selection of 1124 individuals. The list of these individuals (FID and IID) has been saved to `gwas_core_cluster_sd_individuals_to_keep.txt`. Seventy-one (71) individuals were excluded by this criterion.

*Next Steps: Prepare Core Subgroup for GWAS*
The following steps will create a new dataset for these 1124 individuals and re-calculate Principal Components specifically for this subgroup to be used as covariates.

**Input Files for this Stage:**
1.  **Fully QCed Genotype Data (for 1196 individuals, pre-GWAS):** `final_analysis_cohort_hwe_filtered.bed/bim/fam`
2.  **List of Core Individuals to Keep:** `gwas_core_cluster_sd_individuals_to_keep.txt` (contains FIDs/IIDs for 1124 individuals)
3.  **List of High-LD Regions:** `high_ld_regions_exclude.txt`
4.  **Phenotype Data (Subsetted for 1196, or original):** `gwas_pheno_1196.txt` (or `eye_color_updated.txt`)
5.  **FAM file for Sex (for 1196):** `final_analysis_cohort_hwe_filtered.fam`

**Step 2: Create New PLINK Dataset for the Core Subgroup (1124 Individuals)**

plink --bfile final_analysis_cohort_hwe_filtered \
      --keep gwas_core_cluster_sd_individuals_to_keep.txt \
      --make-bed \
      --out gwas_core_1124

**step 3: LD Pruning for PCA on the Core Subgroup**

3a. Exclude known high-LD regions from the core subgroup dataset

plink --bfile gwas_core_1124 \
      --exclude range high_ld_regions_exclude.txt \
      --make-bed \
      --out gwas_core_1124_ldregions_excluded

3b. Perform window-based LD pruning on the core subgroup dataset

plink --bfile gwas_core_1124_ldregions_excluded \
      --indep-pairwise 50 5 0.2 \
      --out gwas_core_1124_pruning_results

plink --bfile gwas_core_1124_ldregions_excluded \
      --extract gwas_core_1124_pruning_results.prune.in \
      --make-bed \
      --out gwas_core_1124_ld_pruned_for_pca

Output: gwas_core_1124_ld_pruned_for_pca.bed/bim/fam. This is the input for the new PCA.
Action: Note the number of SNPs remaining after pruning for this core subgroup.

**Step 4: Calculate New Principal Components for the Core Subgroup**

plink --bfile gwas_core_1124_ld_pruned_for_pca \
      --pca 20 \
      --out gwas_core_1124_pca_20pcs

Outputs: gwas_core_1124_pca_20pcs.eigenvec and gwas_core_1124_pca_20pcs.eigenval.
These are the NEW PCs that will be used as covariates for the GWAS on the core subgroup.

**Step 5: Create Scree Plot for the New Core Subgroup PCs**

Using the screeplot R script, we identified the first 6 PCs to be used in later steps.

**Step 6: Prepare Final Phenotype and Covariate Files for the Core Subgroup GWAS**

Using the covariate_file_creation R script, we created a new .txt file containing the individuals we are moving on with. Again, two of them still have missing sex so they are removed. The txt file is: gwas_covariates_6pcs_sex_1124.txt

**Step 7: Run GWAS on the Core Subgroup**

plink --bfile gwas_core_1124 \
      --pheno gwas_pheno_1196.txt \
      --covar gwas_covariates_6pcs_sex_1124.txt \
      --linear \
      --allow-no-sex \
      --out gwas_eyecolor_CORE_SUBGROUP_6pcs_sex


**Step 8: Evaluate Lambda and QQ Plot for the Core Subgroup GWAS**

Lambda is down to 1.859, but the shapw of the QQ plot persists. The hit on the Manhattan plot is still just one peak on chromsome 15, most likely the HERC2/OCA2 variants.



# Plan A: Go with the results we have now from the 6 and 20 PC's QQ and Man. 

Okay, thank you for sharing the 20 PC QQ plot (Lambda = 1.973).
You are absolutely right in your assessment:
The lambda has improved (down from 6.6 with 7 PCs on this same core, and from 9.7 initially).
However, the shape of the QQ plot is still fundamentally problematic. That early, strong deviation and the "shelf/step-ladder" pattern persist. The "little more hugging in the beginning" is a very minor improvement in the face of the overall distortion.
Your observation "It seems there are just a few outliers that completely ruin the plot right?" is insightful regarding the tail of the distribution (the extremely significant SNPs, likely HERC2/OCA2). But the problem is actually genome-wide: the entire distribution of p-values is shifted upwards from the null line long before you get to that extreme tail. The "shelves" indicate that large numbers of SNPs are more significant than they should be.
Is it realistic to control this population stratification in this current 1124-person "core" subgroup with only PCs and a simple linear model?
At this point, with lambda still ~2.0 and that QQ plot shape after 20 PCs, it's becoming highly unlikely that simply adding more PCs (e.g., 25 or 30) to this specific 1124-person subgroup will magically fix the QQ plot to an ideal state (lambda ~1.0, clean diagonal). The confounding is too strong and too structured.
You are right: you've iterated, and the standard PC approach on this particular "core" subset isn't fully resolving the deep-seated stratification issues for a trait as ancestrally patterned as eye color in a diverse public dataset.
"We need results NOW." - Understood.
Given the time constraints and the need for interpretable results for your project, continuing to chase perfect stratification control in this current 1124-person "core" with more PCs might be a frustrating exercise with diminishing returns.
The Most Pragmatic Path to "Results NOW" that are as interpretable as possible given the data's nature:
We need to accept that perfect control might be very hard with this dataset and model, but we still want to present the most reasonable analysis. Your primary goal for the project is likely to demonstrate the GWAS process, identify major signals, and discuss limitations.
Strategy: Focus on the Clearest Known Biology and Acknowledge Limitations
The HERC2/OCA2 region on Chromosome 15 is the undeniable major player in eye color variation, especially blue vs. brown. Your Manhattan plot (even the one from the 6 PC run) clearly shows this. This is a true biological signal, even if its p-value is currently inflated.
What to do for your report (balancing rigor with project deadlines):
Since Path B (More Stringent Definition of the "Core" Subgroup) is the most methodologically sound next step to try and improve lambda/QQ plot if you had infinite time, but you need results now:
Option A: Proceed with current best results (20 PCs on 1124 core) and HEAVILY CAVEAT.
This is risky because a lambda of ~2 is generally not publishable and indicates unreliable p-values genome-wide. However, if you are forced to produce some result from the current state:
Use the results from your GWAS on the gwas_core_1124 individuals with 20 PCs + Sex (Lambda = 1.973).
Output file: gwas_eyecolor_CORE_SUBGROUP_20pcs_sex.assoc.linear (or similar based on your naming).
In your report:
Show the QQ plot (Lambda = 1.973).
Explicitly state that there is significant residual genomic inflation (Lambda = 1.973), indicating that population stratification is not fully controlled even with 20 PCs in this subgroup. This is a critical acknowledgement.
Show the Manhattan plot. State that due to the inflation, p-values genome-wide are likely anti-conservative.
Focus your "findings" almost exclusively on the Chromosome 15 peak (HERC2/OCA2).
Report the top SNP(s) in this region.
State that this is a well-known major locus for eye color, and your GWAS successfully identifies it.
Acknowledge that the p-value for this peak is likely inflated but the signal itself is consistent with known biology.
Do NOT over-interpret other "suggestive" peaks across the genome. With this lambda, they are highly likely to be false positives.
In your Discussion/Conclusion: Discuss the challenges of controlling for population stratification in diverse cohorts like openSNP, especially for traits like eye color. Mention that further refinement of the cohort (e.g., even stricter ancestral homogeneity) or more advanced models (like mixed models) would be needed for a more robust analysis.


# Plan B: Implement a Quick Version of "More Stringent Core" (If you have a few hours)

This has a chance of improving lambda/QQ plot enough to be more presentable.
Define a "Super-Core":
Go back to your R script for defining the core using new_pca_results_1196ind_30pcs.eigenvec.
Use a stricter SD_MULTIPLIER. Try SD_MULTIPLIER = 1.5 first. WENT TO 0.5
This will give you a smaller list of individuals: gwas_super_core_individuals_to_keep.txt. Note the new N (e.g., N_super_core). Let's say it's ~800-900.
Create the Super-Core PLINK Dataset:
plink --bfile final_analysis_cohort_hwe_filtered \
      --keep gwas_super_core_sd0.5_individuals_to_keep.txt \
      --make-bed \
      --out gwas_super_core_N${N_super_core}

Re-run MAF and HWE on this Super-Core:
plink --bfile gwas_super_core_N${N_super_core} --maf 0.01 --make-bed --out gwas_super_core_N${N_super_core}_maf
plink --bfile gwas_super_core_N${N_super_core}_maf --hwe 1e-6 --make-bed --out gwas_super_core_N${N_super_core}_finalqc

Re-run PCA on this Super-Core LD-pruned data (calculate e.g., 10-15 PCs):
First LD prune gwas_super_core_N${N_super_core}_finalqc (exclude regions, then --indep-pairwise). Let the output be gwas_super_core_N${N_super_core}_ld_pruned_for_pca.

plink --bfile gwas_super_core_N${N_SUPER_CORE_VAR}_finalqc \
      --exclude range high_ld_regions_exclude.txt \
      --make-bed \
      --out gwas_super_core_N${N_SUPER_CORE_VAR}_finalqc_ldregions

plink --bfile gwas_super_core_N${N_SUPER_CORE_VAR}_finalqc_ldregions \
      --indep-pairwise 50 5 0.2 \
      --out gwas_super_core_N${N_SUPER_CORE_VAR}_pruning_results

plink --bfile gwas_super_core_N${N_SUPER_CORE_VAR}_finalqc_ldregions \
      --extract gwas_super_core_N${N_SUPER_CORE_VAR}_pruning_results.prune.in \
      --make-bed \
      --out gwas_super_core_N${N_SUPER_CORE_VAR}_ld_pruned_for_pca

Then PCA:
plink --bfile gwas_super_core_N${N_super_core}_ld_pruned_for_pca \
      --pca 15 \
      --out gwas_super_core_N${N_super_core}_pca_15pcs

Prepare Covariates: Use R to make a covariate file with a small number of these new super-core PCs (e.g., top 4-6 PCs from gwas_super_core_N${N_super_core}_pca_15pcs.eigenvec) + Sex.

Run "r_script_of_death.rmd".

Ensure N_SUPER_CORE_VAR matches N_super_core used in R
Ensure N_SUPER_CORE_PCS_FILENAME matches N_SUPER_CORE_PCS_TO_USE in R
Ensure N_SUPER_CORE_ANALYZED_FILENAME matches num_in_final_covar_file from R output

N_SUPER_CORE_VAR=1048 # Example
N_SUPER_CORE_PCS_FILENAME=6 # Example
N_SUPER_CORE_ANALYZED_FILENAME=1046 # Example (897 - 2 with NA sex)

plink --bfile gwas_super_core_N${N_SUPER_CORE_VAR}_finalqc \
      --pheno gwas_super_core_pheno_N${N_SUPER_CORE_VAR}.txt \
      --covar gwas_super_core_covar_${N_SUPER_CORE_PCS_FILENAME}pcs_sex_N${N_SUPER_CORE_ANALYZED_FILENAME}.txt \
      --linear \
      --allow-no-sex \
      --out gwas_eyecolor_SUPERCORE_${N_SUPER_CORE_PCS_FILENAME}pcs_sex

Run GWAS on Super-Core:
plink --bfile gwas_super_core_N${N_super_core}_finalqc \
      --pheno gwas_super_core_pheno_N${N_super_core}.txt \
      --covar gwas_super_core_covar_pcs_sex_N.txt \
      --linear \
      --allow-no-sex \
      --out gwas_eyecolor_SUPERCORE_FewPCs_sex

Check Lambda and QQ plot. (r_script_of_death_vis.rmd)
If lambda is now much better (e.g., < 1.2) and QQ plot is cleaner: This is your best result to report!


If still problematic: Then you've done due diligence, and you'll have to report results like in Option A, but mention you also tried this more stringent subsetting.
Regarding Removing the "Winner Takes It All" SNPs (Chromosome 15 peak):
Do NOT do this as a primary analysis strategy to "fix" your QQ plot or lambda.
Removing true strong signals because they are "too strong" is not how you control for genome-wide inflation caused by stratification.
This technique (conditional analysis or excluding a region) is an exploratory tool used after you have a well-controlled GWAS to see if secondary, independent signals exist.
If you remove the Chr15 peak now, your lambda calculated on the remaining SNPs will still be inflated if the underlying stratification problem isn't solved, and you won't be able to trust other emerging peaks.
Decision and Exact Steps "NOW":
Given the need for results and the QQ plot shape:
I strongly recommend attempting a quick version of Path B (Option B above: More Stringent "Super-Core").
Exact Steps for "Super-Core" Approach:
Decide on SD_MULTIPLIER for Super-Core:
In your R script that uses new_pca_results_1196ind_30pcs.eigenvec, set SD_MULTIPLIER = 1.5.
Run it to see how many individuals are kept. If it's too few (e.g., < 600-700), try SD_MULTIPLIER = 1.75. Aim for a balance.
Once you have a multiplier that gives you a "super-core" of, say, N=700-900 individuals, save the list as gwas_super_core_individuals_to_keep.txt. Let's say N_super_core is this new number.
Create Super-Core PLINK Dataset:
plink --bfile final_analysis_cohort_hwe_filtered \
      --keep gwas_super_core_individuals_to_keep.txt \
      --make-bed \
      --out gwas_super_core_N${N_super_core}
Use code with caution.
Bash
(Replace ${N_super_core} with the actual number)
MAF Filter for Super-Core:
plink --bfile gwas_super_core_N${N_super_core} \
      --maf 0.01 \
      --make-bed \
      --out gwas_super_core_N${N_super_core}_maf
Use code with caution.
Bash
HWE Filter for Super-Core:
plink --bfile gwas_super_core_N${N_super_core}_maf \
      --hwe 1e-6 \
      --make-bed \
      --out gwas_super_core_N${N_super_core}_finalqc
Use code with caution.
Bash
(This gwas_super_core_N${N_super_core}_finalqc is your genotype file for the super-core GWAS)
LD Pruning for Super-Core PCA:
plink --bfile gwas_super_core_N${N_super_core}_finalqc \
      --exclude range high_ld_regions_exclude.txt \
      --make-bed \
      --out gwas_super_core_N${N_super_core}_finalqc_ldregions

plink --bfile gwas_super_core_N${N_super_core}_finalqc_ldregions \
      --indep-pairwise 50 5 0.2 \
      --out gwas_super_core_N${N_super_core}_pruning_results

plink --bfile gwas_super_core_N${N_super_core}_finalqc_ldregions \
      --extract gwas_super_core_N${N_super_core}_pruning_results.prune.in \
      --make-bed \
      --out gwas_super_core_N${N_super_core}_ld_pruned_for_pca
Use code with caution.
Bash
PCA on Super-Core (Calculate 10-15 PCs):
plink --bfile gwas_super_core_N${N_super_core}_ld_pruned_for_pca \
      --pca 15 \
      --out gwas_super_core_N${N_super_core}_pca_15pcs
Use code with caution.
Bash
Prepare Phenotype and Covariate Files for Super-Core GWAS (R):
Create gwas_super_core_pheno_N${N_super_core}.txt by subsetting eye_color_updated.txt for the N_super_core individuals.
Create gwas_super_core_covar_Xpcs_sex.txt using a small number of the new PCs from gwas_super_core_N${N_super_core}_pca_15pcs.eigenvec (e.g., top 4-6 PCs based on its scree plot) + Sex. Remember this will likely be for N_super_core - 2 individuals if 2 still have undefined sex. Let N_super_core_analyzed = N_super_core - 2.
Run GWAS on Super-Core:
 Example using 4 new PCs from super-core PCA
plink --bfile gwas_super_core_N${N_super_core}_finalqc \
      --pheno gwas_super_core_pheno_N${N_super_core}.txt \
      --covar gwas_super_core_covar_4pcs_sex_N${N_super_core_analyzed}.txt \
      --linear \
      --allow-no-sex \
      --out gwas_eyecolor_SUPERCORE_4pcs_sex
Use code with caution.
Bash
CRITICALLY EVALUATE Lambda and QQ Plot for this gwas_eyecolor_SUPERCORE_4pcs_sex.assoc.linear result.
This "super-core" approach is your most promising bet for getting a reasonably well-controlled result quickly. If this still results in a very bad lambda/QQ plot, then the dataset has fundamental challenges for a simple linear model GWAS of eye color that might require even more advanced techniques or different study designs not feasible for this project. But let's try this path first.
If the "super-core" GWAS results in a lambda around 1.0-1.1 and a much cleaner QQ plot, then THOSE are the results you should focus on for your report. You would then explain the steps taken to arrive at this more homogeneous subgroup.

# Conclusion of GWAS of Eye Color: Results and Interpretation

**Objective:** To identify genetic variants associated with a 4-category self-reported eye color phenotype (Blue, Green/Blue-Green, Hazel, Brown, scored 0-3) in a diverse cohort from openSNP, while controlling for population stratification.

**Summary of Quality Control and Cohort Definition:**
The initial dataset comprised [Original N from gwas_data.fam] individuals. A rigorous quality control (QC) pipeline was implemented. Key steps included:
1.  **Phenotype Harmonization:** Self-reported eye colors were mapped to a 4-category ordered scale (0=Blue, 1=Green/Blue-Green, 2=Hazel, 3=Brown).
2.  **Per-Chip QC:** Individuals were grouped by genotyping chip platform. SNPs with high missingness (>5-20%) or low call rates per individual (>2-5%) were filtered *within each chip group* before merging. The problematic OmniExpress chip data was excluded entirely from the merge.
3.  **Merged Cohort QC:**
    *   Individuals without phenotype data were removed.
    *   Heterozygous haploid genotypes were set to missing.
    *   Sex check was performed, and individuals with ambiguous genetic sex (F-statistic between 0.2-0.8 or NaN for non-OmniExpress chips) were removed. This resulted in a cohort of 1433 individuals.
4.  **Initial PCA-based Outlier Removal:** Principal Component Analysis (PCA) was performed on an LD-pruned subset of the 1433 individuals. Individuals falling beyond +/- 6 standard deviations on PC1 or PC2 were removed as ancestral outliers, leading to a cohort of 1412 individuals.
5.  **Relatedness Filtering:** Identity-by-Descent (IBD) was estimated. A greedy algorithm removed 216 individuals to ensure no pairs remained with PI_HAT > 0.1875, resulting in an analysis cohort of 1196 unrelated individuals.
6.  **SNP QC:** SNPs were filtered based on Minor Allele Frequency (MAF < 0.01) and deviation from Hardy-Weinberg Equilibrium (HWE p < 1e-6) within this 1196-person cohort. This resulted in a final set of 949,579 SNPs for association testing.

**Addressing Population Stratification for GWAS:**
Initial GWAS attempts on the 1196-person cohort, using a linear regression model with the 4-category eye color score as the outcome and correcting with an increasing number of Principal Components (PCs, up to 10 derived from this set) plus sex, revealed significant genomic inflation (initial Lambda ≈ 9.7, improving to Lambda ≈ 2.8 with 10 PCs, but with a persistently poor QQ plot shape). This indicated that substantial population stratification, strongly correlated with the eye color phenotype, was not adequately controlled.

To mitigate this, a more stringent approach was adopted:
1.  **Definition of a "Super-Core" Subgroup:** Based on the PCA of the 1196 individuals (using `new_pca_results_1196ind_30pcs.eigenvec`), individuals falling within +/- 2.0 standard deviations of the mean for both PC1 and PC2 were selected. This defined a more ancestrally homogeneous "super-core" subgroup of 1124 individuals.
2.  **Re-QC of Super-Core:** This super-core subgroup underwent MAF (>0.01) and HWE (p>1e-6) filtering.
3.  **New PCA for Super-Core:** PCs were re-calculated *specifically for this super-core subgroup* using an LD-pruned subset of its SNPs. A scree plot of these new PCs (up to 15 calculated) indicated that the first ~6 PCs captured the most significant remaining variance (see Figure X - your super-core scree plot).
4.  **Final GWAS Model on Super-Core:** A linear regression was performed on the 1124 super-core individuals (effective N=1122 after excluding 2 individuals with undefined sex for the sex covariate). The model included the top **6 PCs derived from the super-core subgroup** and sex as covariates.

**Results from GWAS on the "Super-Core" Subgroup (N=1122, 6 PCs + Sex):**

*   **Genomic Inflation Factor (Lambda):** The lambda value for this analysis was **1.859** (see QQ Plot, Figure Y).
    *   *Interpretation:* While this is a substantial improvement from the initial lambda of ~9.7 and the lambda of ~2.8 on the larger "core" group, a lambda of 1.859 still indicates considerable residual genomic inflation. This suggests that even within this more stringently defined super-core, population stratification highly correlated with eye color remains a significant challenge to fully control with the current approach. P-values across the genome are likely to be anti-conservative (i.e., too small).

*   **QQ Plot (Figure Y - your super-core QQ plot):**
    *   The QQ plot shows a significant early and sustained deviation of observed p-values from the expected null distribution (the red diagonal line). The "shelf" or "step-ladder" pattern, though perhaps slightly less pronounced than in initial analyses, is still evident.
    *   This visual evidence corroborates the high lambda value and confirms that systematic inflation of test statistics persists. The tail of the distribution, representing the most significant SNPs, deviates sharply upwards, as expected if true strong signals are present, but the body of the distribution is also inflated.

*   **Manhattan Plot (Figure Z - your super-core Manhattan plot):**
    *   The Manhattan plot displays association results across all autosomes (and numerically coded X, Y, MT if included).
    *   A single, overwhelmingly significant peak of association is observed on **Chromosome 15**. The top SNPs in this region achieve extremely high -log10(P) values (e.g., >150-200). This peak corresponds to the well-known HERC2/OCA2 locus, a major genetic determinant of human eye color, particularly blue versus brown eyes.
    *   The "background" or "lawn" of SNPs across other chromosomes appears elevated, consistent with the observed genomic inflation (lambda = 1.859). Many SNPs may cross suggestive or even genome-wide significance thresholds due to this inflation rather than true association.
    *   No other distinct, major peaks comparable to the Chromosome 15 signal are clearly evident above the inflated background across the rest of the genome from this analysis.

**Top Locus (Chromosome 15 - HERC2/OCA2):**
The strongest signals of association were found in the HERC2/OCA2 region on Chromosome 15. For example, the SNP [Insert rsID of your top SNP here, e.g., rs12913832 if it's a top hit] showed a P-value of [Insert P-value] and a BETA coefficient of [Insert BETA].
    *   *Interpretation:* This finding strongly replicates the most significant known genetic association with human eye color. The HERC2 enhancer region regulates OCA2 gene expression, which is critical for melanin production in the iris. Variants in this region are major determinants of the blue/brown eye color spectrum. While the p-value from this analysis is likely inflated, the identification of this locus is a robust biological finding.

**Discussion of Results and Limitations:**

This GWAS successfully identified the primary genetic locus (HERC2/OCA2 on Chromosome 15) known to be associated with eye color. However, the analysis was significantly challenged by population stratification inherent in the diverse openSNP cohort.

*   **Persistent Genomic Inflation:** Despite rigorous QC steps, including multiple iterations of PCA-based outlier removal, subsetting to a more ancestrally homogeneous "super-core" subgroup (N=1122), and using up to 6 PCs (derived from this super-core) plus sex as covariates, a substantial genomic inflation factor (Lambda = 1.859) remained. This indicates that the p-values generated are anti-conservative, and associations outside of extremely well-established loci (like HERC2/OCA2) must be interpreted with extreme caution, as they are likely false positives.

*   **Challenges of Stratification Control:** The strong correlation between genetic ancestry and eye color makes this trait particularly susceptible to confounding by population structure. The "super-core" definition, while aiming for homogeneity, may still contain residual fine-scale structure that the PCs do not fully capture, or the relationship between this structure and eye color may not be perfectly linear. The observation that distinct "arms" in the initial PCA plots (largely excluded from the super-core) predominantly comprised individuals with brown/hazel eyes highlights the strength of this phenotype-ancestry correlation.

*   **Model Choice (Linear Regression):** A linear regression model was used, treating the 4-category eye color score as quantitative. While a simplification, this model was sufficient to detect the major HERC2/OCA2 signal. However, this model assumes equal intervals between eye color categories, which may not be biologically accurate and could reduce power or precision for other loci. An ordinal logistic regression model might be statistically more appropriate but was not pursued due to the primary challenge of stratification. It is unlikely that the choice of linear model is the main driver of the observed genome-wide inflation.

*   **Interpretation of Findings:** Given the residual inflation, the primary robust finding from this GWAS is the strong association at the HERC2/OCA2 locus. Other SNPs appearing "significant" on the Manhattan plot are highly likely to be artifacts of the remaining stratification and should not be considered novel findings without further validation in a perfectly controlled study.

**Conclusion of GWAS Component:**
This GWAS project successfully demonstrated the key steps of data QC, association testing, and results interpretation. It highlighted the critical importance and significant challenge of controlling for population stratification in genetically diverse cohorts, especially for traits like eye color that are strongly differentiated across ancestries. While the major known locus for eye color (HERC2/OCA2) was robustly identified, persistent genomic inflation indicates that the current analysis is not sufficiently controlled for confident discovery of novel, weaker associations. Further methodological refinements, such as applying linear mixed models (LMMs) or even more stringent cohort definition based on specific ancestral populations (if such data were available), would be necessary to achieve optimal control and explore the genetic architecture of eye color more comprehensively in this dataset.

The experience underscores the iterative nature of GWAS and the necessity of careful diagnostic checks (like QQ plots and lambda calculation) to ensure the validity of findings.

# Additional Analysis D.1: Phenotype Distribution by Genotype at Top SNP

*Purpose*
To visually examine the relationship between the genotypes of the most significant SNP (identified in the "super-core" GWAS) and the distribution of the 4-category eye color phenotype. This provides a direct look at the SNP's effect and can offer insights into the genotype-phenotype relationship.

*Methodology*

**Step 1: Identify the Most Significant SNP**
From your best GWAS results on the "super-core" subgroup (e.g., from `gwas_eyecolor_CORE_SUBGROUP_6pcs_sex.assoc.linear` which gave Lambda = 1.859), identify the SNP with the smallest p-value. This is likely on Chromosome 15 in the HERC2/OCA2 region.
*   *Action: Note down the rsID of this SNP (e.g., `rs12913832` is a common candidate).* Let's call it `TOP_GWAS_SNP_ID`.

Using R, we can filter by p value. 
	
CHR SNP         BP          A1  TEST    NMISS   BETA    STAT    P
15  rs12913832  28365618    A   ADD     1121    1.3650  38.67   5.261e-208


**Step 2: Extract Genotypes for the Top SNP and Phenotypes (PLINK + R)**

You need a file containing the FID, IID, phenotype score, and the genotype (0, 1, or 2 copies of the minor/effect allele) for `TOP_GWAS_SNP_ID` for your super-core individuals.

1.  **Create a list file for PLINK containing just the top SNP:**
 
    # In your terminal
    # Replace rsXXXXXX with your actual TOP_GWAS_SNP_ID
    echo "rs12913832" > top_snp_for_recode.txt 
  

2.  **Use PLINK to recode genotypes for this SNP into an additive format (0, 1, 2) and output to a text file:**
    The input bfile should be your final QCed super-core dataset.
   
    # Replace ${N_super_core} with the actual number (e.g., 897)
    # Replace rsXXXXXX with your actual TOP_GWAS_SNP_ID
    GWAS_INPUT_BFILE_BASENAME="gwas_core_1124"
    TOP_SNP_ID_VAR="rs12913832" 

plink --bfile gwas_core_1124 \
      --snp rs12913832 \
      --recode A \
      --out top_snp_genotypes_for_supercore_analysis
    ```
    *   `--snp rs12913832`: Tells PLINK to focus on this specific SNP.
    *   `--recode A`: Outputs a `.raw` file where genotypes are coded as counts of one allele (usually the minor allele as A1 in your .bim file). This file will contain FID, IID, PaternalID, MaternalID, Sex, Phenotype (from .fam, can be ignored), and then a column named `rsXXXXXX_A` (or similar, where A is the counted allele). The values will be 0, 1, 2, or NA (missing).
    *   This creates `top_snp_genotypes_supercore.raw`.

3.  **Prepare Data in R:**
    Load the `.raw` file and merge it with your phenotype file for the super-core individuals.


**Step 3: Interpretation**

Biological Background for rs12913832 (HERC2/OCA2 region):
rs12913832 G allele: This is the ancestral allele and is strongly associated with brown eyes. It's linked to higher melanin production.
rs12913832 A allele: This is the derived allele and is strongly associated with blue eyes (and generally lighter eye/hair/skin pigmentation). It's linked to reduced melanin production due to its effect on OCA2 expression.

*   Examine the bar plot (`pheno_dist_plot`). Does it show a clear trend? For a typical HERC2/OCA2 SNP where the minor allele is associated with lighter eyes (e.g., 'A' allele in rs12913832 is linked to blue eyes):
 *  Genotype 0 (Count of 'A' = 0): This means the genotype is GG. These individuals should predominantly have darker eyes (Brown/Hazel).
 *  Genotype 1 (Count of 'A' = 1): This means the genotype is GA. These individuals should have intermediate eye colors or a mix, often Hazel/Green, some Brown, some Blue/Green.    The effect is somewhat co-dominant/additive.
 *  Genotype 2 (Count of 'A' = 2): This means the genotype is AA. These individuals should predominantly have lighter eyes (Blue, Blue/Green, potentially Green).
*   The boxplot (`pheno_box_plot`) will show if the mean/median of your 0-3 eye color score changes across the genotype groups.
*   **In your report:** Present one of these plots. Describe the observed pattern and how it aligns (or doesn't) with the SNP's known role or the GWAS finding. Discuss if the effect appears additive.

This analysis is relatively quick, directly uses your data, provides visual insight into your top hit, and fulfills one of the Section D requirements effectively.

If you want to do the conditional analysis instead/additionally, let me know, and I can detail those steps. But the phenotype distribution plot is a very good, straightforward choice.

**Esther comments**

When looking at the .raw data in R, we see the opposite distribution compared to the known biology of HERC2/OCA2. After chekcing if our data was correctly being read and genotypes counted correctly, which they were, it seems that the problem lays elsewhere. A common cause for this is what is called a Strand Flip. 