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


# GWAS 4 category Quantitative Trait Linear Model visualization

