
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
    
  # How should you group them into a binary phenotype (e.g., light vs. dark)?

  # In the paper on GWAS eyecolor, in their suplementary materials (file abd1239_sm.pdf), they use a 7 category scale for eye color phenotype classification:

  # 1.1.3 Phenotypic definition
  # Iris colour phenotype was split into 7 distinct categories which were then transformed into an ordinal trait for analysis with the following values: 
    0=blue
    1=greenish-blue
    2=green
    3=light hazel
    4=dark hazel
    5=light brown
    6=dark brown.

  # However, our reported eyecolors are a bit hard to classify the exact same wya, so maybe we need to think about if we want the same scale, or if we want to create a more binary scale like Light and Dark?
  # current plan is to stick as much to the paper as possible, with some simplifications. This could be something like:

| Scale | Label      | Included Categories                                                                                                                    |
| ----- | ---------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| 0     | Blue       | `blue`, `dark_blue` , `blue-grey`                                                                                                      |
| 1     | Blue-green | `blue-green`, `blue-green-grey`, `blue-green-gold`, `blue-green_heterochromia`, `blue/green/grey_-_changes_with_lighting_and_clothing` |
| 2     | Green      | `green`, `green-gray`                                                                                                                  |
| 3     | Hazel      | `hazel/brown-green`, `hazel_green`, `green_with_brown_freckles`                                                                        |
| 4     | Brown      | `brown`, `amber-brown`                                                                                                                 |
| 5     | Dark brown | `dark_brown`, `black`, `red/blood`                                                                                                     |


# 4: QC and PLINK info

  # 4.1: What are the different files needed for PLINK?
    .bed = binary genotype matrix (SNPs × individuals) - Binary-encoded genotypes
    .bim = SNP information (chromosome, position, alleles, etc.) - variant metadata
    .fam = sample information (IDs, sex, phenotype placeholders) - individual metadata

  # 4.2: QC filtering steps for samples
    - Missing genotype rate per individual
    - Extreme heterozygosity
    - Related individuals
  
  # 4.3: QC filtering steps for SNPs
    - SNPs with high missingness
    - SNPs with low minor allele frequency (MAF)
    - Hardy-Weinberg equilibrium (HWE) deviation

# 5: QC pipeline

  # 5.1: Initial sample and SNP filtering
    Goal:
     Remove:
      Individuals with too much missing genotype data
      SNPs with high missingness
      SNPs with low minor allele frequency (MAF)
    Why This Step?
      This removes low-quality individuals and unreliable SNPs that could cause false positives or noise in downstream analyses.

    Command to run on GenomeDK (After activating your gwas env):
      plink --bfile gwas_data \
      --allow-no-sex \
      --mind 0.05 \
      --geno 0.05 \
      --maf 0.01 \
      --make-bed \
      --out qc_step1

    | Option           | What it does                               |
    | ---------------- | ------------------------------------------ |
    | `--bfile`        | Uses your `.bed/.bim/.fam` dataset         |
    | `--allow-no-sex` | Allows processing despite missing sex info |
    | `--mind 0.05`    | Removes individuals with >5% missing data  |
    | `--geno 0.05`    | Removes SNPs with >5% missing data         |
    | `--maf 0.01`     | Removes SNPs with minor allele freq <1%    |
    | `--make-bed`     | Outputs new `.bed/.bim/.fam` files         |
    | `--out`          | Sets output filename prefix                |

    Information from the terminal output from PLINK:
      When running the above command, with the 5% missingness threshold, we see that only 4 individuals make it through the filtering step. 
      In the assignment, we are told that the data comes from different companies that use different chips. This can lead to some SNPs only being present on a subset of individuals. That then can make each person have high missingness overall, when we combine all the parkers across platforms. Another reason could be that id the files from the different companies were not harminiced correctly, we would see a lot of missingess. So far, we are looking into the former explination. 

      We start by looking into how bad the missingess actuall is:

      plink --bfile gwas_data --missing --out raw_missing
      awk '{print $6}' raw_missing.imiss | sort -n | tail
      awk '{if(NR>1) print $6}' raw_missing.imiss | awk '{sum+=$1} END {print "Average missingness:", sum/NR}'

      From these commands, we see that the top 6 highest missingness values are all 1, and that the average missingness is 0.521106
        "Due to the openSNP dataset containing data from different direct-to-consumer genotyping platforms, the average individual had ~52% missing genotypes, likely reflecting chip incompatibility or data merging artifacts."

      Now we run the filtering again, but this time with a 0,5 missingess threshold, rather than 0,05
      plink --bfile gwas_data \
      --allow-no-sex \
      --mind 0.5 \
      --geno 0.05 \
      --maf 0.01 \
      --make-bed \
      --out qc_step1_loose

        This sitll filters out over half of our individuals, so our next plan is to try and seperate the chips, do QC on each chip seperately, and then merge the cleaned data.
        We are essentially penalizing individuals for not having SNPs that were never on their chip — that’s unfair and unnecessary, anmd so trying to seperate them by chip would 
        fix this. 
  
  # 5.2: Metadata exploration

    After exploring the metadata file, We found that: 
      The metadata file contains the following columns:
      user, 	build, 	chip, 	chip_version, 	inferred_sex, 	source

      We have 4 different chips:
      HTS iSelect HD (587 entries),    Illumina GSAs (248 entires),      OmniExpress (291 entries),	 OmniExpress plus (484 entries).

      We have a lot of different companies, some who use the same chips, some who dont, and some who never report their chip. In total, we have 1610 chip reports, and 399 that did not report a chip.

    From this info, we can now group our data by chip, and do QC on each individual chip type. 

  # 5.3: Chip grouping

    First, we split our main dataset into 4 chip-based groups:

      HTS iSelect HD

      Illumina GSAs

      OmniExpress

      OmniExpress Plus

    Using the split_by_chip.R on the cluster, and the command Rscript to run it, we create 4 individual chip txt files, that we can then use alongside PLINK to --keep and --remove during the individual QC steps. The .txt files contain the Family ID (FID), and the Individual ID (IID), which in our case are the same. But PLINK needs both to operate the --keep and the --remove. 