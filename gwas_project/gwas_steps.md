
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
| 0     | blue       | `blue`, `dark_blue` , `blue-grey`                                                                                                      |
| 1     | blue_green | `blue-green`, `blue-green-grey`, `blue-green-gold`, `blue-green_heterochromia`, `blue/green/grey_-_changes_with_lighting_and_clothing` |
| 2     | green      | `green`, `green-gray`                                                                                                                  |
| 3     | hazel      | `hazel/brown-green`, `hazel_green`, `green_with_brown_freckles`                                                                        |
| 4     | brown      | `brown`, `amber-brown`                                                                                                                 |
| 5     | dark_brown | `dark_brown`, `black`, `red/blood`                                                                                                     |


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
        We are essentially penalizing individuals for not having SNPs that were never on their chip — that’s unfair and unnecessary, and so trying to seperate them by chip would 
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

    We might also consider doing one more group, which includs the 399 entries that did not report a chip, and do a QC on that group also. This might tell us something. 

  # 5.4: Chip specific genotype data

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

  # 5.5: Prune the SNP set per chip, sp that each .bed/.bim/.fam triplet only includes the SNPs actually genotyped on that chip.
    We do this so that each .bed/.bim/.fam triplet only includes SNPs that were actually genotyped on that specific chip.
    If we don’t do this, Plink will calculate missingness using all SNPs in the merged dataset, including those not present on a given chip, 
    which artificially inflates the missingness per individual. 

    For each chip, filter out SNPs with >5% missing data (per SNP), and save only the remaining SNPs.
    This creates chip_HTS_iSelect_HD_snps_filtered.snplist, a list of good-quality SNPs genotyped on this chip:

      plink --bfile chip_HTS_iSelect_HD \
      --geno 0.05 \
      --write-snplist \
      --out chip_HTS_iSelect_HD_snps_filtered


    Use the list above to make a new, chip-specific cleaned dataset.
    This ensures the new dataset only includes variants truly present on this chip, so downstream missingness estimates are accurate:

      plink --bfile chip_HTS_iSelect_HD \
      --extract chip_HTS_iSelect_HD_snps_filtered.snplist \
      --make-bed \
      --out chip_HTS_iSelect_HD_pruned


    Run the following command to inspect the number of SNPs retained after pruning:

      wc -l chip_HTS_iSelect_HD_pruned.bim


    Calculate the missingness per individual, before then calculating the overall average missingess: 
    
      plink --bfile chip_HTS_iSelect_HD_pruned \
      --missing \
      --out chip_HTS_iSelect_HD_pruned_missing

      awk '{if(NR>1) print $6}' chip_HTS_iSelect_HD_pruned_missing.imiss | \
      awk '{sum+=$1} END {print "Average missingness:", sum/NR}'

  # 5.5.1: unknown chip group prune settings:
   Since this group did not have a reported chip type, we expect greater heterogeneity. These individuals most likely used different chips or came from different companies.  
   As a result, we also expect higher missingness due to low SNP overlap across individuals. Therefore, we need to **adjust our pruning strategy** to account for this.  
    
    We relax the SNP-level missingness threshold from 5% to 20% (`--geno 0.2`) and test whether this retains a reasonable number of variants.

      plink --bfile chip_unknown \
      --geno 0.2 \
      --write-snplist \
      --out chip_unknown_snps_filtered

    Then, we check for SNPs retained, and if they are a reasonable amount, we continue with the PLINK command:

      wc -l chip_unknown_snps_filtered.snplist

      plink --bfile chip_unknown \
      --extract chip_unknown_snps_filtered.snplist \
      --make-bed \
      --out chip_unknown_pruned


    Then, we run the missingness calculation as well as the average missingess calc:

      plink --bfile chip_unknown_pruned \
      --missing \
      --out chip_unknown_pruned_missing

      awk '{if(NR>1) print $6}' chip_unknown_pruned_missing.imiss | \
      awk '{sum+=$1} END {print "Average missingness:", sum/NR}'


  # 5.6: Report average missingness for all chip groups.
    To assess the quality of each chip-specific group, we calculated the average missingness per individual after pruning.
    Run a command that loops through and reports all the average missingness:

    for file in chip_*_pruned_missing.imiss; do
      chip=$(echo $file | sed 's/_pruned_missing.imiss//')
      avg=$(awk '{if(NR>1) print $6}' "$file" | awk '{sum+=$1} END {if(NR>0) print sum/NR; else print "NA"}')
      printf "%-35s  Average missingness: %.4f\n" "$chip" "$avg"
    done

    Our output: 
      chip_HTS_iSelect_HD                  Average missingness: 0.0080
      chip_Illumina_GSAs                   Average missingness: 0.0242
      chip_OmniExpress_plus                Average missingness: 0.0155
      chip_OmniExpress                     Average missingness: 0.0314
      chip_unknown                         Average missingness: 0.1516

    Then, looking into the normal missiness thresholds, which are:

      | Threshold  | Typical Meaning                                                                                                 |
      | ---------- | --------------------------------------------------------------------------------------------------------------- |
      | **< 2–5%** | Excellent — gold standard for most GWAS datasets                                                                |
      | **5–10%**  | Acceptable for exploratory analysis or if justified                                                             |
      | **10–20%** | **Tolerable in special cases** — like our unknown chip — with proper documentation and possibly extra filtering |
      | **>20%**   | Risky — usually excluded unless extremely valuable samples                                                      |

      Our results indicate that all named chip groups meet or exceed best practices, while the unknown group shows elevated missingness (~15%).
      This is expected due to platform heterogeneity, and we currently retain this group for analysis with caution, noting it may require separate treatment or exclusion in stricter analyses.


  # 5.7: Heterozygosity and Missingness Outlier Detection
    This is a multistep process, where we aim to detect sample contamination or failed genotyping. We also want to either flag or remove outliers. 
  
  # 5.7.1: Generate Missingness and Heterozygosity Stats (HTS Group)
    Now we will use two different PLINK commands, one for finding the missingness per individual and another to find hte heterozygosity per individual. 
    The two commands will create .imiss and .het files. 

    Individual missingness:

      plink --bfile chip_HTS_iSelect_HD_pruned \
      --missing \
      --out chip_HTS_iSelect_HD_qc

    Individual Heterozygosity:

      plink --bfile chip_HTS_iSelect_HD_pruned \
      --het \
      --out chip_HTS_iSelect_HD_qc

  # Step 5.7.2: Outlier detection using heterozygosity and missingness
    The goal here is to identify very high or very low heterozygosity, whcih could be due to contamination or inbreeding,
    as well as look for high missingness individuals (poor genotyping). I did this manually for all chips, but in hindsight, a workflow would have been a smart solution to automate this process. 

    The plan is to do this in R. Another possibility is using jupyter notebooks, however I feel like the files are an alright size to work on on my local machine.

    Summary of what we did in R and why we did it: 

    | Step                     | What it does                            | Why it matters                                               |
    | ------------------------ | --------------------------------------- | ------------------------------------------------------------ |
    | Load `.imiss` and `.het` | Get missingness and heterozygosity info | These are two independent indicators of sample quality       |
    | Calculate `HET_RATE`     | % heterozygous genotypes per individual | Outliers here could mean contamination or other issues       |
    | Merge data               | Align individuals across metrics        | Allows joint analysis                                        |
    | Plot                     | Visual sanity check                     | Spot trends, tight clustering, or spread                     |
    | ±3 SD rule               | Statistical outlier filter              | Removes samples that deviate significantly from normal range |
    | Save to file             | Prepares for filtering in Plink         | Enables reproducible sample removal                          |

  # 5.7.3: Remove Ouliers using PLINK.
    Here, we are simply using the --remove command from PLINK to remove the outliers from the chip files. 

      plink --bfile chip_HTS_iSelect_HD_pruned \
      --remove chip_HTS_iSelect_HD_outliers.txt \
      --make-bed \
      --out chip_HTS_iSelect_HD_pruned_clean

      Command explinations: 
      --bfile: Specifies the base file (without extensions) for the chip-specific dataset.
      --remove: Tells Plink to exclude individuals listed in the provided outlier .txt file.
      --make-bed: Plink will output a new cleaned .bed/.bim/.fam dataset, which excludes the outliers.
      --out: Specifies the name of the new output file for the cleaned dataset.

# 6: More data analysis and thoughts
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

# 7:  Population Structure and Confounder Checks
  We will focus on the Principal Component Analysis (PCA), Relatedness Check, and Sex Check before moving on to the GWAS.

  # 7.1: Perform Principal Component Analysis (PCA)
    PCA helps identify population stratification by calculating the first few principal components (PCs). If our dataset contains samples from different populations or ancestry groups, the PCs will help capture that variance.

    We can simply use PLINK for the PCA:

      plink --bfile chip_HTS_iSelect_HD_pruned_clean --pca 10 --mind 0.05 --out chip_HTS_iSelect_HD_pca

        Explanation of the command:
          --bfile chip_HTS_iSelect_HD_pruned_clean: This specifies the cleaned dataset for the HTS iSelect HD chip group.

          --pca 10: This tells Plink to calculate the first 10 principal components.

          --mind 0.05: removes individuals with over 5% missingness

          --out chip_HTS_iSelect_HD_pca: Specifies the output files. It will save the PCA results in files starting with chip_HTS_iSelect_HD_pca.

        This will generate the following outputs:

          chip_HTS_iSelect_HD_pca.eigenvec: The file with the principal components.

          chip_HTS_iSelect_HD_pca.eigenval: The eigenvalues of the components.

    After doing the PLINK commands, we visualize the PCAs in R. When comparing our PCA plots, which are chip specific, to the PCA plots made by Mark Simcoe et al., we see similar V shapes as in the 23AndMe European Cohorts, which were determined to have greater than 97% european ancestry. Since we did this mainly to look for population stratification, and we see the same V shapes as the eauropean cohorts, we feel like it is safe to assume that ours are also mainly of european anceestry. 

      However !! It is safer to potentially do a PCA with referance populations to see if our data mostly falls within specific population groups. Another possibility is to run Admixture analysis. This we can either do with or without a referance data set. 







