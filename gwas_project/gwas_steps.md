
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

    | Scale | Label (Modified)       | Categories from your data                              |
    |-------|------------------------|---------------------------------------------------------|
    | 0     | Blue                   | `blue`, `dark_blue`                                    |
    | 1     | Blue-green             | `blue-green`, `blue-green-grey`, `blue-green-gold`, `blue-green_heterochromia` |
    | 2     | Green                  | `green`, `green-gray`,  `green_with_brown_freckles` |
    | 3     | Hazel (light mix)      | `hazel/brown-green`, `hazel_green`                                   |
    | 4     | Brown (standard)       | `brown`, `amber-brown`                                 |
    | 5     | Dark brown             | `dark_brown`, `black`, `red/blood`                     |  

  # How to create a proper phenotype file for Plink.