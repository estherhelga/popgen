# Load required libraries
library(data.table)
library(dplyr)

# Load metadata and fam file
meta <- fread("metadata.txt", header = TRUE)
# Convert to character to avoid type mismatch in join
fam <- fread("gwas_data.fam", header = FALSE)
colnames(fam)[1:2] <- c("FID", "IID")
fam$FID <- as.character(fam$FID)
fam$IID <- as.character(fam$IID)


# Merge metadata and fam on user ID (assuming user == IID)
merged <- inner_join(fam, meta, by = c("IID" = "user"))

# Count how many per chip (check)
table(merged$chip)

# Define the chip groups you want
chips <- unique(na.omit(merged$chip))

# Write a file for each chip group
for (chip_type in chips) {
  subset_chip <- merged %>% filter(chip == chip_type) %>% select(FID, IID)
  chip_filename <- gsub(" ", "_", chip_type)  # Make filename safe
  fwrite(subset_chip, paste0("chip_group_", chip_filename, ".txt"), col.names = FALSE, sep = "\t")
}
