# Load necessary libraries
library(dplyr)
library(purrr)
library(readr)
library(stringr)

# Define the directory containing the TSV files
straingr_dir <- "cheep.straingr.tsv"

# List all TSV files in the directory
tsv_files <- list.files(path = straingr_dir, pattern = "\\.tsv$", full.names = TRUE)

# Initialize lists to store dataframes and sample names
df_list <- list()
sample_names <- c()

# Read each TSV file, remove 'TOTAL' row, and store in the list
for (f in tsv_files) {
  df <- read_tsv(f, col_names = TRUE) %>%
    filter(!str_detect(rownames(.), "TOTAL")) # Remove TOTAL statistics
  
  df_list[[f]] <- df
  sample_names <- c(sample_names, tools::file_path_sans_ext(basename(f)))
}

# Combine all dataframes into one, adding sample names as a key
straingr_df <- bind_rows(df_list, .id = "sample")
head(straingr_df$sample)
straingr_df$sample <- gsub("Data/cheep.straingr.tsv/|\\.filter.tsv", "", straingr_df$sample)
head(straingr_df$sample)

# Load the straingst_df (assuming it's already available)
# straingst_df <- read_csv("path_to_straingst_df.csv") # Example, adjust accordingly

# Add straingst_present, is_plasmid, and enough_cov columns
straingr_df <- straingr_df %>%
  mutate(#straingst_present = ref %in% straingst_df$index, # Adjust if `straingst_df` column is different
         is_plasmid = length < 4e6,
         enough_cov = coverage > 0.5)

table(straingr_df$enough_cov)

# Filter and re-index the dataframe
straingr_df <- straingr_df %>%
  filter(!is_plasmid & enough_cov) 
table(straingr_df$is_plasmid)


# Print the resulting dataframe
print(straingr_df)
