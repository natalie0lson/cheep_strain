# Load necessary libraries
library(dplyr)

# Define the directory containing the TSV files
directory <- "13C1079T_compare"

# Get a list of all TSV files in the directory
file_list <- list.files(path = directory, pattern = "*summary.tsv", full.names = TRUE)

# Read all TSV files into a list of dataframes
df_list <- lapply(file_list, read.delim)

# Combine all dataframes into a single dataframe
combined_df <- do.call(rbind, df_list)

# View the combined dataframe
print(combined_df)

# Perform some basic analysis
# Example: Summary statistics
summary_stats <- combined_df %>%
  summarise(across(everything(), list(mean = ~mean(.), sd = ~sd(.), median = ~median(.))))

# Print summary statistics
print(summary_stats)
# Load necessary libraries
library(ggplot2)

# Assuming 'compare_df' is your dataframe and contains the columns: gapJaccardSim, singleAgreePct, commonPct
# Adjust the plot aesthetics to match the seaborn plot

ggplot(combined_df, aes(x = gapJaccardSim, y = singleAgreePct, size = commonPct)) +
  geom_point() +
  scale_x_continuous(limits = c(0.970, 1.0), name = "Gap similarity") +
  scale_y_continuous(limits = c(99.9, 100), name = "ACNI", labels = scales::percent_format(accuracy = 1)) +
  labs(size = "Common\nCallable [%]") +
  theme_minimal() +
  theme(legend.position = "right") +
  guides(size = guide_legend(title.position = "top", title.hjust = 0.5))
