# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(epitools)
library(MASS)
library(tidyr)

# Define the directory containing the TSV files
directory <- "Data/ATL_CHEEP_Strains/13C1079T"

# Get a list of all TSV files in the directory
file_list <- list.files(path = directory, pattern = "*summary.tsv", full.names = TRUE)

# Read all TSV files into a list of dataframes
df_list <- lapply(file_list, read.delim)

# Combine all dataframes into a single dataframe
combined_df <- do.call(rbind, df_list)

# Summary statistics
summary_stats <- combined_df %>%
  summarise(across(everything(), list(mean = ~mean(.), sd = ~sd(.), median = ~median(.))))

# Print summary statistics
print(summary_stats)
# Load necessary libraries
library(ggplot2)

# Zoomed-in plot (Limited to threshold values)
zoom_13C1079T <- ggplot() +
  # Plot 'ChEEP x MapSan' first so it appears behind 'ChEEP x FOCAL'
  geom_point(data = combined_df, aes(x = gapJaccardSim, y = singleAgreePct, size = commonPct)) +
  scale_x_continuous(limits = c(0.970, 1.0), name = " ") +
  scale_y_continuous(limits = c(99.95, 100), name = " ") +
  labs(size = "Common\nCallable [%]") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.title= element_blank(),
        text = element_text(size = 24), 
        axis.text.x = element_text(angle = 40, hjust = 1)) +
  guides(size = guide_legend(title.position = "top", title.hjust = 0.5))

zoom_13C1079T

#Full plot (not just threshold)
full_13C1079T <- ggplot() +
 geom_point(data = combined_df, aes(x = gapJaccardSim, y = singleAgreePct, size = commonPct)) +
 # scale_x_continuous(limits = c(0.970, 1.0), name = " ") +
  # scale_y_continuous(limits = c(99.95, 100), name = " ") +
  labs(size = "Common\nCallable [%]") +
  xlab("Gap Similarity") + ylab("ACNI") + 
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 24), axis.text.x = element_text(angle = 40, hjust = 1)  # Rotate x-axis text by 40 degrees
  ) +
  guides(size = guide_legend(title.position = "top", title.hjust = 0.5))

full_13C1079T


fig_13C1079T <- ggarrange(full_13C1079T,zoom_13C1079T,  cols=2)

#save the plots as image files
ggsave("figures/atl_13C1079T_fig.png", plot=fig_13C1079T, width =30, height = 30, units = "cm")

write.csv(combined_df, "atl_13C1079T.csv")