# Load necessary libraries
library(dplyr)

## ChEEP Metadata 
cheep_meta <- read.csv("Data/phy_meta_2024-04-07.csv")

table(cheep_meta$location)
market_ids <- cheep_meta$ID[cheep_meta$location=="Market"]
market_ids

# Define the directory containing the TSV files
#directory <- "Data/Comps/81"
directory <- "/Users/nmolson/Desktop/81"

# Get a list of all TSV files in the directory
file_list <- list.files(path = directory, pattern = "*summary.tsv", full.names = TRUE)
file_list
# Read all TSV files into a list of dataframes
df_list <- lapply(file_list, read.delim)

# Combine all dataframes into a single dataframe
combined_df <- do.call(rbind, df_list)

## Label possible combinations of samples 
combined_df <- combined_df %>%
  mutate(comp = case_when(
    (grepl("^A|^EUP", sample1)& grepl("^SRR", sample2)) |
      (grepl("^A|^EUP", sample2)  & grepl("^SRR", sample1)) ~ "ChEEP x MapSan",
    
    (grepl("^A|^EUP", sample1) & grepl("^[0-9]+$", sample2)) |
      (grepl("^A|^EUP", sample2) & grepl("^[0-9]+$", sample1)) ~ "ChEEP x FOCAL",
    
    TRUE ~ NA_character_
  ))

## Drop all pairs that are not chicken x child 
combined_df <- combined_df %>%
  filter(!is.na(comp))

# Summary statistics
summary_stats <- combined_df %>%
  summarise(across(everything(), list(mean = ~mean(.), sd = ~sd(.), median = ~median(.))))
print(summary_stats)
# Split the data by 'comp'
focal_df <- combined_df %>% filter(comp == "ChEEP x FOCAL")
mapsan_df <- combined_df %>% filter(comp == "ChEEP x MapSan")

# Define custom colors for the 'comp' variable
custom_colors <- c("ChEEP x FOCAL" = "#D81B60", "ChEEP x MapSan" = "#FFC107")

# Zoomed-in plot (Limited to threshold values)
zoom_ST540 <- ggplot() +
  # Plot 'ChEEP x MapSan' first so it appears behind 'ChEEP x FOCAL'
  geom_point(data = mapsan_df, aes(x = gapJaccardSim, y = singleAgreePct, size = commonPct, color = comp)) +
  # Plot 'ChEEP x FOCAL' on top
  geom_point(data = focal_df, aes(x = gapJaccardSim, y = singleAgreePct, size = commonPct, color = comp)) +
  scale_x_continuous(limits = c(0.970, 1.0), name = " ") +
  scale_y_continuous(limits = c(99.95, 100), name = " ") +
  labs(size = "Common\nCallable [%]") +
  theme_minimal() +
  scale_color_manual(values = custom_colors) +
  theme(legend.position = "right",
        legend.title= element_blank(),
        text = element_text(size = 24), 
        axis.text.x = element_text(angle = 40, hjust = 1)) +
  guides(size = guide_legend(title.position = "top", title.hjust = 0.5))

#Full plot (not just threshold)
full_ST540 <- ggplot() +
  # Plot 'ChEEP x MapSan' first so it appears behind 'ChEEP x FOCAL'
  geom_point(data = mapsan_df, aes(x = gapJaccardSim, y = singleAgreePct, size = commonPct, color = comp)) +
  # Plot 'ChEEP x FOCAL' on top
  geom_point(data = focal_df, aes(x = gapJaccardSim, y = singleAgreePct, size = commonPct, color = comp)) +
  # scale_x_continuous(limits = c(0.970, 1.0), name = " ") +
  # scale_y_continuous(limits = c(99.95, 100), name = " ") +
  labs(size = "Common\nCallable [%]") +
  xlab("Gap Similarity") + ylab("ACNI") + 
  theme_minimal() +
  scale_color_manual(values = custom_colors) +
  theme(legend.position = "none",
        text = element_text(size = 24), axis.text.x = element_text(angle = 40, hjust = 1)  # Rotate x-axis text by 40 degrees
  ) +
  guides(size = guide_legend(title.position = "top", title.hjust = 0.5))


#fig_ST540 <- ggarrange(full_ST540,zoom_ST540,  cols=2)

#save the plots as image files
#ggsave("figures/ST540_fig.png", plot=fig_ST540, width =30, height = 30, units = "cm")

write.csv(combined_df, "ST216_2025-01-26.csv")
######################################################
### Only Above Strain-sharing threshold ##############
######################################################

fig3 <- combined_df %>% 
  filter((gapJaccardSim>=0.97) & (singleAgreePct>=99.5))

fig3$ID1 <- paste0(fig3$sample1, "_", fig3$sample2)
obs_pairs <- length(unique(fig3$ID1))

fig3$ID2 <- paste0(fig3$sample2, "_", fig3$sample1)

write.csv(fig3, "ST540_2024-09-18.csv")