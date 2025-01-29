# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(epitools)
library(MASS)
library(tidyr)

## ChEEP Metadata 
cheep_meta <- read.csv("Data/phy_meta_2024-04-07.csv")

market_ids <- cheep_meta$ID[cheep_meta$location=="Market"]

cheep_ids <- cheep_meta$ID

fecal_ids <- cheep_meta$ID[cheep_meta$Feces==1]

carcass_ids <- cheep_meta$ID[cheep_meta$Carcass==1]

water_ids <- cheep_meta$ID[cheep_meta$Water==1]

com_ids <- cheep_meta$ID[cheep_meta$breed=="Commercial"]

ind_ids <- cheep_meta$ID[cheep_meta$breed=="Indigenous"]

# Define the directory containing the TSV files
#directory <- "Data/Comps/13C1079T"
directory <- "/Users/nmolson/Desktop/13C1079T"

# Get a list of all TSV files in the directory
file_list <- list.files(path = directory, pattern = "*summary.tsv", full.names = TRUE)
file_list

# Read all TSV files into a list of dataframes
#df_list <- lapply(file_list, read.delim)

# 1. List files safely
file_list <- list.files(
  path = directory,
  pattern = "*summary.tsv",
  full.names = TRUE
)

# 2. Remove empty files
file_sizes <- file.info(file_list)$size
file_list <- file_list[file_sizes > 0]

# 3. Read files with error handling
df_list <- lapply(file_list, function(file) {
  tryCatch(
    read.delim(file),
    error = function(e) {
      message(sprintf("Skipping %s: %s", file, e$message))
      return(NULL)
    }
  )
})


# Combine all dataframes into a single dataframe
combined_df <- do.call(rbind, df_list)
head(combined_df)

combined_df <- combined_df %>%
  mutate(comp = case_when(
    (grepl("^A|^EUP", sample1)  & grepl("^SRR", sample2)) |
      (grepl("^A|^EUP", sample2)  & grepl("^SRR", sample1)) ~ "ChEEP x MapSan",
    
    (grepl("^A|^EUP", sample1) & grepl("^[0-9]+$", sample2)) |
      (grepl("^A|^EUP", sample2) &  grepl("^[0-9]+$", sample1)) ~ "ChEEP x FOCAL",

    TRUE ~ NA_character_
  ))

table(combined_df$comp)

combined_df <- combined_df %>%
  filter(!is.na(comp))

# Summary statistics
summary_stats <- combined_df %>%
  summarise(across(everything(), list(mean = ~mean(.), sd = ~sd(.), median = ~median(.))))

# Print summary statistics
print(summary_stats)
# Load necessary libraries
library(ggplot2)

focal_df <- combined_df %>% filter(comp == "ChEEP x FOCAL")
mapsan_df <- combined_df %>% filter(comp == "ChEEP x MapSan")
# Define custom colors for the 'comp' variable
custom_colors <- c("ChEEP x FOCAL" = "#af8dc3", "ChEEP x MapSan" = "#7fbf7b")

# Zoomed-in plot (Limited to threshold values)
zoom_13C1079T <- ggplot() +
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

zoom_13C1079T

#Full plot (not just threshold)
full_13C1079T <- ggplot() +
  # Plot 'ChEEP x MapSan' first so it appears behind 'ChEEP x FOCAL'
  geom_point(data = mapsan_df, aes(x = gapJaccardSim, y = singleAgreePct, size = commonPct, color = comp)) +
  # Plot 'ChEEP x FOCAL' on top
  geom_point(data = focal_df, aes(x = gapJaccardSim, y = singleAgreePct, size = commonPct, color = comp), alpha = 0.5)  +
  geom_vline(xintercept = 0.97, linetype = "dotted", color = "goldenrod1", size=2) +
  geom_hline(yintercept = 99.95, linetype = "dotted", color = "goldenrod1", size=2) +

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

full_13C1079T

# Arrange the plots and add labels "A" and "B"
fig_13C1079T <- ggarrange(
  full_13C1079T, 
  zoom_13C1079T, 
  labels = c("A", "B"), # Add panel labels
  ncol = 2
)

# Add a title to the combined plot
fig_13C1079T <- annotate_figure(
  fig_13C1079T,
  top = text_grob("E. coli Strain 13C1079T", face = "bold", 
                  size = 24,  # Increase font size
                  hjust = 2   # Left-justify the title
                  )
)

fig_13C1079T
#save the plots as image files
ggsave("figures/13C1079T_fig.png", plot=fig_13C1079T, width =40, height = 35, units = "cm")


write.csv(combined_df, "Data/13C1079T_2025-01-26.csv")