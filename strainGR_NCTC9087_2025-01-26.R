# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
#install.packages("epitools")
library(epitools)
library(MASS)
library(tidyr)

## ChEEP Metadata 
cheep_meta <- read.csv("Data/phy_meta_2024-04-07.csv")

## Reformat IDs 
cheep_ids <- cheep_meta$ID

market_ids <- cheep_meta$ID[cheep_meta$location=="Market"]

fecal_ids <- cheep_meta$ID[cheep_meta$Feces==1]
carcass_ids <- cheep_meta$ID[cheep_meta$Carcass==1]

water_ids <- cheep_meta$ID[cheep_meta$Water==1]

com_ids <- cheep_meta$ID[cheep_meta$breed=="Commercial"]

ind_ids <- cheep_meta$ID[cheep_meta$breed=="Indigenous"]

# Define the directory containing the TSV files
directory <- "Data/Comps/NCTC9087"

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
custom_colors <- c("ChEEP x FOCAL" = "#5ab4ac", "ChEEP x MapSan" = "#f16433")

# Zoomed-in plot (Limited to threshold values)
zoom_NCTC9087 <- ggplot() +
  # Plot 'ChEEP x MapSan' first so it appears behind 'ChEEP x FOCAL'
  geom_point(data = mapsan_df, aes(x = gapJaccardSim, y = singleAgreePct, size = commonPct, color = comp)) +
  # Plot 'ChEEP x FOCAL' on top
  geom_point(data = focal_df, aes(x = gapJaccardSim, y = singleAgreePct, size = commonPct, color = comp), alpha = 0.5) +
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

zoom_NCTC9087

#Full plot (not just threshold)
full_NCTC9087 <- ggplot() +
  # Plot 'ChEEP x MapSan' first so it appears behind 'ChEEP x FOCAL'
  geom_point(data = mapsan_df, aes(x = gapJaccardSim, y = singleAgreePct, size = commonPct, color = comp)) +
  # Plot 'ChEEP x FOCAL' on top
  geom_point(data = focal_df, aes(x = gapJaccardSim, y = singleAgreePct, size = commonPct, color = comp), alpha = 0.5) +
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

full_NCTC9087
library(ggpubr)

# Arrange the plots and add labels "A" and "B"
fig_NCTC9087 <- ggarrange(
  full_NCTC9087, 
  zoom_NCTC9087, 
  labels = c("A", "B"),  # Add panel labels
  ncol = 2
)

# Add a title to the combined plot
fig_NCTC9087 <- annotate_figure(
  fig_NCTC9087,
  top = text_grob("E. coli Strain 13C1079T", face = "bold", 
                  size = 24,  # Increase font size
                  hjust = 2   # Left-justify the title
  )
)
ggsave("figures/NCTC9087_fig.png", plot=fig_NCTC9087, width =40, height = 35, units = "cm")


write.csv(combined_df, "Data/NCTC9087_2025-01-26.csv")
######################################################
### Only Above Strain-sharing threshold ##############
######################################################

fig3 <- combined_df %>% 
  filter((gapJaccardSim>=0.97) & (singleAgreePct>=99.5))

fig3$ID1 <- paste0(fig3$sample1, "_", fig3$sample2)
obs_pairs <- length(unique(fig3$ID1))

fig3$ID2 <- paste0(fig3$sample2, "_", fig3$sample1)

write.csv(fig3, "NCTC9087_2024-09-18.csv")

### Total Possible Combos = cheep x (focal + mapsan)
possible_pairs <- 101 * (62+177)
obs_pairs
obs_pairs/possible_pairs*100


fig3$human
fig3$human[fig3$comp=="ChEEP x MapSan"] <- "Human (Community)"
fig3$human[fig3$comp=="ChEEP x FOCAL"] <- "Human (Diarrhea)"
table(fig3$human)

fig3$type <- NA
fig3$type[(fig3$sample1 %in% fecal_ids) |(fig3$sample2 %in% fecal_ids) ] <- "Fecal"
fig3$type[(fig3$sample1 %in% carcass_ids) | (fig3$sample2 %in% carcass_ids)] <- "Carcass"
fig3$type[(fig3$sample1 %in% water_ids) | (fig3$sample2 %in% water_ids) ] <- "Rinse Bucket"
table(fig3$type)

fig3$breed <- "Commercial"
fig3$breed[(fig3$sample1 %in% com_ids) |(fig3$sample2 %in% com_ids) ] <- "Commercial"
fig3$breed[(fig3$sample1 %in% ind_ids) | (fig3$sample2 %in% ind_ids)] <- "Local"
table(fig3$breed)

fig3$cheep_id <- ifelse(fig3$sample1 %in% cheep_ids, fig3$sample1, 
                        ifelse(fig3$sample2 %in% cheep_ids, fig3$sample2, NA))

# Generate the frequency table
table(fig3$cheep_id)

# Create the mapsan_id variable 
fig3$mapsan_id <- ifelse(grepl("^SRR", fig3$sample1), fig3$sample1, 
                         ifelse(grepl("^SRR", fig3$sample2), fig3$sample2, NA))
table(fig3$mapsan_id)

# Create the focal_id variable 
fig3$focal_id <- ifelse(grepl("^[0-9]+$", fig3$sample1), fig3$sample1, 
                        ifelse(grepl("^[0-9]+$", fig3$sample2), fig3$sample2, NA))
table(fig3$focal_id)


##### How many mapsan samples per cheep ID 
mapsan <- fig3%>% group_by(mapsan_id) %>% 
  summarise( 
    cheep=n_distinct(cheep_id, na.rm=TRUE))
head(mapsan)
length(mapsan$mapsan_id)

##### How many focal samples per cheep ID 
focal <- fig3%>% group_by(focal_id) %>% 
  summarise( 
    cheep=n_distinct(cheep_id, na.rm=TRUE))
head(focal)
length(focal$focal_id)
# Collapsing data and counting unique focal_id and mapsan_id per cheep_id
cheep <- fig3 %>%
  group_by(cheep_id) %>%
  summarise(
    focal = n_distinct(focal_id, na.rm = TRUE),
    mapsan = n_distinct(mapsan_id, na.rm = TRUE)
  )
head(cheep)
table(cheep$cheep_id)

# Create a new ID variable with ".filter" removed from cheep_id
cheep$ID <- gsub("\\.filter$", "", cheep$cheep_id)

cheep2 <- merge(cheep, cheep_meta, by="ID", all.x=TRUE)
head(cheep2)
cheep2$focal <- ifelse(is.na(cheep2$focal), 0, cheep2$focal)
table(cheep2$focal)
cheep2$mapsan <- ifelse(is.na(cheep2$mapsan), 0, cheep2$mapsan)
table(cheep2$mapsan)

cheep2$human <- cheep2$mapsan + cheep2$focal
summary(cheep2$mapsan)
