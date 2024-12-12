# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(epitools)
library(MASS)
library(tidyr)

## ChEEP Metadata 
cheep_meta <- read.csv("phy_meta_2024-04-07.csv")

## Reformat IDs 
cheep_ids <- paste0(cheep_meta$ID, ".filter")

market_ids <- cheep_meta$ID[cheep_meta$location=="Market"]
market_ids <- paste0(market_ids, ".filter")

fecal_ids <- cheep_meta$ID[cheep_meta$Feces==1]
fecal_ids <- paste0(fecal_ids, ".filter")

carcass_ids <- cheep_meta$ID[cheep_meta$Carcass==1]
carcass_ids <- paste0(carcass_ids, ".filter")

water_ids <- cheep_meta$ID[cheep_meta$Water==1]
water_ids <- paste0(water_ids, ".filter")

com_ids <- cheep_meta$ID[cheep_meta$breed=="Commercial"]
com_ids <- paste0(com_ids, ".filter")

ind_ids <- cheep_meta$ID[cheep_meta$breed=="Indigenous"]
ind_ids <- paste0(ind_ids, ".filter")

# Define the directory containing the TSV files
directory <- "NCTC9087"

# Get a list of all TSV files in the directory
file_list <- list.files(path = directory, pattern = "*summary.tsv", full.names = TRUE)

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
zoom_NCTC9087 <- ggplot() +
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

zoom_NCTC9087

#Full plot (not just threshold)
full_NCTC9087 <- ggplot() +
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

full_NCTC9087


fig_NCTC9087 <- ggarrange(full_NCTC9087,zoom_NCTC9087,  cols=2)

#save the plots as image files
ggsave("figures/NCTC9087_fig.png", plot=fig_NCTC9087, width =30, height = 30, units = "cm")



write.csv(combined_df, "NCTC9087_2024-10-29.csv")
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
possible_pairs <- 75 * (43+171)

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

### Odds Ratio - Shared Strain w/ Human or Not 
cheep2$shared <- 0
cheep2$shared[(cheep2$focal>0) | (cheep2$mapsan>0)] <- 1 
table(cheep2$shared)

table(cheep2$shared[cheep2$Water==1])
#create matrix
shared <- c('No', 'Yes')
breed <- c('Commercial', 'Local')
data <- matrix(c(46, 8, 18, 3), nrow=2, ncol=2, byrow=TRUE)
dimnames(data) <- list( 'Breed'=breed, 'Shared Strain'=shared)


8/(8+46)
3/21
#view matrix
data

#calculate odds ratio
oddsratio(data, method="fisher")

table(cheep2$Market)
model <- lm(focal~Indigenous + Feces + Market , data=cheep2)
summary(model)
model <- glm.nb(human ~ Indigenous, data = cheep2)
summary(model)

sd(cheep2$human[cheep2$Indigenous==1])

cheep_long <- cheep2 %>%
  pivot_longer(cols = c(mapsan, focal), 
               names_to = "study", 
               values_to = "count")

# View the reshaped dataset
table(cheep_long$study)

cheep_long$study2 <- NA
cheep_long$study2[cheep_long$study=="focal"] <- "Human (Diarrhea)"
cheep_long$study2[cheep_long$study=="mapsan"] <- "Human (Community)"
table(cheep_long$study2)

cheep_long$breed2 <- "Commercial"
cheep_long$breed2[cheep_long$breed=="Commercial"] <- "Commercial"
cheep_long$breed2[cheep_long$breed=="Indigenous"] <- "Local"

table(cheep_long$Group)

cheep_long$Group <- factor(cheep_long$Group, levels=c("Production", "Market", "Processing"))

table(cheep_long$loc.type.n)

cheep_long$xnum <- NA
cheep_long$xnum[cheep_long$loc.type.n=="Commercial/Market/Carcass, N=42"] <- "Market \n (Carcass) \n N=5"
cheep_long$xnum[cheep_long$loc.type.n=="Commercial/Market/Fecal, N=8"] <- "Market \n (Fecal) \n N=3"
cheep_long$xnum[cheep_long$loc.type.n=="Indigenous/Market/Fecal, N=7"] <- "  Market \n (Fecal) \n N=3"
cheep_long$xnum[cheep_long$loc.type.n=="Market Rinse Water, N=26"] <- "Market \n (Rinse Bucket) \n N=3"

cheep_long$xnum[cheep_long$loc.type.n=="Indigenous/Household/Fecal, N=7"] <- "Production \n (Fecal) \n N=1"

cheep_long$xnum[cheep_long$loc.type.n=="Commercial/Farm/Fecal, N=4"] <- "Production \n (Fecal) \n N=2"
table(cheep_long$xnum)

table(cheep_long$ID)

cheep_long$xnum <- factor(cheep_long$xnum, levels=c("Production \n (Fecal) \n N=1",
                                                    "Production \n (Fecal) \n N=2",
                                                    "Market \n (Fecal) \n N=3",
                                                    "  Market \n (Fecal) \n N=3" ,
                                                    "Market \n (Carcass) \n N=5",
                                                    "Market \n (Rinse Bucket) \n N=3" ))

table(cheep_long$loc.type.n)
ggplot(cheep_long, aes(x=xnum, y=count, fill=breed2)) + 
  geom_violin(#trim = FALSE, position="dodge",
    alpha=0.5
    ) +
  xlab("") +
  theme_classic() + 
  stat_summary(aes(color=breed2), 
               fun.data=mean_sdl, fun.args = list(mult = 1),
               geom="pointrange",
               shape = 18, size = 0.75,
               position = position_dodge(width = 0.9))+ 
  theme(legend.title = element_blank(), text = element_text(size = 16)) + 
scale_fill_manual(values = c("Commercial" = "#D81B60", "Local" = "#FFC107")) + 
scale_color_manual(values = c("Commercial" = "#D81B60", "Local" = "#FFC107")) + 
  ylab("Human Samples") + 
  ggtitle("E. coli Strain NCTC9087 Shared Between Chicken and Human Samples \n (ACNI > 99.5% & Gap Similarity > 97%)") + 
  facet_wrap(~study2) 
