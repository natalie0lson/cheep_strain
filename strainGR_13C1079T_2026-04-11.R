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
fecal_ids <- cheep_meta$ID[cheep_meta$Feces==1]
carcass_ids <- cheep_meta$ID[cheep_meta$Carcass==1]
water_ids <- cheep_meta$ID[cheep_meta$Water==1]
com_ids <- cheep_meta$ID[cheep_meta$breed=="Commercial"]
ind_ids <- cheep_meta$ID[cheep_meta$breed=="Indigenous"]

focal_ids <- read.csv("Data/focal_ids.csv")
mapsan_ids <- read.csv("Data/mapsan_ids.csv")

cheep_ids <- cheep_meta %>% dplyr::select(ID)
cheep_ids$study <- "ChEEP ChEEP"
write.csv(cheep_ids, "Data/cheep_ids.csv")
# Define the directory containing the TSV files
#directory <- "Data/13C1079T"

# Get a list of all TSV files in the directory
#file_list <- list.files(path = directory, pattern = "*summary.tsv", full.names = TRUE)

# Read all TSV files into a list of dataframes
#df_list <- lapply(file_list, read.delim)

# List all summary TSV files recursively

file_list <- list.files("Data/straingr_2026", pattern = "comparisons.tsv$", 
                        full.names = TRUE, recursive = TRUE)

# file_list <- list.files("/Users/nmolson/Desktop/cheep_strain/cheep_strain/Data/straingr_2026", pattern = "comparisons.tsv$", 
#                         full.names = TRUE, recursive = TRUE)
file_list
# Function to safely read TSV files
safe_read_tsv <- function(f) {
  nlines <- length(readLines(f, warn = FALSE))
  if (nlines > 1) {  # more than just header
    read.delim(f)
  } else {
    message("Skipping empty or header-only file: ", f)
    return(NULL)
  }
}

# Apply to all files
df_list <- lapply(file_list, safe_read_tsv)

# Remove NULL entries
df_list <- df_list[!sapply(df_list, is.null)]

# Combine all dataframes
combined_df <- do.call(rbind, df_list)

combined_df <- combined_df %>%
  mutate(comp = case_when(
    (grepl("^A|^EUP", sample1)  & grepl("^SRR", sample2)) |
      (grepl("^A|^EUP", sample2)  & grepl("^SRR", sample1)) ~ "Chicken x Child Metagenome",
    
    (grepl("^A|^EUP", sample1) & grepl("^[0-9]+$", sample2)) |
      (grepl("^A|^EUP", sample2) &  grepl("^[0-9]+$", sample1)) ~ "Chicken x Child Isolate",
    
    TRUE ~ NA_character_
  ))

combined_df <- combined_df %>%
  filter(!is.na(comp))

combined_df$ref <- gsub("Esch_coli_|\\.fa\\.gz", "", combined_df$ref)

df_13C1079T <- subset(combined_df, combined_df$ref=="13C1079T")

# Summary statistics
summary_stats <- combined_df %>%
  summarise(across(everything(), list(mean = ~mean(.), sd = ~sd(.), median = ~median(.))))

# Print summary statistics
print(summary_stats)

focal_df <- df_13C1079T %>% filter(comp == "Chicken x Child Isolate")
mapsan_df <- df_13C1079T %>% filter(comp == "Chicken x Child Metagenome")


# Define custom colors for the 'comp' variable
custom_colors <- c("Chicken x Child Isolate" = "cornflowerblue", "Chicken x Child Metagenome" = "goldenrod1")


# Zoomed-in plot (Limited to threshold values)
zoom_13C1079T <- ggplot() +
  # Plot 'ChEEP x MapSan' first so it appears behind 'ChEEP x FOCAL'
  geom_point(data = mapsan_df, aes(x = gapJaccardSim, y = singleAgreePct, size = commonPct, color = comp)) +
  # Plot 'ChEEP x FOCAL' on top
  geom_point(data = focal_df, aes(x = gapJaccardSim, y = singleAgreePct, size = commonPct, color = comp)) +
  scale_x_continuous(limits = c(0.970, 1.0), name = "Gap Similarity") +
  scale_y_continuous(limits = c(99.95, 99.98), name = "ACNI") +
  
  scale_size_continuous(
    name = "Common\nCallable [%]",
    range = c(0.5, 4),       
    limits = c(0.5, 80)    
  ) +

  labs(size = "Common\nCallable [%]") + #ggtitle("E. coli 13C1079T") + 
  theme_minimal() +
  scale_color_manual(values = custom_colors) +
  theme(legend.position = "right",
        legend.title= element_blank(),
        text = element_text(size = 20)) +
  guides(size = guide_legend(title.position = "top", title.hjust = 0.5))
zoom_13C1079T

#save the plots as image files
ggsave("figures/13C1079T_fig.pdf", plot=zoom_13C1079T, width =30, height = 30, units = "cm")


#Full plot (not just threshold)
full_13C1079T <- ggplot() +
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

full_13C1079T


fig_13C1079T <- ggarrange(full_13C1079T,zoom_13C1079T,  cols=2,
                          
                          labels = c("A", "B"),  # Label the plots
                          widths = c(1, 1.8)  # make right plot 30% wider
                          )


#save the plots as image files
ggsave("figures/13C1079T_fig.pdf", plot=fig_13C1079T, width =30, height = 30, units = "cm")


write.csv(df_13C1079T, "13C1079T_2026-04-11.csv")
  ####################################################
  ####### Only Above Strain-sharing threshold ######## 
  ####################################################


atl <- read.delim("Data/straingr_2026/cheep_atl_sample_comparisons_ATL.tsv")
atl$ref <- gsub("Esch_coli_|\\.fa\\.gz", "", atl$ref)
table(atl$ref)

# 
# atl <- atl %>%   filter((singleAgreePct>=99.5) & (commonPct>=0.5))
# table(atl$ref)

atl <- atl %>%   filter((gapJaccardSim>=0.97) & (singleAgreePct>=99.5) & (commonPct>=0.5))
table(atl$ref)

atl <- atl %>%   filter((gapJaccardSim>=0.97) & (singleAgreePct>=99.5) & (commonPct>0.5))
table(atl$ref)

  fig3 <- combined_df %>% 
    filter((gapJaccardSim>=0.97) & (singleAgreePct>=99.5))
table(fig3$ref)

table(fig3$ref)


  fig3$ID1 <- paste0(fig3$sample1, "_", fig3$sample2)
 
  obs_pairs <- length(unique(fig3$ID1))
  obs_pairs
  2391/16050
  436/16050
  
  fig3$ID2 <- paste0(fig3$sample2, "_", fig3$sample1)
  
  
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
  
  
  
  
  write.csv(fig3, "Data/combined_df.csv")
  
  
  
  mapsan <- fig3%>% group_by(mapsan_id) %>% 
    summarise( 
      cheep=n_distinct(cheep_id, na.rm=TRUE))
  head(mapsan)
  length(mapsan$mapsan_id)
  
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
  
  model <- lm(human~ Indigenous + Feces, data=cheep2)
  summary(model)
  
  #model <- glm.nb(human ~ Feces, data = cheep2)
  summary(model)
  
  model <- glm.nb(human ~ Indigenous, data = cheep2)
  summary(model)
  
  sd(cheep2$human[cheep2$Indigenous==0])
  
  cheep2$human <- cheep2$mapsan + cheep2$focal
  summary(cheep2$mapsan)
  
  ### Odds Ratio - Shared Strain w/ Human or Not 
  cheep2$shared <- 0
  cheep2$shared[(cheep2$focal>0) | (cheep2$mapsan>0)] <- 1 
  table(cheep2$shared)
  
  table(cheep2$shared[cheep2$breed=="Commercial"])
  #create matrix
  shared <- c('No', 'Yes')
  breed <- c('Commercial', 'Local')
  data <- matrix(c(51, 3, 10, 11), nrow=2, ncol=2, byrow=TRUE)
  dimnames(data) <- list( 'Breed'=breed, 'Shared Strain'=shared)
  
  #calculate odds ratio
  oddsratio(data, method="fisher")
  
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
  cheep_long$xnum[cheep_long$loc.type.n=="Commercial/Market/Carcass, N=42"] <- "Market \n (Carcass) \n N=1"
  cheep_long$xnum[cheep_long$loc.type.n=="Commercial/Farm/Fecal, N=4"] <- "Farm \n (Fecal) \n N=2"
  cheep_long$xnum[cheep_long$loc.type.n=="Indigenous/Household/Fecal, N=7"] <- "Household \n (Fecal) \n N=5"
  cheep_long$xnum[cheep_long$loc.type.n=="Indigenous/Market/Fecal, N=7"] <- "Market \n (Fecal) \n N=6"
  cheep_long$xnum[cheep_long$loc.type.n=="Market Rinse Water, N=26"] <- "Market \n (Rinse Bucket) \n N=2"
  
  table(cheep_long$ID)
  table(cheep_long$xnum)
  
  cheep_long$xnum <- factor(cheep_long$xnum, levels=c("Household \n (Fecal) \n N=5", 
                                                      "Farm \n (Fecal) \n N=2",
                                                      "Market \n (Fecal) \n N=6",
                                                      "Market \n (Carcass) \n N=1", 
                                                      "Market \n (Rinse Bucket) \n N=2" ))
  summary(cheep_long$count)
  
  ggplot(cheep_long, aes(x=xnum, y=count, fill=breed2)) + 
    geom_violin(alpha=0.5
      ) +  # Creating the violin plot with some transparency
    xlab("") +                # Removing the x-axis label
    theme_classic() +         # Applying a classic theme
    stat_summary(aes(color=breed2), 
                 fun.data=mean_sdl, fun.args = list(mult = 1),
                 geom="pointrange", 
                 shape = 18, size = 0.75,
                 position = position_dodge(width = 0.9)) + 
    theme(legend.title = element_blank(), text = element_text(size = 17)) + 
    scale_fill_manual(values = c("Commercial" = "#785EF0", "Local" = "#FE6100")) + 
    scale_color_manual(values = c("Commercial" = "#785EF0", "Local" = "#FE6100")) + 
    ylab("Human Samples") + 
    ggtitle("E. coli Strain 13C1079T Shared Between Chicken and Human Samples \n (ACNI > 99.5% & Gap Similarity > 97%)") + 
    facet_wrap(~study2)
  zoom_13C1079T
  
  figs <- ggarrange(zoom_NCTC9087,zoom_13C1079T,  cols=2, 
                    labels = c("A", "B"),  # Label the plots
                    widths = c(1, 1.8)  # make right plot 30% wider
                    )
  
  #save the plots as image files
  ggsave("figures/13C1079T_NCTC9087_fig.pdf", plot=figs, width =30, height = 30, units = "cm")
  
  