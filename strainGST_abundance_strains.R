library(readr)
library(tidyverse)
library(ggplot2)
library(dplyr)

cheep <- read_tsv("Data/straingst_cheep_results.tsv")
mapsan <- read_tsv("Data/straingst_mapsan_results.tsv")
atl <- read_tsv("Data/straingst_atl_results.tsv")
mapsan_ids <- read.csv("Data/mapsan_ids.csv")
cheep_meta <- read.csv("Data/phy_meta_2024-04-07.csv")
cheep_abundance <- read.csv("Data/ecoli_abundance_.csv")
cheep_vars<- c("ID", "breed", "Group")
cheep_ids <- cheep_meta[cheep_vars]

refs <- read_tsv("Data/refs_concat4.collapsed.tsv")
refs$post_cluster

#ref_13C1079T <- refs$pre_cluster[refs$post_cluster=="Esch_coli_13C1079T"]

#Clean ID Numbers
### Sample Variable
cheep$sample <- gsub("kmerized", "", cheep$sample)
mapsan$sample <- gsub("kmerized", "", mapsan$sample)
atl$sample <- gsub("kmerized", "", atl$sample)

### ID Variable
# Remove ".R1" or ".R2" at the end of each value in cheep$sample
cheep$ID <- sub("\\.R[12]$", "", cheep$sample)

mapsan$ID <- gsub("_\\d$", "",  mapsan$sample)
atl$ID <- gsub("_\\d$", "",  atl$sample)

### Read variable 
mapsan$read <- as.numeric(str_extract(mapsan$sample, "(?<=\\d{3}_)[0-9]{1,3}"))
table(mapsan$read)

atl$read <- as.numeric(str_extract(atl$sample, "(?<=\\d{3}_)[0-9]{1,3}"))
table(atl$read)

cheep$read <- as.numeric(substr(cheep$sample, nchar(cheep$sample), nchar(cheep$sample)))
table(cheep$read)

### Study & Group Variables 
mapsan$Study <- "MapSan"
cheep$Study <- "ChEEP ChEEP"
atl$Study <- "Atlanta"

mapsan$Group <- "Human"
cheep$Group <- "Chicken"
atl$Group <- "Human"

### Merge All Together
all_df <- bind_rows(mapsan, focal, atl, cheep)
head(all_df)
write.csv(all_df, "Data/all_straingst.csv")

## Merge with metadata 
mapsan_meta <- read.csv("Data/mapsan_meta_maputo_only_2024-05-26.csv")
mapsan_meta$ID <- as.numeric(gsub("^Map_", "", mapsan_meta$Sample.ID))

## ChEEP Metadata 
cheep_meta <- read.csv("Data/phy_meta_2024-04-07.csv")
cheep_df <- merge(cheep, cheep_meta, by="ID", all=TRUE)

market_ids <- cheep_df$ID[cheep_df$location=="Market"]

## Clean strain names
all_df$strain2 <- gsub("^Esch_coli_", "", all_df$strain)

refs$pre2 <- gsub("^Esch_coli_", "", refs$pre_cluster)

refs$post2 <- gsub("^Esch_coli_", "", refs$post_cluster)

library(dplyr)
# Perform the left join and create the strain_post variable
all_df <- all_df %>%
  left_join(refs, by = c("strain2" = "pre2")) %>%
  mutate(strain_post = ifelse(!is.na(post2), post2, strain2)) 

## MLST & Phylogroup Data
mlst <- read.csv("Data/strainge_refdb_mlst_phylotype_2024-06-20.csv")
#clean strain names
mlst$strain <- gsub("^Esch_coli_", "", mlst$symlink)
mlst$strain_post <- gsub(".fa.gz", "", mlst$strain)

## Merge MLST & Phylogroup with sample data
all_strains <- merge(mlst, all_df, by="strain_post", all.y=TRUE)

##Subset to high confidence scores
all_conf <- subset(all_strains, all_strains$score >0.2)

# Count unique pre_cluster values per ID
result <- all_conf %>%
  group_by(ID) %>%
  summarise(unique_pre_clusters = n_distinct(pre_cluster), 
            unique_post_clusters = n_distinct(post_cluster))

table(result$ID[result$unique_pre_clusters>1])

result$ID <-gsub('\\.filter', '', result$ID)
cheep_gst <- subset(result, (result$ID %in% cheep_ids$ID))

all_conf$ID <-gsub('\\.filter', '', all_conf$ID)
cheep_conf <- subset(all_conf, all_conf$ID %in% cheep_ids$ID)
length(unique(cheep_conf$pre_cluster))

cheep_gst2 <- merge(cheep_ids, cheep_gst,by="ID", all.x=TRUE)
head(cheep_gst2)
table(cheep_gst$unique_pre_clusters)

#ChEEP ChEEP Strains per Sample vs. E.coli Abundance 
head(cheep_abundance)
vars <- c("ID", "count", "bac_count")
cheep_abund <- cheep_abundance[vars]
cheep_strain_abund <- merge(cheep_gst, cheep_abund, by="ID")
head(cheep_strain_abund)
cheep_strain_abund$study <- "ChEEP ChEEP"

cheep_strain_abund$ecoli_abund <- cheep_strain_abund$count/cheep_strain_abund$bac_count*100
summary(cheep_strain_abund$ecoli_abund)

vars <- c("ID", "unique_pre_clusters", "ecoli_abund", "study")
cheep_strain_abund2 <- cheep_strain_abund[vars]

## Mapsan Strains per sample vs relative abundance 
mapsan_gst <- subset(result, result$ID %in% mapsan_ids$Filename)
head(mapsan_gst)
table(mapsan_gst$unique_pre_clusters)
mapsan_rel <- read.csv("Data/mapsan_ecoli.csv")
head(mapsan_rel)
mapsan_rel$ecoli_abund <- mapsan_rel$count/mapsan_rel$total_count *100 
summary(mapsan_rel$ecoli_abund)
sd(mapsan_rel$ecoli_abund)
mapsan_strain_abund <- merge(mapsan_gst, mapsan_rel, by="ID")

mapsan_strain_abund$study <- "MapSan"
mapsan_strain_abund2 <- mapsan_strain_abund[vars]
head(mapsan_strain_abund2)

## Atlanta Strains per sample vs relative abundance 

atl_gst <- subset(result, result$ID %in% atl$ID)
head(atl_gst)
table(atl_gst$unique_pre_clusters)
atl_rel <- read.csv("Data/atl_ecoli.csv")
head(atl_rel)
atl_rel$ecoli_abund <- atl_rel$count/atl_rel$total_count *100 
summary(atl_rel$ecoli_abund)
sd(atl_rel$ecoli_abund)
atl_strain_abund <- merge(atl_gst, atl_rel, by="ID")

atl_strain_abund$study <- "Atlanta"
atl_strain_abund2 <- atl_strain_abund[vars]
head(atl_strain_abund2)

strain_abund <- rbind(mapsan_strain_abund2, cheep_strain_abund2, atl_strain_abund2)
head(strain_abund)

strain_abund$unique_pre_clusters <- as.character(strain_abund$unique_pre_clusters)

ggplot(strain_abund, aes(x=unique_pre_clusters, y=ecoli_abund, fill=study)) + 
  
  geom_boxplot(alpha=0.5, position = position_dodge2(preserve = "single")) +
  
 xlab("Unique E. coli Strains") +
  theme_classic() +  scale_y_continuous(labels = scales::comma) +
  
  theme(legend.title = element_blank(), text = element_text(size = 16)) + 
  ylab("E. coli Relative Abundance (%)") + 
  ggtitle("E. coli Relative Abundance vs. Unique E. coli Strains Detected") #+ 
 # facet_wrap(~study2) 

lm = lm(unique_pre_clusters~ecoli_abund, data = strain_abund)
summary(lm)
