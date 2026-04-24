library(tidyverse)
library(dplyr)

#Read datasets
mlst <- read.csv("Data/strainge_refdb_mlst_phylotype_2024-06-20.csv")
atl <- read_tsv("Data/straingst_2026/atl_straingst_strains.tsv")
cheep <- read_tsv("Data/straingst_2026/cheep_straingst_strains.tsv")
focal <- read_tsv("Data/straingst_2026/focal_straingst_strains.tsv")
mapsan <- read_tsv("Data/straingst_2026/mapsan_straingst_strains.tsv")
focal
mlst$strain <- mlst$symlink
cheep$study <- "ChEEP ChEEP"
focal$study <- "Focal"
mapsan$study <- "MapSan"
#Build table

### Create dataset with all Maputo samples (exclude ATL samples)
maputo <- bind_rows(focal, cheep, mapsan)
sort(table(atl$strain), decreasing = TRUE)
head(mapsan)

strain_counts <- atl %>%
  group_by(sample) %>%
  summarise(n_strains = n_distinct(strain))

summary_counts <- strain_counts %>%
  count(n_strains)

#What strains were common between cheep and focal
common_focal<-bind_rows(focal, cheep) %>% group_by(study, strain) %>% slice(1) %>% ungroup() %>% group_by(strain) %>% summarize(count=n()) %>% filter(count==2)

#What strains were common between cheep and mapsan
common_mapsan<-bind_rows(mapsan, cheep) %>% group_by(study, strain) %>% slice(1) %>% ungroup() %>% group_by(strain) %>% summarize(count=n()) %>% filter(count==2)

#What strains were detected across at least 2 of three
common<-bind_rows(common_focal, common_mapsan) %>% group_by(strain) %>% slice(1)

#Make count table
counts<-bind_rows(mapsan, cheep, focal) %>% group_by(strain, study) %>%
  summarize(count=n()) %>% filter(strain %in% common$strain) %>% 
  pivot_wider(names_from = "study", values_from = "count") %>%
  mutate(strain=gsub("_new","",strain)) %>% left_join(mlst) %>%
  mutate(strain=gsub("Esch_coli_","",strain), strain=gsub(".fa.gz","",strain))


final_table <- counts %>%
  rename(
    `Reference Genome` = strain,
    `Chicken metagenome (N=101)` = `ChEEP ChEEP`,
    `Child isolate (N=60)` = Focal,
    `Child metagenome (N=177)` = MapSan,
    `Phylotype` = phylotype
  ) %>%
  select(
    `Reference Genome`,
    MLST,
    Phylotype,
    `Chicken metagenome (N=101)`,
    `Child isolate (N=60)`,
    `Child metagenome (N=177)`
  )

final_table

write.csv(final_table, "final_table.csv", row.names = FALSE)

## Create "Group" variable to distinguish human from chicken samples 
mapsan$Group <- "Human"
focal$Group <- "Human"
cheep$Group <- "Chicken"

all_df <- bind_rows(mapsan, focal, cheep)

# Identify common strains between groups
# Create a list of strains for each group
strains_by_group <- all_df %>%
  group_by(Group) %>%
  summarise(Strains = list(unique(strain_post)))

# Find common strains between groups
common_strains <- Reduce(intersect, strains_by_group$Strains)

