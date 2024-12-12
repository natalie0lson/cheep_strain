library(readr)
library(tidyverse)
library(ggplot2)
cheep <- read_tsv("Data/straingst_cheep_results.tsv")
focal <- read_tsv("Data/straingst_focal_results.tsv")
mapsan <- read_tsv("Data/straingst_mapsan_results.tsv")


focal_ids <- read.csv("Data/focal_ids.csv")
mapsan_ids <- read.csv("Data/mapsan_ids.csv")

head(focal_ids)

cheep_meta <- read.csv("Data/phy_meta_2024-04-07.csv")
head(cheep_meta)
cheep_vars<- c("ID", "breed", "Group")
cheep_ids <- cheep_meta[cheep_vars]
head(cheep_ids)

refs <- read_tsv("Data/refs_concat.collapsed.tsv")
ref_13C1079T <- refs$pre_cluster[refs$post_cluster=="Esch_coli_13C1079T"]

#Clean ID Numbers
### Sample Variable
cheep$sample <- gsub("kmerized", "", cheep$sample)
mapsan$sample <- gsub("kmerized", "", mapsan$sample)

### ID Variable
# Remove ".R1" or ".R2" at the end of each value in cheep$sample
cheep$ID <- sub("\\.R[12]$", "", cheep$sample)
table(cheep$ID)

#focal$ID <- substr(focal$sample, 1, 4)
focal$ID <- sub("\\_R[12]$", "", focal$sample)
table(focal$ID)
length(unique(focal$ID))

mapsan$ID <- gsub("_\\d$", "", # gsub("^SRR15209", "", 
                  mapsan$sample)

length(unique(mapsan$ID))

length(unique(cheep$ID))
### Read variable 
mapsan$read <- as.numeric(str_extract(mapsan$sample, "(?<=\\d{3}_)[0-9]{1,3}"))
table(mapsan$read)

focal$read <- as.numeric(substr(focal$sample, nchar(focal$sample), nchar(focal$sample)))
table(focal$read)

cheep$read <- as.numeric(substr(cheep$sample, nchar(cheep$sample), nchar(cheep$sample)))
table(cheep$read)

### Study & Group Variables 
mapsan$Study <- "MapSan"
focal$Study <- "FOCAL"
cheep$Study <- "ChEEP ChEEP"

mapsan$Group <- "Human"
focal$Group <- "Human"
cheep$Group <- "Chicken"

### Merge All Together
all_df <- bind_rows(mapsan, focal, cheep)
head(all_df)
write.csv(all_df, "Data/all_straingst.csv")

## Merge with metadata 
mapsan_meta <- read.csv("Data/mapsan_meta_maputo_only_2024-05-26.csv")
mapsan_meta$ID <- as.numeric(gsub("^Map_", "", mapsan_meta$Sample.ID))
table(mapsan_meta$ID)

## ChEEP Metadata 
cheep_meta <- read.csv("Data/phy_meta_2024-04-07.csv")
cheep_df <- merge(cheep, cheep_meta, by="ID", all=TRUE)

table(cheep_df$location)
market_ids <- cheep_df$ID[cheep_df$location=="Market"]
market_ids

## Clean strain names
all_df$strain2 <- gsub("^Esch_coli_", "", all_df$strain)
table(all_df$strain2)

refs$pre2 <- gsub("^Esch_coli_", "", refs$pre_cluster)
table(refs$pre2)

refs$post2 <- gsub("^Esch_coli_", "", refs$post_cluster)
table(refs$post2)

library(dplyr)
# Perform the left join and create the strain_post variable
all_df <- all_df %>%
  left_join(refs, by = c("strain2" = "pre2")) %>%
  mutate(strain_post = ifelse(!is.na(post2), post2, strain2)) %>%
  select(-post2)  # Optionally remove the post2 column if not needed

# View the updated dataframe
table(all_df$strain_post)
table(all_df$strain2)
## MLST & Phylogroup Data
mlst <- read.csv("Data/strainge_refdb_mlst_phylotype_2024-06-20.csv")
#clean strain names
mlst$strain <- gsub("^Esch_coli_", "", mlst$symlink)
mlst$strain_post <- gsub(".fa.gz", "", mlst$strain)

## Merge MLST & Phylogroup with sample data
all_strains <- merge(mlst, all_df, by="strain_post", all.y=TRUE)

##Subset to high confidence scores
all_conf <- subset(all_strains, all_strains$score >0.2)


head(all_conf)

library(dplyr)


# Count unique pre_cluster values per ID
result <- all_conf %>%
  group_by(ID) %>%
  summarise(unique_pre_clusters = n_distinct(pre_cluster), 
            unique_post_clusters = n_distinct(post_cluster))

# View the result
table(result$unique_pre_clusters)
table(result$unique_post_clusters)


table(result$ID[result$unique_pre_clusters>1])


result$ID <-gsub('\\.filter', '', result$ID)
cheep_gst <- subset(result, (result$ID %in% cheep_ids$ID))
head(cheep_gst)
length(cheep_gst$ID)

cheep_gst2 <- merge(cheep_ids, cheep_gst,by="ID", all.x=TRUE)
head(cheep_gst2)
table(cheep_gst$unique_pre_clusters)
9/33
33/101
101+44+177
178+46+8
17/101
54/177
232/322
#### Strain Summary 
# Load necessary libraries
# Count the number of each strain in each group
strain_count <- all_strains %>%
  group_by(Study, strain_post) %>%
  summarise(Count = round(n()/2)) %>%
  arrange(Study, strain_post)

print(strain_count)

phylo_count <- all_conf %>%
  group_by(phylotype) %>%
  summarise(Count = round(n() / 2)) %>%
  arrange(phylotype) %>%
  mutate(Percent = (Count / sum(Count)) * 100)



### ChEEP ChEEP 
conf_cheep <- subset(all_conf, all_conf$Study=="ChEEP ChEEP")
table(conf_cheep$Study)

phylo_cheep <- conf_cheep %>%
  group_by(phylotype) %>%
  summarise(Count = round(n() / 2)) %>%
  arrange(phylotype) %>%
  mutate(Percent = (Count / sum(Count)) * 100)

print(phylo_cheep)

cheep_strain <- conf_cheep %>%
  group_by(strain2) %>%
  summarise(Count = round(n()/2)) %>%
  
  arrange(desc(Count)) %>%
  mutate(Percent = (Count / sum(Count)) * 100)
print(cheep_strain)

### FOCAL 


conf_focal <- subset(all_conf, all_conf$Study=="FOCAL")
table(conf_focal$Study)

phylo_focal <- conf_focal %>%
  group_by(phylotype) %>%
  summarise(Count = round(n() / 2)) %>%
  arrange(phylotype) %>%
  mutate(Percent = (Count / sum(Count)) * 100)

print(phylo_focal)


focal_strain <- conf_focal %>%
  group_by(strain2) %>%
  summarise(Count = round(n()/2)) %>%
  
  arrange(desc(Count)) %>%
  mutate(Percent = (Count / sum(Count)) * 100)
print(focal_strain)

table(all_conf$MLST[all_conf$strain2=="NCTC9087"])

### MapSan 


conf_mapsan <- subset(all_conf, all_conf$Study=="MapSan")
table(conf_mapsan$Study)

phylo_mapsan <- conf_mapsan %>%
  group_by(phylotype) %>%
  summarise(Count = round(n() / 2)) %>%
  arrange(phylotype) %>%
  mutate(Percent = (Count / sum(Count)) * 100)

print(phylo_mapsan)

mapsan_strain <- conf_mapsan %>%
  group_by(strain2) %>%
  summarise(Count = round(n()/2)) %>%
  
  arrange(desc(Count)) %>%
  mutate(Percent = (Count / sum(Count)) * 100)
print(mapsan_strain)


ecoli <-  strain_count %>% ggplot(aes(x=Study,y=Count)) + 
  geom_bar(aes(fill=strain_post, color=strain_post), stat="identity", position="fill") + 
  ggtitle("Relative Abundance of E. coli Phylotypes") + 
  ylab("Relative Abundance") +xlab("")+ theme_classic() +  theme( legend.title=element_blank(), text=element_text(size=24))
ecoli



# Identify common strains between groups
# Create a list of strains for each group
strains_by_group <- all_conf %>%
  group_by(Group) %>%
  summarise(Strains = list(unique(strain_post)))

# Find common strains between groups
common_strains <- Reduce(intersect, strains_by_group$Strains)

print(common_strains)

all_conf$shared <- NA
all_conf$shared[all_conf$strain_post %in% common_strains] <- 1
all_conf$shared[!(all_conf$strain_post %in% common_strains)]<- 0
table(all_conf$shared)

table(all_conf$Study[all_conf$shared==1])
table(all_conf$Study)

plot_strains1 <- all_conf[all_conf$strain_post %in% common_strains,]
plot_strains1$strainstack <- plot_strains1$strain_post

plot_strains2 <- all_conf[!all_conf$strain_post %in% common_strains,]
plot_strains2$strainstack <- "Other"

plot_strains <- rbind(plot_strains1, plot_strains2)
table(plot_strains$strainstack)

plot_strains$strainstack2 <- NA
plot_strains$strainstack2[plot_strains$strainstack=="13C1079T"] <- "13C1079T (ST 155)"
plot_strains$strainstack2[plot_strains$strainstack=="APEC_102026"] <- "APEC_102026 (ST 117)"
plot_strains$strainstack2[plot_strains$strainstack=="NCTC9087"] <- "NCTC9087 (ST 10)"
plot_strains$strainstack2[plot_strains$strainstack=="ST540"] <- "ST540 (ST 540)"
plot_strains$strainstack2[plot_strains$strainstack=="Other"] <- "Other"
table(plot_strains$strainstack2)


plot_strains$ST155 <- "Other"
plot_strains$ST155[plot_strains$strainstack=="13C1079T"] <-  "ST 155 (13C1079T)"
table(plot_strains$ST155)

plot_strains$ST155 <- factor(plot_strains$ST155, levels=c("ST 155 (13C1079T)", "Other")) 


plot_strains$ST10 <- "Other"
plot_strains$ST10[plot_strains$strainstack=="NCTC9087"] <-  "ST 10 (NCTC9087)"
table(plot_strains$ST10)



plot_strains$ST10 <- factor(plot_strains$ST10, levels=c("ST 10 (NCTC9087)", "Other")) 


strain_count_155 <- plot_strains %>%
  group_by(Study, ST155) %>%
  summarise(Count = round(n()/2)) %>%
  arrange(Study, ST155)



plot_strains$phylotype
phylo_count <- plot_strains %>%
  group_by(phylotype) %>%
  summarise(Count = round(n()/2)) %>%
  arrange(phylotype)

(phylo_count)


strain_count_155$study2 <- NA
strain_count_155$study2[strain_count_155$Study=="ChEEP ChEEP"] <- "Chicken"
strain_count_155$study2[strain_count_155$Study=="FOCAL"] <- "Human (Diarrhea)"
strain_count_155$study2[strain_count_155$Study=="MapSan"] <- "Human (Community)"



head(strain_count_155)

16/(16+22)
17/(17+29)
55/(55+142)

st155 <- strain_count_155 %>% 
  ggplot(aes(x = study2, y = Count #, label = Count
             )) + 
  geom_bar(aes(fill = ST155, color = ST155), stat = "identity", position = "stack") + 
  #geom_text(size = 6, position = position_stack(vjust = 0.8)) + 
  ggtitle("Abundance of E. coli ST 155 (Strain 13C1079T)") + 
  ylab("Number of Samples") + 
  xlab("") + 
  theme_classic() +  
  theme(legend.title = element_blank(), text = element_text(size = 24)) + 
  scale_fill_manual(values = c("ST 155 (13C1079T)" = "#785EF0", "Other" = "#FE6100")) + 
  scale_color_manual(values =c("ST 155 (13C1079T)" = "#785EF0", "Other" = "#FE6100"))

st155


strain_count_10 <- plot_strains %>%
  group_by(Study, ST10) %>%
  summarise(Count = round(n()/2)) %>%
  arrange(Study, ST10)



strain_count_10$study2 <- NA
strain_count_10$study2[strain_count_10$Study=="ChEEP ChEEP"] <- "Chicken"
strain_count_10$study2[strain_count_10$Study=="FOCAL"] <- "Human (Diarrhea)"
strain_count_10$study2[strain_count_10$Study=="MapSan"] <- "Human (Community)"



strain_count_10



16/(16+22)
10/(10+36)
53/(53+144)

st10 <- strain_count_10 %>% 
  ggplot(aes(x = study2, y = Count)) + 
  geom_bar(aes(fill = ST10, color = ST10), stat = "identity", position = "stack") + 
 # geom_text(size = 6, position = position_stack(vjust = 0.85)) + 
  ggtitle("Abundance of E. coli ST 10 (Strain NCTC9087)") + 
  ylab("Number of Samples") + 
  xlab("") + 
  theme_classic() +  
  theme(legend.title = element_blank(), text = element_text(size = 24)) + 
  scale_fill_manual(values = c("ST 10 (NCTC9087)" = "#D81B60", "Other" = "#FFC107")) + 
  scale_color_manual(values =c("ST 10 (NCTC9087)" = "#D81B60", "Other" = "#FFC107"))

st10