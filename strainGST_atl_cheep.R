library(readr)
library(tidyverse)
library(ggplot2)
library(dplyr)

cheep <- read_tsv("Data/straingst_cheep_results.tsv")
atl <- read_tsv("Data/straingst_atl_results.tsv")

cheep_meta <- read.csv("Data/phy_meta_2024-04-07.csv")
cheep_abundance <- read.csv("Data/ecoli_abundance_.csv")
cheep_vars<- c("ID", "breed", "Group")
cheep_ids <- cheep_meta[cheep_vars]

refs <- read_tsv("Data/refs_concat_atl.collapsed.tsv")
length(unique(refs$pre_cluster))
length(unique(refs$post_cluster))
ref_13C1079T <- refs$pre_cluster[refs$post_cluster=="Esch_coli_13C1079T"]

#Clean ID Numbers
### Sample Variable
cheep$sample <- gsub("kmerized", "", cheep$sample)
atl$sample <- gsub("kmerized", "", atl$sample)

### ID Variable
# Remove ".R1" or ".R2" at the end of each value in cheep$sample
cheep$ID <- sub("\\.R[12]$", "", cheep$sample)

#focal$ID <- substr(focal$sample, 1, 4)
atl$ID <- gsub("_\\d$", "",  atl$sample)

length(unique(cheep$ID))
length(unique(atl$ID))
length(unique(atl$strain))
40/2904
### Read variable 
atl$read <- as.numeric(str_extract(atl$sample, "(?<=\\d{3}_)[0-9]{1,3}"))
table(atl$read)

cheep$read <- as.numeric(substr(cheep$sample, nchar(cheep$sample), nchar(cheep$sample)))
table(cheep$read)

### Study & Group Variables 
cheep$Study <- "ChEEP ChEEP"
atl$Study <- "Atlanta"

cheep$Group <- "Chicken"
atl$Group <- "Human"

### Merge All Together
all_df <- bind_rows(atl, cheep)
head(all_df)
write.csv(all_df, "Data/all_straingst_atl_cheep.csv")

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
  mutate(strain_post = ifelse(!is.na(post2), post2, strain2)) 

# View the updated dataframe
table(all_df$strain_post)
table(all_df$strain2)

head(all_df)
length(unique(all_df$ID))
table(all_df$ID)
134/2904
## MLST & Phylogroup Data
mlst <- read.csv("Data/strainge_refdb_mlst_phylotype_2024-06-20.csv")
#clean strain names
mlst$strain <- gsub("^Esch_coli_", "", mlst$symlink)
mlst$strain_post <- gsub(".fa.gz", "", mlst$strain)

## Merge MLST & Phylogroup with sample data
all_strains <- merge(mlst, all_df, by="strain_post", all.y=TRUE)

##Subset to high confidence scores
all_conf <- subset(all_strains, all_strains$score >0.2)
atl_conf <- subset(atl, atl$score>0.2)
length(unique(atl_conf$ID))
33/60
29/33
1/33
3/33
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

head(cheep_gst)
length(cheep_gst$ID)

cheep_gst2 <- merge(cheep_ids, cheep_gst,by="ID", all.x=TRUE)
head(cheep_gst2)
table(cheep_gst$unique_pre_clusters)

atl_gst <- subset(result, result$ID %in% atl$ID)
table(atl_gst$unique_pre_clusters)


#### Strain Summary 
# Load necessary libraries
# Count the number of each strain in each group
strain_count <- all_strains %>%
  group_by(Study, strain_post) %>%
  summarise(Count = round(n()/2)) %>%
  arrange(Study, strain_post)

print(strain_count)
head(all_conf)

all_atl <- subset(all_conf, all_conf$Study=="Atlanta")
phylo_count <- all_atl %>%
  group_by(phylotype) %>%
  summarise(Count = round(n() / 2)) %>%
  arrange(phylotype) %>%
  mutate(Percent = (Count / sum(Count)) * 100)
phylo_count

strain_count <- all_atl %>%
  group_by(strain_post) %>%
  summarise(Count = round(n() / 2)) %>%
  arrange(strain_post) %>%
  mutate(Percent = (Count / sum(Count)) * 100)
strain_count

all_atl$MLST[all_atl$strain_post=="Z0117EC0054"]

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

table(all_conf$MLST[all_conf$strain2=="NCTC9087"])

ecoli <-  strain_count %>% ggplot(aes(x=Study,y=Count)) + 
  geom_bar(aes(fill=strain_post, color=strain_post), stat="identity", position="fill") + 
  ggtitle("Relative Abundance of E. coli Phylotypes") + 
  ylab("Relative Abundance") +xlab("")+ theme_classic() +  theme( legend.title=element_blank(), text=element_text(size=24))
ecoli

table(all_conf$Group)
table(all_conf$Study)
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

table(all_conf$post_cluster[all_conf$shared==1])

table(all_conf$Study[all_conf$shared==1])
table(all_conf$Study)



strain_13C1079T <- subset(all_conf, all_conf$strain_post=="13C1079T")
head(strain_13C1079T)
strain_13C1079T$Study
df_13C1079T <- data.frame(ID = strain_13C1079T$ID, Study = strain_13C1079T$Study)
table(df_13C1079T$Study)
write.csv(df_13C1079T, "Data/df_13C1079T_atl_cheep.csv")


strain_D4 <- subset(all_conf, all_conf$strain_post=="D4")
head(strain_D4)
strain_D4$Study
df_D4 <- data.frame(ID = strain_D4$ID, Study = strain_D4$Study)
table(df_D4$Study)
write.csv(df_D4, "Data/df_D4.csv")

strain_KW21 <- subset(all_conf, all_conf$strain_post=="KW21")
head(strain_KW21)
strain_KW21$Study
df_KW21 <- data.frame(ID = strain_KW21$ID, Study = strain_KW21$Study)
table(df_KW21$Study)
write.csv(df_KW21, "Data/df_KW21.csv")

plot_strains1 <- all_conf[all_conf$strain_post %in% common_strains,]
plot_strains1$strainstack <- plot_strains1$strain_post

plot_strains2 <- all_conf[!all_conf$strain_post %in% common_strains,]
plot_strains2$strainstack <- "Other"

plot_strains <- rbind(plot_strains1, plot_strains2)
table(plot_strains$strainstack)

table(plot_strains$MLST[plot_strains$strainstack=="KW21"])

plot_strains$strainstack2 <- NA
plot_strains$strainstack2[plot_strains$strainstack=="13C1079T"] <- "13C1079T (ST 155)"
plot_strains$strainstack2[plot_strains$strainstack=="D4"] <- "D4 (MLST Unknown)"
plot_strains$strainstack2[plot_strains$strainstack=="NCTC9087"] <- "KW21 (ST 224)"
plot_strains$strainstack2[plot_strains$strainstack=="Other"] <- "Other"
table(plot_strains$strainstack2)

plot_strains$ST155 <- "Other"
plot_strains$ST155[plot_strains$strainstack=="13C1079T"] <-  "ST 155 (13C1079T)"
table(plot_strains$ST155)

plot_strains$ST155 <- factor(plot_strains$ST155, levels=c("ST 155 (13C1079T)", "Other")) 

plot_strains$ST224 <- "Other"
plot_strains$ST224[plot_strains$strainstack=="KW21"] <-  "ST 224 (KW21)"
table(plot_strains$ST224)

plot_strains$ST224 <- factor(plot_strains$ST224, levels=c("ST 224 (KW21)", "Other")) 

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
strain_count_155$study2[strain_count_155$Study=="ChEEP ChEEP"] <- "Chicken (Maputo)"
strain_count_155$study2[strain_count_155$Study=="Atlanta"] <- "Human (Atlanta)"

head(strain_count_155)

st155 <- strain_count_155 %>% 
  ggplot(aes(x = study2, y = Count #, label = Count
             )) + 
  geom_bar(aes(fill = ST155, color = ST155), stat = "identity", position = "stack") + 
  #geom_text(size = 6, position = position_stack(vjust = 0.8)) + 
  ggtitle("Frequency of E. coli ST 155 (Strain 13C1079T) Detection") + 
  ylab("Number of Samples") + 
  xlab("") + 
  theme_classic() +  
  theme(legend.title = element_blank(), text = element_text(size = 24)) + 
  scale_fill_manual(values = c("ST 155 (13C1079T)" = "#785EF0", "Other" = "#FE6100")) + 
  scale_color_manual(values =c("ST 155 (13C1079T)" = "#785EF0", "Other" = "#FE6100"))

st155

