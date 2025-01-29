#load required packages
library(phyloseq)
library(vegan)
library(biomformat)
library(readr)
#Directories
fig_dir <- "C:/Users/nmolson/OneDrive - Emory University/cheep_strain/figures/"
data_dir <- "C:/Users/nmolson/OneDrive - Emory University/cheep_strain/Data/"

#Import biom file into phyloseq object
today <- format(Sys.Date(), "%Y-%m-%d")

biom_file <- "Data/bracken.species.biom"
data <- import_biom(biom_file, parseFunction=parse_taxonomy_default)

# pull OTU table from imported biom file
OTU <- otu_table(data@otu_table, taxa_are_rows = TRUE)
# Get the current column names
current_names <- colnames(OTU)
head(current_names)

# Use gsub to remove everything after the underscore in each column name
new_names <- gsub("_.*", "", current_names)
head(new_names)

# Assign the new column names to the OTU table
colnames(OTU) <- new_names
head(OTU)

# pull taxonomy table from imported biom file
TAX <- tax_table(data@tax_table)

# Extract characters after underscore
TAX@.Data <- gsub(".*_", "", TAX@.Data)

# Print the modified dataframe
head(TAX@.Data)
head(OTU@.Data)

# use OTU table sample names to label sample data object
#sample_names(sampledata)  <- sample_names(OTU)

## sampledata
sampledata <- sample_data(data@sam_data)

# Rename the values in the sample data
sampledata@row.names <- gsub("_.*", "", sampledata@row.names)
sampledata$Id <- gsub("_.*", "", sampledata$Id)

########### make phyloseq object
physeq <- phyloseq(TAX, sampledata, OTU)
physeq@sam_data
class(physeq@sam_data$Id)

# label taxonomy table columns
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
physeq@tax_table

########### pre-process/clean phyloseq object
# remove taxa with unassigned sequences
physeq <- subset_taxa(physeq, Genus != "-1")

physeq2 <- subset_taxa(physeq, Kingdom!="Bacteria")
physeq@sam_data$nonbac_count <- sample_sums(physeq2)
table(physeq@sam_data$nonbac_count)

physeq <- subset_taxa(physeq, Kingdom=="Bacteria")
physeq@sam_data$bac_count <- sample_sums(physeq)

physeq <- filter_taxa(physeq, function(x) sum(x > 2) > (0.2*length(x)), TRUE)

otu_table_df <- physeq@otu_table
head(otu_table_df)

head(physeq@sam_data)
head(physeq@tax_table)
nsamples(physeq)
sample_names(physeq)
ntaxa(physeq)

# check out our taxonomy table
rank_names(physeq)

################################
################################
#### Save Phyloseq Object ###### 
################################
################################ 

write_rds(physeq,  paste0("mapsan_phyloseq_obj_", today, ".rds"))

phyloseq_obj <- paste0("mapsan_phyloseq_obj_", today, ".rds")
read_phy <- readRDS(phyloseq_obj)
(read_phy@sam_data$ID)

### write dataframes from phyloseq object(s) 
bac.tax<-as.data.frame(physeq@tax_table)
bac.otu <- as.data.frame(physeq@otu_table)

row_names <- row.names(bac.otu)
bac.otu$tax_id <- row_names

row_names <- row.names(bac.tax)
bac.tax$tax_id <- row_names

merged <- merge(bac.tax, bac.otu, by="tax_id", all.y=T)
head(merged)

library(tidyr)
df_long <- merged%>% 
  gather(key="ID", value="count", -c("tax_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
head(df_long)


# Genus & species name
df_long$g.species <- paste0(substr(df_long$Genus, 1, 1), ".", " ", df_long$Species)
table(df_long$g.species)
head(df_long)
summary(df_long$count[df_long$g.species=="E. coli"])


library(dplyr)
# Summarize total count for each ID
summed_data <- df_long %>%
  group_by(ID) %>%
  summarise(total_count = sum(count, na.rm = TRUE))
head(summed_data)

ecoli_count <- subset(df_long, df_long$g.species=="E. coli")
head(ecoli_count)

relabund <- merge(ecoli_count, summed_data, by="ID", all=TRUE)
head(relabund)
relabund$relabund <- relabund$count/relabund$total_count 
summary(relabund$relabund)

vars <- c("ID", "count", "total_count")
mapsan_ecoli <- relabund[vars]
head(mapsan_ecoli)
write.csv(mapsan_ecoli, "Data/mapsan_ecoli.csv" )
