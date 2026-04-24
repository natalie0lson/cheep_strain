library(ggplot2)
library(scales)
###############################################################
######## Virulence Factors & AMR  in FOCAL Samples ############
###############################################################
shared_focal <- subset(shared, shared$focal==1)
length(unique(shared_focal$human_id[shared_focal$threshold==1]))

focal_ids$human_id <- focal_ids$Filename 
mapsan_ids$human_id <- mapsan_ids$Filename

all_focal <- merge(focal_ids, shared_focal, by="human_id", all.x=TRUE)

########## ARGs ###############

amr <- read.delim("Data/amr.summary.tab")

amr_long <- gather(amr, arg, value, aac.3..IId:tet.D., factor_key=TRUE)
amr_long$human_id <- gsub("\\.tab", "", amr_long$X.FILE)

amr_long$value <- as.numeric(amr_long$value)
amr_long <- amr_long[!is.na(amr_long$value), ]

amr_shared <- merge(shared_focal, amr_long, by="human_id", all.x=TRUE)

# Remove periods from the 'varia ble' column
amr_shared$arg <- gsub("\\.", "", amr_shared$arg)

hr_args <- read.csv("Data/clinically_relevant_amr_zhang_2022.csv")

hr_args_q1 <- subset(hr_args, hr_args$Rank=="Q1")

hr_args_q1$HR_ARG <- 1 

# Remove hyphens and parentheses from the ARG name 
hr_args_q1$arg  <- gsub("[-()']", "", hr_args_q1$ARO.Term)

hr_args_shared <- merge(amr_shared,hr_args_q1,  by="arg", all.x=TRUE)

hr_args2 <- hr_args_shared[!is.na(hr_args_shared$ARG.Class), ]

hr_args2$shared[hr_args2$threshold==0] <- "Not Shared"
hr_args2$shared[hr_args2$threshold==1] <- "Shared"

strains <- c("13C1079T", "NCTC9087")
hr_args3 <- subset(hr_args2, hr_args2$strain %in% strains)

hr_args3$Class <- hr_args3$ARG.Class

arg_chart <- hr_args3 %>%
  group_by(arg, strain) %>% 
  summarise(count = n_distinct(human_id)) 


# Stacked
ggplot(arg_chart, aes(fill=arg, y=count, x=strain)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme_classic() +   
  scale_y_continuous(labels = comma) +  # Add comma separators for the y-axis
  ylab("Samples") + 
  xlab("E. coli Strain") + 
#  ggtitle("High Risk ARGs Detected Among E. coli Strains Shared Between Chickens and Children") + 
  guides(fill = guide_legend(title = NULL))  # Hide legend title

########## Virulence Factors ###############

vir <- read.delim("Data/vir.summary.tab")

vf_long <- gather(vir, vf, value, APECO1_3696:yjaa, factor_key=TRUE)

vf_long$value <- as.numeric(vf_long$value)
vf_long <- vf_long[!is.na(vf_long$value), ]

vf_long$human_id <- gsub("\\.tab$", "", vf_long$X.FILE)
vf_long$ID <- vf_long$human_id

####### PATHOTYPE VF COMBINATIONS - JESSER & LEVY 

EAEC_aggR <- vf_long$ID[vf_long$vf=="aggR"]
EAEC_aaiC <- vf_long$ID[vf_long$vf=="aaiC"]

### Pathotype Combinations - 
EAEC_agg <- vf_long$ID[grepl("^agg", vf_long$vf)]
EAEC_aai <- vf_long$ID[grepl("^aai", vf_long$vf)]
EAEC_aat <- vf_long$ID[grepl("^aat", vf_long$vf)]


EPEC_eae <- vf_long$ID[grepl("^eae", vf_long$vf)]
EPEC_bfp <- vf_long$ID[grepl("^bfp", vf_long$vf)]

ETEC_elt <- vf_long$ID[grepl("^elt", vf_long$vf)]
ETEC_est <- vf_long$ID[grepl("^est", vf_long$vf)]

## NO ETEC VFs DETECTED 

EHEC_stx <- vf_long$ID[grepl("^stx", vf_long$vf)]

## NO EHEC VF DETECTED 

DAEC_afa<- vf_long$ID[grepl("^afa", vf_long$vf)]
DAEC_drg <- vf_long$ID[grepl("^drg", vf_long$vf)]

### 1/2 DAEC VFs not DETECTED 

EIEC_ipa <- vf_long$ID[grepl("^ipa", vf_long$vf)]


vf_long$EAEC
vf_long$EAEC[(vf_long$ID %in% EAEC_agg) & (vf_long$ID %in% EAEC_aai)] <- 1 
vf_long$EAEC[!(vf_long$ID %in% EAEC_agg) | !(vf_long$ID %in% EAEC_aai)] <- 0


vf_long$EPEC
vf_long$EPEC[(vf_long$ID %in% EPEC_eae) & (vf_long$ID %in% EPEC_bfp)] <- 1 
vf_long$EPEC[!(vf_long$ID %in% EPEC_eae) | !(vf_long$ID %in% EPEC_bfp)] <- 0

merged <- merge(shared_focal, vf_long, by="human_id", all.x=TRUE) 


sd(merged$NUM_FOUND[df_unique$threshold==1])


merged$shared <- NA
merged$shared[merged$threshold==0] <- "Not Shared"
merged$shared[merged$threshold==1] <- "Shared"


ggplot(merged, aes(x=shared, y=NUM_FOUND)) + 
  geom_boxplot() +
  
  xlab("Strain Shared with Any Chicken Sample") +
  theme_classic() +  scale_y_continuous(labels = scales::comma) +
  
  theme(legend.title = element_blank(), text = element_text(size = 16)) + 
  ylab("Number of Virulence Genes Detected") + 
  ggtitle("Virulence Genes Detected Among Strains Shared Between Chicken and Children") #+ 
# facet_wrap(~study2) 

strains <- c("13C1079T", "NCTC9087", "Not Shared")

merged$strain[merged$shared=="Not Shared"] <- "Not Shared"

merged2 <- subset(merged, merged$strain %in% strains)

#shared_vf_path <- subset(merged2, merged2$threshold==1)

merged2$Pathotype <- NA
merged2$Pathotype[merged2$EPEC==1] <- "EPEC"
merged2$Pathotype[merged2$EAEC==1] <- "EAEC"
shared_vf_path <- subset(merged2, !is.na(shared_vf_path$Pathotype))


patho_chart <- shared_vf_path %>%
  group_by(Pathotype, strain, threshold) %>% 
  summarise(count = n_distinct(human_id)) 


# Stacked
library(scales)
ggplot(patho_chart, aes(fill=Pathotype, y=count, x=strain)) + 
  geom_bar(position="dodge", stat="identity") + theme_classic() +   scale_y_continuous(labels = comma) + # Add comma separators for the y-axis
  ylab("Samples") + xlab("E. coli Strain") + ggtitle("Pathotypes Detected Among E. coli Strains Shared Between Chickens and Children") 
 # ylim(0, 15)

####################
### Validating MLST 
######################
mlst <- read.csv("Data/strainge_refdb_mlst_phylotype_2024-06-20.csv")
#clean strain names
mlst$strain <- gsub("^Esch_coli_", "", mlst$symlink)
mlst$strain <- gsub(".fa.gz", "", mlst$strain)

shared_focal$strain <-shared_focal$ref
## Merge MLST & Phylogroup with sample data
all_strains <- merge(mlst, shared_focal, by="strain", all.y=TRUE)

mlst_focal <- read.csv("Data/mlst.csv")

mlst_focal$focal_id <- gsub(".contigs.fasta", "", mlst_focal$ID)


mlst_focal2 <- merge(mlst_focal, all_strains, by="focal_id", all.y=TRUE)

mlst_focal2$mlst_valid <- NA
mlst_focal2$mlst_valid[mlst_focal2$mlst_focal==mlst_focal2$MLST] <- 1
mlst_focal2$mlst_valid[mlst_focal2$mlst_focal!=mlst_focal2$MLST] <- 0
mlst_focal3 <- subset(mlst_focal2, mlst_focal2$threshold==1)

vars <- c("focal_id", "mlst_focal", "MLST", "strain", "threshold", "mlst_valid")
mlst2 <- mlst_focal3[vars]


unique_mlst <- mlst2 %>%
  distinct()
      
####################
### Validating MLST vs. StrainGST Results
######################

mlst_focal <- read.csv("Data/mlst.csv")
phylo_focal <- read.table("Data/clermont_results.txt", header = FALSE, sep = "\t")

phylo_focal$ID <- gsub(".contigs", "", phylo_focal$V1) 
mlst_focal$ID <- gsub(".contigs.fasta", "", mlst_focal$ID)

phylo_focal$phylotype <-  phylo_focal$V2

mlst_phylo_focal <- merge(mlst_focal, phylo_focal, by="ID")
mlst_phylo <- merge(mlst_phylo_focal, all2, by="ID", all.y=TRUE)

# Assuming your dataset is named df
unique_mp <- mlst_phylo %>%
  distinct()

unique_mp$phylo_valid<- NA
unique_mp$phylo_valid[unique_mp$phylotype==unique_mp$phylo_post] <- 1
unique_mp$phylo_valid[unique_mp$phylotype!=unique_mp$phylo_post] <- 0

invalid_phylo <- subset(unique_mp, unique_mp$phylo_valid==0)

unique_mp$phylo_pre_valid<- NA
unique_mp$phylo_pre_valid[unique_mp$phylotype==unique_mp$phylo_pre] <- 1
unique_mp$phylo_pre_valid[unique_mp$phylotype!=unique_mp$phylo_pre] <- 0

invalid_phylo_pre <- subset(unique_mp, unique_mp$phylo_pre_valid==0)


shared_ids <- c("31271121", "31331121", "31581121", "33271121", "33431121", "33481121")
unique_mp$shared[unique_mp$ID %in% shared_ids] <- 1
unique_mp$shared[!(unique_mp$ID %in% shared_ids)] <- 0

mlst_focal2$mlst_valid <- NA
mlst_focal2$mlst_valid[mlst_focal2$mlst_focal==mlst_focal2$MLST] <- 1
mlst_focal2$mlst_valid[mlst_focal2$mlst_focal!=mlst_focal2$MLST] <- 0

vars <- c("ID", "mlst_focal", "MLST", "strain_pre", "mlst_valid")
mlst2 <- mlst_focal2[vars]

# Assuming your dataset is named df
unique_mlst <- mlst2 %>%
  distinct()
