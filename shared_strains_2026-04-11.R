library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)

shared <- read.csv("Data/combined_df.csv")
focal_ids <- read.csv("Data/focal_ids.csv")
mapsan_ids <- read.csv("Data/mapsan_ids.csv")
cheep_ids <- read.csv("Data/cheep_ids.csv")

focal_ids <- focal_ids %>% 
  rename(Filename = ID, 
         Study = study) %>% 
  dplyr::select(-X)

human_ids <- rbind(focal_ids, mapsan_ids)

# Create a consistent ordering between sample1 and sample2
shared$sample_min <- pmin(shared$sample1, shared$sample2)  # The smaller value goes in sample_min
shared$sample_max <- pmax(shared$sample1, shared$sample2)  # The larger value goes in sample_max

shared$samples <- paste0(shared$sample_min, shared$sample_max)

#drop duplicate rows based on sample_min and sample_max
#shared <- all_strains[!duplicated(all_strains[, c("sample_min", "sample_max")]), ]

# Drop the helper columns if you don't need them
#shared <- shared[, !(names(shared) %in% c("sample_min", "sample_max"))]

## Get human ID numbers 
shared$human_id <- NA
shared$human_id[shared$sample1 %in% human_ids$Filename] <- shared$sample1[shared$sample1 %in% human_ids$Filename]
shared$human_id[shared$sample2 %in% human_ids$Filename] <- shared$sample2[shared$sample2 %in% human_ids$Filename]
length(unique(shared$human_id))
## Get Human Study Name
shared$comp
shared$human <- NA
shared$human[shared$comp=="Chicken x Child Metagenome"] <- "Human (Community)"
shared$human[shared$comp=="Chicken x Child Isolate"] <- "Human (Diarrhea)"
table(shared$human)

## ChEEP Metadata 
cheep_meta <- read.csv("Data/phy_meta_2024-04-07.csv")

market_ids <- cheep_meta$ID[cheep_meta$Group=="Market"]
prod_ids <- cheep_meta$ID[cheep_meta$Group=="Production"]
proc_ids <- cheep_meta$ID[cheep_meta$Group=="Processing"]
cheep_ids <- cheep_meta$ID

shared$chicken_id <- NA
shared$chicken_id[shared$sample1 %in% cheep_ids] <- shared$sample1[shared$sample1 %in% cheep_ids]
shared$chicken_id[shared$sample2 %in% cheep_ids] <- shared$sample2[shared$sample2 %in% cheep_ids]
table(shared$chicken_id)

fecal_ids <- cheep_meta$ID[cheep_meta$Feces==1]
carcass_ids <- cheep_meta$ID[cheep_meta$Carcass==1]
water_ids <- cheep_meta$ID[cheep_meta$Water==1]
com_ids <- cheep_meta$ID[cheep_meta$breed=="Commercial"]
ind_ids <- cheep_meta$ID[cheep_meta$breed=="Indigenous"]

shared$group <- NA
shared$group[(shared$sample1 %in% prod_ids) |(shared$sample2 %in% prod_ids) ] <- "Production"
shared$group[(shared$sample1 %in% market_ids) | (shared$sample2 %in% market_ids)] <- "Market"
shared$group[(shared$sample1 %in% proc_ids) | (shared$sample2 %in% proc_ids) ] <- "Processing"

shared$focal_id <- NA
shared$focal_id[shared$human_id %in% focal_ids$Filename]<- shared$human_id[shared$human_id %in% focal_ids$Filename]

shared$focal <- NA
shared$focal[shared$human_id %in% focal_ids$Filename]<- 1 
shared$focal[shared$human_id %in% mapsan_ids$Filename] <- 0 

shared$mapsan <- NA
mapsan_ids
shared$mapsan[shared$human_id %in% mapsan_ids$Filename]<- shared$human_id[shared$human_id %in% mapsan_ids$Filename]

shared$breed <- "Commercial"
shared$breed[(shared$sample1 %in% com_ids) |(shared$sample2 %in% com_ids) ] <- "Commercial"
shared$breed[(shared$sample1 %in% ind_ids) | (shared$sample2 %in% ind_ids)] <- "Local"

### Threshold values based on ACNI (and/or gap similarity)
shared$threshold <- NA
shared$threshold[#(shared$gapJaccardSim>=0.97) &
  (shared$singleAgreePct>=99.95)] <- 1
shared$threshold[#(shared$gapJaccardSim<0.97) |
  (shared$singleAgreePct<99.95)] <- 0

commonpct <- subset(shared, shared$threshold==1)

focal_ids$human_id <- focal_ids$Filename 
mapsan_ids$human_id <- mapsan_ids$Filename
  
  ###############################################################
  ######## Virulence Factors & AMR  in FOCAL Samples ############
  ###############################################################
shared_focal <- subset(shared, shared$focal==1)

all_focal <- merge(focal_ids, shared_focal, by="human_id", all.x=TRUE)

shared_focal <- all_focal %>%
  # Replace NA in threshold with 0
  mutate(threshold = if_else(is.na(threshold), 0, threshold)) %>%
  group_by(human_id) %>%
  arrange(desc(threshold), .by_group = TRUE) %>%
  # Get the first non-NA strain in the sorted group
  mutate(ref = coalesce(ref, na.omit(ref)[1])) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(human_id, ref, threshold)

length(shared_focal$human_id[shared_focal$threshold==1])

########## ARGs ###############
  amr <- read.delim("Data/amr.summary.tab")
   
  amr_num <- amr %>% dplyr::select(c("X.FILE", "NUM_FOUND"))
  amr_num$human_id <- gsub("\\.tab", "", amr_num$X.FILE)
  
  amr_shared <- merge(shared_focal, amr_num, by="human_id", all.x=TRUE)

  # Conduct two-sample t-test
  t_test_result <- t.test(NUM_FOUND ~ threshold, data = amr_shared)
  # View results
  t_test_result
  
  # Create long version of amr dataset
  amr_long <- gather(amr, arg, value, aac.3..IId:tet.D., factor_key=TRUE)
  amr_long$human_id <- gsub("\\.tab", "", amr_long$X.FILE)
  
  amr_long$value <- as.numeric(amr_long$value)
  amr_long <- amr_long[!is.na(amr_long$value), ]

  amr_shared <- merge(shared_focal, amr_long, by="human_id", all.x=TRUE)

  mean(amr_shared$NUM_FOUND[amr_shared$threshold==1])
  
  mean(amr_shared$NUM_FOUND[amr_shared$threshold==0])
  
  # Remove periods from the 'variable' column
  amr_shared$arg <- gsub("\\.", "", amr_shared$arg)
  
  hr_args <- read.csv("Data/clinically_relevant_amr_zhang_2022.csv")
  hr_args_q1 <- subset(hr_args, hr_args$Rank=="Q1")
  hr_args_q1$HR_ARG <- 1 
 # Remove hyphens and parentheses from the ARG name 
  hr_args_q1$arg  <- gsub("[-()']", "", hr_args_q1$ARO.Term)

hr_args_shared <- merge(amr_shared,hr_args_q1,  by="arg", all.x=TRUE)

table(hr_args_shared$ARO.Term)
length(unique(hr_args_shared$human_id[(hr_args_shared$threshold==1) & (hr_args_shared$HR_ARG==1)])) 
table(hr_args_shared$ref)
length(unique(hr_args_shared$human_id[(hr_args_shared$threshold==1) & (hr_args_shared$HR_ARG==1)& (hr_args_shared$ref=="NCTC9087")])) 

table(hr_args_shared$ARO.Term[(hr_args_shared$threshold==1) & (hr_args_shared$HR_ARG==1)& (hr_args_shared$ref=="13C1079T")])
table(hr_args_shared$ARO.Term[hr_args_shared$threshold==0])
  
hr_args2 <- hr_args_shared[!is.na(hr_args_shared$ARG.Class), ]
hr_args2$shared <- NA 
hr_args2$shared[hr_args2$threshold==0] <- "Not Shared"
hr_args2$shared[hr_args2$threshold==1] <- "Shared"
table(hr_args2$shared)

#strains <- c("13C1079T", "NCTC9087")
#hr_args3 <- subset(hr_args2, hr_args2$ref %in% strains)

length(unique(hr_args2$human_id[hr_args2$shared=="Shared"]))

hr_args2$Class <- hr_args2$ARG.Class

length(unique(hr_args2$human_id[hr_args3$ARO.Term=="tet(B)"]))
  arg_chart <- hr_args2%>%
    group_by(arg, Class, ref) %>% 
    summarise(count = n()) 
  
  # Stacked
  library(scales)
  ggplot(arg_chart[!is.na(arg_chart$ref), ], aes(fill=ref, y=count, x=Class)) + 
    geom_bar(position="stack", stat="identity") + theme_classic() +   ylim(0, 30) + labs(fill = "E. coli Strain") + 
    ylab("Samples") + xlab("Antimicrobial Resistance Class") + ggtitle("HR-ARGs (Grouped by Class) \n Detected Among E. coli Strains Shared Between Chickens and Children")
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


vf_long$EAEC[(vf_long$ID %in% EAEC_agg) & (vf_long$ID %in% EAEC_aai)] <- 1 
vf_long$EAEC[!(vf_long$ID %in% EAEC_agg) | !(vf_long$ID %in% EAEC_aai)] <- 0

  patho <- read.csv("Data/Pakbin2021_pathotype_vf.csv")
  
  patho_det <- subset(patho, tolower(patho$Virulence.Factor) %in% tolower(vf_long$vf))

  patho_det$vf <- tolower(patho_det$Virulence.Factor)
  vf_path <- merge(vf_long, patho_det, by="vf", all.y=TRUE)

  
  merged <- merge(shared_focal, vf_path, by="human_id") 

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

 # strains <- c("13C1079T", "NCTC9087")
  
#merged2 <- subset(merged, merged$ref %in% strains)
shared_vf_path <- subset(merged, merged$threshold==1)
table(shared_vf_path$Pathotype)


patho_chart <- shared_vf_path %>%
  group_by(Pathotype, ref) %>% 
  summarise(count = n()) 


# Stacked
library(scales)
ggplot(patho_chart, aes(fill=ref, y=count, x=Pathotype)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() +  labs(fill = "E. coli Strain") + 
  ylab("Samples") + xlab("Pathotype") + ggtitle("Virulence Factors (Grouped by Pathotype) \n Detected Among E. coli Strains Shared Between Chickens and Children") + 
  ylim(0, 15)
