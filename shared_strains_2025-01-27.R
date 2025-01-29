library(tidyr)

ST155 <- read.csv("Data/13C1079T_2025-01-26.csv")
ST10 <- read.csv("Data/NCTC9087_2025-01-26.csv")
ST117 <- read.csv("Data/APEC_2025-01-26.csv")
ST216 <- read.csv("Data/ST216_2025-01-26.csv")

focal_ids <- read.csv("Data/focal_ids.csv")
length(focal_ids$Filename)
mapsan_ids <- read.csv("Data/mapsan_ids.csv")
length(mapsan_ids$Filename)
human_ids <- rbind(focal_ids, mapsan_ids)
length(human_ids$Filename)

ST155$strain <- "13C1079T"
ST10$strain <- "NCTC9087"
ST117$strain <- "APEC102026"
ST216$strain <- "81"

shared <- rbind(ST155, ST10, ST117, ST216)
length(shared$comp)
head(shared)
# Create a consistent ordering between sample1 and sample2
shared$sample_min <- pmin(shared$sample1, shared$sample2)  # The smaller value goes in sample_min
shared$sample_max <- pmax(shared$sample1, shared$sample2)  # The larger value goes in sample_max

shared$samples <- paste0(shared$sample_min, shared$sample_max)
length(unique(shared$samples))
#drop duplicate rows based on sample_min and sample_max
#shared <- all_strains[!duplicated(all_strains[, c("sample_min", "sample_max")]), ]

# Drop the helper columns if you don't need them
#shared <- shared[, !(names(shared) %in% c("sample_min", "sample_max"))]

## Get human ID numbers 
shared$human_id <- NA
shared$human_id[shared$sample1 %in% human_ids$Filename] <- shared$sample1[shared$sample1 %in% human_ids$Filename]
shared$human_id[shared$sample2 %in% human_ids$Filename] <- shared$sample2[shared$sample2 %in% human_ids$Filename]

## Get Human Study Name
shared$human[shared$comp=="ChEEP x MapSan"] <- "Human (Community)"
shared$human[shared$comp=="ChEEP x FOCAL"] <- "Human (Diarrhea)"
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

length(prod_ids)
11*(177+62)
296/2629

length(market_ids)
15*(177+62)
445/3585

length(proc_ids)
75*(177+62)
560/17925
shared$breed <- "Commercial"
shared$breed[(shared$sample1 %in% com_ids) |(shared$sample2 %in% com_ids) ] <- "Commercial"
shared$breed[(shared$sample1 %in% ind_ids) | (shared$sample2 %in% ind_ids)] <- "Local"

length(com_ids)
80*(177+62)
825/19120

length(ind_ids)
47*(177+62)
46/11233

N_focal <- 44 
N_mapsan <- 177 
N_cheep <- 101 

table(cheep_meta$Group)

17680/(44+177)
75*(62+177)
summary(shared$gapJaccardSim)
summary(shared$singleAgreePct)

shared$threshold <- NA
shared$threshold[(shared$gapJaccardSim>=0.97) & (shared$singleAgreePct>=0.995)] <- 1
shared$threshold[(shared$gapJaccardSim<0.97) | (shared$singleAgreePct<0.995)] <- 0
table(shared$threshold)

length(unique(shared$samples[(shared$threshold==1) & (shared$comp=="ChEEP x FOCAL")]))
length(unique(shared$samples[(shared$threshold==1) & (shared$comp=="ChEEP x MapSan")]))

length(unique(shared$samples[(shared$threshold==1) & (shared$breed=="Commercial")]))

length(unique(shared$samples[(shared$threshold==1) & (shared$breed=="Local")]))

length(unique(shared$samples[(shared$threshold==1) & (shared$group=="Production")]))

length(unique(shared$samples[(shared$threshold==1) & (shared$group=="Market")]))
length(unique(shared$samples[(shared$threshold==1) & (shared$group=="Processing")]))
table(shared$strain)
length(unique(shared$focal_id[(shared$threshold==1) & (shared$strain=="13C1079T")]))
23/62
55/177
18/62
8/24139
1147/24139
1149/24139
764/17925
818/3585
1740/17877
1204/19120
1081/5019
703/2629
length(unique(shared$samples[shared$threshold==0]))
545/6262
head(shared)

2285/24139

focal_ids$human_id <- focal_ids$Filename 
mapsan_ids$human_id <- mapsan_ids$Filename
  
  ###############################################################
  ######## Virulence Factors & AMR  in FOCAL Samples ############
  ###############################################################
shared_focal <- subset(shared, shared$focal==1)
head(shared_focal)
########## ARGs ###############
  library(ggplot2)
  amr <- read.delim("Data/amr.summary.tab")
  head(amr)
  
  amr_long <- gather(amr, arg, value, aac.3..IId:tet.D., factor_key=TRUE)
  amr_long$human_id <- gsub("\\.tab", "", amr_long$X.FILE)
  nrow(amr_long)
  
  amr_long$value <- as.numeric(amr_long$value)
  amr_long <- amr_long[!is.na(amr_long$value), ]
  nrow(amr_long)

  amr_shared <- merge(shared_focal, amr_long, by="human_id", all.x=TRUE)
  head(amr_shared)
  sd(amr$NUM_FOUND)
  sd(amr_shared$NUM_FOUND[amr_shared$threshold==1])
  
  # Remove periods from the 'variable' column
  amr_shared$arg <- gsub("\\.", "", amr_shared$arg)
  
  hr_args <- read.csv("Data/clinically_relevant_amr_zhang_2022.csv")
  head(hr_args)
  hr_args_q1 <- subset(hr_args, hr_args$Rank=="Q1")
  
  hr_args_q1$HR_ARG <- 1 

 # Remove hyphens and parentheses from the ARG name 
  hr_args_q1$arg  <- gsub("[-()']", "", hr_args_q1$ARO.Term)
 
hr_args_shared <- merge(amr_shared,hr_args_q1,  by="arg", all.x=TRUE)

table(hr_args_shared$ARO.Term)
length(unique(hr_args_shared$human_id[(hr_args_shared$threshold==1) & (hr_args_shared$HR_ARG==1)])) 
table(hr_args_shared$strain)
length(unique(hr_args_shared$human_id[(hr_args_shared$threshold==1) & (hr_args_shared$HR_ARG==1)& (hr_args_shared$strain=="NCTC9087")])) 

table(hr_args_shared$ARO.Term[(hr_args_shared$threshold==1) & (hr_args_shared$HR_ARG==1)& (hr_args_shared$strain=="13C1079T")])
table(hr_args_shared$ARO.Term[hr_args_shared$threshold==0])
  

strains <- c("13C1079T", "NCTC9087")
hr_args2 <- subset(hr_args_shared, hr_args_shared$strain %in% strains)

hr_args3 <- hr_args2[!is.na(hr_args2$ARG.Class), ]
hr_args4 <- subset(hr_args3, hr_args3$threshold==1)

hr_args4$Class <- hr_args4$ARG.Class

  arg_chart <- hr_args4 %>%
    group_by(arg, Class, strain) %>% 
    summarise(count = n()) 
  
  head(arg_chart)
  
  # Stacked
  library(scales)
  ggplot(arg_chart, aes(fill=Class, y=count, x=strain)) + 
    geom_bar(position="stack", stat="identity") + theme_classic() +   scale_y_continuous(labels = comma) + # Add comma separators for the y-axis
    ylab("Sample Pairs") + xlab("E. coli Strain") + ggtitle("HR-ARGs (Grouped by Class) \n Detected Among E. coli Strains Shared Between Chickens and Children")
  ########## Virulence Factors ###############
  
  vir <- read.delim("Data/vir.summary.tab")
  head(vir)
  
  table(vir$X.FILE)
  summary(vir$NUM_FOUND)
  sd(vir$NUM_FOUND)
  length(unique(vir$X.FILE))
  

  vf_long1 <- gather(vir, vf, value, APECO1_3696:yjaa, factor_key=TRUE)
  head(vf_long1)
  vf_long1$value <- as.numeric(vf_long1$value)

vf_long <- subset(vf_long1, vf_long1$value>=80)
length(unique(vf_long$vf))
head(vf_long)
vf_long
  table(vf_long$value)
  
  vf_long$human_id <- gsub("\\.tab$", "", vf_long$X.FILE)
  table(vf_long$human_id)
  length(unique(vf_long$human_id))
  vf_long$ID <- vf_long$human_id
  
  patho <- read.csv("Data/Pakbin2021_pathotype_vf.csv")
  
  patho_det <- subset(patho, tolower(patho$Virulence.Factor) %in% tolower(vf_long$vf))
  head(patho_det)
  
  patho_det$vf <- tolower(patho_det$Virulence.Factor)
  table(patho_det$vf)
  vf_path <- merge(vf_long, patho_det, by="vf", all.y=TRUE)
  head(vf_path)
  
  table(vf_path$Class[vf_path$shared=="Shared"])
  
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

  strains <- c("13C1079T", "NCTC9087")
  
merged2 <- subset(merged, merged$strain %in% strains)
shared_vf_path <- subset(merged2, merged2$shared=="Shared")
table(shared_vf_path$Pathotype)

head(shared_vf_path)
table(shared_vf_path$Class)

shared_vf_path2 <- subset(shared_vf_path, shared_vf_path$value>90)

patho_chart <- shared_vf_path2 %>%
  group_by(Pathotype, strain) %>% 
  summarise(count = n()) 

head(patho_chart)

# Stacked
library(scales)
ggplot(patho_chart, aes(fill=Pathotype, y=count, x=strain)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() +   scale_y_continuous(labels = comma) + # Add comma separators for the y-axis
  ylab("Sample Pairs") + xlab("E. coli Strain") + ggtitle("Virulence Factors (Grouped by Pathotype) \n Detected Among E. coli Strains Shared Between Chickens and Children")
