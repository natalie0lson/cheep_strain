library(tidyr)

ST155 <- read.csv("Data/13C1079T_2025-01-26.csv")
ST10 <- read.csv("Data/NCTC9087_2025-01-26.csv")
ST117 <- read.csv("Data/APEC_2025-01-26.csv")
ST216 <- read.csv("Data/ST216_2024-10-29.csv")

focal_ids <- read.csv("Data/focal_ids.csv")
mapsan_ids <- read.csv("Data/mapsan_ids.csv")
human_ids <- rbind(focal_ids, mapsan_ids)

ST155$strain <- "13C1079T"
ST10$strain <- "NCTC9087"
ST117$strain <- "APEC102026"
ST216$strain <- "81"

all_strains <- rbind(ST155, ST10, ST117, ST216)

all_strains$sample1 <-gsub('\\.filter', '', all_strains$sample1)
all_strains$sample2 <-gsub('\\.filter', '', all_strains$sample2)

# Create a consistent ordering between sample1 and sample2
all_strains$sample_min <- pmin(all_strains$sample1, all_strains$sample2)  # The smaller value goes in sample_min
all_strains$sample_max <- pmax(all_strains$sample1, all_strains$sample2)  # The larger value goes in sample_max

#drop duplicate rows based on sample_min and sample_max
shared <- all_strains[!duplicated(all_strains[, c("sample_min", "sample_max")]), ]

# Drop the helper columns if you don't need them
shared <- shared[, !(names(shared) %in% c("sample_min", "sample_max"))]

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
water_ids <- paste0(water_ids, ".filter")

com_ids <- cheep_meta$ID[cheep_meta$breed=="Commercial"]

ind_ids <- cheep_meta$ID[cheep_meta$breed=="Indigenous"]

shared$group <- NA
shared$group[(shared$sample1 %in% prod_ids) |(shared$sample2 %in% prod_ids) ] <- "Production"
shared$group[(shared$sample1 %in% market_ids) | (shared$sample2 %in% market_ids)] <- "Market"
shared$group[(shared$sample1 %in% proc_ids) | (shared$sample2 %in% proc_ids) ] <- "Processing"

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

focal_ids$human_id <- focal_ids$Filename 
mapsan_ids$human_id <- mapsan_ids$Filename

## Include ALL Human sample IDs 
shared_all_ids <- merge(shared, focal_ids, by="human_id", all=TRUE)
head(shared_all_ids)
length(unique(shared_all_ids$human_id))

all_ids <- merge(shared_all_ids, mapsan_ids, by="human_id", all=TRUE)
length(unique(all_ids$human_id))

cheep_df <- data.frame(cheep_ids)
cheep_df$chicken_id <- cheep_df$cheep_ids
all_ids <- merge(all_ids, cheep_df, by="chicken_id", all=TRUE)

length(unique(all_ids$chicken_id))

table(is.na(all_ids$Study.x))

all_ids$Study <-NA
all_ids$Study[is.na(all_ids$Study.x)] <- all_ids$Study.y[is.na(all_ids$Study.x)]
all_ids$Study[is.na(all_ids$Study.y)] <- all_ids$Study.x[is.na(all_ids$Study.y)]
table(all_ids$Study)
head(all_ids)
table(all_ids$strain)

all_155 <- subset(all_ids, all_ids$strain=="13C1079T")
all_10 <- subset(all_ids, all_ids$strain=="NCTC9087")
all_117 <- subset(all_ids, all_ids$strain=="APEC102026")
all_216 <- subset(all_ids, all_ids$strain=="ST216")

## Remove Unnecessary Variables
variables <- c("human_id", "chicken_id", "threshold")
matrix <- all_ids[variables]

matrix_155 <- all_155[variables]
matrix_10 <- all_10[variables]
matrix_117 <- all_117[variables]
matrix_216 <- all_216[variables]

matrix$threshold[is.na(matrix$threshold)] <- 0
matrix_155$threshold[is.na(matrix_155$threshold)] <- 0
matrix_10$threshold[is.na(matrix_10$threshold)] <- 0
matrix_117$threshold[is.na(matrix_117$threshold)] <- 0
matrix_216$threshold[is.na(matrix_216$threshold)] <- 0

## Split Matrices into Comparison Groups 
matrix_focal <- subset(matrix, matrix$human_id %in% focal_ids$Filename)
head(matrix_focal)

matrix_mapsan <- subset(matrix, matrix$human_id %in% mapsan_ids$Filename)

matrix_com <- subset(matrix, matrix$chicken_id %in% com_ids)
head(matrix_com)

matrix_loc <- subset(matrix, matrix$chicken_id %in% ind_ids)

matfoc <-  spread(matrix_focal, key = human_id, value = threshold)
matfoc[is.na(matfoc)] <- 0

matmap <-  spread(matrix_mapsan, key = human_id, value = threshold)
matmap[is.na(matmap)] <- 0 

matloc <-  spread(matrix_loc, key = human_id, value = threshold)
matloc[is.na(matloc)] <- 0
head(matloc)

matcom <-  spread(matrix_com, key = human_id, value = threshold)
matcom[is.na(matcom)] <- 0

########################################################
  ### Updated Function (from Daniel D. Kim)
########################################################

  permutation.test <- function(a, b) {
    set.seed(1992)
    
    # Ensure a and b are numeric and have no NA values
    if (any(is.na(a)) || any(is.na(b))) {
      stop("Input data contains NA values. Please clean your data.")
    }
    
    rep1 <- replicate(10000, sample(c(a, b), replace = FALSE))
    a1 <- rep1[1:length(a), ]
    b1 <- rep1[(length(a) + 1):(length(a) + length(b)), ]
    
    diffs <- apply(a1, 2, mean) - apply(b1, 2, mean)
    if (mean(a, na.rm = TRUE) == 0 & mean(b, na.rm = TRUE) == 0) {
      p <- NA
    } else {
      p <- sum(abs(diffs) >= abs(mean(a) - mean(b)))/1000
    }
    return(round(p, 40))
  }
  
################################################
  # Ensure no NA values in each matrix
  matcom[is.na(matcom)] <- 0
  matloc[is.na(matloc)] <- 0
  
  # Ensure data frames are converted to numeric matrices
  matcom <- as.matrix(matcom)
  matloc <- as.matrix(matloc)
  
  matloc <- apply(matloc, 2, as.numeric)
  matcom <- apply(matcom, 2, as.numeric)
  
  matcom[is.na(matcom)] <- 0
  matloc[is.na(matloc)] <- 0
  
  matloc
  # Apply the permutation test
  p_value_com_vs_loc <- permutation.test(as.vector(matcom), as.vector(matloc))
  print(p_value_com_vs_loc)
  
  ###############################################################
  ##### Permutation Function for multilevel variables ###########
  ###############################################################
  
  permutation.test.multi <- function(..., num_permutations = 10000, seed = 1992) {
    set.seed(seed)
    
    # Combine inputs into a list
    groups <- list(...)
    
    # Ensure all groups are numeric and have no NA values
    if (any(sapply(groups, function(x) any(is.na(x))))) {
      stop("Input data contains NA values. Please clean your data.")
    }
    
    # Flatten all groups into one vector and track group sizes
    all_data <- unlist(groups)
    group_sizes <- sapply(groups, length)
    num_groups <- length(groups)
    
    # Original group means
    group_means <- sapply(groups, mean)
    test_stat_original <- var(group_means) # Variance of means as test statistic
    
    # Permutation testing
    perm_stats <- replicate(num_permutations, {
      shuffled_data <- sample(all_data)
      perm_means <- sapply(1:num_groups, function(i) {
        start <- sum(group_sizes[1:(i-1)]) + 1
        end <- sum(group_sizes[1:i])
        mean(shuffled_data[start:end])
      })
      var(perm_means)
    })
    
    # Calculate p-value
    p_value <- mean(perm_stats >= test_stat_original)
    
    return(p_value)
  }
  
### Site-Level Comparison (Production, Market, Processing)
  
  ## Split Matrices into Comparison Groups 
  matrix_prod <- subset(matrix, matrix$chicken_id %in% prod_ids)
  head(matrix_prod)
  
  matrix_mkt <- subset(matrix, matrix$chicken_id %in% market_ids)
  
  matrix_proc <- subset(matrix, matrix$chicken_id %in% proc_ids)
  head(matrix_proc)
  
  matprod <-  spread(matrix_prod, key = human_id, value = threshold)
  matprod[is.na(matprod)] <- 0
  
  matmkt <-  spread(matrix_mkt, key = human_id, value = threshold)
  matmkt[is.na(matmkt)] <- 0
  
  matproc <-  spread(matrix_proc, key = human_id, value = threshold)
  matproc[is.na(matproc)] <- 0
  
  # Ensure no NA values in each matrix
  matprod[is.na(matprod)] <- 0
  matmkt[is.na(matmkt)] <- 0
  matproc[is.na(matproc)] <- 0
  
  # Ensure data frames are converted to numeric matrices
  matprod <- as.matrix(matprod)
  matmkt <- as.matrix(matmkt)
  matproc <- as.matrix(matproc)
  
  matprod <- apply(matprod, 2, as.numeric)
  matmkt <- apply(matmkt, 2, as.numeric)
  matproc <- apply(matproc, 2, as.numeric)
  
  # Ensure no NA values in each matrix
  matprod[is.na(matprod)] <- 0
  matmkt[is.na(matmkt)] <- 0
  matproc[is.na(matproc)] <- 0

  # Apply the permutation test
  p_value_sites <- permutation.test.multi(as.vector(matprod), as.vector(matmkt), as.vector(matproc))
  print(p_value_sites)
  
  ###########
  ### Strain Comparison
  mat155 <-  spread(matrix_155, key = human_id, value = threshold)
  matprod[is.na(mat155)] <- 0
  
  mat10 <-  spread(matrix_10, key = human_id, value = threshold)
  matmkt[is.na(mat10)] <- 0
  
  mat117 <-  spread(matrix_117, key = human_id, value = threshold)
  matproc[is.na(mat117)] <- 0
  
  mat216 <-  spread(matrix_216, key = human_id, value = threshold)
  matproc[is.na(mat216)] <- 0
  
  # Ensure no NA values in each matrix
  mat155[is.na(mat155)] <- 0
  mat10[is.na(mat10)] <- 0
  mat117[is.na(mat117)] <- 0
  mat216[is.na(mat216)] <- 0
  
  # Ensure data frames are converted to numeric matrices
  mat155 <- as.matrix(mat155)
  mat10 <- as.matrix(mat10)
  mat117 <- as.matrix(mat117)
  mat216 <- as.matrix(mat216)
  
  mat155 <- apply(mat155, 2, as.numeric)
  mat10 <- apply(mat10, 2, as.numeric)
  mat117 <- apply(mat117, 2, as.numeric)
  mat216 <- apply(mat216, 2, as.numeric)
  
  # Ensure no NA values in each matrix
  mat155[is.na(mat155)] <- 0
  mat10[is.na(mat10)] <- 0
  mat117[is.na(mat117)] <- 0
  mat216[is.na(mat216)] <- 0
  
  # Apply the permutation test
  p_value_strains <- permutation.test.multi(as.vector(mat155), as.vector(mat10), as.vector(mat117), as.vector(mat216))
  print(p_value_strains)
  
  ###############################################################
  ######## Virulence Factors & AMR  in FOCAL Samples ############
  ###############################################################
  library(ggplot2)
  amr <- read.delim("Data/amr.summary.tab")
  head(amr)
  
  vir <- read.delim("Data/vir.summary.tab")
  head(vir)
  
  table(vir$X.FILE)
  summary(vir$NUM_FOUND)
  
  # Remove ".tab" from the "ID" column
  amr$human_id <- gsub("\\.tab$", "", amr$X.FILE)
  table(amr$human_id)
  vir$human_id <- gsub("\\.tab$", "", vir$X.FILE)
  table(vir$human_id)
  length(unique(vir$human_id))
vir$ID <- vir$human_id
var <- c("ID")
focal_IDs <- vir[var]
focal_IDs$study <- "FOCAL"
focal_IDs
write.csv(focal_IDs, "Data/focal_ids.csv")
  
  head(matrix_focal)
  table(matrix_focal$threshold)
  
  merged <- merge(matrix_focal, vir, by="human_id") 
  head(merged)
  table(merged$NUM_FOUND[merged$threshold==1])
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
  
patho <- read.csv("Data/Pakbin2021_pathotype_vf.csv")
head(patho)

head(merged)
vf_long <- gather(merged, vf, value, APECO1_3696:yjaa, factor_key=TRUE)
head(vf_long)
table(vf_long$vf)

patho_det <- subset(patho, tolower(patho$Virulence.Factor) %in% tolower(vf_long$vf))
head(patho_det)

patho_det$vf <- tolower(patho_det$Virulence.Factor)
table(patho_det$vf)
vf_path <- merge(vf_long, patho_det, by="vf", all=TRUE)
head(vf_path)

table(vf_path$Class[vf_path$shared=="Shared"])

shared_vf_path <- subset(vf_path, vf_path$shared=="Shared")
head(shared_vf_path)
table(shared_vf_path$Class)


