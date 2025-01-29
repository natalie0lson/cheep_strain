#########################################################
########## Set up Matrices for Permutation Testing ###### 
#########################################################

# Ensure sample1 and sample2 are always in a consistent order for grouping
data <- shared %>%
  rowwise() %>%
  mutate(
    pair_id = paste(sort(c(sample1, sample2)), collapse = "_")
  ) %>%
  ungroup()
# Subset to one row per pair with the desired threshold logic
result <- data %>%
  group_by(pair_id) %>%
  summarise(
    sample1 = first(sample1),  # Retain one representative row for sample1
    sample2 = first(sample2),  # Retain one representative row for sample2
    threshold = max(threshold), # Use `max()` to assign 1 if any row has 1, else 0
    human_id = first(human_id), 
    chicken_id=first(chicken_id)
    
  ) %>%
  ungroup() 

# View the result
print(result)

## Include ALL Human sample IDs 
shared_all_ids <- merge(result, focal_ids, by="human_id", all=TRUE)
head(shared_all_ids)
length(unique(shared_all_ids$human_id))

all_ids <- merge(shared_all_ids, mapsan_ids, by="human_id", all=TRUE)
length(unique(all_ids$human_id))

cheep_df <- data.frame(cheep_ids)
cheep_df$chicken_id <- cheep_df$cheep_ids
all_ids <- merge(all_ids, cheep_df, by="chicken_id", all=TRUE)

## Remove Unnecessary Variables
variables <- c("human_id", "chicken_id", "threshold")
matrix <- all_ids[variables]

matrix$threshold[is.na(matrix$threshold)] <- 0

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
matcom

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

################################
##### Focal vs. Mapsan #########
#################################

# Ensure no NA values in each matrix
matmap[is.na(matmap)] <- 0
matfoc[is.na(matfoc)] <- 0

# Ensure data frames are converted to numeric matrices
matmap <- as.matrix(matmap)
matfoc <- as.matrix(matfoc)

matfoc <- apply(matfoc, 2, as.numeric)
matmap <- apply(matmap, 2, as.numeric)

matmap[is.na(matmap)] <- 0
matfoc[is.na(matfoc)] <- 0

matfoc
# Apply the permutation test
p_value_map_vs_foc <- permutation.test(as.vector(matmap), as.vector(matfoc))
print(p_value_map_vs_foc)


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
