# Function used to find the threshold p-value in Gilbert's method
# Requires the find_BH_threshold function

# Input:
# p_val: the realized p-values for each test
# p_val_supp_min: the minimum achievable p-value for each test
# alpha: nominal significance level

# Output:
# A single numeric value as the threshold p-value in Gilbert's method

source(file = "Helper/find_BH_threshold.R")
find_Gilbert_threshold <- function(p_val, p_val_supp_min, alpha) {
  K <- which(unlist(lapply(1:m, FUN = function(x, supp_min, alpha) sum(supp_min < alpha/x) <= x, supp_min = p_val_supp_min, alpha = alpha)))[1]
  K_idx <- which(p_val_supp_min < alpha/K)
  thres <- find_BH_treshold(p_val[K_idx], alpha)
  return(thres)
}