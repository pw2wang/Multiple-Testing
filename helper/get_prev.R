# Function used to get the largest achievable p-value that is smaller to the realized p-value
# for a single test. Used in MCF method

# Input:
# pval_supp: the support of p-values for the test
# p_org: the realized p-value

# Output:
# If p_org is not the smallest achievable p-value on the support, return the largest p-value
# that is smaller than p_org on the support
# If p_org is the smallest achievable p=value on the support, return 0

get_prev <- function(pval_supp, p_org) {
  # Assuming the pval_supp is already sorted ascendingly
  idx <- which(pval_supp == p_org)[1]
  if (idx > 1) {
    return(pval_supp[idx - 1])
  } else {
    return(0)
  }
}
