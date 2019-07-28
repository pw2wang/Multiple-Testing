# Function used to find the threshold p-value in Benjamini-Hochberg procedure

# Input:
# p_val: a vector of p-values
# alpha: nominal significance level

# Output:
# A single numeric value as the threshold p-value in BH procedure
find_BH_treshold <- function(p_val, alpha) {
  idx <- which(sort(p_val) <= alpha*(1:length(p_val))/length(p_val))
  return(sort(p_val)[idx[length(idx)]])
}