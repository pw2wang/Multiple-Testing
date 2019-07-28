# Function used to calculate the statistical power for a LSMT result

# Input:
# detection: the indices for the detected non-null (alternative) hypotheses
# true_non_null: the indices for the true non-null (alternative) hypotheses

# Output:
# power of the detection result. # of truely detected non-nulls/ # of total true non-nulls

get_power <- function(detection, true_non_null) {
  if (length(detection) == 0) {
    return(0)
  } else {
    return(sum(detection %in% true_non_null)/length(true_non_null)) 
  }
}