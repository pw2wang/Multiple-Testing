# Function used to calculate the false discovery rate (FDR) for a LSMT result

# Input:
# detection: the indices for the detected non-null (alternative) hypotheses
# true_null: the indices for the true null (alternative) hypotheses

# Output:
# FDR of the detection result. # of falsely detected nulls/ # of total true nulls

get_fdr <- function(detection, true_null) {
  if (length(detection) == 0) {
    return(0)
  } else {
    return(sum(detection %in% true_null)/length(detection))
  }
}