# Function used to detect (make rejections) for MCF method

# Input:
# lambda_star: the thresholding p-value from Habiger's method
# p_prev: the vector of the largest p-value that is smaller than p_org on the support of each test
# p_org: the vector of realized p-value
# randp_ecdf: randomized p-values drawn between each pair of p_org and p_prev, used to calculate
#             the expected proportion of rejections in Habiger's method under lambda_star

# Output:
# Indices of detected (rejected) tests

mcf_detect <- function(lambda_star, p_prev, p_org, randp_ecdf) {
  cdf <- sum(randp_ecdf < lambda_star)/length(randp_ecdf)
  r <- unlist(lapply(1:length(p_org), FUN = function(x) (lambda_star > p_org[x])*1 + (lambda_star < p_prev[x])*0 + 
                       (lambda_star >= p_prev[x] && lambda_star <= p_org[x])*(lambda_star - p_prev[x])/(p_org[x] - p_prev[x])))
  detection <- which(r >= quantile(r, probs = 1 - cdf))
  return(detection)
}
