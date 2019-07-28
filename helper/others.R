get.sign <- function(alpha,
                     stepf,
                     c.m.ls,
                     loc) {
  m <- length(stepf)
  sign <- (sum(unlist(lapply(stepf,
                             FUN = function(x)
                               x(c.m.ls[loc])/(1 - x(c.m.ls[loc]))))) <= alpha*m)
  return(sign)
}

find.cm <- function(alpha, stepf, c.m.ls) {
  left <- 1
  right <- length(c.m.ls)
  left.sign <- get.sign(alpha, stepf, c.m.ls, left)
  right.sign <- get.sign(alpha, stepf, c.m.ls, right)
  while (right - left > 1 && left.sign != right.sign) {
    center <- round((right + left)/2)
    center.sign <- get.sign(alpha, stepf, c.m.ls, center)
    if (center.sign != left.sign) {
      right <- center
      right.sign <- center.sign
    } else {
      left <- center
      left.sign <- center.sign
    }
  }
  
  return(c.m.ls[left])  
}