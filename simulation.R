library(expm)

#' create basis for the set of skew-symmetric matrices in R^d
#' 
sk_sm_basis = function (d) {
  res = array(NA, dim = c(d, d, d * (d - 1) / 2))
  
  counter = 1
  for (i in 1:d) {
    for (j in 1:d) {
      if (i == j) {
        next
      } else if (i <= j) {
        next
      }
      
      a = rep(0, d)
      b = rep(0, d)
      a[i] = 1
      b[j] = 1
      res[,,counter] = outer(a, b) - outer(b, a)
      counter = counter + 1
    }
  }
  
  return (res)
}

check_type <- function(x) {
  if (is.matrix(x)) {
    return("matrix")
  } else if (is.atomic(x) && length(x) == 1) {
    return("scalar")
  } else if (is.atomic(x) && is.null(dim(x))) {
    return("vector")
  } else {
    return("other")
  }
}

#' generates a sequence of spherical differences
#'
#' @param type specifies type of the stationary component
#'             (0: iid)
#' @param d dimension of the ambient Euclidean space
#' @param B coefficients in the polynomial trends

SPD_gen = function (n, d, B, sigma = 0.1, type = 0) {
  # generate stationary components
  stationary_components = array(NA, dim = c(n - 1, d * (d - 1) / 2))
  if (type == 0) {
    for (i in 1:(n - 1)) {
      stationary_components[i,] = rnorm(d * (d - 1) / 2, sd = sigma)
    }
  }
  
  B_type = check_type(B)
  if (B_type == "scalar") {
    # No trend case
    
    if (B != 0) {
      stop("spherical difference cannot be a nonzero scalar")
    }
    res = stationary_components
    
  } else if (B_type == "vector") {
    # linear trend case
    
    res = array(0, dim = c(n - 1, d * (d - 1) / 2))
    
    for (i in 1:(n - 1)) {
      res[i,] = B + stationary_components[i,]
    }
    
  } else if (B_type == "matrix") {
    # Polynomial trend case
    p = ncol(B)
    
    res = array(0, dim = c(n - 1, d * (d - 1) / 2))
    
    for (i in 1:(n - 1)) {
      for (j in 1:p) {
        res[i,] = res[i,] + B[,j] * ((i + 1) / n)^(j - 1)
      }
      res[i,] = res[i,] + stationary_components[i,]
    }
  } else {
    stop ("SPD_gen: B is neither matrix, vector, nor scalar.")
  }
  
  return (list("SPD" = res, "stationary_comp" = stationary_components))
}










