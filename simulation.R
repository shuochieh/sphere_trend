library(expm)
library(rgl)
library(RcppEigen)

#' create basis for the set of skew-symmetric matrices in R^d
#' 
SK_basis = function (d) {
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

check_type = function(x) {
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

SK_to_coordinate <- function(R) {
  # R: d × d × n
  d = dim(R)[1]
  n = dim(R)[3]
  
  # Basis: d × d × (d * (d - 1) / 2)
  B = SK_basis(d) 
  k = dim(B)[3]  
  
  # Flatten to (d^2 × n) and (d^2 × k)
  R_mat = matrix(R, nrow = d^2, ncol = n)
  B_mat = matrix(B, nrow = d^2, ncol = k)
  
  # Result: (n × k)
  # Crossprod trick: (B' R) but R is d^2 × n, B is d^2 × k
  res = crossprod(B_mat, R_mat)   # k × n
  
  return (t(res))
}

#' fast linear regression with multivariate Y
#'
#' @param X n by k
#' @param Y n by d
fastLmMulti = function (X, Y) {
  d = ncol(Y)
  k = ncol(X)
  coefs = matrix(NA, nrow = d, ncol = k)
  
  for (j in 1:d) {
    fit = fastLm(X, Y[,j])
    coefs[j,] = coef(fit)
  }
  
  return (coefs)
}

#' computes spherical differences as skew-symmetric matrix
#' 
SPD_SK = function (x1, x2) {
  u1 = x1
  u2 = x2 - c(t(u1) %*% x2) * u1
  u2 = u2 / sqrt(sum(u2^2))
  
  theta = acos(c(t(x1) %*% x2))
  
  return (theta * (outer(u1, u2) - outer(u2, u1)))
}

#' generates a sequence of spherical differences
#'
#' @param type specifies type of the stationary component
#'             (0: iid; 1:MA(1))
#' @param d dimension of the ambient Euclidean space
#' @param B coefficients in the polynomial trends
#' 
SPD_gen = function (n, d, B, sigma = 0.1, type = 1) {
  # generate stationary components
  stationary_components = array(NA, dim = c(n, d * (d - 1) / 2))
  if (type == 0) {
    for (i in 1:n) {
      stationary_components[i,] = rnorm(d * (d - 1) / 2, sd = sigma)
    }
  } else if (type == 1) {
    for (i in 1:n) {
      if (i == 1) {
        stationary_components[i,] = rnorm(d * (d - 1) / 2, sd = sigma)
      } else {
        stationary_components[i,] = rnorm(d * (d - 1) / 2, sd = sigma) - stationary_components[i - 1,]
      }
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
    
    res = array(0, dim = c(n, d * (d - 1) / 2))
    
    for (i in 1:n) {
      res[i,] = B + stationary_components[i,]
    }
    
  } else if (B_type == "matrix") {
    # Polynomial trend case
    p = ncol(B)
    
    res = array(0, dim = c(n, d * (d - 1) / 2))
    
    for (i in 1:n) {
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

#' generates sphere-valued data with spherical polynomial trend
#' 
#' @param x0 initial condition (a point on the sphere)
#' 
SPT_gen = function (n, d, B, x0, sigma = 0.1, type = 1) {
  temp = SPD_gen(n, d, B, sigma, type)
  SPD = temp$SPD
  
  basis = SK_basis(d)
  
  if (dim(basis)[3] != ncol(SPD)) {
    stop("SPT_gen: number of basis != SPD dimension")
  }
  
  x = array(NA, dim = c(n, d))
  
  # turn SPD coordinates into skew-symmetric matrices
  basis_flat = matrix(basis, nrow = d * d, ncol = d * (d - 1) / 2)
  Q_flat = basis_flat %*% t(SPD)
  Q = array(Q_flat, dim = c(d, d, n))
  
  # skew-symmetric matrices --> rotation matrices
  R = array(NA, dim = c(d, d, n))
  for (i in 1:n) {
    R[,,i] = expm(Q[,,i])
  }
  
  for (i in 1:n) {
    if (i == 1) {
      x[i,] = x0 # R[,,1] %*% x0
    } else {
      x[i,] = R[,,i] %*% x[i - 1,]
    }
  }
  
  return (list("x" = x, "basis" = basis, "R" = R, "Q" = Q, "x0" = x0))
}

#' estimates the SPT model
#' 
#' @param order order of the polynomial (must >= 1). (1: linear trend;...)
#' @param x n by d data
#' 
SPT_estimate = function (x, order) {
  n = nrow(x)
  d = ncol(x)

  SPD_in_SK = array(NA, dim = c(d, d, n - 1))
  for (i in 1:(n - 1)) {
    SPD_in_SK[,,i] = SPD_SK(x[i + 1,], x[i,])
  }
  
  SPD_in_basis = SK_to_coordinate(SPD_in_SK)
  
  X = array(NA, dim = c(n - 1, order))
  for (i in 1:order) {
    X[,i] = (c(2:n) / n)^(i - 1)
  }
  
  B_hat = fastLmMulti(X, SPD_in_basis)
  
  return (B_hat)
}






############## Example
set.seed(1)
close3d()

B = cbind(rnorm(3, sd = 0.05), rnorm(3, sd = 0.05))
dta = SPT_gen(100, 3, B[,1], x0 = c(0,0,1), 
              sigma = 0.000)

B_hat = SPT_estimate(dta$x, 1)


spheres3d(0, 0, 0, radius = 1, color = "lightblue", alpha = 0.3)

for (i in 1:1) {
  dta = SPT_gen(100, 3, B[,1], x0 = c(0,0,1), 
                sigma = 0.005)
  x = dta$x
  lines3d(x, col = "steelblue", lwd = 1.5)
  points3d(x, col = "lightblue", size = 5)
}

lines3d(x, col = "red", lwd = 2)
points3d(x, col = "red", size = 5)

dta = SPT_gen(100, 3, B, x0 = c(0,0,1), 
              sigma = 0.02)
y = dta$x
lines3d(y, col = "blue", lwd = 2)
points3d(y, col = "blue", size = 5)

axis_length <- 1.5  
segments3d(rbind(c(-axis_length, 0, 0), c(axis_length, 0, 0)), col = "black",   lwd = 3) # X-axis
segments3d(rbind(c(0, -axis_length, 0), c(0, axis_length, 0)), col = "black", lwd = 3) # Y-axis
segments3d(rbind(c(0, 0, -axis_length), c(0, 0, axis_length)), col = "black",  lwd = 3) # Z-axis

# rgl.postscript("sphere_trajectory.pdf", fmt = "pdf")






