
skew_2_rotation = function (a) {
  res = matrix(c(cos(a), -sin(a), sin(a), cos(a)), ncol = 2)
  
  return (res)
}

circle_perturb = function (x, e) {
  
  res = cos(abs(e)) * x + sin(abs(e)) * c(x[2], - x[1])
  
  return (res)
}

time_basis_core = function (grid, type = "poly", p = NULL, q = NULL, n = NULL) {
  if (type == "poly") {
    time_basis = outer(grid, 1:p, "^")
  } else if (type == "season") {
    time_basis = 2 * pi * outer(grid, 1:q, "*") / n
    time_basis = cbind(sin(time_basis), cos(time_basis))
  }
  
  return (time_basis)
}

#' generate synthetic data
#' @param n sample size
#' @param B coefficient matrix (d by p), or a vector
#' @param grid a vector of time points
#' @param mu0 a (d by 2) matrix or vector of starting point
#' @param sigma noise level
#' @param q number of Fourier terms if type == "season"
#' 
#' @return a list containing dta (d by 2 by n or 2 by n) data array and 
#'                           trend (d by 2 by n or 2 by n) noiseless data and
#'                           B the coefficient matrix
#' 
dta_gen = function (n, B, type = "poly", grid, mu0 = NULL, 
                    sigma = 0, q = NULL) {
  if (is.matrix(B)) {
    d = nrow(B)
    p = ncol(B)
    if (is.null(mu0)) {
      mu0 = matrix(rnorm(2 * d), nrow = d, ncol = 2)
      mu0 = mu0 / apply(mu0, 1, norm, "2")
    }
    B.matrix = TRUE
  } else {
    d = 1
    p = length(B)
    if (is.null(mu0)) {
      mu0 = rnorm(2)
      mu0 = mu0 / norm(mu0, "2")
    }
    B.matrix = FALSE
  }
  
  time_basis = time_basis_core(grid, type, p = p, q = q, n = n)
  
  if (B.matrix) {
    trend = array(NA, dim = c(d, 2, n))
    res = array(NA, dim = c(d, 2, n))
    
    for (j in 1:d) {
      temp = c(time_basis %*% c(B[j,]))
      for (i in 1:n) {
        trend[j,,i] = c(skew_2_rotation(temp[i]) %*% c(mu0[j,]))
        
        if (sigma > 0) {
          res[j,,i] = circle_perturb(trend[j,,i], rnorm(1, sd = sigma))
        }
        
      }
    }
  } else {
    trend = array(NA, dim = c(2, n))
    res = array(NA, dim = c(2, n))
    
    temp = c(time_basis %*% B)
    for (i in 1:n) {
      trend[,i] = c(skew_2_rotation(temp[i]) %*% mu0)
      
      if (sigma > 0) {
        res[j,,i] = circle_perturb(trend[j,,i], rnorm(1, sd = sigma))
      }
      
    }
  }
  

  return (list("dta" = res, "trend" = trend, "B" = B))

}

#' convert to angles mod 2pi
#' 
circ_2_angle = function (X) {
  if (is.vector(X)) {
    res = atan2(X[2], X[1])
    return (res)
  } else {
    res = atan2(X[,2], X[,1])
    return (res)
  }
}

Log_circ = function (m, x) {
  res = x - m - (m %*% t(m)) %*% (x - m)
  res = res / norm(res, "2")
  res = c(acos(t(m) %*% x)) * res
  
  return (res)
}

#' single equation gradient
#' 
#' @param X 2 by n data matrix
#' 
gradient_core = function (time_basis, b, mu0, X) {
  R = time_basis %*% b
  res = array(NA, dim = c(length(b), ncol(X)))
  Q = matrix(c(0, -1, 1, 0), ncol = 2)
  
  for (i in 1:ncol(X)) {
    g = c(skew_2_rotation(R[i]) %*% mu0)
    res[,i] = c(t(g) %*% Q %*% Log_circ(g, X[,i])) * time_basis[i,]
  }
  
  return (res)
}

MST = function (A, lambda) {
  model = svd(A)
  D = diag(pmax(model$d - lambda, 0))
  
  return (model$u %*% D %*% model$v)
}

prox_grad_core = function ()






