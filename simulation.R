library(expm)
library(rgl)
library(RcppEigen)

sphere_diff_core = function (x, y) {
  u1 = x
  u2 = y - c(t(c(y) %*% c(u1))) * u1
  u2 = u2 / sqrt(sum(u2^2))
  
  Q = outer(u1, u2) - outer(u2, u1)
  theta = acos(c(c(x) %*% c(y)))
  R = expm(theta * Q) # diag(1, length(x)) + sin(theta) * Q + (1 - cos(theta)) * (Q %*% Q)
  
  return (list("Q" = Q, "R" = R, "theta" = theta))
}

#' evaluates the position of the spherical polynomial trend at time t
#' 
#' @param t a time index or a grid of time index
#' @param Q an d by d by p array, where each d by d slice is a skew-symmetric matrix
#'          corresponding to the coefficients associated from the lower degree to the 
#'          highest degree
#' @param x reference point
#' 
sphere_trend = function (t, Q, x) {
  if (is.matrix(Q)) {
    Q = array(Q, dim = c(dim(Q), 1))
  } else if (length(dim(Q)) == 2) {
    Q = array(Q, dim = c(dim(Q), 1))
  } else if (length(dim(Q)) != 3) {
    stop("Q must be either an d by d matrix or an d by d by p array")
  }
  
  n = dim(Q)[1]
  p = dim(Q)[3]
  
  if (length(t) == 1) {
    r = 1  # length of output
    
    time_poly = c(t^c(0:(p - 1)))
    M = matrix(Q, n * n, p) %*% time_poly
    M = array(M, dim = c(n, n, 1))
  } else {
    r = length(t)
    
    time_poly = t(outer(t, 0:(p - 1), `^`))
    M = matrix(Q, n * n, p) %*% time_poly
    M = array(M, dim = c(n, n, r))
  }
  
  R = array(0, dim = dim(M))
  res = array(NA, dim = c(n, r))
  for (j in 1:r) {
    R[,,j] = expm(M[,,j])
    res[,j] = R[,,j] %*% x
  }

  return (t(res))
}

#' computes the exponential map of the sphere
#' 
#' @param x an (n by d) array or a vector on the tangent space
#' @param mu the reference point
#' 
Exp_sphere = function (x, mu) {
  # x: an (m by d) matrix or a vector on the tangent space
  # mu: an d-dimensional vector of the reference point
  # tangent space --> sphere
  
  if (is.matrix(x)) {
    
    if (sum(abs(x)) == 0) {
      return (matrix(mu, nrow = nrow(x), ncol = ncol(x), byrow = T))
    }
    
    x_norm = sqrt(rowSums(x^2))
    x_norm[which(x_norm == 0)] = 1
    std_x = x / x_norm
    res = outer(cos(x_norm), mu) + sin(x_norm) * std_x
    
    res = res / sqrt(rowSums(res^2)) # normalize again to avoid numerical instability
  } else {
    
    if (sum(abs(x)) == 0) {
      return (mu)
    }
    
    x_norm = sqrt(sum(x^2))
    res = cos(x_norm) * mu + sin(x_norm) * x / x_norm
    
    res = res / sqrt(sum(res^2))
  }
  
  return (res)
}

#' computes the logarithm map of the sphere 
#' 
#' @param x an (n by d) array or a vector on the sphere
#' @param mu the reference point
#' 
Log_sphere = function (x, mu) {
  # x: an (m by d) matrix or a vector on the sphere
  # mu: an d-dimensional vector of the reference point
  # sphere --> tangent space
  
  if (is.matrix(x)) {
    w = x - matrix(mu, nrow = nrow(x), ncol = ncol(x), byrow = T)
    Proj = w - outer(c(w %*% mu), mu)
    Proj = Proj / sqrt(rowSums(Proj^2))
    
    res = acos(c(x %*% mu)) * Proj
    
    if (any(rowSums(abs(w)) < 1e-7)) {
      res[which(rowSums(abs(w)) < 1e-7),] = 0
    }
    
    
  } else {
    w = x - mu
    Proj = w - mu * c(t(mu) %*% w)
    
    if (sum(abs(w)) < 1e-7) {
      res = rep(0, length(mu))
    } else {
      res = acos(c(x %*% mu)) * Proj / sqrt(sum(Proj^2))
    }
  }
  
  return (res)
}

#' samples truncated normal distribution
t_isonormal_sampler = function (n, p, norm_cut = Inf, sd = 1) {
  counter = 0
  res = array(NA, dim = c(n, p))
  while (counter < n) {
    n_sampler = n * 10
    draw = matrix(rnorm(n_sampler * p, sd = sd), nrow = n_sampler, ncol = p)
    idx_effective = which(sqrt(rowSums(draw^2)) < norm_cut)
    n_eff = length(idx_effective)
    if (n_eff == 0) {
      next
    }
    if (counter + n_eff < n) {
      res[(counter + 1):(counter + n_eff),] = draw[idx_effective,]
      counter = counter + n_eff
    } else {
      res[(counter + 1):n,] = draw[idx_effective[1:(n - counter)],]
      counter = counter + n_eff
    }
  }
  
  if (p > 1) {
    return (res)
  } else {
    return (c(res))
  } 
}

#' samples truncated normal distribution on the tangent space
#' 
#' @param n number of samples
#' @param L magnitude of truncation of the norm
#' @param mu reference point
#' 
trunc_normal_tangent = function (n, L, mu) {
  d = length(mu)
  bas = svd(mu, nu = d)$u[,-1]
  
  temp = t_isonormal_sampler(n, d - 1, norm_cut = L, sd = 1)
  res = bas %*% t(temp)
  
  return (res)
}

#' adds noise to sphere data
#' 
#' @param x a d-dimensional vector or n by d array of data
#' @param type noise distribution ("truncate_normal")
#' @param L norm of the truncated normal
#' 
noise_inject = function (x, type = "truncate_normal", L = NULL) {
  if (type == "truncate_normal") {
    if (is.null(L)) {
      stop("noise_inject: for truncate_normal type, L must be supplied")
    }
    if (is.vector(x)) {
      d = length(x)
      
      res = Exp_sphere(trunc_normal_tangent(1, L = L, mu = x) , mu = x)
    } else if (is.matrix(x)) {
      n = nrow(x)
      d = ncol(x)
      
      res = matrix(0, nrow = n, ncol = d)
      
      for (i in 1:n) {
        res[i,] = Exp_sphere(c(trunc_normal_tangent(1, L = L, mu = x[i,])) , mu = x[i,])
      }
    }
  } else {
    stop("noise_inject: noise type not supported")
  }
  
  return (res)
}

#' computes the gradient with respect to the skew-symmetric coefficients and the reference point
#' 
#' @param U tangent vectors on the spheres (n by d array or vector)
#' @param Q current polynomial coefficient estimate (d by d by p array or matrix)
#' @param x current reference point estimate
#' @param time specific time points associated with U (default is equally spaced grids)
#' 
skew_gradient = function (U, Q, x, time = NULL) {
  if (is.matrix(Q)) {
    Q = array(Q, dim = c(dim(Q), 1))
  } else if (length(dim(Q)) == 2) {
    Q = array(Q, dim = c(dim(Q), 1))
  } else if (length(dim(Q)) != 3) {
    stop("skew_gradient: Q must be either an d by d matrix or an d by d by p array")
  }
  
  if (is.matrix(U)) {
    n = nrow(U)
  } else if (is.vector(U)) {
    n = 1
  } else {
    stop("skew_gradient: U must be a matrix or a vector")
  }
  
  d = dim(Q)[1]
  p = dim(Q)[3]
  if (is.null(time)) {
    time = c(1:n) / n
  } else {
    if (length(time) != n) {
      stop("skew_gradient: time length must equal the first dimension of U")
    }
  }

  if (n == 1) {
    time_poly = time^(0:(p - 1))
    M = matrix(Q, d * d, p) %*% time_poly
    M = array(M, dim = c(d, d, 1))
  } else {
    time_poly = t(outer(time, 0:(p - 1), `^`))
    M = matrix(Q, d * d, p) %*% time_poly
    M = array(M, dim = c(d, d, n))
  }
  
  res_skew = array(0, dim = c(d, d, p))
  res_x = rep(0, d)
  
  for (i in 1:n) {
    res_x = res_x + c(expm(-M[,,i]) %*% U[i,])
    
    temp = expmFrechet(-M[,,i], outer(U[i,], x), expm = FALSE)$Lexpm
    temp = skewpart(temp)
    for (j in 1:p) {
      res_skew[,,j] = res_skew[,,j] + temp * (i / n)^(j - 1)
    }
  }
  

  return (list("grd_skew" = res_skew, "grd_x" = res_x))
  
}

#' estimate the spherical polynomial model via gradient descent
#' 
#' @param y data matrix (n by d)
#' @param p order of the polynomial (linear = 1)
#' @param Q0 initialization of the coefficients (d by d by (p + 1), including intercept term)
#' @param mu0 initialization of the reference point 
#' @param alpha learning rate
#' @param save_iter whether to save the gradient descent path (default is FALSE)
#' 
spt = function (y, p, Q0, mu0, alpha = 0.5, max.iter = 1000,
                tol = 1e-5, save_iter = FALSE, verbose = FALSE) {
  if (is.matrix(Q0)) {
    Q0 = array(Q0, dim = c(dim(Q0), 1))
  } else if (length(dim(Q0)) == 2) {
    Q0 = array(Q0, dim = c(dim(Q0), 1))
  } else if (length(dim(Q0)) != 3) {
    stop("spt: Q0 must be either an d by d matrix or an d by d by (p + 1) array")
  }
  
  # initialize
  if (save_iter) {
    res_skew = array(0, dim = c(dim(Q0), max.iter))
    res_skew[,,,1] = Q0
    res_mu = array(0, dim = c(max.iter, length(mu0)))
    res_mu[1,] = mu0
  } else {
    res_skew = Q0
    res_mu = mu0
  }
  Q = Q0
  mu = mu0
  n = dim(y)[1] ; d = dim(y)[2]
  time = c(1:n) / n
  U = array(0, dim = c(n, d))
  
  
  for (i in 2:max.iter) {
    if (verbose) {
      cat("iteration", i - 1, "\n")
    }
    
    trend = sphere_trend(time, Q, mu)
    
    for (j in 1:n) {
      U[j,] = Log_sphere(y[j,], trend[j,])
    }
    grad = skew_gradient(U, Q, mu, time)
    
    Q_star = Q - alpha * grad$grd_skew
    mu_star = Exp_sphere(mu + alpha * grad$grd_x, mu = mu)
    
    if (save_iter) {
      res_skew[,,,i] = Q_star
      res_mu[i,] = mu_star
    } else {
      res_skew = Q_star
      res_mu = mu_star
    }
    
    if (sqrt(mean((Q_star - Q)^2)) < tol) {
      if (save_iter) {
        res_skew = res_skew[,,,1:i]
        res_mu = res_mu[1:i,]
        cat("spt: early stopping triggered\n")
        break
      } else {
        Q = Q_star
        mu = mu_star
      }
    }
    
  }
  
  return (list("Q" = res_skew, "mu" = res_mu))
}

#' computes the Frechet mean on the sphere
#' 
#' @param x (n by q) array of data
#' 
mean_on_sphere = function (x, tau = 0.1, tol = 1e-8, max.iter = 1000, verbose = FALSE) {
  
  n = nrow(x)
  mu = x[sample(n, 1),]
  for (i in 1:max.iter) {
    grad = colMeans(Log_sphere(x, mu))
    mu_new = Exp_sphere(mu + tau * grad, mu)
    
    temp = c(x %*% mu_new / sqrt(rowSums(x^2)))
    if (any(temp > 1.01) || any(temp < -1.01)) {
      stop("mean_on_sphere: something must be wrong")
    } else {
      temp = pmin(pmax(temp, -1), 1)
    }
    loss = mean(acos(temp))
    if (i > 1 && (loss_old - loss < tol)) {
      mu = mu_new
      break
    }
    if (verbose) {
      cat("mean_on_sphere: iter", i, "; loss", round(loss, 4), "\n")
    }
    mu = mu_new
    loss_old = loss
  }
  
  return (mu)
}

#' add longitude and latitude grids to rgl sphere object
#' 
#' @param n_lat number of latitude lines
#' @param n_lon number of longitude lines
#' @param alpha degree of transparency
#'  
add_sphere_grid = function (n_lat = 9, n_lon = 18, col = "black", alpha = 0.3) {
  phi = seq(-pi / 2, pi / 2, length.out = n_lat)
  theta = seq(0, 2 * pi, length.out = n_lon)
  
  for (lat in phi) {
    x = cos(lat) * cos(theta)
    y = cos(lat) * sin(theta)
    z = sin(lat) * rep(1, n_lon)
    lines3d(x, y, z, col = col, alpha = alpha)
  }
  
  for (lon in theta) {
    x = cos(phi) * cos(lon)
    y = cos(phi) * sin(lon)
    z = sin(phi)
    lines3d(x, y, z, col = col, alpha = alpha)
  }
}








### Example drawing
# x = c(0, 0, 1)
# k = 500  # number of time points
# c = pi   # maximum time
# 
# Q1 = matrix(c(0, 0, 1, 0, 0, 0, -1, 0, 0), ncol = 3)
# res = sphere_trend(c(0:(c * k)) / k, Q1, x)
# 
# spheres3d(0, 0, 0, radius = 1, color = "gray", alpha = 0.6)
# points3d(head(res, 1), col = "lightblue", size = 10)
# points3d(tail(res, 1), col = "salmon", size = 10)
# cols = colorRampPalette(c("lightblue", "salmon"))(c * k + 1)
# points3d(res, col = cols, size = 5)
# 
# Q2 = matrix(c(0, 0, 0, 0, 0, 1, 0, -1, 0), ncol = 3)
# res = sphere_trend(c(0:(c * k)) / k, Q2, x)
# points3d(head(res, 1), col = "lightblue", size = 10)
# points3d(tail(res, 1), col = "darkgreen", size = 10)
# cols = colorRampPalette(c("lightblue", "darkgreen"))(c * k + 1)
# points3d(res, col = cols, size = 5)
# 
# res = sphere_trend(c(0:(c * k)) / k, array(c(Q1, Q2), dim = c(3, 3, 2)), x)
# points3d(head(res, 1), col = "lightblue", size = 10)
# points3d(tail(res, 1), col = "gold", size = 10)
# cols <- colorRampPalette(c("lightblue", "gold"))(c * k + 1)
# points3d(res, col = cols, size = 5)
# 
# add_sphere_grid(alpha = 0.8)
# close3d()















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
  
  if (length(dim(R)) == 2) {
    R = array(R, dim = c(dim(R), 1))
  }
  
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
  res = t(R_mat) %*% B_mat / 2   # k × n
  
  return (res)
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
        res[i,] = res[i,] + B[,j] * (i / n)^(j - 1)
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
  
  return (list("x" = x, "basis" = basis, "SPD" = SPD, "R" = R, "Q" = Q, "x0" = x0))
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






