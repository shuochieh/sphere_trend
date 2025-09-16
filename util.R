library(expm)
library(rgl)
library(RcppEigen)
library(deSolve)

#' computes spherical difference of x - y
#' 
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
#' @param bias whether a constant bias is added to the polynomial trend
#' 
sphere_trend = function (t, Q, x, bias = FALSE) {
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
    
    if (bias) {
      time_poly = c(t^c(0:(p - 1)))
    } else {
      time_poly = c(t^c(1:p))
    }
    M = matrix(Q, n * n, p) %*% time_poly
    M = array(M, dim = c(n, n, 1))
  } else {
    r = length(t)
    
    if (bias) {
      time_poly = t(outer(t, 0:(p - 1), `^`))
    } else {
      time_poly = t(outer(t, 1:p, `^`))
    }
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

#' evaluates the derivative of the spherical polynomial trend at time t
#' (Only internally used)
#' 
#' @param t a time index 
#' @param Q an d by d by p array, where each d by d slice is a skew-symmetric matrix
#'          corresponding to the coefficients associated from the lower degree to the 
#'          highest degree
#' @param x reference point
#' @param bias whether a constant bias is added to the polynomial trend
#' 
diff_sphere_trend = function (t, Q, x, bias = FALSE) {
  if (is.matrix(Q)) {
    Q = array(Q, dim = c(dim(Q), 1))
    if (bias) {
      return (rep(0, length(x)))
    }
  } else if (length(dim(Q)) == 2) {
    Q = array(Q, dim = c(dim(Q), 1))
    if (bias) {
      return (rep(0, length(x)))
    }
  } else if (length(dim(Q)) != 3) {
    stop("Q must be either an d by d matrix or an d by d by p array")
  }
  if (bias) {
    Q_nobias = Q[,,-1]
  } else {
    Q_nobias = Q
  }

  n = dim(Q_nobias)[1]
  p = dim(Q_nobias)[3]
  
  # the derivative of the polynomial part
  time_poly = c(1:p) * t^(c(0:(p - 1)))
  E = matrix(Q_nobias, n * n, p) %*% time_poly
  E = array(E, dim = c(n, n))
  
  # compute the polynomial part
  p = dim(Q)[3]
  if (bias) {
    time_poly = c(t^c(0:(p - 1)))
  } else {
    time_poly = c(t^c(1:p))
  }
  M = matrix(Q, n * n, p) %*% time_poly
  M = array(M, dim = c(n, n))
  
  # compute the derivative of the trend
  res = expmFrechet(M, E)$Lexpm 
  res = res %*% x

  return (c(res))
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

#' turns vector into skew-symmetric matrix
#' 
vec_2_skew = function (v) {
  p = length(v)
  d = (1 + sqrt(1 + 8 * p)) / 2
  
  res = matrix(0, ncol = d, nrow = d)
  res[upper.tri(res)] = v
  res = res - t(res)
  
  return (res)
}

#' helper function for parallel transport on the sphere
#' 
f_rhs = function (t, V, params) {
  Q = params$Q
  x = params$x
  bias = params$bias
  gamma = sphere_trend(t, Q, x, bias = bias)
  gamma_dot = diff_sphere_trend(t, Q, x, bias = bias)
  
  res = -c(t(c(V)) %*% c(gamma_dot)) * c(gamma)
  
  return (list(res))
}

#' parallel transport on the sphere
#' 
#' @param parms list(Q = Q, x = x, bias = bias)
#' 
pt_sphere = function (t0, t1, V0, parms) {
  times = seq(from = t0, to = t1, length.out = 50)
  para_V = ode(y = V0, times = times, func = f_rhs, parms = parms)
  
  return (para_V)
}

#' noisy trend data
#' 
#' @param t a time index or a grid of time index
#' @param Q a d by d by p array
#' @param x reference point
#' @param alpha AR coefficient
#' @param s standard deviation in the spherical AR model
#' @param bias 
#' 
noisy_trend = function (t, Q, x, alpha = 0.8, s = 0.1, bias = FALSE) {
  if (is.matrix(Q)) {
    Q = array(Q, dim = c(dim(Q), 1))
  } else if (length(dim(Q)) == 2) {
    Q = array(Q, dim = c(dim(Q), 1))
  } else if (length(dim(Q)) != 3) {
    stop("Q must be either a d by d matrix or a d by d by p array")
  }
  
  d = length(x)
  n = length(t)
  
  V0 = rnorm(d, sd = s)
  V0 = V0 - c((t(x) %*% V0)) * x
  
  V = matrix(NA, nrow = n, ncol = d)
  res = matrix(0, nrow = n, ncol = d)
  trend = matrix(0, nrow = n, ncol = d)

  for (j in 1:n) {
    # transport V0 from time 0 to time t
    time = t[j]
    if (j == 1) {
      last_time = 0
    } else {
      last_time = t[j - 1]
    }
    para_V = pt_sphere(last_time, time, V0, list(Q = Q, x = x, bias = bias))
    para_V = tail(para_V, 1)[-1]
    trend_value = c(sphere_trend(time, Q, x, bias))
    
    V0 = rnorm(d, sd = s)
    V0 = V0 - c(t(trend_value) %*% V0) * trend_value
    V0 = alpha * para_V + V0
    
    V[j,] = V0
    res[j,] = Exp_sphere(V0, trend_value)
    trend[j,] = trend_value
  }
  
  return (list("dta" = res, "V" = V, "trend" = trend))
}

#' computes the gradient with respect to the skew-symmetric coefficients and the reference point
#' 
#' @param U tangent vectors on the spheres (n by d array or vector)
#' @param Q current polynomial coefficient estimate (d by d by p array or matrix)
#' @param x current reference point estimate
#' @param time specific time points associated with U (default is equally spaced grids)
#' @param bias whether a constant term is added to the polynomial
#' 
skew_gradient = function (U, Q, x, time = NULL, bias = FALSE,
                          x_exact = FALSE, dta = NULL) {
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
  
  if (x_exact) {
    transformed_data = array(0, dim = c(n, d))
  }

  if (n == 1) {
    if (bias) {
      time_poly = time^(0:(p - 1))
    } else {
      time_poly = time^(1:p)
    }
    M = matrix(Q, d * d, p) %*% time_poly
    M = array(M, dim = c(d, d, 1))
  } else {
    if (bias) {
      time_poly = t(outer(time, 0:(p - 1), `^`))
    } else {
      time_poly = t(outer(time, 1:p, `^`))
    }
    M = matrix(Q, d * d, p) %*% time_poly
    M = array(M, dim = c(d, d, n))
  }
  
  res_skew = array(0, dim = c(d, d, p))
  res_x = rep(0, d)
  
  for (i in 1:n) {
    if (!x_exact) {
      res_x = res_x + c(expm(-M[,,i]) %*% U[i,])
    } else {
      transformed_data[i,] = expm(-M[,,i]) %*% dta[i,]
    }

    temp = expmFrechet(-M[,,i], outer(U[i,], x), expm = FALSE)$Lexpm
    temp = skewpart(temp)
    for (j in 1:p) {
      if (bias) {
        res_skew[,,j] = res_skew[,,j] + temp * (i / n)^(j - 1)
      } else {
        res_skew[,,j] = res_skew[,,j] + temp * (i / n)^(j)
      }
    }
  }
  
  if (x_exact) {
    res_x = mean_on_sphere(transformed_data)
  }

  return (list("grd_skew" = res_skew, "grd_x" = res_x))
  
}

#' estimate the spherical polynomial model via gradient descent
#' 
#' @param y data matrix (n by d)
#' @param p order of the polynomial (linear = 1)
#' @param Q0 initialization of the coefficients (d by d by (p + 1), including intercept term)
#' @param mu0 initialization of the reference point 
#' @param bias whether a constant term is in the polynomial
#' @param alpha_Q learning rate for the skew symmetric part
#' @param alpha_mu learning rate for the reference point
#' @param save_iter whether to save the gradient descent path (default is FALSE)
#' 
spt = function (y, p, Q0, mu0, bias = FALSE,
                alpha_Q = 0.5, alpha_mu = 0.01, max.iter = 1000, 
                x_exact = FALSE, tol = 1e-5, save_iter = FALSE, 
                verbose = FALSE) {
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
  
  loss = rep(NA, max.iter)
  
  for (i in 2:max.iter) {
    if (verbose) {
      cat("iteration", i - 1, ": ")
    }
    
    trend = sphere_trend(time, Q, mu)
    temp = sqrt(mean(geod_sphere(trend, y)^2))
    loss[i - 1] = temp
    if (verbose) {
      cat(round(temp, 4), "\n")
    }
    
    for (j in 1:n) {
      U[j,] = Log_sphere(y[j,], trend[j,])
    }
    grad = skew_gradient(U, Q, mu, time, x_exact = x_exact, dta = y, bias = bias)
    
    Q_star = Q + alpha_Q * grad$grd_skew
    if (x_exact) {
      mu_star = grad$grd_x
    } else {
      mu_star = Exp_sphere(alpha_mu * grad$grd_x, mu = mu) # gradient update
    }

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
        loss = loss[1:(i - 1)]
        cat("spt: early stopping triggered\n")
        break
      } else {
        cat("spt: early stopping triggered\n")
        break
      }
    } else {
      Q = Q_star
      mu = mu_star
    }
    
  }
  
  return (list("Q" = res_skew, "mu" = res_mu, "loss" = loss))
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
    mu_new = Exp_sphere(tau * grad, mu) # Exp_sphere(mu + tau * grad, mu)
    
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


#' computes the geodesic distance between two (array) of points on the spheres
#' 
#' If x and y are both arrays, a vector of n distances corresponding to pointwise
#' geodesic distances is returned.
#' 
#' @param x an $(n \times q)$ array of points or a $q$-dimensional vector
#' @param y an $(n \times q)$ array of points or a $q$-dimensional vector
geod_sphere = function (x, y) {
  if (is.vector(x) && is.vector(y)) {
    
    if (is_on_sphere(x)) {
      x = x / norm(x, "2")   # for numerical stability
    } else {
      stop("geod_sphere: x not on sphere")
    }
    
    if (is_on_sphere(y)) {
      y = y / norm(y, "2")   # for numerical stability
    } else {
      stop("geod_sphere: y not on sphere")
    }
    
    temp = c(x %*% y)
    res = acos(temp)
  } else if (is.matrix(x) && is.vector(y)) {
    
    if (prod(apply(x, 1, is_on_sphere)) == 1) {
      x = x / apply(x, 1, is_on_sphere)
    } else {
      stop("geod_sphere: x not on sphere")
    }
    
    if (is_on_sphere(y)) {
      y = y / norm(y, "2")   # for numerical stability
    } else {
      stop("geod_sphere: y not on sphere")
    }
    
    temp = c(x %*% y)
    res = acos(temp)
  } else if (is.vector(x) && is.matrix(y)) {
    
    if (is_on_sphere(x)) {
      x = x / norm(x, "2")   # for numerical stability
    } else {
      stop("geod_sphere: x not on sphere")
    }
    
    if (prod(apply(y, 1, is_on_sphere)) == 1) {
      y = y / apply(y, 1, is_on_sphere)
    } else {
      stop("geod_sphere: y not on sphere")
    }
    
    temp = c(y %*% x)
    res = acos(temp)
  } else if (is.matrix(x) && is.matrix(y)) {
    
    if (prod(apply(x, 1, is_on_sphere)) == 1) {
      x = x / apply(x, 1, is_on_sphere)
    } else {
      stop("geod_sphere: x not on sphere")
    }
    
    if (prod(apply(y, 1, is_on_sphere)) == 1) {
      y = y / apply(y, 1, is_on_sphere)
    } else {
      stop("geod_sphere: y not on sphere")
    }
    
    temp = c(diag(x %*% t(y)))
    res = acos(temp)
  } else {
    stop("geod_sphere: x and y must be 2-dim arrays or vectors.")
  }
  
  return (res)
}

#' check if x is on the sphere (within some tolerance)
is_on_sphere = function (x, tol = 1e-6) {
  if (abs(norm(x, "2") - 1) > tol) {
    return (FALSE)
  } else {
    return (TRUE)
  }
}

mzero = function (n) {
  return (matrix(0, nrow = n, ncol = n))
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


### Defunct codes
#' adds noise to sphere data
#' 
#' @param x a d-dimensional vector or n by d array of data
#' @param type noise distribution ("truncate_normal", "dependent")
#' @param L norm of the truncated normal
#' @param alpha AR coefficient if type = "dependent"
#' @param s standard deviation in the spherical AR model if type = "dependent"
#' 
# noise_inject = function (x, type = "truncate_normal", L = NULL, alpha = NULL,
#                          s = NULL) {
#   if (type == "truncate_normal") {
#     if (is.null(L)) {
#       stop("noise_inject: for truncate_normal type, L must be supplied")
#     }
#     if (is.vector(x)) {
#       res = Exp_sphere(c(trunc_normal_tangent(1, L = L, mu = x)) , mu = x)
#     } else if (is.matrix(x)) {
#       n = nrow(x)
#       d = ncol(x)
#       
#       res = matrix(0, nrow = n, ncol = d)
#       
#       for (i in 1:n) {
#         res[i,] = Exp_sphere(c(trunc_normal_tangent(1, L = L, mu = x[i,])) , mu = x[i,])
#       }
#     } else {
#       stop("noise_inject: input x must either vector or matrix")
#     }
#   } else if (type == "dependent") {
#     if (is.null(alpha)) {
#       stop("noise_inject: for dependent type, alpha must be supplied")
#     }
#     if (is.null(s)) {
#       stop("noise_inject: for dependent type, s must be supplied")
#     }
#     if (is.vector(x)) {
#       res = Exp_sphere(trunc_normal_tangent(1, L = s, mu = x), mu = x)
#     } else if (is.matrix(x)) {
#       n = nrow(x)
#       d = ncol(x)
#       
#       res = matrix(0, nrow = n, ncol = d)
#       
#       res[1,] = Exp_sphere(c(trunc_normal_tangent(1, L = s, mu = x[1,])), mu = x[1,])
#       
#       for (j in 2:n) {
#         v = Log_sphere(res[1,], )
#       }
#       
#       return (list("dta" = res, "Q" = Q_save))
#       
#     } else {
#       stop("noise_inject: input x must either vector or matrix")
#     }
#   } else {
#     stop("noise_inject: noise type not supported")
#   }
#   
#   return (res)
# }




