library(abind)
source("./sphere_generic_util.R")

skew_2_rotation = function (a) {
  res = matrix(c(cos(a), -sin(a), sin(a), cos(a)), ncol = 2)
  
  return (res)
}

circle_perturb = function (x, e) {
  
  res = cos(abs(e)) * x + sin(abs(e)) * sign(e) * c(x[2], - x[1])
  
  return (res)
}

time_basis_core = function (grid, type = "poly", p = NULL, n = NULL) {
  if (type == "poly") {
    time_basis = outer(grid, 1:p, "^")
  } else if (type == "season") {
    time_basis = 2 * pi * outer(grid, 1:p, "*") / n
    time_basis = sin(time_basis)      # Using only sines due to f(0) = 0
  }
  
  return (time_basis)
}

#' generate synthetic data
#' @param n sample size
#' @param B coefficient matrix (d by p), or a vector
#' @param grid a vector of time points
#' @param mu0 a (d by 2) matrix or vector of starting point
#' @param sigma noise level
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
  
  time_basis = time_basis_core(grid, type, p = p, n = n)
  
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
        res[,i] = circle_perturb(trend[,i], rnorm(1, sd = sigma))
      }
      
    }
  }
  

  return (list("dta" = res, "trend" = trend, "B" = B, "mu" = mu0))

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

#' single equation gradient
#' 
#' @param X 2 by n data matrix
#' 
#' @returns an (n by p) matrix of gradients
#' 
gradient_core = function (time_basis, b, mu0, X) {
  R = time_basis %*% b
  res = array(NA, dim = c(length(b), ncol(X)))
  Q = matrix(c(0, -1, 1, 0), ncol = 2)
  
  for (i in 1:ncol(X)) {
    g = c(skew_2_rotation(R[i]) %*% mu0)
    res[,i] = c(t(g) %*% Q %*% Log_sphere(X[,i], g)) * time_basis[i,]
    if (abs(c(X[,i] %*% g) + 1) < 1e-10) {
      res[,i] = 0 # handle anti-podal point breakdown
    }
  }
  
  return (res)
}

gradient_wrap = function (time_basis, B, mu, X) {
  if (is.matrix(B)) {
    B.matrix = TRUE
    d = nrow(B)
    p = ncol(B)
  } else {
    B.matrix = FALSE
    d = 1
    p = length(B)
  }
  
  if (B.matrix) {
    B_grad = array(NA, dim = c(d, p))
    for (j in 1:d) {
      temp = gradient_core(time_basis, b = B[j,], mu0 = mu[j,], X[j,,])
      temp = rowMeans(temp)
      B_grad[j,] = temp
    }
    
  } else {
    B_grad = rep(0, p)
    temp = gradient_core(time_basis, b = B, mu0 = mu, X = X)
    temp = rowMeans(temp)
    B_grad = temp
  }
  
  return (B_grad)
}

MST = function (A, lambda) {
  model = svd(A)
  D = diag(pmax(model$d - lambda, 0))
  
  return (model$u %*% D %*% t(model$v))
}

#' construct proximal gradient update by matrix soft thresholding
#' 
#' @param B a (d by p) matrix or a length-p vector of current parameters
#' @param alpha step size
#' @param mu a (d by 2) matrix or a length-2 vector of current estimated mu's
#' @param X a (d by 2 by n) or (2 by n) array consisting of data
#' @param time_basis an (n by p) matrix of time basis functions
#' @param lambda regularization parameter
#'
#' @returns a (d by p) matrix or a length-p vector of one proximal gradient update
#'
prox_grad_core = function (B, alpha, mu, X, time_basis, lambda) {
  if (is.matrix(B)) {
    B.matrix = TRUE
    d = nrow(B)
    p = ncol(B)
  } else {
    B.matrix = FALSE
    d = 1
    p = length(B)
  }
  
  if (B.matrix) {
    B_grad = array(NA, dim = c(d, p))
    for (j in 1:d) {
      temp = gradient_core(time_basis, b = B[j,], mu0 = mu[j,], X[j,,])
      temp = rowMeans(temp)
      B_grad[j,] = B[j,] - alpha * temp
    }
    
    res = MST(B_grad, alpha * lambda)
  } else {
    B_grad = rep(0, p)
    temp = gradient_core(time_basis, b = B, mu0 = mu, X = X)
    temp = rowMeans(temp)
    B_grad = B - alpha * temp
    
    # soft thresholding becomes some sort of l2-regularization for vectors
    res = max((1 - alpha * lambda / norm(B_grad, "2")), 0) * B_grad 
  }
  
  return (res)
}

#' runs the Frank-Wolfe algorithm to estimate B
#' 
#' @param X a (d by 2 by n) or (2 by n) array consisting of data
#' @param grid a vector of time series
#' @param type type of time basis ("poly" or "season") Default is "poly"
#' @param B_init initialization for B (0 by default)
#' @param mu_init initialization for mu's (x1 by default)
#' @param rule update rules for Frank-Wolfe 
#'             ("predefine" uses 2 / (t + 2)),
#'              "exact" uses exact line search, (not yet supported)
#'              "backtrack" uses backtracking) (not yet supported)
#' @param L hyperparameter for FW
#' @param p number of time basis
#' 
FW_model_B = function (X, grid, type = "poly", B_init = NULL, mu_init = NULL,
                       rule = "predefine", L = 10, max_iter = 500, p = 1,
                       iter_start = 0,
                       verbose = TRUE) {
  if (length(dim(X)) == 3) {
    d = dim(X)[1]
    n = dim(X)[3]
    
    if (dim(X)[2] != 2) {
      stop("incorrect X dimensions")
    }
    if (d == 1) {
      stop("first dimension of X is 1. Should be dropped")
    }
    
  } else if (length(dim(X)) == 2) {
    d = 1
    n = dim(X)[2]
    
    if (dim(X)[1] != 2) {
      stop("incorrect X dimensions")
    }
  } else {
    stop("incorrect X dimensions")
  }
  
  # Initializations
  time_basis = time_basis_core(grid, type, p = p, n = n)
  if (is.null(B_init)) {
    if (d > 1) {
      B = array(0, dim = c(d, p))
    } 
    if (d == 1) {
      B = rep(0, p)
    }
  } else {
    B = B_init
  }
  if (is.null(mu_init)) {
    if (d > 1) {
      mu = array(NA, dim = c(d, 2))
      for (j in 1:d) {
        mu[j,] = mean_on_sphere(t(X[j,,1:5]), max.iter = 500)
      }
    } 
    if (d == 1) {
      mu = X[,1]
    }
  } else {
    mu = mu_init
  }
  
  B_history = array(NA, dim = c(d, p, max_iter))

  for (z in 1:max_iter) {
    if (d > 1) {
      grad_mu = array(0, dim = c(d, 2))
    } else if (d == 1) {
      grad_mu = c(0, 0)
    }
    
    B_grad = gradient_wrap(time_basis, B, mu, X)
    temp = svd(-B_grad, nu = 1, nv = 1)
    uv = c(temp$u) %*% t(c(temp$v))
    
    alpha = 2 / (2 + z)
    
    B = (1 - alpha) * B + alpha * L * uv
    B_history[,,z] = B
    
    if ((z + iter_start) %% 10 == 0) {
      loss = 0
      fitted_c = fit_trend(B, grid, mu, type)
      if (d > 1) {
        for (j in 1:d) {
          loss = loss + loss_compute(X[j,,], fitted_c[j,,]) 
        }
      } 
      
      if (d == 1) {
        loss = loss_compute(X, fitted_c) 
      }
      if (verbose) {
        cat("iteration", z + iter_start, ": loss", round(loss, 4), "\n")
      }
    }
    
  }
  
  return (list("B" = B, "B_history" = B_history, "mu" = mu))
}

mu_updates = function (X, grid, B, mu, type = "poly", p = 1, alpha = 1e-3,
                       max_iter = 500, 
                       iter_start = 0,
                       verbose = TRUE) {
  
  if (length(dim(X)) == 3) {
    d = dim(X)[1]
    n = dim(X)[3]
    
    if (dim(X)[2] != 2) {
      stop("incorrect X dimensions")
    }
    if (d == 1) {
      stop("first dimension of X is 1. Should be dropped")
    }
    
  } else if (length(dim(X)) == 2) {
    d = 1
    n = dim(X)[2]
    
    if (dim(X)[1] != 2) {
      stop("incorrect X dimensions")
    }
  } else {
    stop("incorrect X dimensions")
  }
  
  mu_history = array(NA, dim = c(d, 2, max_iter))
  time_basis = time_basis_core(grid, type, p = p, n = n)
  
  # update mu 
  for (zz in 1:max_iter) {
    if (d > 1) {
      grad_mu = array(0, dim = c(d, 2))
    } else if (d == 1) {
      grad_mu = c(0, 0)
    }
    
    for (j in 1:d) {
      if (d > 1) {
        temp = c(time_basis %*% B[j,])
      } else {
        temp = c(time_basis %*% B)
      }
      
      # (compute Riemannian gradients)
      for (i in 1:n) {
        R = skew_2_rotation(-temp[i])
        if (d > 1) {
          rot_x = R %*% X[j,,i]
          grad_mu[j,] = grad_mu[j,] - 2 * Log_sphere(c(rot_x), c(mu[j,])) / n
        }
        if (d == 1) {
          rot_x = c(R %*% X[,i])
          grad_mu = grad_mu - 2 * Log_sphere(rot_x, mu) / n
        }
      }
      
      # (Exponential map)
      if (d > 1) {
        for (j in 1:d) {
          if (alpha == 0) {
            mu_history[j,,zz] = mu[j,]
          } else {
            mu[j,] = Exp_sphere(-alpha * grad_mu[j,], mu[j,])
            mu_history[j,,zz] = mu[j,]
          }
        }
      } else {
        mu = Exp_sphere(-alpha * grad_mu, mu)
        mu_history[1,,zz] = mu
      }
    }
    
    if ((zz + iter_start) %% 10 == 0) {
      loss = 0
      fitted_c = fit_trend(B, grid, mu, type)
      if (d > 1) {
        for (j in 1:d) {
          loss = loss + loss_compute(X[j,,], fitted_c[j,,]) 
        }
      } 
      
      if (d == 1) {
        loss = loss_compute(X, fitted_c) 
      }
      if (verbose) {
        cat("  (mu) iteration", zz + iter_start, ": loss", round(loss, 4), "\n")
      }
    }
  }
  
  
  return (list("mu" = mu, "mu_history" = mu_history))
}

  
FW_model_alternate = function (X, grid, type = "poly", B_init = NULL, mu_init = NULL,
                               rule = "predefine", L = 10, alpha = 1e-3,
                               p = 1, max_iter = 5, FW_iter = 100, mu_iter = 30,
                               verbose = TRUE) {
  B_history = NULL
  mu_history = NULL

  for (zz in 1:max_iter) {
    model_B = FW_model_B(X, grid, type, B_init, mu_init, rule, L, max_iter = FW_iter,
                         p = p, iter_start = ((zz - 1) * (FW_iter + mu_iter)), verbose)
    B_init = model_B$B
    B_history = abind(B_history, model_B$B_history, along = 3)
    
    if (zz == max_iter) {
      break
    }
    model_mu = mu_updates(X, grid, model_B$B, model_B$mu, type, p, alpha,
                          max_iter = mu_iter, 
                          iter_start = zz * FW_iter + (zz - 1) * mu_iter, verbose)
    mu_init = model_mu$mu
    mu_history = abind(mu_history, model_mu$mu_history, along = 3)
  }
  
  return (list("B" = model_B$B, "mu" = model_mu$mu, "B_history" = B_history,
               "mu_history" = mu_history))
}


#' estimate the multivariate torus trend via proximal gradient descent
#'
#' @param X a (d by 2 by n) or (2 by n) array consisting of data
#' @param grid a vector of time points
#' @param type type of time basis ("poly" or "season") Default is "poly"
#' @param B_init initialization for B (0 by default)
#' @param mu_init initialization for mu's (x1 by default)
#' @param alpha1 step size for updating B
#' @param alpha2 step size for updating mu's
#' @param lambda regularization parameter for the nuclear norm
#' @param max_iter maximum iteration of gradient steps (500 by default)
#' @param p number of polynomial terms (1 by fault)
#' @param q number of sine/cosine pairs in seasonal model (p = 2q). Default is 1.
#' 
#'
prox_model = function (X, grid, type = "poly", B_init = NULL, mu_init = NULL,
                       alpha1 = 0.1, alpha2 = 0.01, lambda = 0, max_iter = 500,
                       p = 1, q = 1, 
                       verbose = TRUE,
                       iter_start = 0) {
  
  # construct time basis functions and initializations
  if (length(dim(X)) == 3) {
    d = dim(X)[1]
    n = dim(X)[3]
    
    if (dim(X)[2] != 2) {
      stop("incorrect X dimensions")
    }
    if (d == 1) {
      stop("first dimension of X is 1. Should be dropped")
    }
    
  } else if (length(dim(X)) == 2) {
    d = 1
    n = dim(X)[2]
    
    if (dim(X)[1] != 2) {
      stop("incorrect X dimensions")
    }
  } else {
    stop("incorrect X dimensions")
  }
  
  time_basis = time_basis_core(grid, type, p = p, q = q, n = n)
  if (type == "season") {
    # should fix this in the future
    p = q
  }
  
  if (is.null(B_init)) {
    if (d > 1) {
      B = array(0, dim = c(d, p))
    } 
    if (d == 1) {
      B = rep(0, p)
    }
  } else {
    B = B_init
  }
  
  if (is.null(mu_init)) {
    if (d > 1) {
      mu = array(NA, dim = c(d, 2))
      for (j in 1:d) {
        mu[j,] = X[j,,1]
      }
    } 
    if (d == 1) {
      mu = X[,1]
    }
  } else {
    mu = mu_init
  }
  
  # Alternating update of gradients
  B_history = array(NA, dim = c(d, p, max_iter))
  mu_history = array(NA, dim = c(d, 2, max_iter))
  loss_history = rep(NA, max_iter)
  
  for (z in 1:max_iter) {
    if (d > 1) {
      grad_mu = array(0, dim = c(d, 2))
    } else if (d == 1) {
      grad_mu = c(0, 0)
    }
    
    # update B
    B = prox_grad_core(B, alpha1, mu, X, time_basis, lambda)
    B_history[,,z] = B
    
    # update mu 
    for (j in 1:d) {
      if (d > 1) {
        temp = c(time_basis %*% B[j,])
      } else {
        temp = c(time_basis %*% B)
      }
      
      # (compute Riemannian gradients)
      for (i in 1:n) {
        R = skew_2_rotation(-temp[i])
        if (d > 1) {
          rot_x = R %*% X[j,,i]
          grad_mu[j,] = grad_mu[j,] - 2 * Log_sphere(c(rot_x), c(mu[j,])) / n
        }
        if (d == 1) {
          rot_x = c(R %*% X[,i])
          grad_mu = grad_mu - 2 * Log_sphere(rot_x, mu) / n
        }
      }
      
      # (Exponential map)
      if (d > 1) {
        for (j in 1:d) {
          if (alpha2 == 0) {
            mu_history[j,,z] = mu[j,]
          } else {
            mu[j,] = Exp_sphere(-alpha2 * grad_mu[j,], mu[j,])
            mu_history[j,,z] = mu[j,]
          }
        }
      } else {
        mu = Exp_sphere(-alpha2 * grad_mu, mu)
        mu_history[1,,z] = mu
      }
    }
    
    if ((z + iter_start) %% 10 == 0) {
      loss = 0
      fitted_c = fit_trend(B, grid, mu, type)
      if (d > 1) {
        for (j in 1:d) {
          loss = loss + loss_compute(X[j,,], fitted_c[j,,]) 
        }
        loss = loss + lambda * sum(svd(B)$d)
      } 
      
      if (d == 1) {
        loss = loss_compute(X, fitted_c) + lambda * sum(svd(B)$d)
      }
      if (verbose) {
        cat("iteration", z, ": loss", round(loss, 4), "\n")
      }
    }
    
  }
    
  return (list("B" = B, "mu" = mu, "loss" = loss_history,
               "B_history" = B_history,
               "mu_history" = mu_history))  
}



#' wrapper for diminishing step size
#' apply decay1 and decay2 to alpha1 and alpha2 respectively every 50 iterations
#' 
prox_model_decay = function(X, grid, type = "poly", B_init = NULL, mu_init = NULL,
                            alpha1 = 0.1, alpha2 = 0.01, decay1 = 1, decay2 = 1,
                            lambda = 0, max_iter = 500,
                            p = 1, q = 1, verbose = TRUE) {
  
  outer_iter = max_iter %/% 50
  B_history = NULL
  mu_history = NULL
  loss = NULL
  
  for (z in 1:outer_iter) {
    model = prox_model(X = X, grid = grid, type = type, B_init = B_init, mu_init = mu_init,
                       alpha1 = alpha1, alpha2 = alpha2, lambda = lambda, max_iter = 50,
                       p = p , q = q, verbose = verbose)
    B_init = model$B
    mu_init = model$mu
    
    B_history = abind(B_history, model$B_history, along = 3)
    mu_history = abind(mu_history, model$mu_history, along = 3)
    loss = c(loss, model$loss)
    
    alpha1 = alpha1 * decay1
    alpha2 = alpha2 * decay2
    
  }
  
  return (list("B" = model$B, "mu" = model$mu, "loss" = loss, "B_history" = B_history,
               "mu_history" = mu_history))
}

prox_model_alternate = function (X, grid, type = "poly", B_init = NULL, mu_init = NULL,
                                 alpha1 = 0.1, alpha2 = 0.01, 
                                 max_iter = 10, 
                                 max_iter_B = 100,
                                 max_iter_mu = 10,
                                 lambda = 0, p = 1, q = 1,
                                 verbose = TRUE) {
  B_history = NULL
  mu_history = NULL
  loss = NULL
  
  for (zz in 1:max_iter) {
    model = prox_model(X = X, grid = grid, type = type, B_init = B_init, mu_init = mu_init,
                       alpha1 = alpha1, alpha2 = 0, lambda = lambda, 
                       max_iter = max_iter_B,
                       p = p, q = q,
                       verbose = verbose,
                       iter_start = (zz - 1) * (max_iter_B + max_iter_mu))
    B_init = model$B
    mu_init = model$mu
    
    B_history = abind(B_history, model$B_history, along = 3)
    mu_history = abind(mu_history, model$mu_history, along = 3)
    loss = c(loss, model$loss)
    
    model = prox_model(X = X, grid = grid, type = type, B_init = B_init, mu_init = mu_init,
                       alpha1 = 0, alpha2 = alpha2, lambda = lambda, 
                       max_iter = max_iter_mu,
                       p = p, q = q,
                       verbose = verbose)
    B_init = model$B
    mu_init = model$mu
    
    B_history = abind(B_history, model$B_history, along = 3)
    mu_history = abind(mu_history, model$mu_history, along = 3)
    loss = c(loss, model$loss)
  }
  
  return (list("B" = model$B, "mu" = model$mu, "loss" = loss, "B_history" = B_history,
               "mu_history" = mu_history))
}

fit_trend = function (B, grid, mu, type = "poly") {
  if (is.matrix(B)) {
    B.matrix = TRUE
    d = dim(B)[1]
    p = dim(B)[2]
  } else {
    B.matrix = FALSE
    d = 1
    p = length(B)
  }
  
  n = length(grid)
  
  # if (type == "season") {
  #   q = p / 2
  # }
  
  time_basis = time_basis_core(grid, type, p = p, n = n)
  res = array(NA, dim = c(d, 2, n))
  
  for (j in 1:d) {
    if (d > 1) {
      temp = c(time_basis %*% B[j,])
    } else {
      temp = c(time_basis %*% B)
    }
    
    for (i in 1:n) {
      R = skew_2_rotation(temp[i])
      if (d > 1) {
        res_temp = R %*% mu[j,]
        res[j,,i] = R %*% mu[j,]
      }
      if (d == 1) {
        res_temp = R %*% mu
        res[1,,i] = R %*% mu
      }
    }
  }
  
  
  return (res)
}

#' compute loss
#' 
#' @param X (d by 2 by n) or (2 by n)
#' @param Y (d by 2 by n) or (2 by n)
#' 
loss_compute = function (X, Y) {
  res = 0
  if (length(dim(X)) == 3) {
    d = dim(X)[1]
    n = dim(X)[3]
    for (j in 1:d) {
      for (i in 1:n) {
        res = res + geod_sphere(c(X[j,,i]), c(Y[j,,i]))^2
      }
    }
  } 
  if (length(dim(X)) == 2) {
    n = dim(X)[2]
    for (i in 1:n) {
      res = res + geod_sphere(c(X[,i]), c(Y[,i]))^2
    }
  }
  
  res = res / n
  
  return (res)
}



