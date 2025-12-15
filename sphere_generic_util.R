library(Matrix)

#' check if x is on the sphere (within some tolerance)
is_on_sphere = function (x, tol = 1e-6) {
  if (abs(norm(x, "2") - 1) > tol) {
    return (FALSE)
  } else {
    return (TRUE)
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
    
    temp = round(c(x %*% y), 10) # resolve numerical issue (not ideal, please revisit in the future)
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

#' exponential map on the sphere
#' 
#' @param x an array of tangent vectors
#' @param mu base point 
Exp_sphere = function (x, mu) {
  
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

#' logarithmic map for the sphere
#' 
#' @param x an array of points on the sphere
#' @param mu a base point
Log_sphere = function (x, mu) {
  
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

#' parallel transport via geodesics for the sphere
#' 
#' @param x starting point
#' @param y end point
#' @param v tangent vector at x
#' 
pt_sphere = function (x, y, v) {
  
  if (abs(c(x %*% v)) > 1e-6) {
    stop("pt_sphere: v is not tangent at x")
  }
  
  e1 = x
  e2 = y - c(x %*% y) * x
  e2 = e2 / norm(e2, "2")
  
  v_perp = v - c(v %*% e2) * e2
  
  a = c(v %*% e2)
  theta = acos(c(y %*% x))
  
  res = a * (cos(theta) * e2 - sin(theta) * e1) + v_perp
  
  return(res)
}







