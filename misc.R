### local constant regression for circular data

gaussian_kernel = function (x, y, h) {
  diff_mat = outer(x, y, "-")
  u = diff_mat / h
  exp(-0.5 * u^2) / (sqrt(2 * pi) * h)
}

#' local constant regression for circular data with Gaussian kernels
#' 
#' @param x a vector of prediction points
#' @param X a vector of observed time points (usually grid)
#' @param theta a vector of principal angles
#' @param h bandwidth
#' 
nw_gaussian_circ = function (x, X, theta, h) {
  m1_hat = gaussian_kernel(x, X, h) %*% sin(theta) / length(X)
  m2_hat = gaussian_kernel(x, X, h) %*% cos(theta) / length(X)
  
  return (atan2(c(m1_hat), c(m2_hat)))
}
