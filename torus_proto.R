source("torus_util.R")
source("misc.R")

low_rank_gen = function (m, n, r, s = 1) {
  A = matrix(rnorm(m * r), ncol = r)
  A = qr.Q(qr(A))
  B = matrix(rnorm(n * r), ncol = r)
  B = qr.Q(qr(B))
  
  if (r > 1) {
    D = diag(c(1:r)^(-1))
    return (A %*% D %*% t(B))
  } else {
    return (A %*% t(B))
  }
}

# Case 1 

set.seed(1)
n = 30
d = 25
q = 3
grid = c(1:n) 
lambda = 0.5

B = 8 * low_rank_gen(d, 2 * q, 1)

dta = dta_gen(n, B, "season", grid, sigma = 0.75, q = q)

### Option 1: Infeasible perturbed truth
# B_init = B + 1.0 * matrix(rnorm(d * 2 * q), nrow = d, ncol = 2 * q)
################

### Option 2: Initialization warm-up
model_unc = prox_model_alternate(dta$dta, grid, "season", 
                                 alpha1 = 1e-2, alpha2 = 1e-3,
                                 max_iter = 20,
                                 max_iter_B = 50,
                                 max_iter_mu = 20,
                                 mu_init = NULL,
                                 B_init = NULL,
                                 q = 2)
B_init = cbind(model_unc$B, matrix(0, nrow = d, ncol = (q - 2) * 2))

model_unc = prox_model_alternate(dta$dta, grid, "season",
                                 alpha1 = 1e-2, alpha2 = 1e-3,
                                 max_iter = 20,
                                 max_iter_B = 50,
                                 max_iter_mu = 20,
                                 mu_init = NULL,
                                 B_init = B_init,
                                 q = q)
fitted_unc = fit_trend(model_unc$B, grid, model_unc$mu, "season")


# model_unc = prox_model(dta$dta, grid, "season", 
#                    alpha1 = 1e-2, alpha2 = 1e-3, 
#                    max_iter = 1000,
#                    q = q, 
#                    lambda = 0.0, 
#                    mu_init = NULL,
#                    B_init = B_init)
# fitted_unc = fit_trend(model_unc$B, grid, model_unc$mu, "season")



lambda_grid = seq(from = 0.5, to = 1.5, length.out = 5)
temp = rep(NA, 5)
for (zz in 1:5) {
  model_c = prox_model_alternate(dta$dta, grid, "season",
                                 alpha1 = 1e-2, alpha2 = 1e-3,
                                 lambda = lambda_grid[zz],
                                 max_iter = 20,
                                 max_iter_B = 50,
                                 max_iter_mu = 20,
                                 mu_init = NULL,
                                 B_init = model_unc$B,
                                 q = q,
                                 verbose = FALSE)
  temp[zz] = sqrt(loss_compute(dta$trend, fitted_c) / d)
  cat("lambda search...", zz, "\n")
}
lambda_opt = lambda_grid[which.min(temp)]

model_c = prox_model_alternate(dta$dta, grid, "season",
                               alpha1 = 1e-2, alpha2 = 1e-3,
                               lambda = lambda_opt,
                               max_iter = 20,
                               max_iter_B = 50,
                               max_iter_mu = 20,
                               mu_init = NULL,
                               B_init = model_unc$B,
                               q = q)
fitted_c = fit_trend(model_c$B, grid, model_c$mu, "season")

# fitted_init = fit_trend(B_init, grid, dta$dta[,,1], "season")

h_grid = seq(from = 0.01, to = 0.1, length.out = 10)
oracle_loss = rep(0, 10)
for (z in 1:10) {
  h = h_grid[z]
  temp = array(NA, dim = c(d, 2, n))
  for (j in 1:d) {
    temp[j,1,] = cos(nw_gaussian_circ(grid/n, grid/n, circ_2_angle(t(dta$dta[j,,])), h = h))
    temp[j,2,] = sin(nw_gaussian_circ(grid/n, grid/n, circ_2_angle(t(dta$dta[j,,])), h = h))
  }
  oracle_loss[z] = sqrt(loss_compute(dta$trend, temp) / d)
}
h_opt = h_grid[which.min(oracle_loss)]


par(mfrow = c(5, 5), mar = c(2, 2, 1, 1))
for (j in 1:d) {
  plot(circ_2_angle(t(dta$dta[j,,])), ylab = "", xlab = "", ylim = c(-pi, pi), pch = 19,
       col = "gray70",
       main = paste("Series", j))
  lines(circ_2_angle(t(dta$trend[j,,])), col = "lightblue", lwd = 2)
  lines(circ_2_angle(t(fitted_unc[j,,])), col = "salmon2", lty = 2, lwd = 2.5)
  lines(circ_2_angle(t(fitted_c[j,,])), col = "darkolivegreen4", lty = 2, lwd = 2.5)
  lines(nw_gaussian_circ(grid/n, grid/n, circ_2_angle(t(dta$dta[j,,])), h = h_opt), 
        col = "purple", lty = 4, lwd = 2)
}

par(mfrow = c(1, 1))
plot(log((model_unc$loss)), type = "b", pch = 19, col = "salmon", 
     ylim = c(0, 5), main = "training error (log scale)")
lines(log((model_c$loss)), type = "b", pch = 19, col = "darkolivegreen4")


training_err_unc = array(0, dim = c(d, 1400))
for (iter in 1:1400) {
  B_hat = model_unc$B_history[,,iter]
  mu_hat = model_unc$mu_history[,,iter]
  fitted_unc = fit_trend(B_hat, grid, mu_hat, "season")
  for (j in 1:d) {
    training_err_unc[j,iter] = sqrt(loss_compute(dta$dta[j,,], fitted_unc[j,,]) / d)
  }
  
  cat("    iter", iter, "\n")
}

training_err_c = array(0, dim = c(d, 1400))
obj_val = rep(0, 1400) 
par(mfrow = c(5, 5), mar = c(2, 2, 1, 1))
for (iter in 1:1400) {
  B_hat = model_c$B_history[,,iter]
  mu_hat = model_c$mu_history[,,iter]
  fitted_c = fit_trend(B_hat, grid, mu_hat, "season")
  for (j in 1:d) {
    training_err_c[j,iter] = sqrt(loss_compute(dta$dta[j,,], fitted_c[j,,]) / d)
    obj_val[iter] = obj_val[iter] + training_err_c[j,iter]^2 
  }
  obj_val[iter] = obj_val[iter] + lambda * sum(svd(B_hat)$d)
  
  cat("    iter", iter, "\n")
}

par(mfrow = c(5, 5), mar = c(2, 2, 1, 1))
for (j in 1:d) {
  plot(training_err_unc[j,], type = "l", col = "salmon", ylim = c(0, max(training_err_c)),
       main = paste("Series", j),
       ylab = "")
  lines(training_err_c[j,], type = "l", col = "darkolivegreen4")
}

par(mfrow = c(1, 1))
plot(obj_val, type = "l", col = "darkolivegreen4",
       main = "Objective function for the regularized problem",
       ylab = "")

sqrt(loss_compute(dta$trend, fitted_unc) / d)
sqrt(loss_compute(dta$trend, fitted_c) / d)
sqrt(loss_compute(dta$trend, dta$dta) / d)
temp = array(NA, dim = c(d, 2, n))
for (j in 1:d) {
  temp[j,1,] = cos(nw_gaussian_circ(grid/n, grid/n, circ_2_angle(t(dta$dta[j,,])), h = h_opt))
  temp[j,2,] = sin(nw_gaussian_circ(grid/n, grid/n, circ_2_angle(t(dta$dta[j,,])), h = h_opt))
}
sqrt(loss_compute(dta$trend, temp) / d)

# Case 2

n = 100
d = 25
p = 10
grid = c(1:n) / (n)

B = 25 * low_rank_gen(d, p, 1)
dta = dta_gen(n, B, "poly", grid, sigma = 1.0, q = q)

model = prox_model(dta$dta, grid, "poly", 
                   alpha1 = 1e-2, alpha2 = 1e-3, 
                   max_iter = 500,
                   p = p, 
                   lambda = 0.0, 
                   mu_init = NULL,
                   B_init = B + matrix(rnorm(d * p), nrow = d, ncol = p))
fitted_unc = fit_trend(model$B, grid, model$mu, "poly")
model = prox_model(dta$dta, grid, "poly", 
                   alpha1 = 1e-2, alpha2 = 1e-3, 
                   max_iter = 500,
                   p = p, 
                   lambda = 0.5, 
                   mu_init = NULL,
                   B_init = B + matrix(rnorm(d * p), nrow = d, ncol = p))
fitted_c = fit_trend(model$B, grid, model$mu, "poly")


par(mfrow = c(5, 5), mar = c(2, 2, 1, 1))
for (j in 1:d) {
  plot(circ_2_angle(t(dta$dta[j,,])), ylab = "", xlab = "", ylim = c(-pi, pi), pch = 19,
       col = "lightblue")
  lines(circ_2_angle(t(dta$trend[j,,])), col = "steelblue", lwd = 2)
  lines(circ_2_angle(t(fitted_unc[j,,])), col = "salmon2", lty = 2, lwd = 1.5)
  lines(circ_2_angle(t(fitted_c[j,,])), col = "firebrick", lty = 3, lwd = 1.5)
}

sqrt(loss_compute(dta$trend, fitted_unc) / d)
sqrt(loss_compute(dta$trend, fitted_c) / d)
sqrt(loss_compute(dta$trend, dta$dta) / d)
