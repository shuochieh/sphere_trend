source("torus_util.R")

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

n = 50
d = 25
q = 5
grid = c(1:n) 

B = 8 * low_rank_gen(d, 2 * q, 1)
dta = dta_gen(n, B, "season", grid, sigma = 0.5, q = q)

B_init = B + 0.0 * matrix(rnorm(d * 2 * q), nrow = d, ncol = 2 * q)

model_unc = prox_model(dta$dta, grid, "season", 
                   alpha1 = 1e-2, alpha2 = 1e-3, 
                   max_iter = 1000,
                   q = q, 
                   lambda = 0.0, 
                   mu_init = dta$mu,
                   B_init = B_init)
fitted_unc = fit_trend(model_unc$B, grid, model_unc$mu, "season")

# model_c = prox_model(dta$dta, grid, "season", 
#                    alpha1 = 1e-2, alpha2 = 1e-3, 
#                    max_iter = 700,
#                    q = q, 
#                    lambda = 1.0, 
#                    mu_init = dta$mu,
#                    B_init = B_init)
# fitted_c = fit_trend(model_c$B, grid, model_c$mu, "season")

model_c = prox_model_decay(dta$dta, grid, "season", 
                           alpha1 = 1e-2, alpha2 = 1e-3, 
                           decay1 = 0.5, decay2 = 1,
                           max_iter = 1000,
                           q = q, 
                           lambda = 1.0, 
                           mu_init = dta$mu,
                           B_init = B_init)
fitted_c = fit_trend(model_c$B, grid, model_c$mu, "season")

par(mfrow = c(5, 5), mar = c(2, 2, 1, 1))
for (j in 1:d) {
  plot(circ_2_angle(t(dta$dta[j,,])), ylab = "", xlab = "", ylim = c(-pi, pi), pch = 19,
       col = "lightblue")
  lines(circ_2_angle(t(dta$trend[j,,])), col = "steelblue", lwd = 2)
  lines(circ_2_angle(t(fitted_unc[j,,])), col = "salmon2", lty = 2, lwd = 2)
  lines(circ_2_angle(t(fitted_c[j,,])), col = "darkolivegreen4", lty = 3, lwd = 2)
}

sqrt(loss_compute(dta$trend, fitted_unc) / d)
sqrt(loss_compute(dta$trend, fitted_c) / d)

par(mfrow = c(1, 1))
plot(log((model_unc$loss)), type = "b", pch = 19, col = "salmon", ylim = c(0, 5))
lines(log((model_c$loss)), type = "b", pch = 19, col = "darkolivegreen4")


training_err_unc = array(0, dim = c(d, 1000))
for (iter in 1:1000) {
  B_hat = model_unc$B_history[,,iter]
  mu_hat = model_unc$mu_history[,,iter]
  fitted_unc = fit_trend(B_hat, grid, mu_hat, "season")
  for (j in 1:d) {
    training_err_unc[j,iter] = sqrt(loss_compute(dta$dta[j,,], fitted_unc[j,,]) / d)
  }
  
  cat("    iter", iter, "\n")
}

training_err_c = array(0, dim = c(d, 1000))
obj_val = array(0, dim = c(d, 1000))
par(mfrow = c(5, 5), mar = c(2, 2, 1, 1))
for (iter in 1:1000) {
  B_hat = model_c$B_history[,,iter]
  mu_hat = model_c$mu_history[,,iter]
  fitted_c = fit_trend(B_hat, grid, mu_hat, "season")
  for (j in 1:d) {
    training_err_c[j,iter] = sqrt(loss_compute(dta$dta[j,,], fitted_c[j,,]) / d)
    obj_val[j, iter] = training_err_c[j,iter]^2 + 1.0 * sum(svd(B_hat)$d)
  }
  
  cat("    iter", iter, "\n")
}

par(mfrow = c(5, 5), mar = c(2, 2, 1, 1))
for (j in 1:d) {
  plot(training_err_unc[j,], type = "l", col = "salmon", ylim = c(0, max(training_err_c)),
       main = paste("Series", j),
       ylab = "")
  lines(training_err_c[j,], type = "l", col = "darkolivegreen4")
}

par(mfrow = c(5, 5), mar = c(2, 2, 1, 1))
for (j in 1:d) {
  plot(obj_val[j,], type = "l", col = "darkolivegreen4",
       main = paste("Series", j),
       ylab = "Obj Func")
}


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
