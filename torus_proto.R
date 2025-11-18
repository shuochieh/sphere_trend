low_rank_gen = function (m, n, r, s = 1) {
  A = matrix(rnorm(m * r), ncol = r)
  B = matrix(rnorm(n * r), ncol = r)
  
  return (A %*% t(B))
}

# Case 1 

n = 50
d = 6
q = 5

B = low_rank_gen(d, 2 * q, 1)
dta = dta_gen(n, B, "season", 1:n, sigma = 1.25, q = q)

model = prox_model(dta$dta, 1:n, "season", alpha1 = 1e-3, alpha2 = 1e-3, 
                   max_iter = 1000,
                   q = q, 
                   lambda = 0.000, 
                   mu_init = NULL,
                   B_init = B + 0 * matrix(rnorm(d * 2 * q, sd = 2), nrow = d, ncol = 2 * q))
fitted = fit_trend(model$B, 1:n, model$mu, "season")

par(mfrow = c(1, 1))
plot(log(model$loss), type = "b")
par(mfrow = c(2, 3), mar = c(2, 2, 1, 1))
for (j in 1:d) {
  plot(circ_2_angle(t(dta$dta[j,,])), ylab = "", xlab = "", ylim = c(-pi, pi), pch = 19,
       col = "lightblue")
  lines(circ_2_angle(t(dta$trend[j,,])), col = "steelblue", lwd = 1.5)
  lines(circ_2_angle(t(fitted[j,,])), col = "salmon2", lty = 2)
}


# Case 2

n = 30
d = 6
p = 5

B = low_rank_gen(d, p, 1)
dta = dta_gen(n, B, "poly", c(1:n)/n, sigma = 0.5, q = q)

model = prox_model(dta$dta, c(1:n)/n, "poly", alpha1 = 1e-2, alpha2 = 1e-2, 
                   max_iter = 1000,
                   p = p, lambda = 0.000, 
                   mu_init = NULL,
                   B_init = B + 1 * matrix(rnorm(d * p, sd = 1), nrow = d, ncol = p))
fitted = fit_trend(model$B, c(1:n)/n, model$mu, "poly")
par(mfrow = c(1, 1))
plot(log(model$loss), type = "b")

par(mfrow = c(2, 3), mar = c(2, 2, 1, 1))
for (j in 1:d) {
  plot(circ_2_angle(t(dta$dta[j,,])), ylab = "", xlab = "", ylim = c(-pi, pi), pch = 19,
       col = "lightblue")
  lines(circ_2_angle(t(dta$trend[j,,])), col = "steelblue", lwd = 1.5)
  lines(circ_2_angle(t(fitted[j,,])), col = "salmon2", lty = 2)
}