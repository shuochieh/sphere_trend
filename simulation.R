source("util.R")
num_sim = 100
max.iter = 1000
type = 1
alpha = 0.8
n = 100

Q_loss = rep(NA, num_sim)
Q0_loss = rep(NA, num_sim)
mu_loss = rep(NA, num_sim)
mu0_loss = rep(NA, num_sim)
training_loss = array(NA, dim = c(num_sim, max.iter))

R2_0 = rep(NA, num_sim)
R2 = rep(NA, num_sim)
alpha_hat = rep(NA, num_sim)

for (sim in 1:num_sim) {
  if (type == 1) {
    # generate data
    mu = c(0, -sin(0.5), cos(0.5))
    Q1 = 1 * matrix(c(0, -1, 0, 1, 0, 0, 0, 0, 0), ncol = 3)
    Q_star = array(Q1, dim = c(3, 3, 1))
    trend_model = noisy_trend(c(1:n) / n, Q_star, mu, s = 0.01, alpha = alpha)
    dta = trend_model$dta
    trend = trend_model$trend
    V = trend_model$V
    
    # initialization
    Q0 = array(sphere_diff_core(dta[10,], dta[1,])$Q, dim = c(3, 3, 1))
    mu0 = dta[1,]
    
    Q0_loss[sim] = sqrt(0.5 * sum((Q_star - Q0)^2))
    mu0_loss[sim] = geod_sphere(mu0, mu)
    
    R2_0[sim] = sqrt(mean(geod_sphere(dta, trend)^2))
  }
  
  # estimate the spherical polynomial trend
  model = spt(dta, 1, Q0, mu0, alpha_Q = 0.01, alpha_mu = 0.001, x_exact = FALSE,
              max.iter = max.iter, verbose = FALSE)
  Q_loss[sim] = sqrt(0.5 * sum((Q_star - model$Q)^2))
  mu_loss[sim] = geod_sphere(model$mu, mu)
  
  if (length(model$loss) != max.iter) {
    temp = model$loss
    training_loss[sim,] = c(temp, rep(tail(temp, 1), max.iter - length(temp)))
  } else {
    temp = model$loss[-max.iter]
    training_loss[sim,] = c(temp, tail(temp, 1))
  }
  
  trend_hat = sphere_trend(c(1:n) / n, model$Q, model$mu)
  R2[sim] = mean(geod_sphere(trend_hat, trend)^2)
  
  V_hat = matrix(NA, nrow = n - 1, ncol = ncol(dta)) 
  V_hat_lag = matrix(NA, nrow = n - 1, ncol = ncol(dta)) 
  for (j in 2:n) {
    V_hat[j - 1,] = Log_sphere(dta[j,], trend_hat[j,])
    V_hat_lag[j - 1,] = tail(pt_sphere((j - 1) / n, j / n, Log_sphere(dta[j - 1,], trend_hat[j - 1,]),
                                       list(Q = model$Q, x = model$mu, bias = FALSE)), 1)[-1]
  }
  
  alpha_hat[sim] = sum(V_hat * V_hat_lag) / sum(V_hat_lag^2)
  
  cat("simulation", sim, "\n")
}