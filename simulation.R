source("util.R")
num_sim = 100
max.iter = 10000
alpha = 0.8
n = 100

for (type in c(1, 3, 4, 5, 6)) {
  Q_loss = rep(NA, num_sim)
  Q0_loss = rep(NA, num_sim)
  mu_loss = rep(NA, num_sim)
  mu0_loss = rep(NA, num_sim)
  training_loss = array(NA, dim = c(num_sim, max.iter))
  
  R2_0 = rep(NA, num_sim)
  R2_var = rep(NA, num_sim)
  R2 = rep(NA, num_sim)
  alpha_hat = rep(NA, num_sim)
  
  start_time <- proc.time()
  for (sim in 1:num_sim) {
    if (type == 1) {
      # generate data
      mu = c(0, -sin(0.5), cos(0.5))
      Q1 = 1 * matrix(c(0, -1, 0, 1, 0, 0, 0, 0, 0), ncol = 3)
      Q_star = array(Q1, dim = c(3, 3, 1))
      trend_model = noisy_trend(c(1:n) / n, Q_star, mu, s = 0.01, alpha = alpha)
      dta = trend_model$dta
      trend = trend_model$trend
      
      # initialization
      Q0 = array(sphere_diff_core(dta[10,], dta[1,])$Q, dim = c(3, 3, 1))
      mu0 = dta[1,]
      
      trend_init = sphere_trend(c(1:n) / n, Q0, mu0)
      
      Q0_loss[sim] = sqrt(0.5 * sum((Q_star - Q0)^2))
      mu0_loss[sim] = geod_sphere(mu0, mu)
      
      R2_var[sim] = sqrt(mean(geod_sphere(dta, trend)^2))
      R2_0[sim] = sqrt(mean(geod_sphere(trend_init, trend)^2))
    }
    
    if (type == 2) {
      # generate data 
      mu = c(0, -sin(0.5), cos(0.5))
      Q1 = 10 * matrix(c(0, -1, 0, 1, 0, 0, 0, 0, 0), ncol = 3)
      Q_star = array(Q1, dim = c(3, 3, 1))
      trend_model = noisy_trend(c(1:n) / n, Q_star, mu, alpha = alpha, s = 0.01)
      dta = trend_model$dta
      trend = trend_model$trend
      
      # (bad) initialization
      Q0 = array(sphere_diff_core(dta[10,], dta[1,])$Q, dim = c(3, 3, 1))
      mu0 = dta[1,]
      
      trend_init = sphere_trend(c(1:n) / n, Q0, mu0)
      
      Q0_loss[sim] = sqrt(0.5 * sum((Q_star - Q0)^2))
      mu0_loss[sim] = geod_sphere(mu0, mu)
      
      R2_var[sim] = sqrt(mean(geod_sphere(dta, trend)^2))
      R2_0[sim] = sqrt(mean(geod_sphere(trend_init, trend)^2))
    }
    
    if (type == 3) {
      # generate data 
      mu = c(0, -sin(0.5), cos(0.5))
      Q1 = 10 * matrix(c(0, -1, 0, 1, 0, 0, 0, 0, 0), ncol = 3)
      Q_star = array(Q1, dim = c(3, 3, 1))
      trend_model = noisy_trend(c(1:n) / n, Q_star, mu, alpha = alpha, s = 0.01)
      dta = trend_model$dta
      trend = trend_model$trend
      
      # (good) initialization
      Q0 = 5 * array(sphere_diff_core(dta[10,], dta[1,])$Q, dim = c(3, 3, 1))
      mu0 = dta[1,]
      
      trend_init = sphere_trend(c(1:n) / n, Q0, mu0)
      
      Q0_loss[sim] = sqrt(0.5 * sum((Q_star - Q0)^2))
      mu0_loss[sim] = geod_sphere(mu0, mu)
      
      R2_var[sim] = sqrt(mean(geod_sphere(dta, trend)^2))
      R2_0[sim] = sqrt(mean(geod_sphere(trend_init, trend)^2))
    }
    
    if (type == 4) {
      # generate data 
      mu = c(0, -sin(0.5), cos(0.5))
      Q1 = 1 * matrix(c(0, 0, 0, 0, 0, 1, 0, -1, 0), ncol = 3)
      Q2 = 2 * matrix(c(0, 1, 0, -1, 0, 0, 0, 0, 0), ncol = 3)
      Q_star = array(c(Q1, Q2), dim = c(3, 3, 2))
      trend_model = noisy_trend(c(1:n) / n, Q_star, mu, alpha = alpha, s = 0.01)
      dta = trend_model$dta
      trend = trend_model$trend
      
      # initialization
      Q0 = 1 * sphere_diff_core(dta[10,], dta[1,])$Q
      Q0 = array(c(Q0, rep(0, 9)), dim = c(3, 3, 2))
      mu0 = dta[1,]
      
      trend_init = sphere_trend(c(1:n) / n, Q0, mu0)
      
      Q0_loss[sim] = sqrt(0.5 * sum((Q_star - Q0)^2))
      mu0_loss[sim] = geod_sphere(mu0, mu)
      
      R2_var[sim] = sqrt(mean(geod_sphere(dta, trend)^2))
      R2_0[sim] = sqrt(mean(geod_sphere(trend_init, trend)^2))
    }
    
    if (type == 5) {
      # generate data
      mu = rep(0.5, 4)
      Q1 = matrix(c(0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, sqrt(2), 0,0, -sqrt(2), 0), ncol = 4)
      Q_star = array(Q1, dim = c(4, 4, 1))
      trend_model = noisy_trend(c(1:n) / n, Q_star, mu, s = 0.01, alpha = alpha)
      dta = trend_model$dta
      trend = trend_model$trend
      
      # initialization
      Q0 = array(sphere_diff_core(dta[10,], dta[1,])$Q, dim = c(4, 4, 1))
      mu0 = dta[1,]
      
      trend_init = sphere_trend(c(1:n) / n, Q0, mu0)
      
      Q0_loss[sim] = sqrt(0.5 * sum((Q_star - Q0)^2))
      mu0_loss[sim] = geod_sphere(mu0, mu)
      
      R2_var[sim] = sqrt(mean(geod_sphere(dta, trend)^2))
      R2_0[sim] = sqrt(mean(geod_sphere(trend_init, trend)^2))
      
    }
    
    if (type == 6) {
      # generate data 
      mu = c(0, -sin(0.5), cos(0.5))
      Q1 = 1 * matrix(c(0, 0, 0, 0, 0, 1, 0, -1, 0), ncol = 3)
      Q2 = 0 * matrix(c(0, 1, 0, -1, 0, 0, 0, 0, 0), ncol = 3)
      Q_star = array(c(Q1, Q2), dim = c(3, 3, 2))
      trend_model = noisy_trend(c(1:n) / n, Q_star, mu, alpha = alpha, s = 0.01)
      dta = trend_model$dta
      trend = trend_model$trend
      
      # initialization
      Q0 = 1 * sphere_diff_core(dta[10,], dta[1,])$Q
      Q0 = array(c(Q0, rep(0, 9)), dim = c(3, 3, 2))
      mu0 = dta[1,]
      
      trend_init = sphere_trend(c(1:n) / n, Q0, mu0)
      
      Q0_loss[sim] = sqrt(0.5 * sum((Q_star - Q0)^2))
      mu0_loss[sim] = geod_sphere(mu0, mu)
      
      R2_var[sim] = sqrt(mean(geod_sphere(dta, trend)^2))
      R2_0[sim] = sqrt(mean(geod_sphere(trend_init, trend)^2))
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
  
  total_time = proc.time() - start_time
  print(round(total_time["elapsed"] / 60, 2))
  
  pdf(paste0("Type", type, "_estimation.pdf"), width = 12, height = 6)
  par(mfrow = c(1, 3))
  boxplot(cbind(Q0_loss, Q_loss), names = c("initialization", "gradient descent"),
          main = "estimation error of the skew-symmetric coefficients",
          ylim = c(0, max(c(Q0_loss, Q_loss))))
  
  boxplot(cbind(mu0_loss, mu_loss), names = c("initialization", "gradient descent"),
          main = "estimation error of the mu_0",
          ylim = c(0, max(c(mu0_loss, mu_loss))))
  
  boxplot(cbind(R2_var, R2), names = c("Noise level", "gradient descent"),
          main = "estimation error of the trend",
          ylim = c(0, max(c(R2_var, R2))))
  dev.off()
  
  pdf(paste0("Type", type, "_AR.pdf"), width = 12, height = 6)
  par(mfrow = c(1, 1))
  hist(alpha_hat, xlim = c(0, 1), main = "Histogram of the estimated AR coefficient",
       xlab = "", ylab = "", probability = T)
  abline(v = 0.8, col = "red")
  dev.off()
}
