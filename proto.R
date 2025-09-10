## Example 1
close3d()
x = c(0, sin(0.5), cos(0.5))
k = 100  # number of time points
c = 1   # maximum time

Q1 = 1 * matrix(c(0, -1, 0, 1, 0, 0, 0, 0, 0), ncol = 3)
trend = sphere_trend(c(0:(c * k)) / k, array(Q1, dim = c(3, 3, 1)), x)
dta = noise_inject(trend, L = 0.1)

spheres3d(0, 0, 0, radius = 1, color = "gray", alpha = 0.6)
add_sphere_grid(alpha = 0.8)
lines3d(trend, col = "red", lwd = 3)
cols = colorRampPalette(c("lightblue", "firebrick"))(c * k + 1)
points3d(dta, col = cols, size = 5)

mu0 = mean_on_sphere(dta[1:5,])
Q0 = sphere_diff_core(dta[10,], trend[1,])$Q

trend_pred = sphere_trend(c(0:(c * k)) / k, Q0, mu0)
cols = colorRampPalette(c("lightblue", "lightgreen"))(c * k + 1)
lines3d(trend_pred, col = cols, lwd = 5)


model = spt(dta, 1, Q0, x, alpha_Q = 0.01, alpha_mu = 0.001, x_exact = FALSE,
            max.iter = 5000, verbose = TRUE)

trend_pred = sphere_trend(c(0:(c * k)) / k, model$Q, model$mu)
cols = colorRampPalette(c("lightgreen", "darkgreen"))(c * k + 1)
lines3d(trend_pred, col = cols, lwd = 5)
points3d(model$mu, col = "steelblue", size = 10)

## Example 2
close3d()
x = c(0, sin(0.5), cos(0.5))
k = 100  # number of time points
c = 1   # maximum time

Q1 = 10 * matrix(c(0, -1, 0, 1, 0, 0, 0, 0, 0), ncol = 3)
trend = sphere_trend(c(0:(c * k)) / k, array(Q1, dim = c(3, 3, 1)), x)
dta = noise_inject(trend, L = 0.1)

spheres3d(0, 0, 0, radius = 1, color = "gray", alpha = 0.6)
add_sphere_grid(alpha = 0.8)
lines3d(trend, col = "red", lwd = 3)
cols = colorRampPalette(c("lightblue", "firebrick"))(c * k + 1)
points3d(dta, col = cols, size = 5)

mu0 = mean_on_sphere(dta[1:5,])
Q0 = 5 * sphere_diff_core(dta[10,], trend[1,])$Q

trend_pred = sphere_trend(c(0:(c * k)) / k, Q0, mu0)
cols = colorRampPalette(c("lightblue", "lightgreen"))(c * k + 1)
lines3d(trend_pred, col = cols, lwd = 5)


model = spt(dta, 1, Q0, x, alpha_Q = 0.01, alpha_mu = 0.001, x_exact = FALSE,
            max.iter = 3000, verbose = TRUE)

trend_pred = sphere_trend(c(0:(c * k)) / k, model$Q, model$mu)
cols = colorRampPalette(c("lightgreen", "darkgreen"))(c * k + 1)
lines3d(trend_pred, col = cols, lwd = 5)
points3d(model$mu, col = "steelblue", size = 10)

## Example 3
close3d()
x = c(0, -sin(0.5), cos(0.5))
k = 100  # number of time points
c = 1   # maximum time

Q1 = -2 * matrix(c(0, -1, 0, 1, 0, 0, 0, 0, 0), ncol = 3)
Q2 = -1 * matrix(c(0, 0, 0, 0, 0, -1, 0, 1, 0), ncol = 3)
trend = sphere_trend(c(0:(c * k)) / k, array(c(Q2, Q1), dim = c(3, 3, 2)), x)
dta = noise_inject(trend, L = 0.05)

spheres3d(0, 0, 0, radius = 1, color = "gray", alpha = 0.6)
add_sphere_grid(alpha = 0.8)
lines3d(trend, col = "red", lwd = 3)
cols = colorRampPalette(c("lightblue", "firebrick"))(c * k + 1)
points3d(dta, col = cols, size = 5)

mu0 = mean_on_sphere(dta[1:5,])
Q0 = 1 * sphere_diff_core(dta[10,], trend[1,])$Q
Q0 = array(c(Q0, mzero(3)), dim = c(3, 3, 2))

trend_pred = sphere_trend(c(0:(c * k)) / k, Q0, mu0)
cols = colorRampPalette(c("lightblue", "lightgreen"))(c * k + 1)
lines3d(trend_pred, col = cols, lwd = 5)

model = spt(dta, 2, Q0, x, alpha_Q = 0.01, alpha_mu = 0.001, x_exact = FALSE,
            max.iter = 5000, verbose = TRUE)

trend_pred = sphere_trend(c(0:(c * k)) / k, model$Q, model$mu)
cols = colorRampPalette(c("lightgreen", "darkgreen"))(c * k + 1)
lines3d(trend_pred, col = cols, lwd = 5)
points3d(model$mu, col = "steelblue", size = 10)
