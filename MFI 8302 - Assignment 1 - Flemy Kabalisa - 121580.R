#-------- Assignment 1 - Computational Methods in Finance ---------#

#------- Flemy Kabalisa - 121580 --------#

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

#Setting up the simulator

simulate_sv <- function(nSim=10000, nSteps=1000, T=1, mu=0.15, vbar=0.04, delta_v=0.2,
                        kappa_v = 1, S0=0.5, v0= NULL, rho=0, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  if(is.null(v0)) v0 <- vbar
  
  dt <- T/nSteps
  sqrt_dt <- sqrt(dt)
  
  #setting up the vectors
  logS <- rep(log(S0), nSim)
  v <- rep(v0, nSim)
  
  for(i in seq_len(nSteps)){
    #correlated normals for increments
    Zv <- rnorm(nSim)
    eps <- rnorm(nSim)
    Zs <- rho * Zv + sqrt(1 - rho^2) * eps
    
    dW_v <- Zv * sqrt_dt
    dW_s <- Zs * sqrt_dt
    
    v_old <- v
    sqrtv_old <- sqrt(pmax(v_old, 0))
    
    #We use v_old to update log S
    logS <- logS + (mu - 0.5 * v_old) * dt + sqrtv_old * dW_s
    
    #update v
    v <- v_old + kappa_v * (vbar - v_old) *dt + delta_v * sqrtv_old * dW_v
    
    #enforce non-negativity
    v <- pmax(v, 0)
  }
  
  S_T <- exp(logS)
  return(S_T)
}

#We now derive the lognormal parameters

lognormal_params_from_meanvar <- function(M,V) {
  if (M <= 0 || V<0 ) stop("Mean must be > 0 and variance >= 0")
  sigma2 <- log(1 + V / (M^2))
  meanlog <- log(M) - 0.05 * sigma2
  sdlog <- sqrt(sigma2)
  return(list(meanlog = meanlog, sdlog = sdlog))
}

#We now bring back the parameters for our experiement

nSim <- 10000
nSteps <- 1000 #We can manipulate this to check for convergence
T <- 1

mu <- 0.15
vbar <- 0.04
delta_v_default <- 0.2
kappa_v <- 1
S0 <- 0.5
v0 <- vbar #It is the assumed initial variance
seed <- 12345

rhos <- c(-0.5,0.5)

#We now try running simulation for different rho

set.seed(seed)
results_list <- list()
summary_stats <- data.frame()

for (rho in rhos){
  cat("Simulating rho =", rho, "...\n")
  S_T <- simulate_sv(nSim = nSim, nSteps = nSteps, T=T, mu = mu, vbar = vbar,
                     delta_v = delta_v_default, kappa_v = kappa_v, S0= S0, v0=v0,
                     rho = rho, seed = NULL) #We set the seed externally
  df <- data.frame(S = S_T, rho = factor(rho, levels = rhos))
  results_list[[as.character(rho)]] <- df
  
  M <- mean(S_T)
  V <- var(S_T)
  lp <- lognormal_params_from_meanvar(M,V)
  
  summary_stats <- rbind(summary_stats,
                         data.frame(rho = rho,
                                    mean_S = M,
                                    var_S = V,
                                    meanlog = lp$meanlog,
                                    sdlog = lp$sdlog,
                                    median_S = median(S_T),
                                    q025 = quantile(S_T, 0.025),
                                    q975 = quantile(S_T, 0.975)))
}

df_all <- bind_rows(results_list)

#We now prepare the lognormal curves

x_min <- max(1e-6, min(df_all$S))
x_max <- quantile(df_all$S, probs = 0.999)
x_grid <- seq(x_min, x_max, length.out = 1000)

ln_curves <- lapply(1:nrow(summary_stats), function(i){
  row <- summary_stats[i, ]
  dens <- dlnorm(x_grid, meanlog = row$meanlog, sdlog = row$sdlog)
  data.frame(x = x_grid, density = dens, rho = factor(row$rho, levels = rho))
}) %>% bind_rows()

#We now plot the densities with an overlay of fitted lognormal

p_density <- ggplot(df_all, aes(x = S, color = rho, fill = rho)) +
  geom_density(alpha = 0.15, adjust=1)+
  geom_line(data = ln_curves, aes(x=x, y=density, color = rho), linetype = "dashed", linewidth = 1) +
  labs(title = "Density of S(1) for different rho (solid = simulated kernel, dashed = lognormal fit)",
       x = "S(1)", y = "Density", color= expression(rho), fill=expression(rho))+
  theme_minimal() +
  theme(legend.position = "right")


#Same plot but x in log-scale to inspect the tails

p_density_logx <- p_density + scale_x_log10()+
  labs(title = "Density of S(1) on log-x scale (tails visible)")

print(p_density)
print(p_density_logx)

#we save the plots

ggsave("sv_density_linear.png", p_density, width = 10, height = 6, dpi = 200)
ggsave("sv_density_logx.png", p_density_logx, width = 10, height = 6, dpi = 200)

#Print summary stats
print(summary_stats)


#---------
#We now further experiment with rho = 0, and various deltas
#---------
delta_vals <- c(0.5,1,1.5)
res_delta <- list()
summary_delta <- data.frame()

cat("\nDelta_v experiment (rho=0):\n")
for (dv in delta_vals){
  cat(" Simulating delta_v =", dv, "...\n")
  S_T <- simulate_sv(nSim = nSim, nSteps = nSteps, T=T, mu = mu, vbar = vbar,
                     delta_v = dv, kappa_v = kappa_v, S0=S0, v0 = v0, rho = 0, seed = NULL)
  df <- data.frame(S = S_T, delta = factor(dv, levels = delta_vals))
  res_delta[[as.character(dv)]] <- df
  
  M <- mean(S_T)
  V <- var(S_T)
  lp <- lognormal_params_from_meanvar(M,V)
  summary_delta <- rbind(summary_delta, data.frame(delta = dv,
                                                   mean_S = M, var_S = V,
                                                   meanlog= lp$meanlog, sdlog= lp$sdlog,
                                                   median_S = median(S_T),
                                                   q025 = quantile(S_T, 0.025),
                                                   q975 = quantile(S_T, 0.975)))
}

df_delta_all <- bind_rows(res_delta)

#Lognormal curves for our experiment

x_min_d <- max(1e-6, min(df_delta_all$S))
x_max_d <- quantile(df_delta_all$S, 0.999)
x_grid_d <- seq(x_min_d, x_max_d, length.out = 1000)
ln_curves_delta <- lapply(1:nrow(summary_delta), function(i){
  row <- summary_delta[i, ]
  dens <- dlnorm(x_grid_d, meanlog = row$meanlog, sdlog = row$sdlog)
  data.frame(x = x_grid_d, density = dens, delta = factor(row$delta, levels = delta_vals))
}) %>% bind_rows()

p_delta <- ggplot(df_delta_all, aes(x = S, color = delta, fill = delta)) +
  geom_density(alpha = 0.15, adjust=1) +
  geom_line(data = ln_curves_delta, aes(x=x, y=density, color=delta), linetype = "dashed",
            size = 1) +
  labs(title = expression("rho = 0 : Density of S(1) for various " * delta[v]),
       x = "S(1)", y="Density", color = expression(delta[v]), fill = expression(delta[v]))+
  theme_minimal()

p_delta_logx <- p_delta + scale_x_log10() +
  labs(title = expression("rho = 0 : Density (log-x) for various " * delta[v]))

print(p_delta)
print(p_delta_logx)

ggsave("sv_delta_linear.png", p_delta, width = 10, height = 6, dpi = 200)
ggsave("sv_delta_logx.png", p_delta_logx, width = 10, height = 6, dpi = 200)

print(summary_delta)


#---------------------------#

#-- Question 2 --#

#---------------------------#

#--- Part a

#We reset our environment

rm(list = ls())
set.seed(123456)

#Parameters 

S0 <- 100
r <- 0.05
sigma <- 0.2
T <- 1
N <- 50
M <- 100000

cat("Running a Monte Carlo simulation for an Asian arithmetic-average-strike call option\n")
cat(sprintf("Parameters: S0=%.1f, r=%.3f, sigma=%.3f, T=%.2f, N=%d, M=%d\n",
            S0, r, sigma, T, N, M))

dt <- T/N
sqrt_dt <- sqrt(dt)
discount_factor <- exp(-r *T)

# We setup the matrix for prices

S_mat <- matrix(0, nrow = M, ncol = N+1)
S_mat[, 1] <- S0

cat("Simulating paths using Euler (log-update) ...\n")
t0 <- proc.time()

for (i in 1:N){
  Z <- rnorm(M)
  increment <- (r - 0.05 * sigma^2) * dt + sigma * sqrt_dt * Z
  # update using: S_{t+dt} = S_t * exp(increment)
  S_mat[, i+1] <- S_mat[, i] * exp(increment)
  #Progress output (optional)
  if (i %% 10 == 0 || i==N){
    cat(sprintf(" Completed step %2d / %2d (t = %.3f)\n", i, N, i*dt))
    flush.console()
  }
}

time_elapsed <- proc.time() - t0
cat(sprintf("Simulation finished. Elapsed (sec): %.2f\n", time_elapsed["elapsed"]))


#-- We now compute the averages, payoffs, price, etc...

#Arithmetic average
S_bar <- rowMeans(S_mat[, 2:(N+1)])
S_T <- S_mat[, N+1]

#Payoffs for each path
payoff <- pmax(S_T - S_bar, 0)

#Monte Carlo estimator
mean_payoff <- mean(payoff)
sd_payoff <- sd(payoff)
se_payoff <- sd_payoff / sqrt(M)

price_mc <- discount_factor * mean_payoff
ci_lower <- discount_factor * (mean_payoff - 1.96 * se_payoff)
ci_upper <- discount_factor * (mean_payoff + 1.96 * se_payoff)

#Alternatively, we can do it by discounting the payoffs first.

discounted_payoffs <- discount_factor * payoff
price_mc_alt <- mean(discounted_payoffs)
sd_disc <- sd(discounted_payoffs)
se_disc <- sd_disc / sqrt(M)

#Confidence intervals using discounted payoffs:
ci_disc_lower <- price_mc_alt - 1.96 * se_disc
ci_disc_upper <- price_mc_alt + 1.96 * se_disc

#We now print the resdiscounted_payoffs#We now print the results & the diagnostics

cat("\nMonte Carlo Results:\n")
cat(sprintf("  Estimated price (discounted mean): %.6f\n", price_mc))
cat(sprintf("  95%% CI (discount-factor applied to mean+/-1.96*SE): [%.6f, %.6f]\n",
            ci_lower, ci_upper))
cat(sprintf("  (Alternative CI computed on discounted payoffs): [%.6f, %.6f]\n",
            ci_disc_lower, ci_disc_upper))
cat(sprintf("  Mean undisc. payoff = %.6f, sd = %.6f, SE = %.6f\n",
            mean_payoff, sd_payoff, se_payoff))
cat(sprintf("  Mean discounted payoff = %.6f, sd = %.6f, SE = %.6f\n",
            price_mc_alt, sd_disc, se_disc))


# A few percentiles and median to observe the distribution

median_disc <- median(discounted_payoffs)
q025 <- quantile(discounted_payoffs, 0.025)
q975 <- quantile(discounted_payoffs, 0.975)
cat(sprintf("  Discounted payoff median = %.6f\n", median_disc))
cat(sprintf("  Discounted payoff 2.5%% = %.6f, 97.5%% = %.6f\n", q025, q975))


#We finish by plotting the histogram and empirical density

hist(discounted_payoffs, breaks = 100, prob = TRUE,
     main = "Histogram of discounted payoffs (Asian arithmetic-strike call)",
     xlab = "Discounted payoff", col = "wheat", border= "darkgray")
lines(density(discounted_payoffs), lwd=2)


#We can look at some sample paths for a visual check

n_show <- 100
matplot(t(S_mat[1:n_show, ]), type = "l", lty = 1, col = 1:n_show,
        main = sprintf("First %d simulated price paths (sample)", n_show),
        xlab = "monitoring index (1..N)", ylab = "S_t")
legend("topright", legend = paste0("path", 1:n_show), col = 1:n_show, lty = 1, cex = 0.7)

#------
# We now do a quick check to compare our results to exact-log discretization in closed form per step
#The log-update we have is exact for GBM, so the exact simulation is in the results above
#-------

cat("\n Using the log-update (S_{t+dt}= S_t * exp(...)) is exact for GBM increments.\n")
cat("If you decided to use the standard Euler on S (S_{t+dt} = S_t + r*S_t*dt + sigma*S_t*sqrt(dt)*Z)\n",
    " you would get a slightly diffferent (less stable) result for coarse dt.\n")

#Finally, we will save our results

saveRDS(list(price = price_mc, ci = c(ci_lower, ci_upper),
        params = list(S0 = S0, r=r, sigma=sigma, T=T, N=N, M=M),
        summary = list(mean_payoff=mean_payoff, sd_payoff=sd_payoff)),
        file = "asian_mc_result.rds")
cat("\nResults saved to 'asian_mc_result.rds'\n")        


#--------------------------------#

#------ Part b: Control variate using geometric-average floating strike ------#

#--------------------------------#

#We first compute the geometric average G for each path

logS_mat <- log(S_mat)
G_vec <- exp(rowMeans(logS_mat[, 2:(N+1)]))
Y_payoff <- pmax(S_T -G_vec, 0)

#We now compute the analytic E[Y] for Y= max(S_T - G, 0)

muY_geo_float <- function(S0, r, sigma, T, N){
  #monitoring times
  t_i <- (1:N) * (T/N)
  #the moments for ln(S_T) and ln(G)
  muX <- log(S0) + (r - 0.5 * sigma^2) *T #E[ln (S_T)]
  varX <- sigma^2 * T
  muZ <- log(S0) + (r- 0.5*sigma^2) * mean(t_i)
  #for var(ln G) we need two things: sum(i,j) and min(t_i, t_j)
  tij <- outer(t_i, t_i, pmin)
  varZ <- sigma^2 * (1/(N^2)) * sum(tij)
  covXZ <- sigma^2 * (1/N) * sum(t_i)
  #We clean up to avoid degenerate cases
  if (varZ <= 0 || varX <= 0) {
    stop("Degenerate variance in muY calculation (varZ or varX <= 0).")
  }
  #correlation
  rho <- covXZ / sqrt(varX * varZ)
  
  # E1= E[e^{X} 1_{X > Z}]
  a <- 1 - covXZ/varX
  b <- (covXZ/varX) * muX - muZ
  s <- sqrt(varZ * (1- rho^2))
  #we add a step to be safe
  if (s <= 0) s<- sqrt(.Machine$double.eps)
  alpha <- a/s
  beta <- b/s
  arg1 <- ( alpha * (muX + varX) + beta) / sqrt(1+ alpha^2 * varX)
  E1 <- exp(muX + 0.5 * varX) * pnorm(arg1)
  
  # E2 = E[e*Z 1_{X > z}]
  a2 <- (covXZ / varZ) - 1
  b2 <- muX - (covXZ/varZ) * muZ
  s2 <- sqrt(varX * (1-rho^2))
  if (s2 <= 0 ) s2<- sqrt(.Machine$double.eps)
  alpha2 <- a2/s2
  beta2 <-  b2/s2
  arg2 <- ( alpha2 * (muZ + varZ) + beta2) / sqrt(1 + alpha2^2 * varZ)
  E2 <- exp(muZ + 0.5* varZ) * pnorm(arg2)
  
  muY <- E1 - E2
  return(list(muY = muY, 
              details  = list(muX = muX, varX=varX, muZ = muZ, varZ = varZ,
                              covXZ = covXZ, rho= rho, E1 = E1, E2 = E2,
                              arg1 = arg1, arg2 = arg2)))
}

#compute analytic mu_Y for current parameters 
muY_info <- muY_geo_float(S0 = S0, r=r, sigma = sigma, T = T, N = N)
muY <- muY_info$muY

cat("\nPart (b) - control variate setup:\n")
cat(sprintf("Analytic mu_Y (undiscounted) = %.8f\n", muY))

# empirical theta estimate
X_payoff <- payoff
if (var(Y_payoff) <= 0) {
  theta_hat <- 0
  warning("Var(Y) is zero or numerically negligible; setting theta_hat = 0.")
} else {
  theta_hat <- cov(X_payoff, Y_payoff) / var(Y_payoff)
}

cat(sprintf("Empirical theta_hat = %.8f\n", theta_hat))

#Adjusted estimator and discounted CV price
Z_adj <- X_payoff + theta_hat * (muY - Y_payoff)
disc_Z_adj <- discount_factor * Z_adj

price_CV <- mean(disc_Z_adj)
sd_CV <- sd(disc_Z_adj)
se_CV <- sd_CV / sqrt(M)
ci_CV <- price_CV + c(-1.96,1.96) * se_CV

cat(sprintf("Control variate price (discounted mean) = %.8f\n", price_CV))
cat(sprintf("95%% CI (CV) = [%.8f, %.8f]\n", ci_CV[1], ci_CV[2]))
cat(sprintf("sd(disc adj payoff) = %.8f, SE = %.8f\n", sd_CV, se_CV))


#----------------------------------#

#----- Part c: Compare plain MC vs control variate ------#

#----------------------------------#
price_plain <- price_mc_alt
sd_plain <- sd_disc
se_plain <- se_disc
ci_plain <- c(ci_disc_lower, ci_disc_upper)

cat("\nPart (c) - comparison:\n")
cat(sprintf("Plain Monte Carlo price = %.8f\n", price_plain))
cat(sprintf("95%% CI (Plain) = [%.8f, %.8f]\n", ci_plain[1], ci_plain[2]))
cat(sprintf("Control variate price = %.8f\n", price_CV))
cat(sprintf("95%% CI (CV) = [%.8f, %.8f]\n", ci_CV[1], ci_CV[2]))

#We notice that the Control Variates method produces a much tighter Confidence interval

#We calculate the percentage variance reduction using undiscounted payoffs

var_plain_undisc <- var(X_payoff)
var_cv_undisc <- var(Z_adj)
pct_var_reduction <- 100 * (var_plain_undisc - var_cv_undisc) / var_plain_undisc

cat(sprintf("Variance plain (undiscounted) - %.8f\n", var_plain_undisc))
cat(sprintf("Variance CV (undiscounted) = %.8f\n", var_cv_undisc))
cat(sprintf("Percentage varince reduction (undiscounted) = %.2f %%\n", pct_var_reduction))

#We now show the narrowness of CI

ci_width_plain <- diff(ci_plain)
ci_width_cv <- diff(ci_CV)
cat(sprintf("Ci width plain = %.8f, CI width CV =%.8f (reduction %.2f%%)\n",
            ci_width_plain, ci_width_cv, 100 * (ci_width_plain - ci_width_cv) / ci_width_plain))

#We overlay the plot of densities

df_plot <- data.frame(plain = discounted_payoffs, cv = disc_Z_adj)
df_long <- data.frame(value = c(df_plot$plain, df_plot$cv),
                      method = rep(c("Plain", "ControlVar"), each = length(df_plot$plain)))
ggplot(df_long, aes(x = value, color = method, fill = method)) +
  geom_density(alpha = 0.15) +
  labs(title = "Discounted payoff densities: Plain vs Control Variate",
        x= "Discounted payoff", y="Density") +
  theme_minimal()

#-----------------------------------------------#

#Part d: Experiments

#-----------------------------------------------#

#CI half-width vs M

M_seq <- c(500, 1000, 2000, 5000, 10000, 20000, 50000, M)
df_CI <- data.frame(M = integer(), method = character(), halfwidth = numeric(),
                    price = numeric(), stringsAsFactors = FALSE)

for (m in M_seq) {
  Xm <- payoff[1:m]
  Gm <- G_vec[1:m]
  Ym <- pmax(S_T[1:m] - Gm, 0)
  disc_Xm <- discount_factor* Xm
  hw_plain_m <- 1.96 * sd(disc_Xm) / sqrt(m)
  price_plain_m <- mean(disc_Xm)
  df_CI <- rbind(df_CI, data.frame(M=m, method = "Plain", halfwidth = hw_plain_m,
                                   price = price_plain_m))
  if (var(Ym) <= 0){
    theta_m <- 0
  } else {
    theta_m <- cov(Xm, Ym) / var(Ym)
  }
  Zm <- Xm + theta_m * (muY - Ym)
  disc_Zm <- discount_factor * Xm
}

ggplot(df_CI, aes(x = M, y= halfwidth, color = method)) +
  geom_line() + geom_point() + scale_x_log10() + scale_y_log10()+
  labs(title = "95% CI half-width vs M (Plain vs Control Variate)", 
       x = "Number of simulations M (log schale)", y= "CI half-width (log scale)")+
  theme_minimal() -> p_ci_vs_M
print(p_ci_vs_M)

#We can tell from our output that increasing M consistently decrease the width of the CI.
#This results in a consistently increasing accuracy at a proportionally heavy cost of computing power.


# CV efficiency depending on N and sigma
#We start with multiple N with a moderate M_rep for each

N_values <- c(10,50,100)
sigma_values <- c(0.1,0.2,0.4)
M_rep <- 20000
set.seed(2025)

exp_results <- data.frame(N= integer(), sigma = numeric(), method = character(), 
                          price = numeric(), sd_disc = numeric(), ci_halfwidth = numeric(),
                          pct_var_reduction = numeric(), stringsAsFactors = FALSE)

simulate_paths <- function(S0, r, sigma, T, N, M, seed= NULL){
  if (!is.null(seed)) set.seed(seed)
  dt <- T/N
  sqrt_dt <- sqrt(dt)
  Z <- matrix(rnorm(M *N), nrow = M, ncol = N)
  logS <- matrix(0, nrow= M, ncol = N+1)
  logS[, 1] <- log(S0)
  for (i in 1:N) {
    logS[, i+1] <- logS[, i] + (r-0.5 *sigma^2) *dt + sigma * sqrt_dt * Z[,i]
  }
  S <- exp(logS)
  S_T_sim <- S[, N+1]
  S_bar_arith_sim <- rowMeans(S[, 2:(N+1)])
  S_bar_geom_sim <- exp(rowMeans(logS[, 2:(N+1)]))
  Xs <- pmax(S_T_sim - S_bar_arith_sim, 0)
  Ys <- pmax(S_T_sim - S_bar_geom_sim, 0)
  return(list(X=Xs, Y=Ys, S_T = S_T_sim, S_arith = S_bar_arith_sim, S_geom = S_bar_geom_sim))
}

cat("\nPart(d) - exploring dependence on N and sigma (M_rep =", M_rep, ")\n")
for (Nv in N_values) {
  for (sigv in sigma_values){
    cat(sprintf("Simulating N=%d, sigma=%.2f ...\n", Nv, sigv))
    sim <- simulate_paths(S0 = S0, r=r, sigma = sigv, T=T, N=Nv, M=M_rep, seed = NULL)
    Xv <- sim$X
    Yv <- sim$Y
    disc_Xv <- discount_factor * Xv
    price_plain_v <- mean(disc_Xv)
    sd_plain_v <- sd(disc_Xv)
    hw_plain_v <- 1.96 * sd_plain_v / sqrt(M_rep)
    muYv <- muY_geo_float(S0 = S0, r=r, sigma = sigv, T=T, N=Nv)$muY
    if (var(Yv) <= 0) theta_v <- 0 else theta_v <- cov(Xv, Yv) / var(Yv)
    Zv <-  Xv + theta_v * (muYv -Yv)
    disc_Zv <- discount_factor * Zv
    price_cv_v <- mean(disc_Zv)
    sd_cv_v <- sd(disc_Zv)
    hw_cv_v <- 1.96 * sd_cv_v / sqrt(M_rep)
    pct_red_v <- 100 * (var(Xv) - var(Zv)) /var(Xv)
    exp_results <- rbind(exp_results,
                         data.frame(N = Nv, sigma = sigv, method = "Plain",
                                    price = price_plain_v, sd_disc = sd_plain_v, ci_halfwidth = hw_plain_v,
                                    pct_var_reduction = NA),
                         data.frame(N = Nv, sigma = sigv, method = "ControlVar",
                                    price = price_cv_v, sd_disc = sd_cv_v, ci_halfwidth = hw_cv_v,
                                    pct_var_reduction = pct_red_v))
    cat(sprintf(" -> pct variance reduction (CV) = %.2f %%\n", pct_red_v))
  }
}

print(exp_results)

df_red <- subset(exp_results, method == "ControlVar")
ggplot(df_red, aes(x= factor(sigma), y = pct_var_reduction, group = factor(N), color = factor(N))) +
  geom_point() + geom_line(aes(group = factor(N))) +
  labs(title = "Percent variance reduction by Control Variate",
       x = expression(sigma), y= "Variance reduction(%)", color = "N (monitor points)") +
  theme_minimal() -> p_red 
print(p_red)

#We notice that the number of time steps N does not seem to create a big difference,
#while on the other hand, and increasing sigma reduces the Variance reduction percentage.

#Final summary prints
cat("\nSummary (final):\n")
cat(sprintf("Plain Monte Carlo price: %.8f (95%% CI width = %.8f)\n", price_plain, ci_width_plain))
cat(sprintf("Control variate price: %.8f (95%% CI width = %.8f)\n", price_CV, ci_width_cv))
cat(sprintf("Emprical theta_hat = %.8f, analytic mu_Y (undiscounted) = %.8f\n", theta_hat, muY))
cat(sprintf("Overall percentage variance reduction (using entire sample) = %.2f %%\n", pct_var_reduction))
cat("Experiment table (N, sigma, method, price, sd_disc, ci_halfwidth, pct_var_reduction):\n")
print(exp_results)


#We now decide to save the outputs
saveRDS(list(plain_price = price_plain, plain_CI = ci_plain,
             cv_price = price_CV, cv_CI = ci_CV,
             theta_hat = theta_hat, muY = muY,
             df_CI = df_CI, exp_results = exp_results),
        file = "asian_mc_controls_results.rds")
cat("\nSaved results to 'asian_mc_controls_results.rds'\n")
