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
