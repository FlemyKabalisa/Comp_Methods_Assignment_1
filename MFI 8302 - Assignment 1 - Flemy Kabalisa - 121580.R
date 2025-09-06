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

