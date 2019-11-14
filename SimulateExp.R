library(jagsUI)
library(survival)

#### fit survival  Weibull trap effort ####
nind_sim <- 150
lam <- 2
tmax <- qexp(p = 0.9, lam)
t_sim <- rexp(nind_sim, lam)
summary(t_sim)
mean_t_sim <- mean(t_sim)
t.cen_sim <- rep(tmax, nind_sim)

cens_sim <- as.integer(t_sim >= tmax)
t_sim[cens_sim==1] <- NA
summary(t_sim)
t.start_sim <- rep(NA, nind_sim)
t.start_sim[cens_sim==1] <- tmax + 0.01
t.comb <- t_sim
t.comb[is.na(t_sim)] <- tmax
sim_df <- data.frame(t_sim, t.cen_sim, t.start_sim, 
                     t.comb=t.comb,
                     status=1 - cens_sim, cens=cens_sim)
plot(survfit(Surv(sim_df$t.comb, sim_df$status)~1))
expmod_eg <- survreg(Surv(sim_df$t.comb, sim_df$status)~1, sim_df, dist = "exponential")
summary(expmod_eg)
mean(predict(expmod_eg, type = "response"))
exp(expmod_eg$coefficients)

# JAGS
dat=list(censored=cens_sim, t=t_sim, t.cen=t.cen_sim, nind=nind_sim)

# arguments 
params = c("mu.b0", "sigma.b0", "b0")
#params = c("mu.b0", "sigma.b0", "t_est", "t")

inits <- function() {
  list(mu.b0=rnorm(1, 0.1, 0.5),
       sigma.b0=runif(1, 0.2, 0.75),
       t=t.start_sim)
}

# call to JAGS
ni<- 5000
nb<- 2000
nt<- 1
nc<- 3
np <- 8 # Number of CPUs

fitExpSim = jags(dat, inits, params, model.file="./Models/SurvModel_ExpTrap_eff_ind.txt", 
                                n.chains=nc, n.iter=ni, n.burnin=nb, 
                                n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                                n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitExpSim, digits=3)
fitExpSim$mean$mu.b0
lam_est <- exp(fitExpSim$mean$mu.b0)
lam_est
1/lam_est # the eman of expo distribution

summary(fitExpSim)
fitExpSim$mean$t

