library(jagsUI)
library(survival)

#### fit survival  Weibull trap effort ####
nobserved <- 150
nmiss <- 2
nind_sim <- nobserved + nmiss
lam <- 1
add_tmax <- 0.01
t_sim <- c(rexp(nobserved, lam), rep(NA, nmiss))
t.cen_sim <- rep(max(t_sim,na.rm = T) + add_tmax, nobserved + nmiss)

cens_sim <- c(rep(0, nobserved), rep(1, nmiss))
t.start_sim <- c(rep(NA, nobserved), rep(max(t_sim,na.rm = T) + 2*add_tmax, nmiss))
sim_df <- data.frame(t_sim, t.cen_sim, t.start_sim, 
                     t.comb=c(t_sim[seq_len(nobserved)], t.cen_sim[-seq_len(nobserved)]),
                     status=1 - cens_sim, cens=cens_sim)
plot(survfit(Surv(sim_df$t.comb, sim_df$status)~1))
expmod_eg <- survreg(Surv(sim_df$t.comb, sim_df$status)~1, sim_df, dist = "exponential")
summary(expmod_eg)
predict(expmod_eg, type = "linear")
predict(expmod_eg, type = "response")
mean(t_sim, na.rm = T)
exp(expmod_eg$coefficients)

# JAGS
dat=list(censored=cens_sim, t=t_sim, t.cen=t.cen_sim, nind=nind_sim)

# arguments 
params = c("mu.b0", "sigma.b0", "b0")
params = c("mu.b0", "sigma.b0", "t_est", "t")

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
exp(lam_est)

summary(fitExpSim)
fitExpSim$mean$t

