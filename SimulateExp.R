library(jagsUI)
library(survival)

#### fit exp survival ####
nind_sim <- 1500
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
params = c("b0")
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

fitExpSim = jags(dat, inits, params, model.file="./Models/SurvModel_ExpTrap_eff.txt", 
                                n.chains=nc, n.iter=ni, n.burnin=nb, 
                                n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                                n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitExpSim, digits=3)
fitExpSim$mean$b0
lam_est <- exp(fitExpSim$mean$b0)
lam_est
1/lam_est # the mean of expo distribution

summary(fitExpSim)
fitExpSim$mean$t[cens_sim==1]

# Try to fit exp with Weibull model (v should be = 1)
params = c("v", "b0")


inits <- function() {
  list(b0=rnorm(1, 0.1, 0.5),
       v=runif(1, 0.2, 1.75),
       t=t.start_sim)
}

fitExpWIthWeiSim = jags(dat, inits, params, model.file="./Models/SurvModel_WeibTrap_eff.txt", 
                 n.chains=nc, n.iter=ni, n.burnin=nb, 
                 n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                 n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitExpWIthWeiSim, digits=3)
fitExpWIthWeiSim$mean$b0
lam_est <- exp(fitExpWIthWeiSim$mean$b0)
lam_est
1/lam_est # the mean of expo distribution

#### fit Weibull survival ####
library(simsurv)
library(flexsurv)

nind_sim <- 1500
lam <- 2
v <- 1.5

# Create a data frame with the subject IDs and treatment covariate
cov <- data.frame(id = 1:nind_sim)

# Simulate the event times
dat <- simsurv(lambdas = lam, 
               gammas = v, 
               x = cov, 
               maxt = tmax
               )
summary(dat$eventtime)
table(dat$status)

# Merge the simulated event times onto covariate data frame
dat <- merge(cov, dat)

# Fit a Weibull proportional hazards model with flexsurv::flexsurvspline
mod <- flexsurv::flexsurvspline(Surv(eventtime, status) ~ 1, data = dat)
mod$res

coef(mod)
exp(coef(mod)[1]) # Lambda

# Jags expect different formatting
#            "censored" 1 if censored, 0 if trapped
#            "TimeTrap" time when trapped (if not trapped NA) calculated as ntrapdays/time_period
#            "t.cen"  censoring time (max time). If trapped NA
#            "t.start"  starting value for censored rows
nind_sim - sum(dat$status)
cens_sim_w <- abs(dat$status - 1)
sum(cens_sim_w)
t_sim_w <- dat$eventtime
t_sim_w[cens_sim_w==1] <- NA
summary(t_sim_w)
t.start_sim_w <- rep(NA, nind_sim)
t.start_sim_w[cens_sim_w==1] <- tmax + 0.01
t.cen_sim_w <- rep(tmax, nind_sim)

# JAGS
dat=list(censored=cens_sim_w, t=t_sim_w, t.cen=t.cen_sim_w, nind=nind_sim)

# arguments 
params = c("v", "b0")


inits <- function() {
  list(b0=rnorm(1, 0.1, 0.5),
       v=runif(1, 0.2, 1.75),
       t=t.start_sim_w)
}

# call to JAGS
ni<- 5000
nb<- 2000
nt<- 1
nc<- 3
np <- 8 # Number of CPUs

fitWeiSim = jags(dat, inits, params, model.file="./Models/SurvModel_WeibTrap_eff.txt", 
                 n.chains=nc, n.iter=ni, n.burnin=nb, 
                 n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                 n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitWeiSim, digits=3)
fitWeiSim$mean$b0
fitWeiSim$mean$v
fitWeiSim$summary
lam_w_est <- exp(fitWeiSim$mean$b0)
lam_w_est
#fitWeiSim$mean$t[cens_sim_w==1]

