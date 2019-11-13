library(jagsUI)


#### fit survival  Weibull trap effort ####
nobserved <- 150
nmiss <- 2
nind_sim <- nobserved + nmiss
lam <- 1

t_sim <- c(rexp(nobserved, lam), rep(NA, nmiss))
t.cen_sim <- rep(max(t_sim,na.rm = T) + 1, nobserved + nmiss)

cens_sim <- c(rep(0, nobserved), rep(1, nmiss))
t.start_sim <- c(rep(NA, nobserved), rep(max(t_sim,na.rm = T) + 2, nmiss))

dat=list(censored=cens_sim, t=t_sim, t.cen=t.cen_sim, nind=nind_sim)

# arguments 
params = c("mu.b0", "sigma.b0", "b0")
params = c("mu.b0", "sigma.b0")

inits <- function() {
  list(mu.b0=rnorm(1, 0.1, 0.5),
       sigma.b0=runif(1, 0.2, 0.75),
       t=t.start_sim)
}

# call to JAGS
ni<- 50000
nb<- 20000
nt<- 1
nc<- 3
np <- 8 # Number of CPUs

fitExpSim = jags(dat, inits, params, model.file="./Models/SurvModel_ExpTrap_eff_ind.txt", 
                                n.chains=nc, n.iter=ni, n.burnin=nb, 
                                n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                                n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitExpSim, digits=3)
exp(fitExpSim$mean$mu.b0)
summary(fitExpSim)

