library(jagsUI)
source("Functions.R")
source("Trapping_DataPrep.R")

#### Bernoulli model ####
#### fit intercept only ####
dat <- list(y=caphist_m, nind=length(IDs), nobs=nobs)
# arguments 
params = c("p", "b0", "y")

inits <- function() {
  list(b0=rnorm(1, -3, 0.5))
}

# call to JAGS
ni<- 2000
nb<- 100
nt<- 1
nc<- 3
np <- 8 # Number of CPUs

fit = jags(dat, inits, params, model.file="./Models/SimpleBernModel.txt", 
              n.chains=nc, n.iter=ni, n.burnin=nb, 
              n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
              n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fit$mean$b0, digits=3)
inv.logit(fit$mean$b0)
1-(1-inv.logit(fit$mean$b0))^(15*28) # Prob for 15 traps deployed 28 days
fit$n.eff$b0
fit$Rhat$b0

ggplot(data.frame(b0=fit$sims.list$b0, Iter=seq_along(fit$sims.list$b0))) + 
  geom_line(aes(Iter, b0))


#----- Test ----#
nc<- 1
fit = jags(dat, inits = list(p=0.0001), params, model.file="./Models/SimpleBernModelv2.txt", 
           n.chains=nc, n.iter=ni, n.burnin=nb, 
           n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
           n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fit$mean$p, digits=3)

fit$n.eff$p
fit$Rhat$p

ggplot(data.frame(p=fit$sims.list$p, Iter=seq_along(fit$sims.list$p))) + 
  geom_line(aes(Iter, p))
#----------#

#### fit with covariates (ndays and distance) ####
dat <- list(y=caphist_m, nind=length(IDs), nobs=nobs, d=dist, t=time_diff)
# arguments 
params = c("b0", "bd", "bt")

inits <- function() {
  list(b0=rnorm(1, -3, 0.5),
       bd=rnorm(1, 0, 0.5),
       bt=rnorm(1, 0, 0.5))
}

# call to JAGS
ni<- 7000
nb<- 5000
nt<- 1
nc<- 3
np <- 8 # Number of CPUs

fitCov = jags(dat, inits, params, model.file="./Models/BernModel_time_dist.txt", 
           n.chains=nc, n.iter=ni, n.burnin=nb, 
           n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
           n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitCov, digits=3)
sapply(fitCov$n.eff, summary)
fitCov$Rhat
plot(fitCov)

#### fit with distance ####
dat <- list(y=caphist_m, nind=length(IDs), nobs=nobs, d=dist)
# arguments 
params = c("b0", "bd")

inits <- function() {
  list(b0=rnorm(1, -3, 0.5),
       bd=rnorm(1, -1, 0.5)
       )
}

# call to JAGS
ni<- 12000
nb<- 8000
nt<- 1
nc<- 3
np <- 8 # Number of CPUs

fitDist = jags(dat, inits, params, model.file="./Models/BernModel_dist.txt", 
              n.chains=nc, n.iter=ni, n.burnin=nb, 
              n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
              n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitDist, digits=3)
sapply(fitDist$n.eff, summary)
fitDist$Rhat
plot(fitDist)

#### Survivival ####

#### fit simple survival ####
dat <- list(y=caphist_m, nind=length(IDs), nobs=nobs, t=time_diff)
# arguments 
params = c("v", "b0")

inits <- function() {
  list(b0=rnorm(1, 0, 0.5),
       v=rgamma(1, 2.5, 3))
}

# call to JAGS
ni<- 2000
nb<- 100
nt<- 1
nc<- 2
np <- 8 # Number of CPUs

fitSurv = jags(dat, inits, params, model.file="./Models/SurvModel.txt", 
           n.chains=nc, n.iter=ni, n.burnin=nb, 
           n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
           n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitSurv, digits=3)
plot(fitSurv)
save(fitSurv, file="fitSurv.rda")
load(file="fitSurv.rda")

ndays <- 28
ntraps <- 15
# S_t = exp(-lambda*t*v)
S_t <- exp(-fitSurv$mean$b0 * seq_len(ndays) * fitSurv$mean$v)
p <- inv.cloglog(fitSurv$mean$b0 + log(seq_len(ndays)) + log(fitSurv$mean$v))
df <- data.frame(p, days=1:58, p_nTraps=1-(1-p)^n)
ggplot(df) + geom_line(aes(days, p)) + geom_line(aes(days, Traps_15), col="red")

#### fit survival with distance ####
dat <- list(y=caphist_m, nind=length(IDs), nobs=nobs, t=time_diff, d=dist)
# arguments 
params = c("v", "b0", "bd")

inits <- function() {
  list(b0=rnorm(1, 0, 0.5),
       bd=rnorm(1, 0, 0.5),
       v=rgamma(1, 2.5, 3))
}

# call to JAGS
ni<- 200
nb<- 100
nt<- 1
nc<- 1
np <- 8 # Number of CPUs

fitSurvDist = jags(dat, inits, params, model.file="./Models/SurvModel_dist.txt", 
               n.chains=nc, n.iter=ni, n.burnin=nb, 
               n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
               n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitSurvDist, digits=3)
sapply(fitSurvDist$n.eff, summary)
fitSurvDist$Rhat
plot(fitSurvDist)

#### fit survival without distance Weibull ####
dat=list(censored=surv_df$censored, t=surv_df$TimeTrap, 
         t.cen=surv_df$t.cen, nobs=nrow(surv_df))

# arguments 
params = c("v", "b0")

inits <- function() {
  list(b0=rnorm(1, 0, 0.5),
       v=rgamma(1, 2.5, 3),
       t=surv_df$t.start)
}

# call to JAGS
ni<- 50000
nb<- 4000
nt<- 1
nc<- 2
np <- 8 # Number of CPUs

fitSurvWeib = jags(dat, inits, params, model.file="./Models/SurvModel_WeibT.txt", 
                   n.chains=nc, n.iter=ni, n.burnin=nb, 
                   n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                   n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitSurvWeib, digits=3)
plot(fitSurvWeib)
save(fitSurvWeib, file="fitSurvWeib.rda")
load(file="fitSurvWeib.rda")

#### fit survival with distance Weibull ####

dat=list(censored=1 - surv_df$NonCensored, t=surv_df$TimeTrap, 
         t.cen=surv_df$TimeTrap, nobs=nrow(surv_df), d=dist_vect)

# arguments 
params = c("v", "b0", "bd")

inits <- function() {
  list(b0=rnorm(1, 0, 0.5),
       bd=rnorm(1, 0, 0.5),
       v=rgamma(1, 2.5, 3))
}

# call to JAGS
ni<- 200
nb<- 100
nt<- 1
nc<- 1
np <- 8 # Number of CPUs

fitSurvWeibDist = jags(dat, inits, params, model.file="./Models/SurvModel_WeibT_dist.txt", 
                   n.chains=nc, n.iter=ni, n.burnin=nb, 
                   n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                   n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitSurvWeibDist, digits=3)
sapply(fitSurvWeibDist$n.eff, summary)
fitSurvWeibDist$Rhat
plot(fitSurvWeibDist)



