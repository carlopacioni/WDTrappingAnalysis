library(jagsUI)
library(ggplot2)
library(readxl)
library(rgdal)
# time_period is used to rescale time periods to a similar range as other covariate (distance) and model parameters. 
# This is calculated as number of days since the begininning of trapping / time_period
# This is modtly because otherwise ndays would reange from 1 to 58.
time_period <- 7

ndays <- 28 # Usual number of days traps are deployed for
ntraps <- 15 # Usual number of traps deployed
fittedMods <- "./FittedModels" # where to save fitted models
data_path <- "./Data" # where the data are stored

# Load helper functions
source("Functions.R")
# Process data to create input objects
source("Trapping_DataPrep.R")



#------------------------------------------------------------------------------#
#   Description of objects used to fit models                                  #

#   IDs char vector with animal IDs identified from scats
#   caphist_m matrix with 0 or 1 indicating capture of individuals with traps
#        rows are individuals, cols each day across traps (e.g. trap01_day01, trap02_day01, ...)
#   time_diff matrix as caphist_m but the the time difference (ndays/time_period) 
#        between the first detection of an animal and the day being considered
#   dist matrix as caphist_m but the the distance (in km)  
#        between the location first detection of an animal and the trap in the matrix
#   nobs is the number of observation in the matrix that should be considered (i.e. before NA)
#   surv_df data frame with one row for each animal and trap.
#      cols: "ID"  char with animal IDs     
#            "ind"  animal ID as factor    
#            "Trap"   char trap ID
#            "censored" 0 if censored, 1 if trapped
#            "TimeTrap" time when trapped (if not rtrrapped NA) calculated as ndays/time_period
#            "t.cen"  censoring time. If trapped NA
#            "t.start"  starting value for censored rows
#            "t.comb" censoring or trapping time combined in one column


#### Bernoulli model ####
#### fit intercept only ####
dat <- list(y=caphist_m, nind=length(IDs), nobs=nobs)
# arguments 
params = c("mu.b0", "sigma.b0", "b0")

inits <- function() {
  list(mu.b0=rnorm(1, -3, 0.5),
       sigma.b0=runif(1, 0.2, 0.75))
}

# call to JAGS
ni<- 20000
nb<- 1000
nt<- 10
nc<- 5
np <- 8 # Number of CPUs

fitBern_indb0 = jags(dat, inits, params, model.file="./Models/SimpleBernModel_indb0.txt", 
              n.chains=nc, n.iter=ni, n.burnin=nb, 
              n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
              n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitBern_indb0, digits=3)
inv.logit(fitBern_indb0$mean$mu.b0)
1-(1-inv.logit(fitBern_indb0$mean$mu.b0))^(ndays*ntraps) # Prob for ntraps deployed ndays
plot(fitBern_indb0)
save(fitBern_indb0, file = file.path(fittedMods, "fitBern_indb0.rda"))
load(file.path(fittedMods, "fitBern_indb0.rda"))

ggplot(data.frame(mu.b0=fitBern_indb0$sims.list$mu.b0, 
                  Iter=seq_along(fitBern_indb0$sims.list$mu.b0))) + 
  geom_line(aes(Iter, mu.b0))

#### fit with covariates (ndays and distance) TO-DO #### 
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

#### fit with distance TO-DO ####
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
dat <- list(y=surv_df$censored, nind=length(IDs), t=surv_df$t.comb, 
            nobs=nrow(surv_df), ind=surv_df$ind)
# arguments 
params = c("v", "b0", "mu.b0", "sigma.b0")

inits <- function() {
  list(mu.b0=rnorm(1, 0, 0.5),
       sigma.b0=runif(1, 0.2, 0.75),
       v=rgamma(1, 2.5, 3))
}

# call to JAGS
ni<- 50000
nb<- 4000
nt<- 10
nc<- 4
np <- 8 # Number of CPUs

fitSurv_indb0 = jags(dat, inits, params, model.file="./Models/SurvModel_indb0.txt", 
           n.chains=nc, n.iter=ni, n.burnin=nb, 
           n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
           n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitSurv_indb0, digits=3)
sapply(fitSurv_indb0$n.eff, summary)
plot(fitSurv_indb0, parameters=c("v", "mu.b0", "sigma.b0"))

save(fitSurv_indb0, file=file.path(fittedMods, "fitSurv_indb0.rda"))
load(file=file.path(fittedMods, "fitSurv_indb0.rda"))


# S_t = exp(-lambda*t^v)
S_t <- exp(-exp(fitSurv_indb0$mean$mu.b0) * (seq_len(ndays)/time_period)^fitSurv_indb0$mean$v)
p <- inv.cloglog(fitSurv_indb0$mean$mu.b0 + log((seq_len(ndays)/time_period)) * fitSurv_indb0$mean$v)
df <- data.frame(p, days=seq_len(ndays), p_nTraps=1-(1-p)^ntraps)
ggplot(df) + geom_line(aes(days, p)) + geom_line(aes(days, p_nTraps), col="red")

#### fit survival with distance TO-DO ####
dat <- list(y=caphist_m, nind=length(IDs), nobs=nobs, t=time_diff, d=dist)
# arguments 
params = c("v", "b0", "bd")

inits <- function() {
  list(b0=rnorm(1, 0, 0.5),
       v=rgamma(1, 2.5, 3))
}

# call to JAGS
ni<- 200
nb<- 100
nt<- 1
nc<- 1
np <- 8 # Number of CPUs

fitSurvDist_indb0 = jags(dat, inits, params, model.file="./Models/SurvModel_dist.txt", 
               n.chains=nc, n.iter=ni, n.burnin=nb, 
               n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
               n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitSurvDist_indb0, digits=3)
sapply(fitSurvDist_indb0$n.eff, summary)
fitSurvDist_indb0$Rhat
plot(fitSurvDist_indb0)
save(fitSurvDist_indb0, file = file.path(fittedMods, "fitSurvDist_indb0.rda"))
load(file.path(fittedMods, "fitSurvDist_indb0.rda"))

#### fit survival without distance Weibull ####
dat=list(censored=surv_df$censored, t=surv_df$TimeTrap, 
         t.cen=surv_df$t.cen, nobs=nrow(surv_df), ind=surv_df$ind, nind=length(IDs))

# arguments 
params = c("v", "mu.b0", "sigma.b0", "b0")

inits <- function() {
  list(mu.b0=rnorm(1, -3, 0.5),
       sigma.b0=runif(1, 0.2, 0.75),
       v=rgamma(1, 2.5, 3),
       t=surv_df$t.start)
}

# call to JAGS
ni<- 50000
nb<- 4000
nt<- 1
nc<- 4
np <- 8 # Number of CPUs

fitSurvWeib_indb0 = jags(dat, inits, params, model.file="./Models/SurvModel_WeibT_ind.txt", 
                   n.chains=nc, n.iter=ni, n.burnin=nb, 
                   n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                   n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitSurvWeib_indb0, digits=3)
fitSurvWeib_indb0$n.eff
plot(fitSurvWeib_indb0)
save(fitSurvWeib_indb0, file=file.path(fittedMods, "fitSurvWeib_indb0.rda"))
load(file=file.path(fittedMods, "fitSurvWeib_indb0.rda"))

# S_t = exp(-lambda*t^v)
S_t_Weib <- exp(-exp(fitSurvWeib_indb0$mean$mu.b0) * (seq_len(ndays)/time_period)^fitSurvWeib_indb0$mean$v)
p_Weib <- inv.cloglog(fitSurvWeib_indb0$mean$mu.b0 + log((seq_len(ndays)/time_period)) * fitSurvWeib_indb0$mean$v)
df_Weib <- data.frame(p, days=seq_len(ndays), p_nTraps=1-(1-p)^ntraps)
ggplot(df_Weib) + geom_line(aes(days, p)) + geom_line(aes(days, p_nTraps), col="red")

#### fit survival with distance Weibull TO-DO ####

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



