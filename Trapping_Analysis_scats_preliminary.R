library(jagsUI)
library(ggplot2)
library(readxl)
library(rgdal)
library(data.table)
# time_period is used to rescale time periods to a similar range as other covariate (distance) and model parameters. 
# This is calculated as number of days since the begininning of trapping / time_period
# This is modtly because otherwise ndays would reange from 1 to 58.
time_period <- 7

ndays <- 28 # Usual number of days traps are deployed for
ntraps <- 15 # Usual number of traps deployed
fittedMods <- "./FittedModels" # where to save fitted models
data_path <- "./Data" # where the data are stored
results <- "./Results" # where the results are stored
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
#            "censored" 1 if censored, 0 if trapped
#            "TimeTrap" time when trapped (if not trapped NA) calculated as ndays/time_period
#            "t.cen"  censoring time. If trapped NA
#            "t.start"  starting value for censored rows
#            "t.comb" censoring or trapping time combined in one column
#   surv_eff data frame with one row for each animal 
#      cols: "ID"  char with animal IDs     
#            "censored" 1 if censored, 0 if trapped
#            "TimeTrap" time when trapped (if not trapped NA) calculated as ntrapdays/time_period
#            "t.cen"  censoring time (max time). If trapped NA
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

fitBern_indb0 = jags(dat, inits, params, model.file="./Models/BernModel_indb0.txt", 
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


#### fit with distance ####
dat <- list(y=caphist_m, nind=length(IDs), nobs=nobs, d=dist)
# arguments 
params = c("mu.b0", "sigma.b0", "mu.bd", "sigma.bd")

inits <- function() {
  list(mu.b0=rnorm(1, -3, 0.5),
       sigma.b0=runif(1, 0.2, 0.75),
       mu.bd=rnorm(1, -3, 0.5),
       sigma.bd=runif(1, 0.2, 0.75)
       )
}

# call to JAGS
ni<- 12000
nb<- 8000
nt<- 1
nc<- 3
np <- 8 # Number of CPUs

fitBernDist_ind = jags(dat, inits, params, model.file="./Models/BernModel_dist_ind.txt", 
              n.chains=nc, n.iter=ni, n.burnin=nb, 
              n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
              n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitBernDist_ind, digits=3)
sapply(fitBernDist_ind$n.eff, summary)
fitBernDist_ind$Rhat
plot(fitBernDist_ind)
inv.logit(fitBernDist_ind$mean$mu.b0 + fitBernDist_ind$mean$mu.bd * mean(dist))
1-(1-inv.logit(fitBernDist_ind$mean$mu.b0 + fitBernDist_ind$mean$mu.bd * mean(dist)))^(ndays*ntraps) # Prob for ntraps deployed ndays

#### Survivival ####



#### fit survival Weibull Individual traps ####
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
df_Weib <- data.frame(p=p_Weib, days=seq_len(ndays), p_nTraps=1-(1-p_Weib)^ntraps)
ggplot(df_Weib) + geom_line(aes(days, p)) + geom_line(aes(days, p_nTraps), col="red")

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

dat=list(censored=1 - surv_df$censored, t=surv_df$TimeTrap, 
         t.cen=surv_df$t.cen, nobs=nrow(surv_df), d=dist_vect)

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


#### fit survival  Weibull trap effort ####
dat=list(censored=surv_eff$censored, t=surv_eff$TimeTrap, 
         t.cen=surv_eff$t.cen, nind=length(IDs))

# arguments 
params = c("v", "b0")

inits <- function() {
  list(b0=rnorm(1, 0.1, 0.5),
       v=rgamma(1, 2.5, 3),
       t=surv_eff$t.start)
}

# call to JAGS
ni<- 60000
nb<- 40000
nt<- 10
nc<- 3
np <- 8 # Number of CPUs

fitSurvWeibTrapEff = jags(dat, inits, params, model.file="./Models/SurvModel_WeibTrap_eff.txt", 
                         n.chains=nc, n.iter=ni, n.burnin=nb, 
                         n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                         n.cores=ifelse(floor(nc/np) < np, nc, np))

save(fitSurvWeibTrapEff, file=file.path(fittedMods, "fitSurvWeibTrapEff.rda"))
load(file=file.path(fittedMods, "fitSurvWeibTrapEff.rda"))

print(fitSurvWeibTrapEff, digits=3)
fitSurvWeibTrapEff$mean$v
fitSurvWeibTrapEff$mean$b0
exp(fitSurvWeibTrapEff$mean$b0) # lambda
#plot(fitSurvWeibTrapEff)

p_WeibTrapEffort <- plot_probs(lam=exp(fitSurvWeibTrapEff$mean$b0), max.ndays=28, 
                               v=fitSurvWeibTrapEff$mean$v, 
                              time.period=time_period, n.traps=ntraps)
p_WeibTrapEffort

ggsave(file.path(results, "plot_WeibTrapEffort.pdf"), plot = p_ExpTrapEffort)

OneDayOnetrap <- 1/(time_period * ntraps)
pOneTrapOneDay <- 1 - S_t(lam = exp(fitSurvWeibTrapEff$mean$b0), t = OneDayOnetrap, v = fitSurvWeibTrapEff$mean$v)
pOneTrapOneDay

#### fit survival  Weibull trap effort individual dog coefs ####
params = c("v", "b0", "mu.b0", "sigma.b0")

inits <- function() {
  list(mu.b0=rnorm(1, 0.1, 0.5),
       sigma.b0=runif(1, 0.2, 0.75),
       v=rgamma(1, 2.5, 3),
       t=surv_eff$t.start)
}

# call to JAGS
ni<- 600000
nb<- 40000
nt<- 10
nc<- 3
np <- 8 # Number of CPUs

fitSurvWeibTrapEff_indDogCoefs = jags(dat, inits, params, model.file="./Models/SurvModel_WeibTrap_eff_indDogCoef.txt", 
                          n.chains=nc, n.iter=ni, n.burnin=nb, 
                          n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                          n.cores=ifelse(floor(nc/np) < np, nc, np))

save(fitSurvWeibTrapEff_indDogCoefs, file=file.path(fittedMods, "fitSurvWeibTrapEff.rda"))
load(file=file.path(fittedMods, "fitSurvWeibTrapEff_indDogCoefs.rda"))

print(fitSurvWeibTrapEff_indDogCoefs, digits=3)
plot(fitSurvWeibTrapEff_indDogCoefs)
fitSurvWeibTrapEff_indDogCoefs$mean$v
fitSurvWeibTrapEff_indDogCoefs$mean$b0
exp(fitSurvWeibTrapEff_indDogCoefs$mean$b0) # lambda
#plot(fitSurvWeibTrapEff)

p_fitSurvWeibTrapEff_indDogCoefs <- plot_probs(lam=exp(fitSurvWeibTrapEff_indDogCoefs$mean$b0), max.ndays=28, 
                               v=fitSurvWeibTrapEff_indDogCoefs$mean$v, 
                               time.period=time_period, n.traps=ntraps)
p_fitSurvWeibTrapEff_indDogCoefs


#### fit survival  Exp trap effort ####
dat=list(censored=surv_eff$censored, t=surv_eff$TimeTrap, 
         t.cen=surv_eff$t.cen, nind=length(IDs))

# arguments 
params = c("b0")

inits <- function() {
  list(b0=rnorm(1, 0.1, 0.5),
       t=surv_eff$t.start)
}

# call to JAGS
ni<- 50000
nb<- 40000
nt<- 10
nc<- 3
np <- 8 # Number of CPUs

fitSurvExpTrapEff = jags(dat, inits, params, model.file="./Models/SurvModel_ExpTrap_eff.txt", 
                                n.chains=nc, n.iter=ni, n.burnin=nb, 
                                n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                                n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitSurvExpTrapEff, digits=3)
fitSurvExpTrapEff$mean$b0
exp(fitSurvExpTrapEff$mean$b0) # lambda
#plot(fitSurvExpTrapEff)
save(fitSurvExpTrapEff, file=file.path(fittedMods, "fitSurvExpTrapEff.rda"))
load(file=file.path(fittedMods, "fitSurvExpTrapEff.rda"))

p_ExpTrapEffort <- plot_probs(lam=exp(fitSurvExpTrapEff$mean$b0), max.ndays=28, v=1, 
           time.period=time_period, n.traps=ntraps)
p_ExpTrapEffort

ggsave(file.path(results,"plot_ExpTrapEffort.pdf"), plot = p_ExpTrapEffort)

fitSurvExpTrapEff$DIC - fitSurvWeibTrapEff$DIC

