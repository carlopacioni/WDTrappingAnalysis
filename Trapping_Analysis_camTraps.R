library(jagsUI)
library(ggplot2)
library(readxl)
library(rgdal)
library(data.table)

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

# The two data.frame generated below (surv_eff_sck_camtraps, surv_eff_sck_camtraps) 
# ara data frames with one row for each animal and have the following columns:
#   
#      cols: "DogIDs" Progressing number as animal IDs     
#            "censored" 1 if censored, 0 if trapped
#            "TimeTrap" time when trapped (if not trapped NA) calculated as ndays/time_period
#            "t.cen"  censoring time. If trapped NA
#            "t.start"  starting value for censored rows

#### Mansfield ####
mfd_nind <- 54 # number of individuals estimated from camtrap data

# generate input data to fit model
surv_eff_mfd_camtraps <- make.censor.table(nind=mfd_nind, tissue.table=tissue_mfd, 
                                           max.trap_eff = max.mfd.trap_eff,
                                           mtrap=trap_matrix)
surv_eff_mfd_camtraps

#### fit survival  Weibull trap effort ####
dat=list(censored=surv_eff_mfd_camtraps$censored, t=surv_eff_mfd_camtraps$TimeTrap, 
         t.cen=surv_eff_mfd_camtraps$t.cen, nind=mfd_nind)

# arguments 
params=c("v", "b0")

inits <- function() {
  list(b0=rnorm(1, 0.1, 0.5),
       v=rgamma(1, 2.5, 3),
       t=surv_eff_mfd_camtraps$t.start)
}

# call to JAGS
ni<- 60000
nb<- 20000
nt<- 10
nc<- 3
np <- 8 # Number of CPUs

fitSurvWeibTrapEffCamTrapMfd=jags(dat, inits, params, model.file="./Models/SurvModel_WeibTrap_eff.txt", 
                          n.chains=nc, n.iter=ni, n.burnin=nb, 
                          n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                          n.cores=ifelse(floor(nc/np) < np, nc, np))

save(fitSurvWeibTrapEffCamTrapMfd, file=file.path(fittedMods, "fitSurvWeibTrapEffCamTrap.rda"))
load(file=file.path(fittedMods, "fitSurvWeibTrapEffCamTrapMfd.rda"))

print(fitSurvWeibTrapEffCamTrapMfd, digits=3)
fitSurvWeibTrapEffCamTrapMfd$mean$v
fitSurvWeibTrapEffCamTrapMfd$mean$b0
exp(fitSurvWeibTrapEffCamTrapMfd$mean$b0) # lambda
#plot(fitSurvWeibTrapEffCamTrap)

p_WeibTrapEffortCamTrapMfd <- plot_probs(lam=exp(fitSurvWeibTrapEffCamTrapMfd$mean$b0), 
                                      max.ndays=28, 
                               v=fitSurvWeibTrapEffCamTrapMfd$mean$v, 
                               time.period=time_period, n.traps=ntraps)
p_WeibTrapEffortCamTrapMfdHPD <- plot_probsHPD(fittedMod = fitSurvWeibTrapEffCamTrapMfd, 
                                         max.ndays=28, 
                                         v="mod", 
                                         time.period=time_period, n.traps=ntraps)
p_WeibTrapEffortCamTrapMfd
p_WeibTrapEffortCamTrapMfdHPD

ggsave(file.path(results, "plot_WeibTrapEffortCamTrapMfd.pdf"), 
       plot=p_WeibTrapEffortCamTrapMfd, width = 18, height = 15, units = "cm")
ggsave(file.path(results, "plot_WeibTrapEffortCamTrapMfd.png"), 
       plot=p_WeibTrapEffortCamTrapMfd, width = 18, height = 15, units = "cm")

OneDayOnetrap <- 1/(time_period * ntraps)
pOneTrapOneDay <- 1 - S_t(lam=exp(fitSurvWeibTrapEffCamTrapMfd$mean$b0), 
                          t=OneDayOnetrap, v=fitSurvWeibTrapEffCamTrapMfd$mean$v)
pOneTrapOneDay

#### fit survival Exp trap effort ####
#Note that this could be fitted with Weinbull and fixing v=1
dat=list(censored=surv_eff_mfd_camtraps$censored, t=surv_eff_mfd_camtraps$TimeTrap, 
         t.cen=surv_eff_mfd_camtraps$t.cen, nind=mfd_nind)

# arguments 
params = c("b0")

inits <- function() {
  list(b0=rnorm(1, 0.1, 0.5),
       t=surv_eff_mfd_camtraps$t.start)
}

# call to JAGS

fitSurvExpTrapEffCamTrapsMfd = jags(dat, inits, params, model.file="./Models/SurvModel_ExpTrap_eff.txt", 
                         n.chains=nc, n.iter=ni, n.burnin=nb, 
                         n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                         n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitSurvExpTrapEffCamTrapsMfd, digits=3)
fitSurvExpTrapEffCamTrapsMfd$mean$b0
exp(fitSurvExpTrapEffCamTrapsMfd$mean$b0) # lambda
#plot(fitSurvExpTrapEff)
save(fitSurvExpTrapEffCamTrapsMfd, file=file.path(fittedMods, "fitSurvExpTrapEffCamTrapMfd.rda"))
load(file=file.path(fittedMods, "fitSurvExpTrapEffCamTrapMfd.rda"))

p_ExpTrapEffortCamTrapMfd <- plot_probs(lam=exp(fitSurvExpTrapEffCamTrapsMfd$mean$b0), 
                                        max.ndays=28, v=1, 
                              time.period=time_period, n.traps=ntraps) + ylim(c(0, 0.18))
p_ExpTrapEffortCamTrapMfd 

ggsave(file.path(results,"plot_ExpTrapEffortCamTrapMfd.pdf"), plot = p_ExpTrapEffortCamTrapMfd)
ggsave(file.path(results,"plot_ExpTrapEffortCamTrapMfd.png"), 
       plot = p_ExpTrapEffortCamTrapMfd, width = 18, height = 15, units = "cm")

fitSurvExpTrapEffCamTrapsMfd$DIC - fitSurvWeibTrapEffCamTrapMfd$DIC

#### Swifts creek ####
sck_nind <- 27 # number of individuals estimated from camtrap data

# generate input data to fit model
surv_eff_sck_camtraps <- make.censor.table(nind=sck_nind, tissue.table=tissue_sck, 
                                           max.trap_eff = max.sck.trap_eff,
                                           mtrap=trap_matrix_sck)
surv_eff_sck_camtraps

#### fit survival  Weibull trap effort ####
dat=list(censored=surv_eff_sck_camtraps$censored, t=surv_eff_sck_camtraps$TimeTrap, 
         t.cen=surv_eff_sck_camtraps$t.cen, nind=sck_nind)

# arguments 
params=c("v", "b0")

inits <- function() {
  list(b0=rnorm(1, 0.1, 0.5),
       v=rgamma(1, 2.5, 3),
       t=surv_eff_sck_camtraps$t.start)
}

# call to JAGS
ni<- 60000
nb<- 40000
nt<- 10
nc<- 3
np <- 8 # Number of CPUs

fitSurvWeibTrapEffCamTrapSck=jags(dat, inits, params, model.file="./Models/SurvModel_WeibTrap_eff.txt", 
                                  n.chains=nc, n.iter=ni, n.burnin=nb, 
                                  n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                                  n.cores=ifelse(floor(nc/np) < np, nc, np))

save(fitSurvWeibTrapEffCamTrapSck, file=file.path(fittedMods, "fitSurvWeibTrapEffCamTrapSck.rda"))
load(file=file.path(fittedMods, "fitSurvWeibTrapEffCamTrapSck.rda"))

print(fitSurvWeibTrapEffCamTrapSck, digits=3)
fitSurvWeibTrapEffCamTrapSck$mean$v
fitSurvWeibTrapEffCamTrapSck$mean$b0
exp(fitSurvWeibTrapEffCamTrapSck$mean$b0) # lambda
#plot(fitSurvWeibTrapEffCamTrapSck)

p_fitSurvWeibTrapEffCamTrapSck <- plot_probs(lam=exp(fitSurvWeibTrapEffCamTrapSck$mean$b0), 
                                         max.ndays=28, 
                                         v=fitSurvWeibTrapEffCamTrapSck$mean$v, 
                                         time.period=time_period, n.traps=ntraps)
p_fitSurvWeibTrapEffCamTrapSck

ggsave(file.path(results, "plot_WeibTrapEffortCamTrapSck.pdf"), 
       plot=p_fitSurvWeibTrapEffCamTrapSck, width = 18, height = 15, units = "cm")
ggsave(file.path(results, "plot_WeibTrapEffortCamTrapSck.png"), 
       plot=p_fitSurvWeibTrapEffCamTrapSck, width = 18, height = 15, units = "cm")

OneDayOnetrap <- 1/(time_period * ntraps)
pOneTrapOneDay <- 1 - S_t(lam=exp(fitSurvWeibTrapEffCamTrapSck$mean$b0), 
                          t=OneDayOnetrap, v=fitSurvWeibTrapEffCamTrapSck$mean$v)
pOneTrapOneDay

#### fit survival Exp trap effort ####
#Note that this could be fitted with Weinbull and fixing v=1
dat=list(censored=surv_eff_sck_camtraps$censored, t=surv_eff_sck_camtraps$TimeTrap, 
         t.cen=surv_eff_sck_camtraps$t.cen, nind=sck_nind)

# arguments 
params = c("b0")

inits <- function() {
  list(b0=rnorm(1, 0.1, 0.5),
       t=surv_eff_sck_camtraps$t.start)
}

fitSurvExpTrapEffCamTrapsSck = jags(dat, inits, params, model.file="./Models/SurvModel_ExpTrap_eff.txt", 
                                    n.chains=nc, n.iter=ni, n.burnin=nb, 
                                    n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                                    n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitSurvExpTrapEffCamTrapsSck, digits=3)
fitSurvExpTrapEffCamTrapsSck$mean$b0
exp(fitSurvExpTrapEffCamTrapsSck$mean$b0) # lambda
#plot(fitSurvExpTrapEffCamTrapsSck)
save(fitSurvExpTrapEffCamTrapsSck, file=file.path(fittedMods, "fitSurvExpTrapEffCamTrapsSck.rda"))
load(file=file.path(fittedMods, "fitSurvExpTrapEffCamTrapsSck.rda"))

p_ExpTrapEffCamTrapsSck <- plot_probs(lam=exp(fitSurvExpTrapEffCamTrapsSck$mean$b0), 
                                        max.ndays=28, v=1, 
                                        time.period=time_period, n.traps=ntraps) #+ ylim(c(0, 0.18))
p_ExpTrapEffCamTrapsSck 

ggsave(file.path(results,"plot_ExpTrapEffortCamTrapSck.pdf"), plot = p_ExpTrapEffCamTrapsSck)
ggsave(file.path(results,"plot_ExpTrapEffortCamTrapSck.png"), 
       plot = p_ExpTrapEffCamTrapsSck, width = 18, height = 15, units = "cm")

fitSurvExpTrapEffCamTrapsSck$DIC - fitSurvWeibTrapEffCamTrapSck$DIC



#### Combined data ####

surv_eff_comb_camtraps <- rbind(surv_eff_mfd_camtraps, surv_eff_sck_camtraps)
surv_eff_comb_camtraps$t.cen <- max(surv_eff_comb_camtraps$t.cen)
surv_eff_comb_camtraps$t.start <- max(surv_eff_comb_camtraps$t.cen) + 1
dat=list(censored=surv_eff_comb_camtraps$censored, t=surv_eff_comb_camtraps$TimeTrap, 
         t.cen=surv_eff_comb_camtraps$t.cen, nind=mfd_nind + sck_nind)

# arguments 
params = c("b0")

inits <- function() {
  list(b0=rnorm(1, 0.1, 0.5),
       t=surv_eff_comb_camtraps$t.start)
}

fitSurvExpTrapEffCamTrapsCbn = jags(dat, inits, params, model.file="./Models/SurvModel_ExpTrap_eff.txt", 
                                    n.chains=nc, n.iter=ni, n.burnin=nb, 
                                    n.thin=nt, parallel=ifelse(nc>1, TRUE, FALSE), 
                                    n.cores=ifelse(floor(nc/np) < np, nc, np))

print(fitSurvExpTrapEffCamTrapsCbn, digits=3)
fitSurvExpTrapEffCamTrapsCbn$mean$b0
exp(fitSurvExpTrapEffCamTrapsCbn$mean$b0) # lambda
#plot(fitSurvExpTrapEffCamTrapsCbn)
save(fitSurvExpTrapEffCamTrapsCbn, file=file.path(fittedMods, "fitSurvExpTrapEffCamTrapsCbn.rda"))
load(file=file.path(fittedMods, "fitSurvExpTrapEffCamTrapsCbn.rda"))

p_ExpTrapEffCamTrapsCbn <- plot_probs(lam=exp(fitSurvExpTrapEffCamTrapsCbn$mean$b0), 
                                      max.ndays=28, v=1, 
                                      time.period=time_period, n.traps=ntraps) #+ ylim(c(0, 0.18))
p_ExpTrapEffCamTrapsCbn 

ggsave(file.path(results,"plot_ExpTrapEffortCamTrapCbn.pdf"), plot = p_ExpTrapEffCamTrapsCbn)
ggsave(file.path(results,"plot_ExpTrapEffortCamTrapCbn.png"), 
       plot = p_ExpTrapEffCamTrapsCbn, width = 18, height = 15, units = "cm")

p_WeibTrapEffortCamTrapCbnHPD <- plot_probsHPD(fittedMod = fitSurvExpTrapEffCamTrapsCbn, 
                                               max.ndays=28, 
                                               v=1, 
                                               time.period=time_period, n.traps=ntraps)
p_WeibTrapEffortCamTrapCbnHPD

ggsave(file.path(results,"plot_ExpTrapEffortCamTrapCbnHPD.pdf"), plot = p_WeibTrapEffortCamTrapCbnHPD)
ggsave(file.path(results,"plot_ExpTrapEffortCamTrapCbnHPD.png"), 
       plot = p_WeibTrapEffortCamTrapCbnHPD, width = 18, height = 15, units = "cm")
