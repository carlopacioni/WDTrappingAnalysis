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

xls <- file.path(data_path,"TissueSamples.xlsx")
Tissue_sheet <- "Tissue"

# IDs
tissue_dt <- read_excel(xls, sheet=Tissue_sheet)


#### Mansfield ####
mfd_start_date <- strptime("2018-11-14", format="%Y-%m-%d")
mfd_nind <- 54

# subset to Mandfield
tissue_mfd <- tissue_dt[tissue_dt$Site == "Mansfield", ]

tissue_mfd$Date <- as.character(format(tissue_mfd$Date, format="%d-%m-%Y"))
tissue_mfd <- tissue_mfd[tissue_mfd$Date %in% names(trap_matrix),]
class(tissue_mfd$Date)
max.mfd.trap_eff <- check.trap_eff(trap_matrix) / (time_period * ntraps)

make.censor.table <- function(nind, tissue.table, mtrap,
                              std_trap_effort=time_period * ntraps) {
  ntrapped <- nrow(tissue.table)
  df <- data.frame(DogIDs=seq_len(nind), censored=rep(1, nind), 
                   TimeTrap=rep(NA, nind), t.cen=rep(max.mfd.trap_eff, nind),
                   t.start=rep(max.mfd.trap_eff, nind) + 1)
  nms_traps <- names(mtrap)
  for(rn in seq_len(ntrapped)) {
    ncol.keep <- which(nms_traps == tissue.table$Date[rn])
    trap_sub <- mtrap[, seq_len(ncol.keep)]
    nTrapDays <- check.trap_eff(trap_sub) / std_trap_effort
    df$censored[rn] <- 0
    df$t.start[rn] <- NA
    df$TimeTrap[rn] <- nTrapDays
  }
  return(df)
}

surv_eff_mfd_camtraps <- make.censor.table(nind=mfd_nind, tissue.table=tissue_mfd, 
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
nb<- 40000
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
p_WeibTrapEffortCamTrapMfd

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
ni<- 50000
nb<- 40000
nt<- 10
nc<- 3
np <- 8 # Number of CPUs

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

