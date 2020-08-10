calc.latlong.dist<- function(xy1,xy2) {
  # uses spherical law of cosines to calculate distance between two lat/long
  # coordinates in decimal degrees
  R <- 6371 # Earths radius
  xy1 <- (pi * xy1)/180 # radians
  xy2 <- (pi * xy2)/180 
  D <- acos(sin(xy1[,1])*sin(xy2[,1]) + cos(xy1[,1])*cos(xy2[,1])*cos(xy2[,2]-xy1[,2]))  
  return(R*D)
}

inv.logit <- function(x) { exp(x) / (1+exp(x))}

inv.cloglog <- function(x) {1-exp(-exp(x))}

check.trap_eff <- function(x) {
  check.zeros <- function(x) {x == 0}
  suppressWarnings(
  m <- apply(x, 2, as.numeric)
  )
  m <- apply(m, 2, check.zeros)
  return(sum(m, na.rm = TRUE))
}

make.censor.table <- function(nind, tissue.table, mtrap, max.trap_eff,
                              std_trap_effort=time_period * ntraps) {
  ntrapped <- nrow(tissue.table)
  df <- data.frame(DogIDs=seq_len(nind), censored=rep(1, nind), 
                   TimeTrap=rep(NA, nind), t.cen=rep(max.trap_eff, nind),
                   t.start=rep(max.trap_eff, nind) + 1)
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

S_t <- function(lam, t, v) {
  return(exp(-lam*t^v))
}


# S_t = exp(-lambda*t^v) | when v=1 is an exponential model
plot_probs <- function(lam, max.ndays=28, v=1, time.period, n.traps) {
  day.1 <- 1/time.period
  day.max <- max.ndays / time.period
  s_t.1trap <- S_t(lam=lam, t=seq(day.1/n.traps, day.max/n.traps, length.out = max.ndays), v=v)
  s_t.ntraps <- S_t(lam=lam, t=seq(day.1, day.max, length.out = max.ndays), v=v)
  s_t.doubletraps <- S_t(lam=lam, t=seq(2*day.1, 2*day.max, length.out = max.ndays), v=v)
  df <- data.frame(Probability=c(1-s_t.1trap, 1-s_t.ntraps, 1-s_t.doubletraps), 
                   Days=rep(1:max.ndays, 3), Effort=factor(rep(c(1, ntraps, 2*ntraps), each=max.ndays)))
  plot <- ggplot(df) + geom_line(aes(Days, Probability, col=Effort))
  return(plot)
}

# v can be either 1 for exponential models or "mod" which uses the parameter estimates from the model
plot_probsHPD <- function(fittedMod, max.ndays=28, v=c(1, "mod"), time.period, n.traps) {
  day.1 <- 1/time.period
  day.max <- max.ndays / time.period
  s_t.1trap <- S_t(lam=exp(fittedMod$mean$b0), 
                   t=seq(day.1/n.traps, day.max/n.traps, 
                   length.out = max.ndays), 
                   v=if(v == 1) 1 else fittedMod$mean$v)
  s_t.1trap_u <- S_t(lam=exp(fittedMod$q2.5$b0), 
                   t=seq(day.1/n.traps, day.max/n.traps, 
                         length.out = max.ndays), 
                   v=if(v == 1) 1 else fittedMod$q2.5$v)
  s_t.1trap_l <- S_t(lam=exp(fittedMod$q97.5$b0), 
                   t=seq(day.1/n.traps, day.max/n.traps, 
                         length.out = max.ndays), 
                   v=if(v == 1) 1 else fittedMod$q97.5$v)
  s_t.ntraps <- S_t(lam=exp(fittedMod$mean$b0), 
                    t=seq(day.1, day.max, length.out = max.ndays), 
                    v=if(v == 1) 1 else fittedMod$mean$v)
  s_t.ntraps_u <- S_t(lam=exp(fittedMod$q2.5$b0), 
                    t=seq(day.1, day.max, length.out = max.ndays), 
                    v=if(v == 1) 1 else fittedMod$q2.5$v)
  s_t.ntraps_l <- S_t(lam=exp(fittedMod$q97.5$b0), 
                    t=seq(day.1, day.max, length.out = max.ndays), 
                    v=if(v == 1) 1 else fittedMod$q97.5$v)
  s_t.doubletraps <- S_t(lam=exp(fittedMod$mean$b0), 
                         t=seq(2*day.1, 2*day.max, length.out = max.ndays), 
                         v=if(v == 1) 1 else fittedMod$mean$v)
  s_t.doubletraps_u <- S_t(lam=exp(fittedMod$q2.5$b0), 
                          t=seq(2*day.1, 2*day.max, length.out = max.ndays), 
                         v=if(v == 1) 1 else fittedMod$q2.5$v)
  s_t.doubletraps_l <- S_t(lam=exp(fittedMod$q97.5$b0), 
                         t=seq(2*day.1, 2*day.max, length.out = max.ndays), 
                         v=if(v == 1) 1 else fittedMod$q97.5$v)
  df <- data.frame(Probability=c(1-s_t.1trap, 1-s_t.ntraps, 1-s_t.doubletraps
                                 ), 
                   p_l=c(1-s_t.1trap_l, 1-s_t.ntraps_l, 1-s_t.doubletraps_l
                         ),
                   p_u=c(1-s_t.1trap_u, 1-s_t.ntraps_u, 1-s_t.doubletraps_u
                         ),
                   Days=rep(1:max.ndays, 3), 
                   Effort=factor(rep(c(1, ntraps, 2*ntraps), each=max.ndays)))
  plot <- ggplot(df) + geom_line(aes(Days, Probability, col=Effort)) +
    geom_ribbon(aes(ymin=p_l, ymax=p_u, x=Days, col=Effort, fill=Effort), alpha=0.2)
  return(plot)
}
