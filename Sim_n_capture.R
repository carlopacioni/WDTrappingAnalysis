library(ggplot2)

n<-1:10
pt<-0.2
nsim <- 1000
pTrapping <- 1-(1-pt)^n
dt <- data.frame(n, pTrapping)
ggplot(dt, aes(n, pTrapping)) + geom_line()
ncapt <- mapply(rbinom, size=n, prob=pt, MoreArgs = list(n=nsim))
Mean <- colMeans(ncapt)
Sd <- apply(ncapt, 2, sd)
Se <- (1.96 * Sd^2)/sqrt(nsim)
dt_ncapture <- data.frame(n, Mean, Lower=Mean - Se, Upper=Mean + Se)
ggplot(dt_ncapture, aes(n, Mean)) + geom_line() + 
  geom_errorbar(aes(ymin = Lower, ymax=Upper))
