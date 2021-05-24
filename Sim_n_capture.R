library(ggplot2)
results <- "./Results" # where the results are stored

n<-1:20
pt <- c(0.0965566668227816, 0.183790143737638)
nsim <- 1000
pTrapping_15 <- 1-(1-pt[1])^n
pTrapping_30 <- 1-(1-pt[2])^n
dt <- data.frame(n=c(n,n), 
                 Probability=c(pTrapping_15, pTrapping_30), 
                 Effort=factor(c(rep(15, length(n)), rep(30, length(n)))))
p_one<-ggplot(dt, aes(n, Probability, col=Effort)) + geom_line() + xlab("Number of wild dogs")
ggsave(file.path(results,"plot_ProbofTrappingOneDog.png"), 
       plot = p_one, width = 13, height = 10, units = "cm")



ncapt <- mapply(rbinom, size=n, prob=pt[2], MoreArgs = list(n=nsim))
Mean <- colMeans(ncapt)
Sd <- apply(ncapt, 2, sd)
Se <- (1.96 * Sd^2)/sqrt(nsim)
dt_ncapture <- data.frame(n, Mean, Lower=Mean - Se, Upper=Mean + Se)
ggplot(dt_ncapture, aes(n, Mean)) + geom_line() + 
  geom_errorbar(aes(ymin = Lower, ymax=Upper))
