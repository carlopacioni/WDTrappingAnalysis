library(survival)
expmod <- survreg(Surv(surv_eff$t.comb, 1-surv_eff$censored) ~ 1, data = surv_eff, dist = "exponential")
summary(expmod)
plot(survfit(Surv(surv_eff$t.comb, 1-surv_eff$censored)~1))
1/exp(expmod$coefficients) # lambda == exp(-coeff)
mean(predict(expmod, type = "response")) # Should be == 1/lambda
