# WDTrappingAnalysis
Script to conduct analysis of simulated trapping exercise based on retrieval and identification of scats.

__*This is not intended for public use*__
## Scripts
- Functions.R This is a collection of helper functions called by the analysis scripts.
- Trapping_DataPrep.R  This collates and formats the data to be ready for the analysis.
- Trapping_Analysis_scats.R is the analysis for the scat data. It calls the initial two scripts and fit the models.
- Trapping_Analysis_scats_preliminary.R is a script with preliminary analyses (dee details below). 
- Trapping_Analysis_camTraps.R is a script to fit the models to the data from cam traps
- SimulateSurvival.R Simulate datasets with exponential and Weibull distribution and fit a model to confirm correct coding of the models.
-Fit_survival_pkg.R fit the scat data with the survival model.

## Details
Before calling `Trapping_DataPrep.R`, the analysis script sets a few variables. Most importantly, `time_period`, `ndays` and `ntraps`. Note that time is expressed as a function of `time_period`. A description of the objects created is at the beginning of the analysis script. *Note* that in `surv_df` time is expressed as ndays/time_period, while in `surv_eff` is expressed as ntrapdays/time_period.
I considered a simple Bernoulli model (p=the probability of being trapped in one day), a survival model with exponential distribution and one with Weibull distribution.
To express occasions|times, I tried two approaches: one where I have considered each trap independently, so here occasions are the traps x dates combinations. Because interaction with traps by the same dog would be expected to be correlated, I indexed the coefficient by dogs, and estimated the mean coefficient using a hyperprior. For example, for the Bernoulli model the linear equation is:
logit(p[i,k]) <- b0[i] # Where i is the dog ID and k is the trap x date
b0[i] ~ dnorm(mu.b0, tau.b0)

All these models are identify by the suffix 'ind' for "individual traps' in the name of the txt files (except for SurvModelWeibT.txt and SurvModelWeibT_dist.txt where I tried to simplify even further the model). However, the estimations with this approach using the Bernoulli model with distance or the survival models are problematic. I think that the reason is because it overestimate the effect of the traps. 

The alternative approach, which seems to work, is that time in the survival models is expressed as trap-days (scaled by `time_period`). I kept this model in the final analysis. To simplify the model and speed up convergence (and because we only have one entry per dog), I estimate directly the coefficient b0 (there is no much difference though in the two models). However, with this approach is not immediately obvious how to include a spatial component (e.g. distance of traps from activity centre).

I explicitly coded the survival model with an exponential distribution, but the Weibull model correctly handle this case by estimating v=1 when the data simulated from an exponential model are passed (indeed the Weibull model, when v=1, reduces to an exponential model).

## Directories
**Models** contains the txt file with the models

**Data** contains the raw data

**FittedModels** contains the saved models
