The provided R code applies very flexible piecewise linear function to estimate the baseline 
hazard function under Box-Cox transformation cure model set-up to estimate parameter and 
associated standard errors in the case where data collected in an interval-censored manner. 
The estimates can be used for predicting cure rates (probabilities) in time-to-event study marked 
by cured proportion, thereby, help users understand treatment efficacy. 

The provided R code contains two parts:

1. Data_Gen.R: 
   This set of R code correspond to the generation of simulated data from parameter settings 
   as described in section 4 of the main article.
    
2. Main_EM_Optim.R 
   This set of R code correspond to the implementation of the EM algorithm and profile likelihood 
   method for $\alpha$, as described in section 4 of the main article, for estimating parameters 
   and standard errors associated with them. 

Note that, both R codes are sufficiently commented to help users understand the intent of each portion.
Users should first generate data using Data_Gen.R and apply Main_EM_Optim.R data on the generated data 
for parameter and cure rate estimation.       