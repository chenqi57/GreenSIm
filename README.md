# GreenSim
Nested Green Simulation Algorithm.

# Description
This package is an implementation of Nested Green Simulation Algorithm. 

Suppose in a nested simulation setting, we simulate scenarios in the outer step and then estimate an objective function's expected value conditional on the corresponding outer scenario in the inner step. Nested Green Simulation Algorithm reuses samples from a single distribution to construct all inner scenarios' estimators via importance sampling.

# Installation
Install `GreenSim` from GitHub:
```r
install.packages("devtools")
library(devtools)
devtools::install_github("chenqi57/GreenSim")
library(GreenSim)
```
or directly from the `GreenSim_0.1.0.tar.gz`:
```r
install.packages("GreenSim_0.1.0.tar.gz", repos = NULL, type = "source")
library(GreenSim)
```
# Function
Fourteen functions are contained in this package. Seven functions are the implementation of nested Green Simulation methods, and the others are helper functions for density calculation, random varaible simulation and weight fitting.

## Nested Green Simulation Methods
* `NGS_MIS`: Calculate Marginal Importance Sampling estimators of an objective function at each outer scenario from nested Green Simulation.
* `NGS_MLR`: Calculate Mixture Likelihood Ratio estimators of an objective function at each outer scenario from nested Green Simulation.
* `NGS_CIS`: Calculate Constant Importance Sampling estimators of an objective function at each outer scenario from nested Green Simulation.
* `NGS_OIS`: Calculate Optimal Importance Sampling estimators of a **bounded** objective function at each outer scenario from nested Green Simulation.
* `NGS_LS.c`: Calculate Least Square estimators of an objective function, mimicing Constant Importance Sampling estimators, at each outer scenario from nested Green Simulation.
* `NGS_LS`: Calculate Least Square estimators of an objective function, mimicing Optimal Importance Sampling estimators, at each outer scenario from nested Green Simulation.
* `NGS_LS.fb`: Calculate Least Square estimators of an objective function with given beta at each outer scenario from nested Green Simulation.

## Density Calculation
* `f.beta`: Calculate the density of mixture distribution, consisting each outer scenario distibution, with given weights.
* `g.c`: Calculate the density of optimal proposal distribution without normalizing constant when the objective function is constant.
* `g.star`: Calculate the density of optimal proposal distribution without normalizing constant for an objective function.

## Random Variable Simulation
* `A_R.c`: Simulate random variables from optimal proposal distribution of a **constant** objective function using Acceptance-Rejection.
* `A_R.b`: Simulate random variables from optimal proposal distribution of a **bounded** objective function using Acceptance-Rejection.

## Weight Fitting
* `ls.c.fitting`: Fit optimal weights of each outer scenario distibution by Least Square, aiming to mimic optimal proposal distibution of constant objective function by a mixture distribution.
* `ls.fitting`: Fit optimal weights of each outer scenario distibution by Least Square, aiming to mimic optimal proposal distibution of the specified objective function by a mixture distribution.

# Documentation
Documentations of the above functions can be accesed by typing `?` before each function's name at the R command. 
For instance, the user can read the function `NGS_MLR`'s argument, output and examples in detail by typing `?NGS_MLR`.

# Example
Suppose you have a financial instrument consisting a zero-coupon bond with face value $10 and time-to-maturity one year, as well as a bonus of $10 if the underlying stock price S1 after one year from time 0 is greater than strike price K = 100. Hence, the objective function is h(x) = exp(-11r/12) [10 + 10I(S1>K)] where r is the risk-free rate. 

We set up the function paramters in the following codes. 
```r setup
rf = 5e-2       # risk-free rate
S0 = 100        # initial stock price
vol = 30e-2     # annual volatility

tau = 1/12      # one month
T = 1           # time-to-maturity from time 0

N_Out = 10      # number of outer scenarios
N_In = 5e2      # number of inner scenarios

T2M = T - tau   # time to maturity from outer scenario time

min = qlnorm(1e-4, meanlog = (rf-0.5*vol^2)*tau + log(S0), sdlog = (vol*sqrt(tau)))
max = qlnorm(1-1e-4, meanlog = (rf-0.5*vol^2)*tau + log(S0), sdlog = (vol*sqrt(tau)))
S_tau = seq(from = min, to = max, length.out = N_Out)  # outer scenario stock price

h = function(x){
  return(exp(-T2M*rf) * (10*as.numeric(x>100)+10))
}              # objective function
M = 10 + 10    # upper bound
sn = FALSE     # without self-normalization
```

Then we use the `GreenSim` package to calculate the MLR, OIS, LS estimates as an example.
```r NGS
library(GreenSim)
set.seed(818)
MLR <- NGS_MLR(S_tau, rf, T2M, vol, N_In, h, sn)
df_MLR <- data.frame(outer = S_tau, est = MLR, Method = rep("MLR", N_Out))

OIS <- NGS_OIS(S_tau, rf, T2M, vol, N_In, h, M, sn)
df_OIS <- data.frame(outer = S_tau, est = OIS, Method = rep("OIS", N_Out))

LS <- NGS_LS(S_tau, rf, T2M, vol, N_In, 0.8, h, sn)
df_LS <- data.frame(outer = S_tau, est = LS, Method = rep("LS", N_Out))
```

Finally we illustrate the estimates in the following plot.
```r plot
df <- rbind(df_MLR, df_OIS, df_LS)
ggplot(df) +
  geom_point(aes(x = outer, y = est, color = Method)) +
  labs(x = "Outer Scenario", y = "Estimated Fair Price", 
       title = "Importance Sampling Estimates without Self-Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
```
![image](https://imgur.com/QhiGtER.png#center)
