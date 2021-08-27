\name{ls.c.fitting}
\alias{ls.c.fitting}
\title{
Least Square Fitting of Optimal Beta for Constant Optimal Proposal Distribution
}
\description{
Fit optimal weights of each outer scenario distibution by Least Square, aiming to mimic optimal proposal distibution of constant objective function by a mixture distribution.
}
\usage{
ls.c.fitting(outer, df, rf, N_In = 100)
}
\arguments{
  \item{outer}{vector of parameters simulated in outer scenario.
}
  \item{df}{density functions for the class of distributions inner simulation random variables follow.
}
  \item{rf}{random generation for the class of distributions inner simulation random variables follow.
}
  \item{N_In}{number of inner replications for each outer scenario.
}
}
\value{
List of both a vector of optimal weights of each outer scenario distibution and a vector of MLR random variables used in fitting these weights.
}
\examples{
library(functional)
r = 5e-2    # risk-free rate
S0 = 100    # initial stock price
vol = 30e-2 # annual volatility

tau = 1/12    # one month
T = 1         # time to maturity (from time 0)

N_Out = 10     # number of outer samples
N_In = 5e2      # number of inner samples

T2M = T - tau

min <- qlnorm(1e-4, meanlog = (r-0.5*vol^2)*tau + log(S0), sdlog = (vol*sqrt(tau)))
max <- qlnorm(1-1e-4, meanlog = (r-0.5*vol^2)*tau + log(S0), sdlog = (vol*sqrt(tau)))
S_tau <- seq(from = min, to = max, length.out = N_Out)
mu <- log(S_tau) + (r-0.5*vol^2)*T2M
sig <- vol * sqrt(T2M)
df <- Curry(dnorm, sd = sig)
rf <- Curry(rnorm, sd = sig)

ls.c.fitting(mu, df, rf, N_In)
}
