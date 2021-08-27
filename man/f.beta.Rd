\name{f.beta}
\alias{f.beta}
\title{
Mixture Distribution of Outer Scenarios
}
\description{
Calculate the density of mixture distribution, consisting each outer scenario distibution, with given weights.
}
\usage{
f.beta(x, beta, outer, df)
}
\arguments{
  \item{x}{vector of random variables.
}
  \item{beta}{vecotr of weights for each outer scenario distribution.
}
  \item{outer}{vector of parameters simulated in outer scenario.
}
  \item{df}{density functions for the class of distributions inner simulation random variables follow.
}
}
\details{
The MLR density is equivalent to having beta as equal wights.
}
\value{
Density of the input random variable vector.
}
\examples{
library(functional)
r = 5e-2    # risk-free rate
S0 = 100    # initial stock price
vol = 30e-2 # annual volatility

tau = 1/12    # one month
T = 1         # time to maturity (from time 0)

N_Out = 10     # number of outer samples

T2M = T - tau

min <- qlnorm(1e-4, meanlog = (r-0.5*vol^2)*tau + log(S0), sdlog = (vol*sqrt(tau)))
max <- qlnorm(1-1e-4, meanlog = (r-0.5*vol^2)*tau + log(S0), sdlog = (vol*sqrt(tau)))
S_tau <- seq(from = min, to = max, length.out = N_Out)
mu <- log(S_tau) + (r-0.5*vol^2)*T2M
sig <- vol * sqrt(T2M)
df <- Curry(dnorm, sd = sig)

x <- seq(from = 4, to = 6, length.out = 10)
beta <- rep(1/N_Out, N_Out)

f.beta(x, beta, mu, df)
}

