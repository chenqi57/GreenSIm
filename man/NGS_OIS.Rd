\name{NGS_OIS}
\alias{NGS_OIS}
\title{
Optimal Importance Sampling Estimators
}
\description{
Calculate Optimal Importance Sampling estimators of the specified bounded objective function at each outer scenario from nested Green Simulation.
}
\usage{
NGS_OIS(outer, df, rf, N_In = 100, h, M, sn = FALSE)
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
  \item{h}{objective function.
}
  \item{M}{upper bound of objective function's absolute value.
}
  \item{sn}{logical; if TRUE, self-normalized the likelihood ratios.
}
}
\value{
A vector of Optimal Importance Sampling estimates for each outer scenario. The length of the result is determined by number of outer scenarios.
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

h <- function(x){return(pmax(exp(x)-90, 0) - pmax(exp(x)-110, 0) + 10)}
M <- 30

NGS_OIS(mu, df, rf, N_In, h, M, TRUE)
NGS_OIS(mu, df, rf, N_In, h, M, FALSE)
}
