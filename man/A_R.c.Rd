\name{A_R.c}
\alias{A_R.c}
\title{
Acceptance-Rejection Algorithm for Constant Objective
}
\description{
Simulate random variables from optimal proposal distribution of a constant objective function using Acceptance-Rejection.
}
\usage{
A_R.c(S_tau, rf = 0.05, T2M = 11/12, vol = 0.3, N_In = 100)
}
\arguments{
  \item{S_tau}{vector of outer scenario stock prices.
}
  \item{rf}{risk-free rate.
}
  \item{T2M}{time-to-maturity from outer scenario time.
}
  \item{vol}{annual volatility of the underlying stock.
}
  \item{N_In}{number of inner replications for each outer scenario.
}
}
\value{
List of both an estimated normalizing constant and a vector of random variables following optimal proposal distribution for constant objective function.
}
\examples{
rf = 5e-2   # risk-free rate
S0 = 100    # initial stock price
vol = 30e-2 # annual volatility

tau = 1/12    # one month
T = 1         # time to maturity (from time 0)

N_Out = 10     # number of outer samples
N_In = 5e2      # number of inner samples

T2M = T - tau

min <- qlnorm(1e-4, meanlog = (rf-0.5*vol^2)*tau + log(S0), sdlog = (vol*sqrt(tau)))
max <- qlnorm(1-1e-4, meanlog = (rf-0.5*vol^2)*tau + log(S0), sdlog = (vol*sqrt(tau)))
S_tau <- seq(from = min, to = max, length.out = N_Out)

A_R.c(S_tau, rf, T2M, vol, N_In)
}
