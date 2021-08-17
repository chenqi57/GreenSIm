\name{A_R.b}
\alias{A_R.b}
\title{
Acceptance-Rejection Algorithm for Bounded Objective
}
\description{
Simulate random variables from optimal proposal distribution of a bounded objective function using Acceptance-Rejection.
}
\usage{
A_R.b(S_tau, rf = 0.05, T2M = 11/12, vol = 0.3, N_In = 100, h, M)
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
  \item{h}{objective function.
}
  \item{M}{upper bound of objective function's absolute value.
}
}
\value{
List of both an estimated normalizing constant and a vector of random variables following optimal proposal distribution for the specified objective function.
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

h <- function(x){return(pmax(x-90, 0) - pmax(x-110, 0) + 10)}
M <- 30

A_R.b(S_tau, rf, T2M, vol, N_In, h, M)
}

