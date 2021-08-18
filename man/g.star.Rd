\name{g.star}
\alias{g.star}
\title{
Optimal Proposal Distribution
}
\description{
Calculate the density of optimal proposal distribution without normalizing constant for the specified objective function.
}
\usage{
g.star(x, S_tau, rf, T2M, vol, h)
}
\arguments{
  \item{x}{vector of random variables.
}
  \item{S_tau}{vector of outer scenario stock prices.
}
  \item{rf}{risk-free rate.
}
  \item{T2M}{time-to-maturity from outer scenario time.
}
  \item{vol}{annual volatility of the underlying stock.
}
  \item{h}{objective function.
}
}
\value{
Density without normalizing constant of the input random variable vector.
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

x <- seq(from = 4, to = 6, length.out = 10)
h <- function(x){return(pmax(x-90, 0))}

g.star(x, S_tau, rf, T2M, vol, h)
}
