\name{NGS_MLR}
\alias{NGS_MLR}
\title{
Mixture Likelihood Ratio Estimators
}
\description{
Calculate Mixture Likelihood Ratio estimators of the specified objective function at each outer scenario from nested Green Simulation.
}
\usage{
NGS_MLR(S_tau, rf = 0.05, T2M = 11/12, vol = 0.3, N_In = 100, h, sn = FALSE)
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
  \item{sn}{logical; if TRUE, self-normalized the likelihood ratios.
}
}
\value{
A vector of Mixture Likelihood Ratio estimates for each outer scenario. The length of the result is determined by number of outer scenarios.
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

NGS_MLR(S_tau, rf, T2M, vol, N_In, h, TRUE)
NGS_MLR(S_tau, rf, T2M, vol, N_In, h, FALSE)
}
