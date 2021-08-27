g.c <- function(x, outer, df){
  g.x <- 0
  for(i in 1:length(outer)){
    g.x <- g.x + df(x, outer[i])^2
  }
  return(sqrt(g.x))
}

g.star <- function(x, outer, df, h){
  return(abs(h(x)) * g.c(x, outer, df))
}

f.beta <- function(x, beta, outer, df){
  f.x <- 0
  for(i in 1:length(outer)){
    f.x <- f.x + beta[i] * df(x, outer[i])
  }
  return(f.x)
}

A_R.b <- function(outer, df, rf, N_In = 1e2, h, M){
  N_Out <- length(outer)

  # Sampling Mixture r.v.
  U <- runif(N_Out*N_In*N_Out)
  U <- ceiling(U*N_Out)
  mixed <- rf(N_Out*N_In*N_Out, outer[U])

  # Computing L1, L2 sum
  beta <- rep(1/N_Out, N_Out)
  f.x <- f.beta(mixed, beta, outer, df) * N_Out
  g.x <- g.c(mixed, outer, df)

  # Acceptance and Rejection
  Z <- runif(N_Out*N_In*N_Out)
  ll <- h(mixed) * g.x
  a <- (Z < (abs(ll)/(M*f.x)))
  rv <- mixed[a]
  sampled <- rv[1:(N_Out*N_In)]
  c <- (N_Out*N_In) / (M * sum(a))

  result <- list("c" = c, "sampled" = sampled, "ll" = c * ll[a][1:(N_Out*N_In)])
  return(result)
}

A_R.c <- function(outer, df, rf, N_In = 1e2){
  N_Out <- length(outer)

  # Sampling Mixture r.v.
  U <- runif(N_Out*N_In*N_Out)
  U <- ceiling(U*N_Out)
  mixed <- rf(N_Out*N_In*N_Out, outer[U])

  # Computing L1, L2 sum
  # Computing L1, L2 sum
  beta <- rep(1/N_Out, N_Out)
  f.x <- f.beta(mixed, beta, outer, df) * N_Out
  g.x <- g.c(mixed, outer, df)

  # Acceptance and Rejection
  Z <- runif(N_Out*N_In*N_Out)
  a <- (Z < (g.x/f.x))
  rv <- mixed[a]
  sampled <- rv[1:(N_Out*N_In)]
  c <- (N_Out*N_In) / sum(a)

  result <- list("c" = c, "sampled" = sampled, "ll" = c * g.x[a][1:(N_Out*N_In)])
  return(result)
}

ls.fitting <- function(outer, df, rf, N_In = 1e2, h){
  N_Out <- length(outer)

  # Sampling Mixture r.v.
  U <- rep((1:N_Out), N_In)
  mixed <- rf(N_Out*N_In, outer[U])

  # Computing L1, L2 sum
  f.x <- 0
  g.x <- 0
  X <- NULL
  for(i in 1:length(outer)){
    l.i <- df(mixed, outer[i])
    f.x <- f.x + l.i
    g.x <- g.x + l.i^2
    X <- cbind(X, l.i)
  }
  g.x <- sqrt(g.x)

  #Fitting beta
  Y = (abs(h(mixed))*g.x)
  soln <- nnls(X, Y)
  beta <- soln$x
  beta <- beta / sum(beta)

  result <- list("beta" = beta, "sampled" = mixed)
  return(result)
}

ls.c.fitting <- function(outer, df, rf, N_In = 1e2){
  N_Out <- length(outer)

  # Sampling Mixture r.v.
  U <- rep((1:N_Out), N_In)
  mixed <- rf(N_Out*N_In, outer[U])

  # Computing L1, L2 sum
  f.x <- 0
  g.x <- 0
  X <- NULL
  for(i in 1:length(outer)){
    l.i <- df(mixed, outer[i])
    f.x <- f.x + l.i
    g.x <- g.x + l.i^2
    X <- cbind(X, l.i)
  }
  g.x <- sqrt(g.x)

  #Fitting beta
  Y = g.x
  soln <- nnls(X, Y)
  beta <- soln$x
  beta <- beta / sum(beta)

  result <- list("beta" = beta, "sampled" = mixed)
  return(result)
}

NGS_MLR <- function(outer, df, rf, N_In = 1e2, h, sn = FALSE){
  N_Out <- length(outer)

  # Sampling Mixture r.v.
  U <- rep((1:N_Out), N_In)
  sampled <- rf(N_Out*N_In, outer[U])
  C_tau_sample = h(sampled)

  # Calculate likelihoods
  beta <- rep(1/N_Out, N_Out)
  llh.g <- f.beta(sampled, beta, outer, df)

  # Calculate option payoffs
  C_tau_MLR <- NULL
  if(sn){
    for (i in 1:N_Out){
      llh.f <- df(sampled, outer[i])
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_MLR <- c(C_tau_MLR, mean(C_tau) / mean(lr))
    }
  }
  else{
    for (i in 1:N_Out){
      llh.f <- df(sampled, outer[i])
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_MLR <- c(C_tau_MLR, mean(C_tau))
    }
  }

  return(C_tau_MLR)
}

NGS_OIS <- function(outer, df, rf, N_In = 1e2, h, M, sn = FALSE){
  N_Out <- length(outer)

  # Sampling OIS r.v.
  A.R <- A_R.b(outer, df, rf, N_In, h, M)
  sampled <- A.R$sampled
  C_tau_sample = h(sampled)

  # Calculate likelihoods
  llh.g <- A.R$ll

  # Calculate option payoffs
  C_tau_OIS <- NULL
  if(sn){
    for (i in 1:N_Out){
      llh.f <- df(sampled, outer[i])
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_OIS <- c(C_tau_OIS, mean(C_tau) / mean(lr))
    }
  }
  else{
    for (i in 1:N_Out){
      llh.f <- df(sampled, outer[i])
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_OIS <- c(C_tau_OIS, mean(C_tau))
    }
  }

  return(C_tau_OIS)
}

NGS_CIS <- function(outer, df, rf, N_In = 1e2, h, sn = FALSE){
  N_Out <- length(outer)

  # Sampling CIS r.v.
  A.R <- A_R.c(outer, df, rf, N_In)
  sampled <- A.R$sampled
  C_tau_sample = h(sampled)

  # Calculate likelihoods
  llh.g <- A.R$ll

  # Calculate option payoffs
  C_tau_CIS <- NULL
  if(sn){
    for (i in 1:N_Out){
      llh.f <- df(sampled, outer[i])
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_CIS <- c(C_tau_CIS, mean(C_tau) / mean(lr))
    }
  }
  else{
    for (i in 1:N_Out){
      llh.f <- df(sampled, outer[i])
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_CIS <- c(C_tau_CIS, mean(C_tau))
    }
  }

  return(C_tau_CIS)
}

NGS_LS <- function(outer, df, rf, N_In = 1e2, p, h, sn = FALSE){
  N_Out <- length(outer)

  # Estimate beta and sample
  result <- ls.fitting(outer, df, rf, N_In, h)
  beta <- result$beta

  sampled <- NULL
  n <- as.integer(beta * N_Out * as.integer(p*N_In))
  for (i in 1:N_Out){
    sampled <- c(sampled, rf(n[i], outer[i]))
  }
  C_tau_sample = h(sampled)

  # Calculate likelihoods
  llh.g <- f.beta(sampled, beta, outer, df)

  # Calculate option payoffs
  C_tau_QP <- NULL
  if(sn){
    for (i in 1:N_Out){
      llh.f <- df(sampled, outer[i])
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_QP <- c(C_tau_QP, mean(C_tau) / mean(lr))
    }
  }
  else{
    for (i in 1:N_Out){
      llh.f <- df(sampled, outer[i])
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_QP <- c(C_tau_QP, mean(C_tau))
    }
  }

  return(C_tau_QP)
}

NGS_LS.c <- function(outer, df, rf, N_In = 1e2, p, h, sn = FALSE){
  N_Out <- length(outer)

  # Estimate beta and sample
  result <- ls.c.fitting(outer, df, rf, N_In)
  beta <- result$beta

  sampled <- NULL
  n <- as.integer(beta * N_Out * as.integer(p*N_In))
  for (i in 1:N_Out){
    sampled <- c(sampled, rf(n[i], outer[i]))
  }
  C_tau_sample = h(sampled)

  # Calculate likelihoods
  llh.g <- f.beta(sampled, beta, outer, df)

  # Calculate option payoffs
  C_tau_QP <- NULL
  if(sn){
    for (i in 1:N_Out){
      llh.f <- df(sampled, outer[i])
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_QP <- c(C_tau_QP, mean(C_tau) / mean(lr))
    }
  }
  else{
    for (i in 1:N_Out){
      llh.f <- df(sampled, outer[i])
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_QP <- c(C_tau_QP, mean(C_tau))
    }
  }

  return(C_tau_QP)
}

NGS_LS.fb <- function(outer, df, rf, N_In = 1e2, beta, p, h, sn = FALSE){
  N_Out <- length(outer)

  # Sample from fixed beta
  sampled <- NULL
  n <- as.integer(beta * N_Out * as.integer(p*N_In))
  for (i in 1:N_Out){
    sampled <- c(sampled, rf(n[i], outer[i]))
  }
  C_tau_sample = h(sampled)

  # Calculate likelihoods
  llh.g <- f.beta(sampled, beta, outer, df)

  # Calculate option payoffs
  C_tau_QP <- NULL
  if(sn){
    for (i in 1:N_Out){
      llh.f <- df(sampled, outer[i])
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_QP <- c(C_tau_QP, mean(C_tau) / mean(lr))
    }
  }
  else{
    for (i in 1:N_Out){
      llh.f <- df(sampled, outer[i])
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_QP <- c(C_tau_QP, mean(C_tau))
    }
  }

  return(C_tau_QP)
}

NGS_MIS <- function(outer, df, rf, dg, rg, N_In = 1e2, h, sn = FALSE){
  N_Out <- length(outer)

  # Sampling Marginal r.v.
  sampled <- rg(N_Out*N_In)
  C_tau_sample = h(sampled)

  # Calculate likelihoods
  llh.g <- dg(sampled)

  # Calculate option payoffs
  C_tau_MIS <- NULL
  if(sn){
    for (i in 1:N_Out){
      llh.f <- df(sampled, outer[i])
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_MIS <- c(C_tau_MIS, mean(C_tau) / mean(lr))
    }
  }
  else{
    for (i in 1:N_Out){
      llh.f <- df(sampled, outer[i])
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_MIS <- c(C_tau_MIS, mean(C_tau))
    }
  }

  return(C_tau_MIS)
}
