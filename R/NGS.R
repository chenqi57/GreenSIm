g.c <- function(x, S_tau, rf, T2M, vol){
  g.x <- 0
  mu <- log(S_tau) + (rf-0.5*vol^2)*T2M
  sig <- vol * sqrt(T2M)
  for(i in 1:length(S_tau)){
    g.x <- g.x + dnorm(x, mu[i], sig)^2
  }
  return(sqrt(g.x))
}

g.star <- function(x, S_tau, rf, T2M, vol, h){
  ST <- exp(x)
  return(abs(h(ST)) * g.c(x, S_tau, rf, T2M, vol))
}

f.beta <- function(x, beta, S_tau, rf, T2M, vol){
  f.x <- 0
  mu <- log(S_tau) + (rf-0.5*vol^2)*T2M
  sig <- vol * sqrt(T2M)
  for(i in 1:length(S_tau)){
    f.x <- f.x + beta[i] * dnorm(x, mu[i], sig)
  }
  return(f.x)
}

A_R.b <- function(S_tau, rf = 5e-2, T2M = 11/12, vol = 30e-2, N_In = 1e2, h, M){
  N_Out <- length(S_tau)

  sig <- vol * sqrt(T2M)
  mu <- log(S_tau) + (rf-0.5*vol^2)*T2M

  # Sampling Mixture Normal r.v.
  U <- runif(N_Out*N_In*N_Out)
  U <- ceiling(U*N_Out)
  mixed <- rnorm(N_Out*N_In*N_Out, mu[U], sig)

  # Computing L1, L2 sum
  beta <- rep(1/N_Out, N_Out)
  f.x <- f.beta(mixed, beta, S_tau, rf, T2M, vol) * N_Out
  g.x <- g.c(mixed, S_tau, rf, T2M, vol)
  ST <- exp(mixed)
  h.x <- h(ST)

  # Acceptance and Rejection
  Z <- runif(N_Out*N_In*N_Out)
  a <- (Z < ((abs(h.x)*g.x) / (M*f.x)))
  rv <- mixed[a]
  sampled <- rv[1:(N_Out*N_In)]

  result <- list("c" = (N_Out*N_In) / (M * sum(a)), "sampled" = sampled)
  return(result)
}

A_R.c <- function(S_tau, rf = 5e-2, T2M = 11/12, vol = 30e-2, N_In = 1e2){
  N_Out <- length(S_tau)

  sig <- vol * sqrt(T2M)
  mu <- log(S_tau) + (rf-0.5*vol^2)*T2M

  # Sampling Mixture Normal r.v.
  U <- runif(N_Out*N_In*N_Out)
  U <- ceiling(U*N_Out)
  mixed <- rnorm(N_Out*N_In*N_Out, mu[U], sig)

  # Computing L1, L2 sum
  # Computing L1, L2 sum
  beta <- rep(1/N_Out, N_Out)
  f.x <- f.beta(mixed, beta, S_tau, rf, T2M, vol) * N_Out
  g.x <- g.c(mixed, S_tau, rf, T2M, vol)

  # Acceptance and Rejection
  Z <- runif(N_Out*N_In*N_Out)
  a <- (Z < (g.x/f.x))
  rv <- mixed[a]
  sampled <- rv[1:(N_Out*N_In)]

  result <- list("c" = (N_Out*N_In) / sum(a), "sampled" = sampled)
  return(result)
}

ls.fitting <- function(S_tau, rf = 5e-2, T2M = 11/12, vol = 30e-2, N_In = 1e2, h){
  N_Out <- length(S_tau)

  sig <- vol * sqrt(T2M)
  mu <- log(S_tau) + (rf-0.5*vol^2)*T2M

  # Sampling Mixture Normal r.v.
  U <- rep((1:N_Out), N_In)
  mixed <- rnorm(N_Out*N_In, mu[U], sig)

  # Computing L1, L2 sum
  f.x <- 0
  g.x <- 0
  X <- NULL
  for(i in 1:length(S_tau)){
    l.i <- dnorm(mixed, mu[i], sig)
    f.x <- f.x + l.i
    g.x <- g.x + l.i^2
    X <- cbind(X, l.i)
  }
  g.x <- sqrt(g.x)
  ST <- exp(mixed)
  h.x <- h(ST)

  #Fitting beta
  Y = (abs(h.x)*g.x)
  soln <- nnls(X, Y)
  beta <- soln$x
  beta <- beta / sum(beta)

  result <- list("beta" = beta, "sampled" = mixed)
  return(result)
}

ls.c.fitting <- function(S_tau, rf = 5e-2, T2M = 11/12, vol = 30e-2, N_In = 1e2){
  N_Out <- length(S_tau)

  sig <- vol * sqrt(T2M)
  mu <- log(S_tau) + (rf-0.5*vol^2)*T2M

  # Sampling Mixture Normal r.v.
  U <- rep((1:N_Out), N_In)
  mixed <- rnorm(N_Out*N_In, mu[U], sig)

  # Computing L1, L2 sum
  f.x <- 0
  g.x <- 0
  X <- NULL
  for(i in 1:length(S_tau)){
    l.i <- dnorm(mixed, mu[i], sig)
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

NGS_MIS <- function(S_tau, rf = 5e-2, T2M = 11/12, vol = 30e-2, N_In = 1e2, S0 = 100, T = 1, h, sn = FALSE){
  N_Out <- length(S_tau)

  mu <- log(S_tau) + (rf-0.5*vol^2)*T2M
  sig <- vol * sqrt(T2M)

  # Sampling Marginal Normal r.v.
  sampled <- rnorm(N_Out*N_In, log(S0) + (rf-0.5*vol^2)*T, vol * sqrt(T))
  ST = exp(sampled)
  C_tau_sample = h(ST)

  # Calculate likelihoods
  llh.g <- dnorm(sampled, log(S0) + (rf-0.5*vol^2)*T, vol * sqrt(T))

  # Calculate option payoffs
  C_tau_MIS <- NULL
  if(sn){
    for (i in 1:N_Out){
      llh.f <- dnorm(sampled, mu[i], sig)
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_MIS <- c(C_tau_MIS, mean(C_tau) / mean(lr))
    }
  }
  else{
    for (i in 1:N_Out){
      llh.f <- dnorm(sampled, mu[i], sig)
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_MIS <- c(C_tau_MIS, mean(C_tau))
    }
  }

  return(C_tau_MIS)
}

NGS_MLR <- function(S_tau, rf = 5e-2, T2M = 11/12, vol = 30e-2, N_In = 1e2, h, sn = FALSE){
  N_Out <- length(S_tau)

  mu <- log(S_tau) + (rf-0.5*vol^2)*T2M
  sig <- vol * sqrt(T2M)

  # Sampling Mixture Normal r.v.
  U <- rep((1:N_Out), N_In)
  sampled <- rnorm(N_Out*N_In, mu[U], sig)
  ST = exp(sampled)
  C_tau_sample = h(ST)

  # Calculate likelihoods
  beta <- rep(1/N_Out, N_Out)
  llh.g <- f.beta(sampled, beta, S_tau, rf, T2M, vol)

  # Calculate option payoffs
  C_tau_MLR <- NULL
  if(sn){
    for (i in 1:N_Out){
      llh.f <- dnorm(sampled, mu[i], sig)
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_MLR <- c(C_tau_MLR, mean(C_tau) / mean(lr))
    }
  }
  else{
    for (i in 1:N_Out){
      llh.f <- dnorm(sampled, mu[i], sig)
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_MLR <- c(C_tau_MLR, mean(C_tau))
    }
  }

  return(C_tau_MLR)
}

NGS_OIS <- function(S_tau, rf = 5e-2, T2M = 11/12, vol = 30e-2, N_In = 1e2, h, M, sn = FALSE){
  N_Out <- length(S_tau)

  # Sampling OIS r.v.
  A.R <- A_R.b(S_tau, rf, T2M, vol, N_In, h, M)
  sampled <- A.R$sampled
  ST = exp(sampled)
  C_tau_sample = h(ST)

  # Calculate likelihoods
  mu <- log(S_tau) + (rf-0.5*vol^2)*T2M
  sig <- vol * sqrt(T2M)
  c <- A.R$c
  llh.g <- c * g.star(sampled, S_tau, rf, T2M, vol, h)

  # Calculate option payoffs
  C_tau_OIS <- NULL
  if(sn){
    for (i in 1:N_Out){
      llh.f <- dnorm(sampled, mu[i], sig)
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_OIS <- c(C_tau_OIS, mean(C_tau) / mean(lr))
    }
  }
  else{
    for (i in 1:N_Out){
      llh.f <- dnorm(sampled, mu[i], sig)
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_OIS <- c(C_tau_OIS, mean(C_tau))
    }
  }

  return(C_tau_OIS)
}

NGS_CIS <- function(S_tau, rf = 5e-2, T2M = 11/12, vol = 30e-2, N_In = 1e2, h, sn = FALSE){
  N_Out <- length(S_tau)

  # Sampling CIS r.v.
  A.R <- A_R.c(S_tau, rf, T2M, vol, N_In)
  sampled <- A.R$sampled
  ST = exp(sampled)
  C_tau_sample = h(ST)

  # Calculate likelihoods
  mu <- log(S_tau) + (rf-0.5*vol^2)*T2M
  sig <- vol * sqrt(T2M)
  c <- A.R$c
  llh.g <- c * g.c(sampled, S_tau, rf, T2M, vol)

  # Calculate option payoffs
  C_tau_CIS <- NULL
  if(sn){
    for (i in 1:N_Out){
      llh.f <- dnorm(sampled, mu[i], sig)
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_CIS <- c(C_tau_CIS, mean(C_tau) / mean(lr))
    }
  }
  else{
    for (i in 1:N_Out){
      llh.f <- dnorm(sampled, mu[i], sig)
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_CIS <- c(C_tau_CIS, mean(C_tau))
    }
  }

  return(C_tau_CIS)
}

NGS_LS <- function(S_tau, rf = 5e-2, T2M = 11/12, vol = 30e-2, N_In = 1e2, p, h, sn = FALSE){
  N_Out <- length(S_tau)

  mu <- log(S_tau) + (rf-0.5*vol^2)*T2M
  sig <- vol * sqrt(T2M)

  # Estimate beta and sample
  result <- ls.fitting(S_tau, rf, T2M, vol, N_In, h)
  beta <- result$beta

  sampled <- NULL
  n <- as.integer(beta * N_Out * as.integer(p*N_In))
  for (i in 1:N_Out){
    sampled <- c(sampled, rnorm(n[i], mu[i], sig))
  }
  ST = exp(sampled)
  C_tau_sample = h(ST)

  # Calculate likelihoods
  llh.g <- f.beta(sampled, beta, S_tau, rf, T2M, vol)

  # Calculate option payoffs
  C_tau_QP <- NULL
  if(sn){
    for (i in 1:N_Out){
      llh.f <- dnorm(sampled, mu[i], sig)
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_QP <- c(C_tau_QP, mean(C_tau) / mean(lr))
    }
  }
  else{
    for (i in 1:N_Out){
      llh.f <- dnorm(sampled, mu[i], sig)
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_QP <- c(C_tau_QP, mean(C_tau))
    }
  }

  return(C_tau_QP)
}

NGS_LS.c <- function(S_tau, rf = 5e-2, T2M = 11/12, vol = 30e-2, N_In = 1e2, p, h, sn = FALSE){
  N_Out <- length(S_tau)

  mu <- log(S_tau) + (rf-0.5*vol^2)*T2M
  sig <- vol * sqrt(T2M)

  # Estimate beta and sample
  result <- ls.c.fitting(S_tau, rf, T2M, vol, N_In)
  beta <- result$beta

  sampled <- NULL
  n <- as.integer(beta * N_Out * as.integer(p*N_In))
  for (i in 1:N_Out){
    sampled <- c(sampled, rnorm(n[i], mu[i], sig))
  }
  ST = exp(sampled)
  C_tau_sample = h(ST)

  # Calculate likelihoods
  llh.g <- f.beta(sampled, beta, S_tau, rf, T2M, vol)

  # Calculate option payoffs
  C_tau_QP <- NULL
  if(sn){
    for (i in 1:N_Out){
      llh.f <- dnorm(sampled, mu[i], sig)
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_QP <- c(C_tau_QP, mean(C_tau) / mean(lr))
    }
  }
  else{
    for (i in 1:N_Out){
      llh.f <- dnorm(sampled, mu[i], sig)
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_QP <- c(C_tau_QP, mean(C_tau))
    }
  }

  return(C_tau_QP)
}

NGS_LS.fb <- function(S_tau, rf = 5e-2, T2M = 11/12, vol = 30e-2, N_In = 1e2, beta, p, h, sn = FALSE){
  N_Out <- length(S_tau)

  mu <- log(S_tau) + (rf-0.5*vol^2)*T2M
  sig <- vol * sqrt(T2M)

  # Sample from fixed beta
  sampled <- NULL
  n <- as.integer(beta * N_Out * as.integer(p*N_In))
  for (i in 1:N_Out){
    sampled <- c(sampled, rnorm(n[i], mu[i], sig))
  }
  ST = exp(sampled)
  C_tau_sample = h(ST)

  # Calculate likelihoods
  llh.g <- f.beta(sampled, beta, S_tau, rf, T2M, vol)

  # Calculate option payoffs
  C_tau_QP <- NULL
  if(sn){
    for (i in 1:N_Out){
      llh.f <- dnorm(sampled, mu[i], sig)
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_QP <- c(C_tau_QP, mean(C_tau) / mean(lr))
    }
  }
  else{
    for (i in 1:N_Out){
      llh.f <- dnorm(sampled, mu[i], sig)
      lr <- llh.f / llh.g
      C_tau <- C_tau_sample * lr
      C_tau_QP <- c(C_tau_QP, mean(C_tau))
    }
  }

  return(C_tau_QP)
}
