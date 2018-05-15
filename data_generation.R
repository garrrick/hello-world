data_generate <- function(size,
                          model.type = 'time',
                          dist.name = 'Weibull',
                          beta.par1,
                          beta.par2,
                          beta.cure,
                          general.dist = 'Weibull',
                          general.par1,
                          general.par2,
                          cen.dist = 'Weibull',
                          cen.par1,
                          cen.par2,
                          c4_from_c234 = 0,
                          covariates.condition = 'independent',
                          mu.par1 = c(1, rep(0, length(beta.par1)-1)),
                          mu.par2 = c(1, rep(0, length(beta.par2)-1)),
                          mu.cure = c(1, rep(0, length(beta.cure)-1)),
                          sigma.par1 = diag(c(0, rep(.25, length(beta.par1) - 1))),
                          sigma.par2 = diag(c(0, rep(1, length(beta.par2) - 1))),
                          sigma.cure = diag(c(0, rep(1, length(beta.cure) - 1)))) {
  # c4_from_c234 = .2; model.type = 'time'; dist.name = 'Weibull'; size = ssize; beta.cure = beta.cure; beta.par1 = beta.par1; beta.par2 = beta.par2; general.par1 = general.par1; general.par2 = general.par2; cen.par1 = cen.par1; cen.par2 = cen.par2; covariates.condition = 'equal'; sigma.par1 = 0; sigma.par2 = 0; sigma.cure = 0; mu.par1 = c(1, rep(0, length(beta.par1)-1)); mu.par2 = c(1, rep(0, length(beta.par2)-1)); mu.cure = c(1, rep(0, length(beta.cure)-1))
  library(MASS)
  if (dist.name == 'Weibull') {
    # model excess hazard with Weibull
    # generate covariates of par2 parameter
    p2x <- mvrnorm(size, mu = mu.par2, Sigma = sigma.par2)
    # par2 parameter
    par2.parameter <- c(exp(crossprod(beta.par2,t(p2x))))
    # par1 parameter p = exp(betaX1), X1 from standard normal
    # generate covariates of par1 parameter
    if ((length(beta.par1) > 1) & (length(beta.par1) == length(beta.par2))) {
      if (covariates.condition == 'equal') {
        p1x <- p2x
      } else if (covariates.condition == 'independent') {
        p1x <- mvrnorm(size, mu = mu.par1, Sigma = sigma.par1)
      } else {
        p1x <- mvrnorm(size, mu = mu.par1, Sigma = sigma.par1)
      }
    } else {
      p1x <- mvrnorm(size, mu = mu.par1, Sigma = sigma.par1)
    }
    # par1 parameter
    par1.parameter <- c(exp(crossprod(beta.par1,t(p1x))))
    # model cure time I(T > c), c=exp(betaX3), X3 from standard normal
    # mean and variance
    # generate covariates of cure model
    if ((length(beta.cure) > 1) & (length(beta.cure) == length(beta.par2))) {
      if (covariates.condition == 'equal') {
        cxx <- p2x
      } else if (covariates.condition == 'independent') {
        cxx <- mvrnorm(size, mu = mu.cure, Sigma = sigma.cure)
      } else {
        cxx <- mvrnorm(size, mu = mu.cure, Sigma = sigma.cure)
      }
    } else {
      cxx <- mvrnorm(size, mu = mu.cure, Sigma = sigma.cure)
    }
    # generate event time from excess hazard
    if (model.type == 'mixture') {
      # cure rate
      cr <- (1+exp(-crossprod(beta.cure,t(cxx))))^(-1)
      mix1 <- 1-rbinom(size, 1, cr)
      e1 <- rweibull(size, shape = par1.parameter, scale = par2.parameter)
    } else if (model.type == 'non-mixture') {
      # cure rate
      cr <- (1+exp(-crossprod(beta.cure,t(cxx))))^(-1)
      N <- rpois(size, -log(cr))
      mix1 <- I(N > 0)
      e1N <- list()
      for(i in 1:size) {
        e1N[[i]] <- rweibull(N[i], shape = par1.parameter, scale = par2.parameter)
      }
      e1 <- sapply(e1N, min)
    } else if (model.type == 'time') {
      # cure time
      cr <- c(exp(crossprod(beta.cure,t(cxx))))
      e1 <- rweibull(size, shape = par1.parameter, scale = par2.parameter)
      mix1 <- e1 <= cr
    }
  } else if (dist.name == 'log-normal') {
    # model excess hazard with log-normal
    # generate covariates of par1 parameter
    p1x <- mvrnorm(size, mu = mu.par1, Sigma = sigma.par1)
    # par1 parameter
    par1.parameter <- c(crossprod(beta.par1,t(p1x)))
    # no covariate in par2 parameter
    # generate covariates of par2 parameter
    p2x <- mvrnorm(size, mu = mu.par2, Sigma = sigma.par2)
    # par2 parameter
    par2.parameter <- c(crossprod(beta.par2,t(p2x)))
    # model cure time I(T > c), c=exp(betaX3), X3 from standard normal
    # mean and variance
    # generate covariates of cure model
    if ((length(beta.cure) > 1) & (length(beta.cure) == length(beta.par1))) {
      if (covariates.condition == 'equal') {
        cxx <- p1x
      } else if (covariates.condition == 'independent') {
        cxx <- mvrnorm(size, mu = mu.cure, Sigma = sigma.cure)
      } else {
        cxx <- mvrnorm(size, mu = mu.cure, Sigma = sigma.cure)
      }
    } else {
      cxx <- mvrnorm(size, mu = mu.cure, Sigma = sigma.cure)
    }
    # generate event time from excess hazard
    if (model.type == 'mixture') {
      # cure rate
      cr <- (1+exp(-crossprod(beta.cure,t(cxx))))^(-1)
      mix1 <- 1-rbinom(size, 1, cr)
      e1 <- rlnorm(size, meanlog = par1.parameter, sdlog = par2.parameter)
    } else if (model.type == 'non-mixture') {
      # cure rate
      cr <- (1+exp(-crossprod(beta.cure,t(cxx))))^(-1)
      N <- rpois(size, -log(cr))
      mix1 <- I(N > 0)
      e1N <- list()
      for(i in 1:size) {
        e1N[[i]] <- rlnorm(N[i], meanlog = par1.parameter, sdlog = par2.parameter)
      }
      e1 <- sapply(e1N, min)
    } else if (model.type == 'time') {
      # cure time
      cr <- c(exp(crossprod(beta.cure,t(cxx))))
      e1 <- rlnorm(size, meanlog = par1.parameter, sdlog = par2.parameter)
      mix1 <- e1 <= cr
    }
  } else if (dist.name == 'log-logistic') {
    # model excess hazard with log-logistic
    # generate covariates of par2 parameter
    p2x <- mvrnorm(size, mu = mu.par2, Sigma = sigma.par2)
    # par2 parameter
    par2.parameter <- c(exp(crossprod(beta.par2,t(p2x))))
    # par1 parameter p = exp(betaX1), X1 from standard normal
    # generate covariates of par1 parameter
    if ((length(beta.par1) > 1) & (length(beta.par1) == length(beta.par2))) {
      if (covariates.condition == 'equal') {
        p1x <- p2x
      } else if (covariates.condition == 'independent') {
        p1x <- mvrnorm(size, mu = mu.par1, Sigma = sigma.par1)
      } else {
        p1x <- mvrnorm(size, mu = mu.par1, Sigma = sigma.par1)
      }
    } else {
      p1x <- mvrnorm(size, mu = mu.par1, Sigma = sigma.par1)
    }
    # par1 parameter
    par1.parameter <- c(exp(crossprod(beta.par1,t(p1x))))
    # model cure time I(T > c), c=exp(betaX3), X3 from standard normal
    # mean and variance
    # generate covariates of cure time model
    if ((length(beta.cure) > 1) & (length(beta.cure) == length(beta.par2))) {
      if (covariates.condition == 'equal') {
        cxx <- p1x
      } else if (covariates.condition == 'independent') {
        cxx <- mvrnorm(size, mu = mu.cure, Sigma = sigma.cure)
      } else {
        cxx <- mvrnorm(size, mu = mu.cure, Sigma = sigma.cure)
      }
    } else {
      cxx <- mvrnorm(size, mu = mu.cure, Sigma = sigma.cure)
    }
    # generate event time from excess hazard
    qloglogistic <- function(x){
      par2.parameter*((1-x)/x)^(-1/par1.parameter)
    }
    if (model.type == 'mixture') {
      # cure rate
      cr <- (1+exp(-crossprod(beta.cure,t(cxx))))^(-1)
      mix1 <- 1-rbinom(size, 1, cr)
      u1 <- runif(size)
      e1 <- qloglogistic(u1)
    } else if (model.type == 'non-mixture') {
      # cure rate
      cr <- (1+exp(-crossprod(beta.cure,t(cxx))))^(-1)
      N <- rpois(size, -log(cr))
      mix1 <- I(N > 0)
      e1N <- list()
      for(i in 1:size) {
        u1 <- runif(N[i])
        e1N[[i]] <- qloglogistic(u1)
      }
      e1 <- sapply(e1N, min)
    } else if (model.type == 'time') {
      # cure time
      cr <- c(exp(crossprod(beta.cure,t(cxx))))
      u1 <- runif(size)
      e1 <- qloglogistic(u1)
      mix1 <- e1 <= cr
    }
  }
  # generate event time from general hazard
  if (general.dist == 'Weibull') {
    e2 <- rweibull(size, shape = general.par1, scale = general.par2)
  } else if (general.dist == 'log-normal') {
    e2 <- rlnorm(size, meanlog = general.par1, sdlog = general.par2)
  } else if (general.dist == 'log-logistic') {
    u2 <- runif(size)
    qloglogistic <- function(x){
      general.par2*((1-x)/x)^(-1/general.par1)
    }
    e2 <- qloglogistic(u2)
  }
  # minimum of e1, e2 for each individual
  e12 <- ifelse(e1 < e2, e1, e2)
  # observed event time (composed of excess hazard and general hazard)
  e3 <- ifelse(mix1, e12, e2)
  # generate censoring time
  if (cen.dist == 'Weibull') {
    c1 <- rweibull(size, shape = exp(cen.par1), scale = exp(cen.par2))
  } else if (cen.dist == 'log-normal') {
    c1 <- rlnorm(size, meanlog = cen.par1, sdlog = cen.par2)
  } else if (cen.dist == 'log-logistic') {
    u3 <- runif(size)
    qloglogistic <- function(x){
      exp(cen.par2)*((1-x)/x)^(-1/exp(cen.par1))
    }
    c1 <- qloglogistic(u3)
  }
  # observed time
  cen1 <- 1*(c1 < e3)
  cen4 <- rbinom(size, 1, c4_from_c234*(e3 < c1))
  cen2 <- (1*(e3 < c1) - cen4)*(e3 == e2)
  cen3 <- (1*(e3 < c1) - cen4)*(e3 == e1)
  obs <- cen1*c1 + cen2*e2 + cen3*e1 + cen4*e3
  #cen1+cen2+cen3+cen4
  # general hazard
  if (general.dist == 'Weibull') {
    gs <- pweibull(obs, shape = general.par1, scale = general.par2, lower.tail = FALSE)
    gf <- dweibull(obs, shape = general.par1, scale = general.par2)
    gh <- gf/gs
  } else if (general.dist == 'log-normal') {
    gs <- plnorm(obs, meanlog = general.par1, sdlog = general.par2, lower.tail = FALSE)
    gf <- dlnorm(obs, meanlog = general.par1, sdlog = general.par2)
    gh <- gf/gs
  } else if (general.dist == 'log-logistic') {
    gs <- 1/(1+(obs/general.par2)^(general.par1))
    gf <- ((general.par1/general.par2)*(obs/general.par2)^(general.par1-1))/(1+(obs/general.par2)^(general.par1))^2
    gh <- gf/gs
  }
  check.gh <- is.na(gh) | is.nan(gh) | is.infinite(gh)
  gh.adjust <- ifelse(check.gh, max(gh), gh)
  # data summary
  sim.data <- list(obs = obs,
                   cen1 = cen1,
                   cen2 = cen2,
                   cen3 = cen3,
                   cen4 = cen4,
                   gh = gh.adjust,
                   p1x = p1x,
                   p2x = p2x,
                   cxx = cxx,
                   excess.event = e1,
                   general.event = e2,
                   event = e3,
                   censor.time = c1)
  return(sim.data)
}
