# Generate data with B-spline model

simulate_data <- function(seed){
  set.seed(seed)
  
  library(MASS)
  library(splines)
  
  # -----------------------------
  # simulation settings
  # -----------------------------
  n <- 1438
  K <- 15
  t.max <- 30
  
  # -----------------------------
  # longitudinal parameters (from Lovblom paper)
  # -----------------------------
  betas <- c(
    "(Intercept)" = 2.4185,
    "B1" = 0.0497,
    "B2" = 0.0653,
    "B3" = 0.4152,
    "B4" = 0.4280,
    "B5" = 0.2390
  )
  
  sigma.y <- 0.5866
  
  # -----------------------------
  # survival parameters
  # -----------------------------
  gammas <- c(
    "(Intercept)" = -6.1133,
    "group1" = -0.2279
  )
  
  alpha <- 0.3928
  phi <- exp(0.0043)
  
  # -----------------------------
  # random effects
  # -----------------------------
  tau_vec <- c(0.5558, 0.4059, 0.7922, 1.7166, 1.8084, 1.2034)
  D <- diag(tau_vec^2)
  
  # -----------------------------
  # spline setup
  # -----------------------------
  Bkn <- c(0, 30)
  kn <- c(6, 10)
  
  # -----------------------------
  # baseline covariate
  # -----------------------------
  Z1 <- rep(0:1, each = n/2)
  
  # longitudinal times
  times <- c(replicate(n, c(0, sort(runif(K - 1, 0, t.max)))))
  
  DF <- data.frame(
    tij = times,
    Z1 = rep(Z1, each = K)
  )
  
  Bs <- bs(DF$tij, knots = kn, Boundary.knots = Bkn, degree = 3)
  
  X <- cbind(1, Bs)
  Z <- cbind(1, Bs)
  
  W <- cbind(1, Z1)
  
  # -----------------------------
  # random effects
  # -----------------------------
  b <- mvrnorm(n, rep(0, nrow(D)), D)
  
  # -----------------------------
  # longitudinal outcomes
  # -----------------------------
  id_full <- rep(1:n, each = K)
  
  eta.y <- as.vector(X %*% betas + rowSums(Z * b[id_full, ]))
  y <- rnorm(n * K, eta.y, sigma.y)
  
  # -----------------------------
  # event times
  # -----------------------------
  eta.t <- as.vector(W %*% gammas)
  adminCens <- 30
  
  invS <- function(t, u, i) {
    h <- function(s) {
      BS <- bs(s, knots = kn, Boundary.knots = Bkn, degree = 3)
      XX <- cbind(1, BS)
      ZZ <- cbind(1, BS)
      f1 <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
      exp(log(phi) + (phi - 1)*log(s) + eta.t[i] + alpha * f1)
    }
    integrate(h, lower = 1e-08, upper = t)$value + log(u)
  }
  
  u <- runif(n)
  trueTimes <- numeric(n)
  
  for (i in 1:n) {
    Root <- try(
      uniroot(invS, interval = c(1e-08, adminCens), u = u[i], i = i)$root,
      silent = TRUE
    )
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else Inf
  }
  
  # censoring
  eventtime <- pmin(trueTimes, adminCens)
  status <- as.integer(trueTimes <= adminCens)
  
  # -----------------------------
  # longitudinal truncation
  # -----------------------------
  ind <- times <= rep(eventtime, each = K)
  
  longdat <- DF[ind, ]
  longdat$id <- id_full[ind]
  longdat$eventtime <- rep(eventtime, each = K)[ind]
  longdat$status <- rep(status, each = K)[ind]
  longdat$Yij_1 <- y[ind]
  
  longdat <- longdat[, c("id", "Z1", "eventtime", "status", "tij", "Yij_1")]
  rownames(longdat) <- NULL
  
  # -----------------------------
  # event dataset
  # -----------------------------
  eventdat <- data.frame(
    id = 1:n,
    Z1 = Z1,
    eventtime = eventtime,
    status = status
  )
  
  # -----------------------------
  # return simjm-style object
  # -----------------------------
  return(list(
    Event = eventdat,
    Long1 = longdat
  ))
}


R <- 100
dat_list <- vector("list", R)

for (r in 1:R) {
  dat_list[[r]] <- simulate_data(seed = r)
  print(r)
}

#saveRDS(dat_list, file = "simulation_riz_100.rds")