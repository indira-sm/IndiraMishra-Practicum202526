# Generate data with linear model

# Longitudinal submodel: intercept/slope
# Survival submodel: Randomization group covariate only

library(simjm)
simulate_data = function(seed){
  set.seed(seed)
  
  # -----------------------------
  # simulate joint data
  # -----------------------------
  dat <- simjm(
    n = 1438, M = 1, # 1441 participants in DCCT, using one longitudinal marker
    fixed_trajectory = "linear",
    random_trajectory = "linear",
    assoc = 'etavalue', #log(AER)
    basehaz = 'weibull',
    
    betaLong_binary = 0,
    betaLong_continuous = 0,
    betaLong_intercept = 2.4103,
    betaLong_linear = 0.023,
    betaLong_aux = 0.6530,
    b_sd  = c(0.6513, 0.0714),
    b_rho = -0.1101,
    
    betaEvent_intercept = -5.9185,
    betaEvent_binary = -0.2343,
    betaEvent_assoc = 0.3315,
    betaEvent_aux = exp(0.0006),
    
    #match cohort:
    prob_Z1 = 0.5, # randomization proportion
    max_fuptime = 30 # maximum followup time
  )
  return(dat)
}

R <- 100
dat_list <- vector("list", R)

for (r in 1:R) {
  dat_list[[r]] <- simulate_data(seed = r)
  print(r)
}

#saveRDS(dat_list, file = "simulation_100.rds")