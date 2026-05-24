# pilot simulation for sample size calculation

library(nlme)
library(survival)
library(JMbayes2)

#--------------------------------------------------
# fit joint model
fit_joint_model = function(dat){
  library(nlme)
  library(splines)
  library(survival)
  library(JMbayes2)
  
  # longitudinal submodel
  lme_fit <- lme(Yij_1 ~ tij, random = ~ tij | id, data = dat$Long1)
  
  # survival submodel
  cox_fit <- coxph(Surv(eventtime, status) ~ Z1, data = dat$Event,
                   x = TRUE)
  # joint model
  jm_fit <- jm(cox_fit,lme_fit,time_var = "tij")
  
  return(jm_fit)
}

# calculate performance metrics for joint model fit
performance_jm = 
  function(fit, dat, s_vals, t_max = 25) {
    
    library(JMbayes2)
    
    # build newdata
    newdata <- merge(
      dat$Long1,
      dat$Event[, c("id", "eventtime", "status")],
      by = "id",
      all.x = TRUE
    )
    
    newdata$eventtime <- newdata$eventtime.y
    newdata$status <- newdata$status.y
    newdata$eventtime.x <- NULL
    newdata$eventtime.y <- NULL
    newdata$status.x <- NULL
    newdata$status.y <- NULL
    
    out <- data.frame()
    
    for (s in s_vals) {
      for (t in (s + 1):t_max) {
        
        Dt <- t - s
        
        # AUC via ROC
        roc_jm <- tvROC(
          fit,
          newdata = newdata,
          Tstart = s,
          Dt = Dt,
          type_weights = "IPCW"
        )
        
        auc_val <- as.numeric(tvAUC(roc_jm)$auc)
        
        out <- rbind(
          out,
          data.frame(
            model   = "JM",
            s       = s,
            t       = t,
            Dt      = Dt,
            AUC     = auc_val
          )
        )
      }
    }
    
    rownames(out) <- NULL
    return(out)
  }


#--------------------------------------------------
# PILOT
#--------------------------------------------------
# run pilot
run_pilot_jm <- function(dat_list, lms, w_set) {
  perf_list <- vector("list", length(dat_list))
  
  for (r in seq_along(dat_list)) {
    message(sprintf("Pilot repetition %d of %d", r, length(dat_list)))
    
    dat <- dat_list[[r]]
    fit <- fit_joint_model(dat)
    perf_list[[r]] <- performance_jm(fit, dat, lms, w_set)
    perf_list[[r]]$rep <- r
  }
  
  perf_all <- do.call(rbind, perf_list)
  rownames(perf_all) <- NULL
  return(perf_all)
}

#--------------------------------------------------
# Monte Carlo sample size calculation
#--------------------------------------------------

# eps = desired MCSE threshold
mc_reps_required <- function(sd_hat, eps = 0.005) {
  ceiling((sd_hat / eps)^2)
}

# pooled SD across all AUC values
calc_R_pooled <- function(perf_all, eps = 0.005) {
  auc_vals <- perf_all$AUC[is.finite(perf_all$AUC)]
  sd_hat <- sd(auc_vals)
  R_req <- mc_reps_required(sd_hat, eps)
  
  list(
    sd_hat = sd_hat,
    eps = eps,
    R_required = R_req
  )
}

# conservative SD: calculate SD at each (landmark, window),
# then use the largest one
calc_R_worst_case <- function(perf_all, eps = 0.005) {
  split_auc <- split(
    perf_all$AUC,
    interaction(perf_all$tLM, perf_all$window, drop = TRUE)
  )
  
  sd_by_setting <- sapply(split_auc, function(x) sd(x[is.finite(x)]))
  sd_by_setting <- sd_by_setting[is.finite(sd_by_setting)]
  
  sd_worst <- max(sd_by_setting)
  R_req <- mc_reps_required(sd_worst, eps)
  
  list(
    sd_by_setting = sd_by_setting,
    sd_worst = sd_worst,
    eps = eps,
    R_required = R_req
  )
}

#--------------------------------------------------
lms <- 0:4
w_set <- c(5, 10)
eps <- 0.005

# Load pilot dataset
dat_list <- readRDS("simulation1.rds")

# Run pilot
perf_all <- run_pilot_jm(dat_list, lms, w_set)

# Calculate required R using pooled SD
pooled_result <- calc_R_pooled(perf_all, eps = eps)
pooled_result

# Calculate required R using worst-case SD
worst_result <- calc_R_worst_case(perf_all, eps = eps)
worst_result

#--------------------------------------------------
# 6. summary table
#--------------------------------------------------
summarize_pilot_auc <- function(perf_all) {
  agg_mean <- aggregate(AUC ~ tLM + window, data = perf_all, FUN = mean, na.rm = TRUE)
  agg_sd   <- aggregate(AUC ~ tLM + window, data = perf_all, FUN = sd, na.rm = TRUE)
  
  out <- merge(agg_mean, agg_sd, by = c("tLM", "window"), suffixes = c("_mean", "_sd"))
  out
}

# Run
summary_auc <- summarize_pilot_auc(perf_all)
print(summary_auc)