############################################################
# Fit joint model and compute metrics
############################################################

# Fit joint model and compute metrics
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

# calculate performance metrics
performance_jm = 
  function(fit, dat, windows) {
    
    library(JMbayes2)
    
    newdata <- merge(
      dat$Long1,
      dat$Event[, c("id", "eventtime", "status")],
      by = "id",
      all.x = TRUE
    )
    
    # remove duplicated variables created in merging
    newdata$eventtime <- newdata$eventtime.y
    newdata$status <- newdata$status.y
    newdata$eventtime.x <- NULL
    newdata$eventtime.y <- NULL
    newdata$status.x <- NULL
    newdata$status.y <- NULL
    
    out <- data.frame()
    
    for (w in windows) {
      if(w == 5){
        lms = c(0,5,10,15,20)
      }else
        if(w == 10){
          lms = c(0,5,10,15)
        }
      for (lm in lms) {
        
        # AUC via ROC
        roc_jm <- tvROC(
          fit,
          newdata = newdata,
          Tstart = lm,
          Dt = w,
          type_weights = "IPCW"
        )
        auc_val <- as.numeric(tvAUC(roc_jm)$auc)
        
        # Brier score
        brier <- tvBrier(
          fit,
          newdata = newdata,
          Tstart = lm,
          Dt = w,
          type_weights = "IPCW"
        )
        
        out <- rbind(
          out,
          data.frame(
            model   = "JM",
            window = w,
            landmark = lm,
            AUC = auc_val,
            Brier = brier$Brier
          )
        )
      }
    }
    rownames(out) <- NULL
    return(out)
  }

# load simulated dataset
dat_list <- readRDS("simulation_riz_100.rds")
windows <- c(5,10)

perf_list_jm_windows <- vector("list", length(dat_list))
fit_list_jm_windows <- vector("list", length(dat_list))

for(r in seq_along(dat_list)) {
  print(r)
  dat <- dat_list[[r]]
  fit <- fit_joint_model(dat)
  fit_list_jm_windows[[r]] <- fit
  perf_list_jm_windows[[r]] <- performance_jm(fit, dat, windows)
}

#saveRDS(perf_list_jm_windows, file = "JMmetrics_simulation_riz_100replaced59.rds")
#saveRDS(fit_list_jm_windows, file = "jmfit_simulation_riz_100replaced59.rds")