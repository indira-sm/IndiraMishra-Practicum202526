############################################################
## fit landmark model and compute metrics
############################################################

library(dplyr)
library(dynamicLM)
library(survival)

# fit a landmark model to simulated data
fit_landmark_model <- function(dat, lm, w) {
  dat_lm <- dat$Long1 %>%
    arrange(id, tij) %>%
    group_by(id) %>%
    mutate(
      tstart = tij,
      tstop  = lead(tij),
      tstop  = if_else(is.na(tstop), eventtime, pmin(tstop, eventtime)),
      event  = as.integer(status == 1 & eventtime <= tstop),
      time   = tij,
      logAER = Yij_1
    ) %>%
    filter(tstart < tstop) %>%
    ungroup() %>%
    dplyr::select(-tij, -Yij_1)
  
  outcome <- list(time = "eventtime", status = "status")
  
  covars <- list(
    fixed   = c("Z1"),
    varying = c("logAER")
  )
  
  # stack data into landmark form
  lmdata <- stack_data(
    dat_lm,
    outcome = outcome,
    lms     = lm,
    w       = w,
    covs = covars,
    format  = "long",
    id      = "id",
    rtime   = "time"
  )
  
  form <- as.formula(
    "Surv(LM, eventtime, status) ~
     Z1 + logAER + cluster(id)"
  )
  
  dynamic_lm(lmdata, form, "coxph", x = TRUE)
}

# load simulated dataset
dat_list <- readRDS("simulation_100.rds")

windows <- c(5,10)
perf_list <- vector("list", length(dat_list))
fit_list_lm_windows <- vector("list", length(dat_list))

for (r in seq_along(dat_list)) {
  print(r)
  dat <- dat_list[[r]]
  res_r <- data.frame()
  
  for(w in windows){
    if(w == 5){
      lms = c(0,5,10,15,20)
    }else
      if(w == 10){
        lms = c(0,5,10,15)
      }
    
    # fit landmark model for this (s, w)
    lm_fit <- fit_landmark_model(dat, lms, w = w)
    fit_list_lm_windows[[r]] <- lm_fit
    
    # predict risk for this window
    preds <- predict(lm_fit)
    
    # score
    sc <- score(
      list(LM = preds),
      times = lms)
    
    briers <- sc$Brier$score %>%
      filter(model == "LM")
    
    # extract AUC / Brier
    auc_val <- sc$AUC$score$AUC
    brier_val <- briers$Brier
    LM <- sc$AUC$score$tLM
    
    res_r <- rbind(
      res_r,
      data.frame(
        model  = "LM",
        landmark  = LM,
        window     = w,
        AUC    = auc_val,
        Brier  = brier_val
      )
    )
  }
  perf_list[[r]] <- res_r
}

perf_listLM_windows = perf_list
#saveRDS(perf_listLM_windows, file = "LM_metrics_simulation_100_corrected.rds")
#saveRDS(fit_list_lm_windows, file = "lmfit_simulation_100_corrected.rds")
