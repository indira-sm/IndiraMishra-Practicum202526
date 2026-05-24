# make plots for calibration metrics


library(survival)
library(JMbayes2)
library(ggplot2)
library(patchwork)
library(JM)
library(dynamicLM)

############################################################
### CREATE CALIBRATION PLOT FOR LANDMARK = 0,5,10,15,20 WINDOW = 5,10
############################################################

############################################################
### DATA AND COORESPONDING JOINT MODEL
############################################################
jm_datas = readRDS("simulation_100 copy.rds") # linear model
jm_datas = readRDS("simulation_riz_100 copy.rds") # spline model
jm_datas_all = jm_datas[[1]]

lme_fit <- lme(Yij_1 ~ tij, random = ~ tij | id, data = jm_datas_all$Long1)

# survival submodel
cox_fit <- coxph(Surv(eventtime, status) ~ Z1, data = jm_datas_all$Event,
                 x = TRUE)
# joint model
jm_fit <- jm(cox_fit,lme_fit,time_var = "tij")

############################################################
### REFIT LANDMARK MODEL TO SINGLE CHOSEN DATA
############################################################
dat_lm <- jm_datas_all$Long1 %>%
  arrange(id, tij) %>%
  group_by(id) %>%
  mutate(
    tstart = tij,
    tstop  = lead(tij),
    tstop  = if_else(is.na(tstop), eventtime, pmin(tstop, eventtime)),
    event  = as.integer(status == 1 & eventtime <= tstop),
    logAER = Yij_1
  ) %>%
  filter(tstart < tstop) %>%
  ungroup() %>%
  dplyr::select(-tij, -Yij_1)

outcome <- list(time = "eventtime", status = "status")

covars <- list(fixed   = c("Z1"),varying = c("logAER"))
w5 <- 5
lms5 <- c(0,5,10,15,20)

lmdata5 <- stack_data(
  dat_lm,
  outcome = outcome,
  lms     = lms5,
  w       = w5,
  covs = covars,
  format  = "long",
  id      = "id",
  rtime   = "tstart"
)
form <- as.formula("Surv(LM, eventtime, status) ~
     Z1 + logAER + cluster(id)")

lm_fit5 <- dynamic_lm(lmdata5, as.formula(form), "coxph", x = TRUE) 

w10 <- 10                  # Predict the 5-year outcome of transplant
lms10 <- c(0,5,10,15)  # Risk assessment time points (every year for 4 years)

lmdata10 <- stack_data(
  dat_lm,
  outcome = outcome,
  lms     = lms10,
  w       = w10,
  covs = covars,
  format  = "long",
  id      = "id",
  rtime   = "tstart"
)
form <- as.formula("Surv(LM, eventtime, status) ~
     Z1 + logAER + cluster(id)")

lm_fit10 <- dynamic_lm(lmdata10, as.formula(form), "coxph", x = TRUE) 

############################################################
### JOINT MODEL
############################################################
plots = list()
windows = c(5,10)
plot = 1
for(w in windows){
  
  if(w == 5){
    lms = c(0,5,10,15,20)
    lm_fit = lm_fit5
  }else
    if(w == 10){
      lms = c(0,5,10,15)
      lm_fit = lm_fit10
    }
  p <- predict(lm_fit)
  
  for(lm in lms){
    jm_cal = calibration_plot(
      jm_fit,
      newdata = jm_datas_all$Long1,
      Tstart = lm, 
      Thoriz = lm+w,
      plot = FALSE)
    
    N_bins = 10
    
    cuts <- quantile(
      jm_cal$predicted,
      probs = seq(0, 1, length.out = N_bins + 1),
      na.rm = TRUE)
    
    bins <- cut(
      jm_cal$predicted,
      breaks = unique(cuts),
      include.lowest = TRUE)
    
    jm_pred <- tapply(jm_cal$predicted, bins, mean, na.rm = TRUE)
    jm_obs  <- tapply(jm_cal$observed, bins, mean, na.rm = TRUE)
    
    jm_abs_err <- abs(jm_pred - jm_obs)
    jm_met <- list(
      ICI = mean(jm_abs_err),
      E50 = median(jm_abs_err),
      E90 = quantile(jm_abs_err, 0.90, names = FALSE))
    
    ############################################################
    ### LANDMARK
    ############################################################
    lm_cal <- calplot(
      list("LM" = p),
      times = lm,
      main = "",
      method = 'quantile', q = 10
    )
    
    lm_pred_surv = lm_cal[[1]]$plotFrames$LM$Pred
    lm_obs_surv = lm_cal[[1]]$plotFrames$LM$Obs
    
    lm_pred = 1-lm_pred_surv 
    lm_obs = 1-lm_obs_surv
    
    lm_abs_err <- abs(lm_pred - lm_obs)
    lm_met <- list(
      ICI = mean(lm_abs_err),
      E50 = median(lm_abs_err),
      E90 = quantile(lm_abs_err, 0.90, names = FALSE))
    
    ############################################################
    ### PLOT
    ############################################################
    
    df <- rbind(
      data.frame(
        pred = jm_pred,
        obs  = jm_obs,
        model = "Joint model"),
      data.frame(
        pred = lm_pred,
        obs  = lm_obs,
        model = "Landmark model"))
    
    title =paste0("Landmark = ", lm, ", Window = ", w)
    jm_label = sprintf("Joint model\n    ICI = %.3f\n    E50 = %.3f\n    E90 = %.3f",
                       jm_met$ICI, jm_met$E50, jm_met$E90)
    lm_label = sprintf("Landmark model\n    ICI = %.3f\n    E50 = %.3f\n    E90 = %.3f",
                       lm_met$ICI, lm_met$E50, lm_met$E90)
    
    
    pl <- ggplot(df, aes(x = pred, y = obs,
                         color = model, shape = model)) +
      geom_point(size = 2) +
      geom_abline(intercept = 0, slope = 1,linetype = "dashed", color = "gray") +
      scale_color_manual(
        values = c("Joint model" = "red","Landmark model" = "blue"),
        labels = c("Joint model" = jm_label,"Landmark model" = lm_label)) +
      scale_shape_manual(
        values = c("Joint model" = 16,"Landmark model" = 1),
        labels = c("Joint model" = jm_label,"Landmark model" = lm_label)) +
      coord_cartesian(xlim = c(0, 0.2),ylim = c(0, 0.2)) +
      labs(
        x = "Predicted Probabilities",
        y = "Observed Probabilities",
        color = NULL,
        shape = NULL,
        title = title) +
      theme_minimal() + 
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = c(0.9, 0.05),
            legend.justification = c("right", "bottom"),
      )
    plots[[plot]] = pl
    print(plot)
    plot = plot + 1
  }
}

library(patchwork)
(plots[[1]] | plots[[2]] | plots[[3]])/
  (plots[[4]] | plots[[5]] | plots[[6]])/
  (plots[[7]] | plots[[8]] | plots[[9]])