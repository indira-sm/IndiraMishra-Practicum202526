# Checking Assumptions

# PH Assumption check
ph_test <- cox.zph(cox_fit)
print(ph_test)
plot(ph_test)

# LMM Assumptions check
plot(lme_fit, resid(., type = "p") ~ fitted(.))

jm_fit <- jm(cox_fit, lme_fit, time_var = "tij")

qqnorm(resid(lme_fit))
qqline(resid(lme_fit))

ranef_vals <- ranef(lme_fit)

qqnorm(ranef_vals[,1])  # random intercept
qqline(ranef_vals[,1])

qqnorm(ranef_vals[,2])  # random slope
qqline(ranef_vals[,2])