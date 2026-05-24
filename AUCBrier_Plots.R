############################################################
## FINAL PLOTTING CODE FOR AUC/BRIER
############################################################

############################################################
## LOAD PACKAGES
############################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

############################################################
## LOAD RESULTS
############################################################
lm_metrics  <- readRDS("LM_metrics_simulation_100.rds")
jm_metrics  <- readRDS("JMmetrics_simulation_100.rds")

lm_fits  <- readRDS("lmfit_simulation_riz_100.rds")
lm_fit1  <- lm_fits[[1]]

############################################################
## COMBINE REPLICATES
############################################################
lm_df  <- bind_rows(lm_metrics)
jm_df  <- bind_rows(jm_metrics)
all_df <- bind_rows(lm_df, jm_df)

############################################################
## AVERAGE OVER REPLICATES
############################################################
avg_df <- all_df %>%
  group_by(model, window, landmark) %>%
  summarise(
    AUC    = mean(AUC, na.rm = TRUE),
    Brier = mean(Brier, na.rm = TRUE),
    .groups = "drop"
  )

############################################################
## Table of statistics
############################################################
stats_df <- all_df %>%
  group_by(model, window, landmark) %>%
  summarise(
    AUC_mean    = mean(AUC, na.rm = TRUE),
    AUC_sd    = sd(AUC, na.rm = TRUE),
    AUC_MCSE = AUC_sd/ sqrt(100),
    Brier_mean = mean(Brier, na.rm = TRUE),
    Brier_sd = sd(Brier, na.rm = TRUE),
    Brier_MCSE = Brier_sd/ sqrt(100),
    .groups = "drop"
  )

LM_fixed <- all_df %>%
  filter(model == "LM") %>%
  mutate(row_id = row_number()) %>%
  mutate(
    cycle18 = ceiling(row_id / 18),                # one replicate per 18 rows
    pos18   = row_id - 18 * (cycle18 - 1)           # position within cycle
  ) %>%
  filter(
    (pos18 >= 6  & pos18 <= 10) |                   # keep second 5
      (pos18 >= 15 & pos18 <= 18)                     # keep second 4
  ) %>%
  mutate(
    replicate = cycle18
  ) %>%
  dplyr::select(-row_id, -cycle18, -pos18)

JM_fixed <- all_df %>%
  filter(model == "JM") %>%
  mutate(
    replicate = ceiling(row_number() / 9)
  )
all_df_fixed <- bind_rows(LM_fixed, JM_fixed)

paired_df <- all_df_fixed %>%
  pivot_wider(
    names_from  = model,
    values_from = c(AUC, Brier)
  )

paired_df <- paired_df %>%
  mutate(
    delta_AUC    = AUC_JM - AUC_LM,
    delta_Brier = Brier_JM - Brier_LM
  )

diff_summary <- paired_df %>%
  group_by(window, landmark) %>%
  summarise(
    delta_AUC_mean   = mean(delta_AUC),
    delta_AUC_sd     = sd(delta_AUC),
    delta_AUC_MCSE   = delta_AUC_sd / sqrt(n()),
    
    delta_Brier_mean = mean(delta_Brier),
    delta_Brier_sd   = sd(delta_Brier),
    delta_Brier_MCSE = delta_Brier_sd / sqrt(n()),
    
    n_reps = n(),
    .groups = "drop"
  )
############################################################
## Make data for tables
############################################################
table_data  <- readRDS("simulation_riz_100.rds")
table_data <- table_data[[1]]
table_data <- table_data$Long1 %>%
  group_by(id) %>%
  summarise(
    eventtime = first(eventtime),
    status = first(status),
    .groups = "drop"
  )

############################################################
## PLOT AUC — s = 5
############################################################
p_auc_s5 <- ggplot(
  avg_df %>% filter(window == 5),
  aes(x = landmark, y = AUC, linetype = model,color = model)
) +
  scale_x_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 5)) +
  scale_y_continuous(
    limits = c(0.5, 0.65),
    breaks = seq(0.5, 0.65, by = .05)
  ) +
  geom_line() +
  scale_color_manual(
    values = c(
      "JM" = "red",
      "LM" = "blue"
    )
  ) +
  labs(
    title = "AUC (window = 5)",
    x = "Landmark time",
    y = "AUC",
    linetype = "Model",
    color = "Model"
  ) +
  theme_bw()

landmarks <- c(0, 5, 10, 15,20)
window <- 5
landmark_table <- map_df(landmarks, function(s) {
  risk_set <- table_data %>%
    filter(eventtime > s)
  events <- sum(risk_set$eventtime <= s + window & risk_set$status == 1)
  person_years <- sum(pmin(risk_set$eventtime, s + window) - s)
  n_participants <- nrow(risk_set)
  tibble(
    landmark = s,
    window = paste0("[", s, ", ", s + window, "]"),
    participants = n_participants,
    events = events,
    person_years = person_years
  )
})
t_auc_s5 <- ggplot(data=landmark_table, aes(x = landmark)) +
  geom_text(aes(y="Participants", label=participants),size = 3.5) +
  geom_text(aes(y="Events", label=events),size = 3.5) +
  scale_y_discrete(limits=rev) +
  scale_x_continuous(breaks = c(0, 5, 10, 15,20),limits = c(0, 20)) +
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
        panel.grid = element_blank(), strip.text = element_blank())

############################################################
## PLOT AUC — s = 10
############################################################
p_auc_s10 <- ggplot(
  avg_df %>% filter(window == 10),
  aes(x = landmark, y = AUC, linetype = model,color = model)
) +
  scale_x_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 5)) +
  scale_y_continuous(
    limits = c(0.5, 0.65),
    breaks = seq(0.5, 0.65, by = .05)
  ) +
  geom_line() +
  scale_color_manual(
    values = c(
      "JM" = "red",
      "LM" = "blue"
    )
  ) +
  labs(
    title = "AUC (window = 10)",
    x = "Landmark time",
    y = "AUC",
    linetype = "Model",
    color = "Model"
  ) +
  theme_bw()

landmarks <- c(0, 5, 10, 15)
window <- 10
landmark_table <- map_df(landmarks, function(s) {
  risk_set <- table_data %>%
    filter(eventtime > s)
  events <- sum(risk_set$eventtime <= s + window & risk_set$status == 1)
  person_years <- sum(pmin(risk_set$eventtime, s + window) - s)
  n_participants <- nrow(risk_set)
  tibble(
    landmark = s,
    window = paste0("[", s, ", ", s + window, "]"),
    participants = n_participants,
    events = events,
    person_years = person_years,
  )
})
landmark_table <- bind_rows(
  landmark_table,
  tibble(
    landmark = Inf,
    window = "",
    participants = NA_real_,
    events = NA_real_,
    person_years = NA_real_
  )
)

t_auc_s10 <- ggplot(data=landmark_table, aes(x = landmark)) +
  geom_text(aes(y="Events", label=events),size = 3.5) +
  geom_text(aes(y="Participants", label=participants),size = 3.5) +
  scale_y_discrete(limits=rev) +
  scale_x_continuous(breaks = c(0, 5, 10, 15),limits = c(0, 20)) +
  labs(y = NULL, x = NULL) +
  theme_minimal() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
        panel.grid = element_blank(), strip.text = element_blank())


############################################################
## PLOT BRIER — s = 5
############################################################
p_brier_s5 <- ggplot(
  avg_df %>% filter(window == 5),
  aes(x = landmark, y = Brier, linetype = model,color = model)
) +
  scale_x_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 5)) +
  scale_y_continuous(
    limits = c(0.02, 0.08),
    breaks = seq(0.02, 0.08, by = .02)
  ) +
  geom_line() +
  scale_color_manual(
    values = c(
      "JM" = "red",
      "LM" = "blue"
    )
  ) +
  labs(
    title = "Brier Score (window = 5)",
    x = "Landmark time",
    y = "Brier score",
    linetype = "Model",
    color = "Model"
  ) +
  theme_bw()

############################################################
## PLOT BRIER — s = 10
############################################################
p_brier_s10 <- ggplot(
  avg_df %>% filter(window == 10),
  aes(x = landmark, y = Brier, linetype = model, colour = model)
) +
  scale_x_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 5)) +
  scale_y_continuous(
    limits = c(0.02, 0.08),
    breaks = seq(0.02, 0.08, by = .02)
  ) +
  geom_line() +
  scale_color_manual(
    values = c(
      "JM" = "red",
      "LM" = "blue"
    )
  ) +
  labs(
    title = "Brier Score (window = 10)",
    x = "Landmark time",
    y = "Brier score",
    linetype = "Model",
    color = "Model"
  ) +
  theme_bw()

############################################################
## DISPLAY PLOTS 
############################################################
library(patchwork)
((p_auc_s5/t_auc_s5+ plot_layout(heights = c(7, 2))) | (p_auc_s10/t_auc_s10+ plot_layout(heights = c(7, 2)))) /
  ((p_brier_s5/t_auc_s5+ plot_layout(heights = c(7, 2))) | (p_brier_s10/t_auc_s10+ plot_layout(heights = c(7, 2))))

p_auc_s5  <- p_auc_s5  + theme(axis.title.x = element_blank())
p_auc_s10 <- p_auc_s10 + theme(axis.title.x = element_blank())
p_auc_s10    <- p_auc_s10    + theme(axis.title.y = element_blank())
p_brier_s10  <- p_brier_s10  + theme(axis.title.y = element_blank())
t_auc_s10    <- t_auc_s10    + theme(axis.title.y = element_blank())
p_brier_s5  <- p_brier_s5  + theme(axis.title.x = element_blank())
p_brier_s10 <- p_brier_s10 + theme(axis.title.x = element_blank())
t_auc_s10 <- t_auc_s10 +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

combined_plot <-
  (
    (p_auc_s5 |p_auc_s10)
    /
      ((p_brier_s5 / t_auc_s5 + plot_layout(heights = c(7, 2))) |
         (p_brier_s10 / t_auc_s10 + plot_layout(heights = c(7, 2))))
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_plot &
  theme(
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
