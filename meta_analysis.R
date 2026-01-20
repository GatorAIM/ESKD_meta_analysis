################################################################################
######    Systematic Review and Meta Analysis for ESRD Health Outcomes    ######
################################################################################


### ============================================================================
### PART I: Descriptive Analysis 
### ============================================================================
library(here)
library(readxl)
library(dplyr)
library(tidyr)
library(forestplot)
library(metamisc)
library(meta)
library(metafor)
library(stringr)
library(ggplot2)
library(forcats)

### ----------------------------------------------------------------------------
### Step 1: Generate the grouped histogram for types of prediction outcomes
### ----------------------------------------------------------------------------

df_hist <- read_excel(
  here("data", "extracted_data.xlsx"),
  sheet = "histogram"
)

df_hist$broad <- factor(df_hist$broad, levels = unique(df_hist$broad), ordered = TRUE)
df_hist <- df_hist %>% 
  group_by(broad) %>%                       
  mutate(granular = fct_reorder(granular,   
                                count,
                                .desc = TRUE)) %>% 
  ungroup()


ggplot(df_hist, aes(x = granular, y = count, fill = granular)) +
  geom_col(width = 0.4) +
  geom_text(aes(label = count),
            vjust =0.11,  size = 3.5) +  # adjust text position and size
  facet_wrap(~ broad, scales = "free_x") +
  labs(
    x = NULL,
    y = "Number of Studies"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text.x  = element_text(face = "bold", size = 12)   
  ) # 850 x 600


### ============================================================================
### PART II: Meta Analysis of Survival Outcomes 
### ============================================================================

### ----------------------------------------------------------------------------
### Step 1: Load meta analysis data with survival outcomes and preprocess data 
### ----------------------------------------------------------------------------

data <- read_excel(
  here("data", "extracted_data.xlsx"),
  sheet = "meta_mortality"
)

data$author_year <- sub("^(\\S+)\\s+(\\d{4})$", "\\1 et al. (\\2)", data$author)
df <- data

df <- df %>%
  mutate(
    se = if_else((!is.na(lb)) & (!is.na(ub)) & is.na(se),
                 (ub - lb)/(2* qnorm(0.975)), 
                 se)
  )

df <- df %>%
  mutate(
    lb = if_else(is.na(lb) & !is.na(se),
                 mean_metric - qnorm(0.975) * se,
                 lb),
    ub = if_else(is.na(ub) & !is.na(se),
                 mean_metric + qnorm(0.975) * se,
                 ub)
  )

df <- df %>%
  mutate(
    prediction_time = factor(
      prediction_time,
      levels = c("30d", "60d", "90d", "180d", "270d", "1yr", 
                 "15m",  "18m",  "21m", "2yr", 
                 "3yr", "4yr", "5yr", "6yr",  
                 "7yr",  "8yr",  "9yr",  "10yr")
    )
  )

df <- df %>%
  mutate(
    prediction_time_days = case_when(
      prediction_time == "30d" ~ 30,
      prediction_time == "60d" ~ 60,
      prediction_time == "90d" ~ 90,
      prediction_time == "180d" ~ 180,
      prediction_time == "270d" ~ 270,
      prediction_time == "1yr" ~ 365,
      prediction_time == "15m" ~ 30 * 15,
      prediction_time == "18m" ~ 30 * 18,
      prediction_time == "21m" ~ 30 * 21,
      prediction_time == "2yr" ~ 365 * 2,
      prediction_time == "3yr" ~ 365 * 3,
      prediction_time == "4yr" ~ 365 * 4,
      prediction_time == "5yr" ~ 365 * 5,
      prediction_time == "6yr" ~ 365 * 6,
      prediction_time == "7yr" ~ 365 * 7,
      prediction_time == "8yr" ~ 365 * 8,
      prediction_time == "9yr" ~ 365 * 9,
      prediction_time == "10yr" ~ 365 * 10,
      TRUE ~ NA_real_
    )
  )

df <- df %>%
  mutate(
    values = if_else(
      !is.na(lb) & !is.na(ub),
      # If lb and ub exist:
      paste0(
        sprintf("%.3f", mean_metric), " (",
        sprintf("%.3f", lb), ", ",
        sprintf("%.3f", ub), ")"
      ),
      # Otherwise just the mean:
      sprintf("%.3f", mean_metric)
    )
  )

# Define outcome groups 
`%notin%` <- Negate(`%in%`)
df <- df %>% 
  mutate(
    outcome_group = factor(
      case_when(
        (outcome_type %in% c("All-cause Mortality", "Survival")) 
        & (prediction_time %in% c("30d", "60d", "90d")) ~ "short_term",
        (outcome_type %in% c("All-cause Mortality", "Survival")) 
        & (prediction_time %in% c("180d", "270d", "1yr")) ~ "mid_term",
        (outcome_type %in% c("All-cause Mortality", "Survival")) 
        & (prediction_time %notin% c("30d", "60d", "90d","180d","270d","1yr")) ~ "long_term",
        (outcome_type == "Survival Time") ~ "surv_time"
      ), 
      levels = c("short_term", "mid_term" ,"long_term", "surv_time"),
      ordered = TRUE
    )
  )

df <- df %>%
  mutate(
    is_short = as.integer(outcome_group == "short_term"),
    is_medium = as.integer(outcome_group == "mid_term"),
    is_long = as.integer(outcome_group == "long_term"),
    is_survtime = as.integer(outcome_group == "surv_time")
  )

# Filter those with reported standard error for the following meta-analysis
df_filtered <- df %>%
  filter((!is.na(se)) & (metric != 'AUPRC')) %>%
  arrange(outcome_group, prediction_time, author_year)


df_filtered$index <- rownames(df_filtered)

df_filtered <- df_filtered %>%
  mutate(
    logit_metric = qlogis(mean_metric),                     # transform estimate
    se_logit  = se / (mean_metric * (1 - mean_metric))      # delta-method SE
  )

df_filtered <- df_filtered %>%
  group_by(author_year) %>%
  mutate(index_within_author = row_number()) %>%
  ungroup()


### ----------------------------------------------------------------------------
### Step 2: Generate the forest plot for meta analysis of survival outcomes (Vcov = I)
### ----------------------------------------------------------------------------
# Fit mixed effect model with clustering
res_surv <- metafor::rma.mv(
  yi     = logit_metric, #mean_metric,          
  V      = se_logit^2,                  
  mods   = ~ 0 + outcome_group,  
  random = ~ 1| author_year,  
  data   = df_filtered, 
  method = "REML",     
  slab   = paste0(author_year, strrep("\u00A0", index_within_author)) # invisible NBSPs
)

res_overall <- metafor::rma.mv(
  yi     = logit_metric,  
  V      = se_logit^2,                  
  random = ~ 1| author_year,  
  data   = df_filtered, 
  method = "REML",     
  slab   = paste0(author_year, strrep("\u00A0", index_within_author)) # invisible NBSPs
)


# Fit mixed effect model without clustering
res_surv_nc <- metafor::rma.mv(
  yi     = logit_metric,          
  V      = se_logit^2,                  
  mods   = ~ 0 + outcome_group, 
  random = ~ 1| index,  
  data   = df_filtered,  
  method = "REML",   
  slab   = paste0(author_year, strrep("\u00A0", index_within_author)) # invisible NBSPs
)

res_overall_nc <- metafor::rma.mv(
  yi     = logit_metric,     
  V      = se_logit^2,
  random = ~ 1| index,  
  data   = df_filtered,  
  method = "REML",     
  slab   = paste0(author_year, strrep("\u00A0", index_within_author))   # invisible NBSPs 
)

# Fit fixed effect model (no random effect)
res_surv_fixed <- metafor::rma.mv(
  yi     = logit_metric,          
  V      = se_logit^2,                  
  mods   = ~ 0 + outcome_group, 
  data   = df_filtered,  
  method = "REML",     
  slab   = paste0(author_year, strrep("\u00A0", index_within_author)) # invisible NBSPs
)


res_overall_fixed <- metafor::rma.mv(
  yi     = logit_metric,     
  V      = se_logit^2,                  
  data   = df_filtered,  
  method = "REML",     
  slab   = paste0(author_year, strrep("\u00A0", index_within_author))   # invisible NBSPs 
)


# Predict subgroup mean 
levs <- levels(df_filtered$outcome_group)   
newmods <- diag(length(levs))               
colnames(newmods) <- paste0("outcome_group", levs) 

# raw 
preds_group_raw <- predict(res_surv, newmods = newmods)  
preds_group_raw_nc <- predict(res_surv_nc, newmods = newmods) 
preds_group_raw_fixed <- predict(res_surv_fixed, newmods = newmods) 

preds_overall_raw <- predict(res_overall)
preds_overall_raw_nc <- predict(res_overall_nc)
preds_overall_raw_fixed <- predict(res_overall_fixed)

# expit transformed
preds_group <- predict(res_surv, newmods = newmods, transf = transf.ilogit)  
preds_group_nc <- predict(res_surv_nc, newmods = newmods, transf = transf.ilogit)
preds_group_fixed <- predict(res_surv_fixed, newmods = newmods, transf = transf.ilogit) 

preds_overall <- predict(res_overall, transf = transf.ilogit)
preds_overall_nc <- predict(res_overall_nc, transf = transf.ilogit)
preds_overall_fixed <- predict(res_overall_fixed, transf = transf.ilogit)

# Generate the forest plot
metafor::forest(
  res_surv$yi,res_surv$vi,
  xlim = c(0.1, 1.2),
  ylim = c(-21.3, 35),
  transf = transf.ilogit,
  ilab = res_surv$data$outcome,
  ilab.xpos = 0.5,
  at = seq(0.65, 1, by = 0.05),
  xlab = expression(bold("AUROC / C-statistics")),
  header = "Study", 
  digits = 3
)

# add summary effects
addpoly(preds_group_raw,
        transf = transf.ilogit,
        row = -1,
        digits = 3,
        mlab = expression(
          bold("Subgroup: ≤90 Days (M.E., cl.)"),
          bold("Subgroup: >90 Days & ≤1 Year (M.E., cl.)"),
          bold("Subgroup: >1 Year (M.E., cl.)"),
          bold("Subgroup: Survival Time (M.E., cl.)")
        ), 
        col = 'steelblue',
        border = 'steelblue',
        cex = 1)

addpoly(preds_overall_raw,
        transf = transf.ilogit,
        row = -6,
        digits = 3,
        mlab = expression(
          bold("Overall Effect (M.E., cl.)")
        ), 
        col = 'steelblue',
        border = 'steelblue',
        cex = 1)

addpoly(preds_group_raw_nc,
        transf = transf.ilogit,
        row = -8,
        digits = 3,
        mlab = expression(
          bold("Subgroup: ≤90 Days (M.E.)"),
          bold("Subgroup: >90 Days & ≤1 Year (M.E.)"),
          bold("Subgroup: >1 Year (M.E.)"),
          bold("Subgroup: Survival Time (M.E.)")
        ),
        col= 'darkred',
        border = 'darkred',
        cex = 1)

addpoly(preds_overall_raw_nc,
        transf = transf.ilogit,
        row = -13,
        digits = 3,
        mlab = expression(
          bold("Overall Effect (M.E.)")
        ), 
        col = 'darkred',
        border = 'darkred',
        cex = 1)

addpoly(preds_group_raw_fixed,
        transf = transf.ilogit,
        row = -15,
        digits = 3,
        mlab = expression(
          bold("Subgroup: ≤90 Days (F.E.)"),
          bold("Subgroup: >90 Days & ≤1 Year (F.E.)"),
          bold("Subgroup: >1 Year (F.E.)"),
          bold("Subgroup: Survival Time (F.E.)")
        ),
        col= 'olivedrab3',
        border = 'olivedrab3',
        cex = 1)

addpoly(preds_overall_raw_fixed,
        transf = transf.ilogit,
        row = -20,
        digits = 3,
        mlab = expression(
          bold("Overall Effect (F.E.)")
        ), 
        col = 'olivedrab3',
        border = 'olivedrab3',
        cex = 1)

# add subgroup separator
abline(h = 0)
abline(h = 7.5, lty = "dashed")
abline(h = 11.5, lty = "dashed")
abline(h = 23.5, lty = "dashed")
abline(h = -7)
abline(h = -14)

# add dashed lines for subgroup effect (clustered, mixed effect)
segments(x0 = preds_group$pred[1], y0 = 23.5, 
         x1 = preds_group$pred[1], y1 = 33, 
         lty = "dashed", col = 'steelblue')
segments(x0 = preds_group$pred[2], y0 = 11.5, 
         x1 = preds_group$pred[2], y1 = 23.5, 
         lty = "dashed", col = 'steelblue')
segments(x0 = preds_group$pred[3], y0 = 7.5, 
         x1 = preds_group$pred[3], y1 = 11.5, 
         lty = "dashed", col = 'steelblue')
segments(x0 = preds_group$pred[4], y0 = 0, 
         x1 = preds_group$pred[4], y1 = 7.5, 
         lty = "dashed", col = 'steelblue')
# add dotted line for overall effect 
segments(x0 = preds_overall$pred, y0 = -7, 
         x1 = preds_overall$pred, y1 = 0, 
         lty = "dotted", col = 'steelblue')

# add dashed lines for subgroup effect (unclustered, mixed effect)
segments(x0 = preds_group_nc$pred[1], y0 = 23.5, 
         x1 = preds_group_nc$pred[1], y1 = 33, 
         lty = "dashed", col = 'darkred')
segments(x0 = preds_group_nc$pred[2], y0 = 11.5, 
         x1 = preds_group_nc$pred[2], y1 = 23.5, 
         lty = "dashed", col = 'darkred')
segments(x0 = preds_group_nc$pred[3], y0 = 7.5, 
         x1 = preds_group_nc$pred[3], y1 = 11.5, 
         lty = "dashed", col = 'darkred')
segments(x0 = preds_group_nc$pred[4], y0 = 0, 
         x1 = preds_group_nc$pred[4], y1 = 7.5, 
         lty = "dashed", col = 'darkred')
# add dotted line for overall effect 
segments(x0 = preds_overall_nc$pred, y0 = -14, 
         x1 = preds_overall_nc$pred, y1 = -7, 
         lty = "dotted", col = 'darkred') 

# add dashed lines for subgroup effect (fixed effect)
segments(x0 = preds_group_fixed$pred[1], y0 = 23.5, 
         x1 = preds_group_fixed$pred[1], y1 = 33, 
         lty = "dashed", col = 'olivedrab3')
segments(x0 = preds_group_fixed$pred[2], y0 = 11.5, 
         x1 = preds_group_fixed$pred[2], y1 = 23.5, 
         lty = "dashed", col = 'olivedrab3')
segments(x0 = preds_group_fixed$pred[3], y0 = 7.5, 
         x1 = preds_group_fixed$pred[3], y1 = 11.5, 
         lty = "dashed", col = 'olivedrab3')
segments(x0 = preds_group_fixed$pred[4], y0 = 0, 
         x1 = preds_group_fixed$pred[4], y1 = 7.5, 
         lty = "dashed", col = 'olivedrab3')
# add dotted line for overall effect 
segments(x0 = preds_overall_fixed$pred, y0 = -21, 
         x1 = preds_overall_fixed$pred, y1 = -14, 
         lty = "dotted", col = 'olivedrab3') 

# add text for outcomes
text(x = 0.38, y = 33.7, expression(bold("Outcome (Prediction Horizon)")), 
     pos = 4, font = 2)

# add legend for meta-regression setting 
legend(x =  0.506,                   
       y = -14.5,
       inset  = 0.01,                 
       legend = c("Mixed-effects, clustered",
                  "Mixed-effects, un-clustered", 
                  "Fixed-effects"),
       col    = c("steelblue", "darkred", "olivedrab3"),
       pch    = 18,                   
       lty    = 2,                    
       cex = 0.8)                 

# add text on test of heterogeneity and fixed effect

text(x = 0.1, y = -5,
     expression(bold(
     Q[E] * "(28) = 3522.17, " * p * "< 0.0001" *
      ", " * Q[M] * "(4) = 5706.95, " * p * "< 0.0001")),
  pos = 4, cex = 0.9)

text(x = 0.1, y = -12,
     expression(bold(
       Q[E] * "(28) = 3522.17, " * p * "< 0.0001" *
         ", " * Q[M] * "(4) = 459.08, " * p * "< 0.0001")),
  pos = 4, cex = 0.9)

text(x = 0.1, y = -19,
     expression(bold(
      Q[E] * "(28) = 3522.17, " * p * "< 0.0001" *
         ", " * Q[M] * "(4) = 379618.60, " * p * "< 0.0001")),
  pos = 4, cex = 0.9)



## ------------------------------------------------------------------ ##
##  1)  FUNCTION: run once for any chosen ρ                           ##
## ------------------------------------------------------------------ ##
run_sensitivity_rho <- function(rho, dat){
  
  ## ---- build working V matrix ------------------------------------
  V <- vcalc(vi      = dat$se_logit^2,
             cluster = dat$author_year,
             obs     = dat$index,          # row-id within study
             data    = dat,
             rho     = rho)
  
  ## ---- mixed-effects meta-regression -----------------------------
  mod_sub <- rma.mv(yi     = dat$logit_metric,
                    V      = V,
                    mods   = ~ 0 + outcome_group,      # ‘-1’ type coding
                    random = ~ 1 | author_year,
                    data   = dat,
                    method = "REML")
  
  ## ---- overall model (same V) ------------------------------------
  mod_ov  <- rma.mv(yi     = dat$logit_metric,
                    V      = V,
                    random = ~ 1 | author_year,
                    data   = dat,
                    method = "REML")
  
  ## ---- grab coefficients (β ± SE) --------------------------------
  cf <- coef(summary(mod_sub))
  ## put into two tidy vectors
  beta  <- cf[,"estimate"];    se <- cf[,"se"]
  
  ## ---- predict subgroup and overall AUROCs -----------------------
  levs       <- levels(dat$outcome_group)
  newmods    <- diag(length(levs)); colnames(newmods) <- paste0("outcome_group", levs)
  pred_sub   <- predict(mod_sub, newmods = newmods, transf = transf.ilogit)
  pred_over  <- predict(mod_ov,  transf = transf.ilogit)
  
  ## ---- build single-row output -----------------------------------
  out <- tibble(
    rho                   = rho,
    beta_short            = beta[1],   se_short   = se[1],
    beta_mid              = beta[2],     se_mid   = se[2],
    beta_long             = beta[3],    se_long   = se[3],
    beta_survtime         = beta[4],    se_surv   = se[4],
    ## predicted means
    auc_short   = pred_sub$pred[levs == "short_term"],
    lb_short    = pred_sub$ci.lb[levs == "short_term"],
    ub_short    = pred_sub$ci.ub[levs == "short_term"],
    auc_mid     = pred_sub$pred[levs == "mid_term"],
    lb_mid      = pred_sub$ci.lb[levs == "mid_term"],
    ub_mid      = pred_sub$ci.ub[levs == "mid_term"],
    auc_long    = pred_sub$pred[levs == "long_term"],
    lb_long     = pred_sub$ci.lb[levs == "long_term"],
    ub_long     = pred_sub$ci.ub[levs == "long_term"],
    auc_surv    = pred_sub$pred[levs == "surv_time"],
    lb_surv     = pred_sub$ci.lb[levs == "surv_time"],
    ub_surv     = pred_sub$ci.ub[levs == "surv_time"],
    auc_overall = pred_over$pred,
    lb_overall  = pred_over$ci.lb,
    ub_overall  = pred_over$ci.ub
  )
  out
}

## ------------------------------------------------------------------ ##
##  2)  RUN FOR A GRID OF ρ VALUES                                   ##
## ------------------------------------------------------------------ ##
rhos       <- c(0, 0.2, 0.5, 0.8)
sens_tbl   <- bind_rows(lapply(rhos, run_sensitivity_rho, dat = df_filtered))

## ------------------------------------------------------------------ ##
##  3)  PRETTIFY & PRINT / SAVE                                      ##
## ------------------------------------------------------------------ ##
pretty_tbl <- sens_tbl %>%                         # round first
  mutate(across(where(is.numeric) & !matches("^rho$"), ~ round(., 3))) %>%
  transmute(
    `ρ`                          = rho,
    `β_short-term`               = str_glue("{beta_short} ({se_short})"),
    `β_mid-term`                 = str_glue("{beta_mid} ({se_mid})"),
    `β_long-term`                = str_glue("{beta_long} ({se_long})"),
    `β_survival-time`            = str_glue("{beta_survtime} ({se_surv})"),
    `Mean AUC_short`             = str_glue("{auc_short} ({lb_short}, {ub_short})"),
    `Mean AUC_mid`               = str_glue("{auc_mid} ({lb_mid}, {ub_mid})"),
    `Mean AUC_long`              = str_glue("{auc_long} ({lb_long}, {ub_long})"),
    `Mean AUC_survival-time`     = str_glue("{auc_surv} ({lb_surv}, {ub_surv})"),
    `Mean AUC_overall`           = str_glue("{auc_overall} ({lb_overall}, {ub_overall})")
  )

readr::write_csv(
  pretty_tbl,
  here("output", "tables", "sensitivity_rho_mortality_summary.csv")
)





### ----------------------------------------------------------------------------
### Step 3: Meta Regression on Prediction Horizon Measured in Days
### ----------------------------------------------------------------------------
df_reg<-df_filtered[(!is.na(df_filtered$prediction_time_days)),]

# Fit mixed effect model clustered at study level
res_survtime <- metafor::rma.mv(
  yi     = logit_metric,       
  V      = se_logit^2,                  
  mods   = ~ splines::ns(log(prediction_time_days), df = 3), 
  random = ~ 1| author_year, 
  data   = df_reg, 
  method = "REML",     
  slab   = paste0(author_year, strrep("\u00A0", index_within_author)) # invisible NBSPs 
)

# Fit un-clustered mixed-effect model
res_survtime_nc <- metafor::rma.mv(
  yi     = logit_metric,        
  V      = se_logit^2,                  
  mods   = ~ splines::ns(log(prediction_time_days), df = 3),
  random = ~ 1| index, 
  data   = df_reg, 
  method = "REML",     
  slab   = paste0(author_year, strrep("\u00A0", index_within_author)) # invisible NBSPs
)

# Fit fixed-effect model
res_survtime_fixed <- metafor::rma.mv(
  yi     = logit_metric,        
  V      = se_logit^2,                  
  mods   = ~ splines::ns(log(prediction_time_days), df = 3),
  data   = df_reg, 
  method = "REML",     
  slab   = paste0(author_year, strrep("\u00A0", index_within_author)) # invisible NBSPs
)



# Predict subgroup mean 
levs <- sort(unique(df_reg$prediction_time_days))  
newmods <- splines::ns(log(levs), df = 3)
#colnames(newmods) <- names(coef(res_survtime))  
colnames(newmods) <- c("ns(log(prediction_time_days), df = 3)1",
                       "ns(log(prediction_time_days), df = 3)2",
                       "ns(log(prediction_time_days), df = 3)3")


preds_reg_raw <- predict(res_survtime, newmods = newmods)  
preds_reg_raw_nc <- predict(res_survtime_nc, newmods = newmods) 
preds_reg_raw_fixed <- predict(res_survtime_fixed, newmods = newmods) 

preds_reg <- predict(res_survtime, newmods = newmods, transf = transf.ilogit) 
preds_reg_nc <- predict(res_survtime_nc, newmods = newmods, transf = transf.ilogit) 
preds_reg_fixed <- predict(res_survtime_fixed, newmods = newmods, transf = transf.ilogit) 

# Plot meta regression curves
# Levels: turn days into factors
day_levels <- c(30, 60, 90, 180, 270, 365, 1095, 1825)
day_labels <- as.character(day_levels)

df_points <- df_reg %>%
  mutate(weight = 1 / sqrt(se^2),
         point_size = 0.5 + 3 * (weight - min(weight)) / (max(weight) - min(weight)),
         days = log(prediction_time_days))

df_pred <- data.frame(
  days              = log(day_levels),
  pred_clustered    = preds_reg$pred,
  ci.lb_clustered   = preds_reg$ci.lb,
  ci.ub_clustered   = preds_reg$ci.ub,
  pred_nonclustered = preds_reg_nc$pred,
  ci.lb_nonclustered = preds_reg_nc$ci.lb,
  ci.ub_nonclustered = preds_reg_nc$ci.ub,
  pred_fixed = preds_reg_fixed$pred,
  ci.lb_fixed = preds_reg_fixed$ci.lb,
  ci.ub_fixed = preds_reg_fixed$ci.ub
)

df_long <- df_pred %>%
  pivot_longer(-days, names_to = c("metric", "model"), names_sep = "_") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(model = factor(model,
                        levels = c("clustered", "nonclustered", "fixed"),
                        labels = c("Mixed-effects, clustered", 
                                   "Mixed-effects, un-clustered",
                                   "Fixed-effects")))
color_dict = c("Mixed-effects, clustered"  = "steelblue",
               "Mixed-effects, un-clustered"  = "darkred",
               "Fixed-effects"  = "olivedrab3"
)
ggplot() +
  geom_point(data = df_points,
             aes(x = days, y = mean_metric, size = point_size),
             shape = 19, alpha = 0.65, colour = "grey30") +
  
  geom_ribbon(data = df_long,
              aes(x = days, ymin = ci.lb, ymax = ci.ub, fill = model),
              alpha = 0.2, show.legend = FALSE) +
  
  geom_line(data = df_long,
            aes(x = days, y = pred, colour = model),
            size = 1) +
  # Add horizontal overall lines
  geom_hline(yintercept = preds_overall$pred, colour = "steelblue", linetype = "dotted", size = 0.9) +
  geom_hline(yintercept = preds_overall_nc$pred, colour = "darkred", linetype = "dotted", size = 0.9) +
  geom_hline(yintercept = preds_overall_fixed$pred, colour = "olivedrab3", linetype = "dotted", size = 0.9) +
  # Annotate the overall lines
  annotate("text", x = 3.95, y = preds_overall$pred + 0.008,
           label = paste0("Overall Effect: ", round(preds_overall$pred, 3)),
           color = "steelblue", hjust = 1, size = 4) +
  annotate("text", x = 3.95, y = preds_overall_nc$pred - 0.008,
           label = paste0("Overall Effect: ", round(preds_overall_nc$pred, 3)),
           color = "darkred", hjust = 1, size = 4) +
  annotate("text", x = 3.95, y = preds_overall_fixed$pred - 0.008,
           label = paste0("Overall Effect: ", round(preds_overall_fixed$pred, 3)),
           color = "olivedrab3", hjust = 1, size = 4) +  
  scale_color_manual(values = color_dict) +
  scale_fill_manual(values  = color_dict) +
  
  #scale_linetype_manual(values = color_dict) +
  
  scale_y_continuous(name = "AUROC / C-statistics",
                     breaks = seq(0.65, 1, 0.05),
                     limits = c(0.65, 1)) +
  scale_x_continuous(breaks = seq(3.5, 7.5, 0.5),
                     limits = c(3.3, 7.6)) +  
  xlab("Prediction Horizons (Log(Days))") +
  guides(size = "none") +
  theme_minimal(base_size = 14) +
  theme(
    legend.title       = element_blank(),
    legend.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.5),
    legend.position    = c(0.15, 0.86),
    panel.grid.minor   = element_blank()
  ) # 1000 x 600



### ----------------------------------------------------------------------------
### Step 4: Generate the funnel plot for publication bias
### ----------------------------------------------------------------------------

#### Option 1: No trimming and filling

# Fit mixed effect model with clustering
res_surv <- metafor::rma.mv(
  yi     = logit_metric, #mean_metric,          
  V      = se_logit^2,                  
  mods   = ~ 0 + outcome_group,  
  random = ~ 1| author_year,  
  data   = df_filtered, 
  method = "REML",     
  slab   = paste0(author_year, strrep("\u00A0", index_within_author)) # invisible NBSPs
)

res_overall <- metafor::rma.mv(
  yi     = logit_metric,  
  V      = se_logit^2,                  
  random = ~ 1| author_year,  
  data   = df_filtered, 
  method = "REML",     
  slab   = paste0(author_year, strrep("\u00A0", index_within_author)) # invisible NBSPs
)


# Fit mixed effect model without clustering
res_surv_nc <- metafor::rma.mv(
  yi     = logit_metric,          
  V      = se_logit^2,                  
  mods   = ~ 0 + outcome_group, 
  random = ~ 1| index,  
  data   = df_filtered,  
  method = "REML",   
  slab   = paste0(author_year, strrep("\u00A0", index_within_author)) # invisible NBSPs
)

res_overall_nc <- metafor::rma.mv(
  yi     = logit_metric,     
  V      = se_logit^2,
  random = ~ 1| index,  
  data   = df_filtered,  
  method = "REML",     
  slab   = paste0(author_year, strrep("\u00A0", index_within_author))   # invisible NBSPs 
)

# Fit fixed effect model (no random effect)
res_surv_fixed <- metafor::rma.mv(
  yi     = logit_metric,          
  V      = se_logit^2,                  
  mods   = ~ 0 + outcome_group, 
  data   = df_filtered,  
  method = "REML",     
  slab   = paste0(author_year, strrep("\u00A0", index_within_author)) # invisible NBSPs
)


res_overall_fixed <- metafor::rma.mv(
  yi     = logit_metric,     
  V      = se_logit^2,                  
  data   = df_filtered,  
  method = "REML",     
  slab   = paste0(author_year, strrep("\u00A0", index_within_author))   # invisible NBSPs 
)


## ------------------------------------------------------------
## 1.  Run Egger tests on the *simple* rma() objects
## ------------------------------------------------------------

res_overall_nc2 <- rma(yi = logit_metric, vi = se_logit^2, data = df_filtered, method = "REML")
res_overall_fixed2 <- rma(yi = logit_metric, vi = se_logit^2, data = df_filtered, method = "FE")

egger_nc    <- regtest(res_overall_nc2,    model = "rma")  # panel E
egger_fixed <- regtest(res_overall_fixed2, model = "rma")  # panel F

## helper to format a line of text ----------------------------
fmt_egger <- function(obj){
  z   <- obj$zval
  p   <- obj$pval
  b0  <- plogis(obj$est)                 # intercept
  #se  <- obj$se[1]
  cil  <- plogis(obj$ci.lb)
  ciu  <- plogis(obj$ci.ub)
  sprintf(" Egger: z = %.2f, p = %.3g\n Bias intercept (95%% CI): %.3f (%.3f, %.3f)",
          z, p, b0, cil, ciu)
}

txt_E <- fmt_egger(egger_nc)
txt_F <- fmt_egger(egger_fixed)

## ------------------------------------------------------------
## 2.  Set up a 2×3 plotting grid
## ------------------------------------------------------------
oldpar <- par(no.readonly = TRUE)   # save user settings
par(mfrow = c(2, 3),  mar = c(4,4,3,1),
    cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.1)

## ------------------------------------------------------------
## 3.  Draw the six funnel plots
## ------------------------------------------------------------
funnel(res_surv,          xlab = "Residual Logit(AUROC/C-stats.)", main = "(A) Mixed-Effects (clust.)")
funnel(res_surv_nc,       xlab = "Residual Logit(AUROC/C-stats.)", main = "(B) Mixed-Effects (un-clust.)")
funnel(res_surv_fixed,    xlab = "Residual Logit(AUROC/C-stats.)", main = "(C) Fixed-Effects")
funnel(res_overall,       xlab = "Residual Logit(AUROC/C-stats.)", main = "(D) Random-Effects (clust.)")
funnel(res_overall_nc,    xlab = "Residual Logit(AUROC/C-stats.)", main = "(E) Random-Effects (un-clust.)")
mtext(txt_E, side = 3, line = -23, adj = 0, cex = 0.75)   # annotate panel E
funnel(res_overall_fixed, xlab = "Residual Logit(AUROC/C-stats.)", main = "(F) Common-Effect")
mtext(txt_F, side = 3, line = -23, adj = 0, cex = 0.75)   # annotate panel F

## ------------------------------------------------------------
par(oldpar)    # restore previous graphic parameters




### ============================================================================
### PART III: Meta Analysis of Non-Survival Outcomes 
### ============================================================================

# Load data with non-survival outcomes and preprocess data 
data <- read_excel(
  here("data", "extracted_data.xlsx"),
  sheet = "meta_nonsurv"
)

data$author_year <- sub("^(\\S+)\\s+(\\d{4})$", "\\1 et al. (\\2)", data$author)
df <- data

df <- df %>%
  mutate(
    se = if_else((!is.na(lb)) & (!is.na(ub)) & is.na(se),
                 (ub - lb)/(2* qnorm(0.975)), 
                 se)
  )

df <- df %>%
  mutate(
    lb = if_else(is.na(lb) & !is.na(se),
                 mean_metric - qnorm(0.975) * se,
                 lb),
    ub = if_else(is.na(ub) & !is.na(se),
                 mean_metric + qnorm(0.975) * se,
                 ub),
    lb = pmax(lb, 0),
    ub = pmin(ub, 1)
  )

df <- df %>%
  mutate(
    values = if_else(
      !is.na(lb) & !is.na(ub),
      # If lb and ub exist:
      paste0(
        sprintf("%.3f", mean_metric), " (",
        sprintf("%.3f", lb), ", ",
        sprintf("%.3f", ub), ")"
      ),
      # Otherwise just the mean:
      sprintf("%.3f", mean_metric)
    )
  )

df_ind <- df %>%
  arrange(outcome, author_year) %>%
  mutate(
    lb = if_else(is.na(lb), mean_metric , lb),
    ub = if_else(is.na(ub), mean_metric , ub)
  )

df_ind_filtered<- df_ind %>%
  filter(metric != "AUPRC") %>%
  rename(
    'Study' = author_year,
    'Outcome' = outcome,
    'AUROC (95% CI)' = values
  )


# Generate forest plot for individual studies with non-survival outcomes
df_ind_filtered$` ` <- paste(rep(" ", 35), collapse = " ")
forestploter::forest(df_ind_filtered[,c(3,9, 11, 10)],
                     est = df_ind_filtered$mean_metric,
                     lower = df_ind_filtered$lb, 
                     upper = df_ind_filtered$ub,
                     ci_column = 3,
                     ref_line = 0.5,
                     xlim = c(0.35, 1),
                     ticks_at = c(0.4,  0.5, 0.6,  0.7, 0.8, 0.9,  1)
)



### ============================================================================
### PART IV: Risk of Bias and Applicability Analysis
### ============================================================================

# Step 1: Import the Excel file
df <- read_excel(
  here("data", "extracted_data.xlsx"),
  sheet = "risk_assess"
)

# Step 2: Reshape the data into long format for plotting
df_long <- df %>%
  pivot_longer(
    cols = -Study,
    names_to = "Domain",
    values_to = "Rating"
  )

# Step 3: Calculate percentage of Low / High / Unclear per Domain
df_summary <- df_long %>%
  group_by(Domain, Rating) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(Domain) %>%
  mutate(Percentage = 100 * N / sum(N))

# Step 4: Define colors
rating_colors <- c("Low" = "#0072B2", "High" = "#D55E00", "Unclear" = "#999999")

# Optional: Reorder the domains to group Risk of Bias and Applicability
df_summary <- df_summary %>%
  mutate(Domain = fct_relevel(
    Domain,
    "RoB.: Participants Selection",
    "RoB.: Predictor",
    "RoB.: Outcome",
    "RoB.: Analysis",
    "Overall RoB.",
    "App.: Participants Selection",
    "App.: Predictor",
    "App.: Outcome",
    "Overall App."
  ))

df_summary <- df_summary %>%
  mutate(Domain = fct_relevel(
    Domain,
    "RoB.: Participants Selection",
    "RoB.: Predictor",
    "RoB.: Outcome",
    "RoB.: Analysis",
    "Overall RoB.",
    "App.: Participants Selection",
    "App.: Predictor",
    "App.: Outcome",
    "Overall App."
  ))


bold_domains <- c("Overall RoB.", "Overall App.")

# Step 2: Custom labeller for y-axis with bold labels for "overall" rows
custom_labels <- function(labels) {
  sapply(labels, function(label) {
    if (label %in% bold_domains) {
      bquote(bold(.(label)))
    } else {
      label
    }
  })
}

# Step 3: Plot without bar borders, with bolded overall labels
ggplot(df_summary, aes(x = Percentage, y = fct_rev(Domain), fill = Rating)) +
  geom_col(width = 0.7) +  # removed color = "black" to eliminate bar edges
  scale_fill_manual(values = rating_colors) +
  geom_text(
    aes(label = paste0(round(Percentage), "%")),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "white"
  ) +
  theme_minimal(base_size = 12) +
  labs(
    #title = "Risk of Bias and Applicability in ML Studies",
    x = "Percentage (%)",
    y = NULL,
    fill = NULL
  ) +
  scale_y_discrete(labels = custom_labels) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    axis.text.y   = element_text(size = 12),  # y-axis label size
    legend.text   = element_text(size = 12)   # legend label size
  ) # 1100 x 550



