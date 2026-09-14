# ============================================================
# DeTox dataset: discrimination types & engagement analysis
# Dataset: Demus et al. (2025)
# The code was sorted and optimized by AI
# ============================================================

# install.packages(c("DBI", "RSQLite", "dplyr", "MASS", "lmtest", "car",
#     "sandwich", "ggplot2"))

library(DBI)
library(RSQLite)
library(dplyr) 
library(MASS)  
library(lmtest)
library(car)
library(sandwich)
library(ggplot2)  

# ------------------------------------------------------------
# 0. Helpers
# ------------------------------------------------------------

# Poisson dispersion statistic (close to 1 => Poisson adequate)
dispersion_stat <- function(model) model$deviance / model$df.residual

# Descriptive engagement statistics by group
summarise_engagement <- function(df, group_col) {
  df %>%
    group_by({{ group_col }}) %>%
    summarise(
      n_comments      = n(),
      mean_likes      = mean(like_count, na.rm = TRUE),
      median_likes    = median(like_count, na.rm = TRUE),
      mean_retweets   = mean(retweet_count, na.rm = TRUE),
      median_retweets = median(retweet_count, na.rm = TRUE),
      mean_replies    = mean(reply_count, na.rm = TRUE),
      median_replies  = median(reply_count, na.rm = TRUE),
      mean_quotes     = mean(quote_count, na.rm = TRUE),
      median_quotes   = median(quote_count, na.rm = TRUE),
      .groups = "drop"
    )
}

# Build a formula from a character vector of predictors
make_formula <- function(outcome, predictors) reformulate(predictors, outcome)

# Fit the same NB model for all four engagement outcomes.
# extra: optional additional RHS term(s), e.g. "log_age" or "offset(log_age)"
fit_nb_all <- function(predictors, data, extra = NULL) {
  outcomes <- c("like_count", "retweet_count", "reply_count", "quote_count")
  models <- lapply(outcomes, function(oc) {
    glm.nb(reformulate(c(predictors, extra), oc),
           data = data, control = glm.control(maxit = 100))
  })
  setNames(models, c("likes", "retweets", "replies", "quotes"))
}

# Same for Poisson baselines
fit_poisson_all <- function(predictors, data, extra = NULL) {
  outcomes <- c("like_count", "retweet_count", "reply_count", "quote_count")
  models <- lapply(outcomes, function(oc) {
    glm(reformulate(c(predictors, extra), oc),
        family = poisson, data = data)
  })
  setNames(models, c("likes", "retweets", "replies", "quotes"))
}

# IRRs and 95% CIs for a list of models
irr_table <- function(models) {
  list(
    IRR = lapply(models, function(m) exp(coef(m))),
    CI  = lapply(models, function(m) exp(confint(m)))
  )
}

# LR tests between two matched lists of models
lrtest_all <- function(restricted, full) {
  for (nm in names(full)) {
    cat("\n---", nm, "---\n")
    print(lrtest(restricted[[nm]], full[[nm]]))
  }
  invisible(NULL)
}

# Complete-case filter for the four engagement metrics
filter_complete_engagement <- function(df) {
  df %>%
    filter(!is.na(like_count), !is.na(retweet_count),
           !is.na(reply_count), !is.na(quote_count))
}

# Parse SQLite TIMESTAMP robustly: ISO text or epoch seconds
parse_ts <- function(x) {
  if (is.numeric(x)) as.POSIXct(x, origin = "1970-01-01", tz = "UTC")
  else               as.POSIXct(x, tz = "UTC")
}

# Add posting time, exposure (comment age) and conversation cluster id.
# cluster_id is only used in the appendix robustness check.
add_exposure <- function(df, collection_time) {
  df %>%
    mutate(
      posted_at  = parse_ts(date),
      age_days   = as.numeric(difftime(collection_time, posted_at,
                                       units = "days")),
      age_days   = pmax(age_days, 1/24),  # floor at 1 hour -> no log(0)
      log_age    = log(age_days),
      # comments without a conversation become singleton clusters
      cluster_id = ifelse(is.na(conv_id),
                          paste0("single_", c_id),
                          as.character(conv_id))
    )
}

# Tidy results table from an NB model (profile-likelihood CIs)
tidy_nb <- function(model, hypothesis, outcome) {
  co <- summary(model)$coefficients
  ci <- suppressMessages(confint(model))
  data.frame(
    hypothesis  = hypothesis,
    outcome     = outcome,
    term        = rownames(co),
    estimate    = co[, "Estimate"],
    std_error   = co[, "Std. Error"],
    z_value     = co[, "z value"],
    p_value     = co[, "Pr(>|z|)"],
    irr         = exp(co[, "Estimate"]),
    irr_ci_low  = exp(ci[, 1]),
    irr_ci_high = exp(ci[, 2]),
    row.names   = NULL,
    stringsAsFactors = FALSE
  )
}

# Cluster-robust variant (appendix only): sandwich::vcovCL + coeftest.
# CIs use the normal approximation (estimate +/- 1.96 * robust SE).
tidy_cluster <- function(model, cluster, hypothesis, outcome) {
  ct <- lmtest::coeftest(
    model, vcov. = sandwich::vcovCL(model, cluster = cluster)
  )
  data.frame(
    hypothesis  = hypothesis,
    outcome     = outcome,
    term        = rownames(ct),
    estimate    = ct[, "Estimate"],
    std_error   = ct[, "Std. Error"],
    z_value     = ct[, "z value"],
    p_value     = ct[, "Pr(>|z|)"],
    irr         = exp(ct[, "Estimate"]),
    irr_ci_low  = exp(ct[, "Estimate"] - 1.96 * ct[, "Std. Error"]),
    irr_ci_high = exp(ct[, "Estimate"] + 1.96 * ct[, "Std. Error"]),
    row.names   = NULL,
    stringsAsFactors = FALSE
  )
}

# Combine a named model list into one results table.
# Holm adjustment is applied within the hypothesis family across the
# focal terms only (i.e., excluding intercept and the log_age control).
build_results <- function(models, hypothesis, cluster = NULL) {
  tab <- do.call(rbind, Map(
    function(m, oc) {
      if (is.null(cluster)) tidy_nb(m, hypothesis, oc)
      else tidy_cluster(m, cluster, hypothesis, oc)
    },
    models, names(models)
  ))
  tab$is_control <- tab$term %in% c("(Intercept)", "log_age")
  focal <- !tab$is_control
  tab$p_holm <- NA_real_
  tab$p_holm[focal] <- p.adjust(tab$p_value[focal], method = "holm")
  tab
}

# Report which focal terms change significance after clustering
compare_sig <- function(conv, cl, label) {
  m <- merge(
    conv[!conv$is_control, c("outcome", "term", "p_value")],
    cl[!cl$is_control,   c("outcome", "term", "p_value")],
    by = c("outcome", "term"), suffixes = c("_conv", "_cluster")
  )
  m$changed <- (m$p_value_conv < 0.05) != (m$p_value_cluster < 0.05)
  cat("\n", label, ": focal terms changing significance after clustering\n")
  print(m[m$changed, , drop = FALSE])
  invisible(m)
}

# ------------------------------------------------------------
# 1. Load database
# ------------------------------------------------------------

database_path <- "DeTox-Dataset_complete.sqlite3"

if (!file.exists(database_path)) {
  stop(
    "Database file not found: ", database_path,
    "\nPlace it in the working directory or set the full path."
  )
}

dbconnect <- tryCatch(
  dbConnect(RSQLite::SQLite(), database_path),
  error = function(e) stop("Database connection failed: ", conditionMessage(e))
)

table_names <- dbListTables(dbconnect)
print(table_names)

db_tables <- setNames(
  lapply(table_names, function(tbl) {
    dbGetQuery(dbconnect, sprintf("SELECT * FROM %s;", tbl))
  }),
  table_names
)
dbDisconnect(dbconnect)

Annotations  <- db_tables$Annotations
Goldstandard <- db_tables$Goldstandard
Comments     <- db_tables$Comments

lapply(db_tables, head)

# Validate tables and columns the analysis depends on
discrim_all <- c(
  "discrim_job", "discrim_attitude", "discrim_engagement",
  "discrim_sexIdent", "discrim_characteristics", "discrim_nation",
  "discrim_religion", "discrim_socialStatus", "discrim_worldview",
  "discrim_Ethnicity"
)
engagement_cols <- c("like_count", "quote_count", "retweet_count", "reply_count")

stopifnot(all(discrim_all %in% colnames(Annotations)))
stopifnot(all(c("c_id", "hate_speech", discrim_all) %in% colnames(Goldstandard)))
stopifnot(all(c("c_id", "date", "conv_id", engagement_cols) %in% colnames(Comments)))
stopifnot(sum(duplicated(Comments$c_id)) == 0)

# ------------------------------------------------------------
# 1b. Exposure: reference time and comment age
# ------------------------------------------------------------

collection_time <- max(parse_ts(Comments$date), na.rm = TRUE)
print(collection_time)

Comments <- Comments %>%
  mutate(posted_at = parse_ts(date))

sum(is.na(Comments$posted_at))  # parsing failures: must be 0
summary(as.numeric(difftime(collection_time, Comments$posted_at,
                            units = "days")))

# ------------------------------------------------------------
# 2. Summary of discrimination variables
# ------------------------------------------------------------

discrim_counts_annotations <- sapply(
  Annotations[discrim_all], function(x) sum(x > 0.5, na.rm = TRUE)
)
discrim_counts_goldstandard <- sapply(
  Goldstandard[discrim_all], function(x) sum(x > 0.5, na.rm = TRUE)
)
print(discrim_counts_annotations)
print(discrim_counts_goldstandard)

# ------------------------------------------------------------
# 3. Regression dataset
# ------------------------------------------------------------

discrim_included <- setdiff(
  discrim_all,
  c("discrim_sexIdent", "discrim_socialStatus", "discrim_worldview")
)

mod_Goldstandard <- Goldstandard %>%
  filter(hate_speech > 0.5) %>%
  dplyr::select(c_id, hate_speech, all_of(discrim_included)) %>%
  mutate(across(all_of(discrim_included), ~ as.integer(. > 0.5))) %>%
  left_join(
    Comments %>% dplyr::select(c_id, date, conv_id, all_of(engagement_cols)),
    by = "c_id"
  ) %>%
  add_exposure(collection_time)

sum(is.na(mod_Goldstandard$like_count))
sum(is.na(mod_Goldstandard$posted_at))  # must be 0
mod_Goldstandard <- filter_complete_engagement(mod_Goldstandard)
nrow(mod_Goldstandard)

# ------------------------------------------------------------
# 4. Poisson vs. negative binomial, with exposure
# ------------------------------------------------------------

poisson_models <- fit_poisson_all(discrim_included, mod_Goldstandard,
                                  extra = "log_age")
sapply(poisson_models, dispersion_stat)

nb_models <- fit_nb_all(discrim_included, mod_Goldstandard,
                        extra = "log_age")
lapply(nb_models, summary)

# Offset variant: coefficient on log(exposure) fixed to 1 (rate model)
nb_models_offset <- fit_nb_all(discrim_included, mod_Goldstandard,
                               extra = "offset(log_age)")

# LR test per outcome, H0: coef(log_age) = 1 (rejected -> keep free term)
lrtest_all(nb_models_offset, nb_models)
sapply(nb_models, function(m) coef(m)["log_age"])

# NB vs. Poisson (both with free log_age)
lrtest_all(poisson_models, nb_models)
AIC(poisson_models$likes,    nb_models$likes)
AIC(poisson_models$retweets, nb_models$retweets)
AIC(poisson_models$replies,  nb_models$replies)
AIC(poisson_models$quotes,   nb_models$quotes)

irr_table(nb_models)

# ------------------------------------------------------------
# 5. Baseline and excluded-type checks
# ------------------------------------------------------------

baseline_comments <- mod_Goldstandard %>%
  filter(if_all(all_of(discrim_included), ~ . == 0))
nrow(baseline_comments)
nrow(mod_Goldstandard)

excluded_check <- Goldstandard %>%
  filter(c_id %in% baseline_comments$c_id) %>%
  dplyr::select(c_id, discrim_sexIdent, discrim_socialStatus, discrim_worldview) %>%
  mutate(any_excluded_type = discrim_sexIdent > 0.5 |
           discrim_socialStatus > 0.5 |
           discrim_worldview > 0.5)
sum(excluded_check$any_excluded_type, na.rm = TRUE)

mod_Goldstandard_clean <- mod_Goldstandard %>%
  filter(!(c_id %in% excluded_check$c_id[excluded_check$any_excluded_type]))
nrow(mod_Goldstandard_clean)

nb_likes_clean <- glm.nb(
  make_formula("like_count", c(discrim_included, "log_age")),
  data = mod_Goldstandard_clean, control = glm.control(maxit = 100)
)
summary(nb_likes_clean)
nb_likes_clean$converged

vif(nb_likes_clean)
table(mod_Goldstandard_clean$discrim_nation,
      mod_Goldstandard_clean$discrim_Ethnicity)
cor(mod_Goldstandard_clean[, c("discrim_nation",
                               "discrim_religion",
                               "discrim_Ethnicity")])

# ------------------------------------------------------------
# 6. Single vs. multiple discrimination types (H2 / H4)
# ------------------------------------------------------------

mod_Goldstandard <- mod_Goldstandard %>%
  mutate(
    n_discrim_types = rowSums(across(all_of(discrim_included)), na.rm = TRUE),
    discrim_category = case_when(
      n_discrim_types == 0 ~ "none",
      n_discrim_types == 1 ~ "single",
      n_discrim_types >= 2 ~ "multiple"
    ),
    discrim_category = factor(discrim_category,
                              levels = c("none", "single", "multiple"))
  )

table(mod_Goldstandard$n_discrim_types)
table(mod_Goldstandard$discrim_category)

engagement_summary <- summarise_engagement(mod_Goldstandard, discrim_category)
print(engagement_summary)

# Age balance check
mod_Goldstandard %>%
  group_by(discrim_category) %>%
  summarise(median_age_days = median(age_days),
            mean_age_days   = mean(age_days), .groups = "drop")

# Reference = "none" (with exposure adjustment)
category_models <- fit_nb_all("discrim_category", mod_Goldstandard,
                              extra = "log_age")
lapply(category_models, summary)
irr_table(category_models)

# Reference = "single" (direct multiple-vs-single test)
mod_Goldstandard <- mod_Goldstandard %>%
  mutate(discrim_category_single_ref = relevel(discrim_category, ref = "single"))

single_ref_models <- fit_nb_all("discrim_category_single_ref",
                                mod_Goldstandard, extra = "log_age")
lapply(single_ref_models, summary)
irr_table(single_ref_models)

# ------------------------------------------------------------
# 7. Hate speech vs. non-hate speech (H1)
# ------------------------------------------------------------

hate_comparison <- Goldstandard %>%
  mutate(hate_speech_binary = as.integer(hate_speech > 0.5)) %>%
  dplyr::select(c_id, hate_speech_binary) %>%
  left_join(
    Comments %>% dplyr::select(c_id, date, conv_id, all_of(engagement_cols)),
    by = "c_id"
  ) %>%
  add_exposure(collection_time) %>%
  filter_complete_engagement()

table(hate_comparison$hate_speech_binary)

engagement_by_hate <- summarise_engagement(hate_comparison, hate_speech_binary)
print(engagement_by_hate)

t.test(age_days ~ hate_speech_binary, data = hate_comparison)

poisson_hate <- fit_poisson_all("hate_speech_binary", hate_comparison,
                                extra = "log_age")
sapply(poisson_hate, dispersion_stat)

hate_models <- fit_nb_all("hate_speech_binary", hate_comparison,
                          extra = "log_age")
lapply(hate_models, summary)

hate_models_offset <- fit_nb_all("hate_speech_binary", hate_comparison,
                                 extra = "offset(log_age)")
lrtest_all(hate_models_offset, hate_models)

lrtest_all(poisson_hate, hate_models)
irr_table(hate_models)

# ------------------------------------------------------------
# 8. Results tables with Holm adjustment
# ------------------------------------------------------------

results_h1 <- build_results(hate_models,       "H1")
results_h2 <- build_results(category_models,   "H2")
results_h3 <- build_results(nb_models,         "H3")
results_h4 <- build_results(single_ref_models, "H4")

results_all <- rbind(results_h1, results_h2, results_h3, results_h4)

# Focal results (no intercept / control), Holm-adjusted within hypothesis
results_focal <- results_all %>% filter(!is_control)

# ------------------------------------------------------------
# 9. ALL HYPOTHESIS OUTPUTS IN ONE PLACE 
# ------------------------------------------------------------

# ---- H1: hate speech vs. non-hate speech --------------------
cat("\n########## H1: hate speech vs. non-hate speech ##########\n")
lapply(hate_models, summary)                 # coefficients, z, p, theta, AIC
irr_table(hate_models)                       # IRRs + 95% profile CIs
print(engagement_by_hate)                    # descriptives
table(hate_comparison$hate_speech_binary)    # group sizes

# ---- H2: discrimination breadth, reference = "none" ---------
cat("\n########## H2: tagged vs. non-tagged (ref = none) ##########\n")
lapply(category_models, summary)
irr_table(category_models)
print(engagement_summary)
table(mod_Goldstandard$discrim_category)

# ---- H3: discrimination types across the four metrics -------
cat("\n########## H3: profiles across metrics ##########\n")
lapply(nb_models, summary)
irr_table(nb_models)

# IRR matrix: rows = discrimination types, columns = metrics
round(cbind(
  likes    = exp(coef(nb_models$likes))[-1],
  retweets = exp(coef(nb_models$retweets))[-1],
  replies  = exp(coef(nb_models$replies))[-1],
  quotes   = exp(coef(nb_models$quotes))[-1]
), 3)

vif(nb_models$likes)                         # collinearity check

# ---- H4: multiple vs. single, reference = "single" ----------
cat("\n########## H4: multiple vs. single (ref = single) ##########\n")
lapply(single_ref_models, summary)
irr_table(single_ref_models)

# ---- Holm-adjusted focal p-values (all hypotheses) ----------
cat("\n########## Holm-adjusted focal results ##########\n")
print(results_focal)

# ------------------------------------------------------------
# 10. Figures (saved as PNGs)
# ------------------------------------------------------------

outcome_levels <- c("likes", "retweets", "replies", "quotes")

# Fig. 1: H1 - hate-speech IRRs across the four metrics
fig_h1 <- results_h1 %>%
  filter(term == "hate_speech_binary") %>%
  mutate(outcome = factor(outcome, levels = outcome_levels)) %>%
  ggplot(aes(outcome, irr)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = irr_ci_low, ymax = irr_ci_high), width = .15) +
  scale_y_log10() +
  labs(x = NULL, y = "IRR (log scale)",
       title = "H1: Hate speech vs. non-hate speech",
       subtitle = "Incidence rate ratios with 95% profile-likelihood CIs; adjusted for comment age")
print(fig_h1)
ggsave("fig_h1_hate_irr.png", fig_h1, width = 6, height = 4, dpi = 300)

# Fig. 2: H3 - engagement profiles of the seven discrimination types
fig_h3 <- results_h3 %>%
  filter(!is_control) %>%
  mutate(
    outcome = factor(outcome, levels = outcome_levels),
    term    = factor(term, levels = discrim_included),
    sig     = ifelse(!is.na(p_holm) & p_holm < .05, "*", "")
  ) %>%
  ggplot(aes(outcome, term, fill = irr)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = paste0(sprintf("%.2f", irr), sig)), size = 3.4) +
  scale_fill_gradient2(low = "#b2182b", mid = "grey92", high = "#2166ac",
                       midpoint = 0, trans = "log10", name = "IRR") +
  labs(x = NULL, y = NULL,
       title = "H3: Discrimination-type engagement profiles (IRRs)",
       caption = "* Holm-adjusted p < .05 within the H3 family") +
  theme_minimal(base_size = 11)
print(fig_h3)
ggsave("fig_h3_profiles_heatmap.png", fig_h3, width = 7, height = 5, dpi = 300)

# Fig. 3: H2 + H4 - discrimination breadth contrasts across metrics
breadth_df <- bind_rows(
  results_h2 %>% filter(term == "discrim_categorysingle") %>%
    mutate(contrast = "Single vs. none"),
  results_h2 %>% filter(term == "discrim_categorymultiple") %>%
    mutate(contrast = "Multiple vs. none"),
  results_h4 %>% filter(term == "discrim_category_single_refmultiple") %>%
    mutate(contrast = "Multiple vs. single")
) %>%
  mutate(outcome = factor(outcome, levels = outcome_levels))

fig_breadth <- ggplot(breadth_df, aes(outcome, irr, colour = contrast)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_point(size = 2.6, position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = irr_ci_low, ymax = irr_ci_high), width = .15,
                position = position_dodge(width = .5)) +
  scale_y_log10() +
  labs(x = NULL, y = "IRR (log scale)", colour = "Contrast",
       title = "H2/H4: Discrimination breadth and engagement",
       subtitle = "Hate-speech subsample; adjusted for comment age")
print(fig_breadth)
ggsave("fig_h2h4_breadth_irr.png", fig_breadth, width = 7, height = 4.5, dpi = 300)

# ------------------------------------------------------------
# 11. Sensitivity: comments old enough that engagement plateaued
# ------------------------------------------------------------

plateau_cutoff_days <- 30
mod_plateau <- mod_Goldstandard %>% filter(age_days >= plateau_cutoff_days)
nrow(mod_plateau)

nb_models_plateau <- fit_nb_all(discrim_included, mod_plateau,
                                extra = "log_age")
lapply(nb_models_plateau, summary)

# ------------------------------------------------------------
# 12. Appendix: exploratory conversation-clustered robustness
# ------------------------------------------------------------
# Reported in the paper's Limitations section only.
# Estimates are unchanged by construction; only SEs differ.

length(unique(mod_Goldstandard$cluster_id))   # conversations, hate subsample
length(unique(hate_comparison$cluster_id))    # conversations, full sample
table(table(mod_Goldstandard$cluster_id))     # cluster-size distribution

results_h1_cl <- build_results(hate_models,       "H1", hate_comparison$cluster_id)
results_h2_cl <- build_results(category_models,   "H2", mod_Goldstandard$cluster_id)
results_h3_cl <- build_results(nb_models,         "H3", mod_Goldstandard$cluster_id)
results_h4_cl <- build_results(single_ref_models, "H4", mod_Goldstandard$cluster_id)

results_all_cluster <- rbind(results_h1_cl, results_h2_cl,
                             results_h3_cl, results_h4_cl)

print(results_all_cluster %>% filter(!is_control))

# Which focal conclusions change under clustered SEs?
compare_sig(results_h1, results_h1_cl, "H1")
compare_sig(results_h2, results_h2_cl, "H2")
compare_sig(results_h3, results_h3_cl, "H3")
compare_sig(results_h4, results_h4_cl, "H4")

# Reproducibility record
sessionInfo()


