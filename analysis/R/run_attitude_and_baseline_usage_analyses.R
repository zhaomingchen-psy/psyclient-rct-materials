options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(foreign)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
  project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
} else {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

analysis_dir <- file.path(project_root, "output", "analysis")
output_dir <- file.path(project_root, "output", "r_attitudes")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

normalize_id <- function(x) {
  x <- as.character(x)
  x <- sub("\\.0$", "", x)
  x[x %in% c("NA", "NaN", "", "NULL")] <- NA_character_
  x
}

group_factor_zh <- function(x) {
  factor(x, levels = c("空白组", "无反馈组", "实验组"))
}

group_factor_en <- function(x) {
  factor(
    x,
    levels = c("空白组", "无反馈组", "实验组"),
    labels = c("blank", "no_feedback", "feedback")
  )
}

group_factor_sav <- function(x) {
  factor(
    x,
    levels = c(3, 2, 1),
    labels = c("空白组", "无反馈组", "实验组")
  )
}

zscore <- function(x) {
  as.numeric(scale(x))
}

mean_sd <- function(x) {
  c(n = sum(!is.na(x)), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
}

extract_anova_term <- function(fit, term_name, model_name, subset_name, outcome_name, moderator_name) {
  tab <- as.data.frame(anova(fit))
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  names(tab) <- make.names(names(tab))
  row <- subset(tab, term == term_name)
  if (nrow(row) == 0) {
    return(NULL)
  }
  row$model <- model_name
  row$subset <- subset_name
  row$outcome <- outcome_name
  row$moderator <- moderator_name
  row
}

extract_lm_terms <- function(fit, model_name, subset_name, outcome_name, moderator_name) {
  tab <- as.data.frame(summary(fit)$coefficients)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  names(tab)[1:4] <- c("estimate", "std_error", "t_value", "p_value")
  tab$model <- model_name
  tab$subset <- subset_name
  tab$outcome <- outcome_name
  tab$moderator <- moderator_name
  tab$r_squared <- summary(fit)$r.squared
  tab$adj_r_squared <- summary(fit)$adj.r.squared
  tab$n <- nobs(fit)
  tab
}

predict_group_means <- function(fit, group_levels, moderator_mean, moderator_sd, moderator_name, subset_name, outcome_name) {
  at_values <- c(low = moderator_mean - moderator_sd, mean = moderator_mean, high = moderator_mean + moderator_sd)
  rows <- list()
  for (level_name in names(at_values)) {
    newdat <- data.frame(
      group = factor(group_levels, levels = group_levels),
      moderator_value = at_values[[level_name]]
    )
    names(newdat)[2] <- moderator_name
    preds <- predict(fit, newdata = newdat)
    rows[[level_name]] <- data.frame(
      subset = subset_name,
      outcome = outcome_name,
      moderator = moderator_name,
      moderator_level = level_name,
      moderator_value = at_values[[level_name]],
      group = group_levels,
      predicted = preds
    )
  }
  do.call(rbind, rows)
}

pm <- read.csv(file.path(analysis_dir, "participant_master.csv"), check.names = FALSE)
wide <- read.csv(file.path(analysis_dir, "rct_scale_totals_wide.csv"), check.names = FALSE)
sav <- read.spss(
  file.path(project_root, "Psyclient_data.sav"),
  to.data.frame = TRUE,
  use.value.labels = FALSE
)
names(sav) <- trimws(names(sav))

pm$full_id <- normalize_id(pm$full_id)
pm$group_label <- as.character(pm$group_label)

wide$full_id <- normalize_id(wide$full_id)
wide$group_label <- as.character(wide$group_label)

sav$full_id <- normalize_id(sav$ID)
sav$group_label <- as.character(group_factor_sav(sav$Group))

# 1. Baseline confidence predicting feedback usage / copying
feedback_dat <- merge(
  subset(wide, group_label == "实验组"),
  subset(
    pm,
    select = c(
      full_id,
      practice_session_count,
      practice_feedback_count,
      practice_overall_mean_score_weighted,
      practice_total_score
    )
  ),
  by = "full_id",
  all.x = TRUE
)

feedback_dat <- feedback_dat[!is.na(feedback_dat$practice_session_count) & feedback_dat$practice_session_count > 0, ]
feedback_dat$feedback_rate <- feedback_dat$practice_feedback_count / feedback_dat$practice_session_count

baseline_predictors <- c("hs_mean_t1", "sm_mean_t1", "cc_mean_t1", "self_efficacy_mean_t1")
usage_outcomes <- c(
  "practice_feedback_count",
  "feedback_rate",
  "practice_session_count",
  "practice_overall_mean_score_weighted"
)

baseline_model_rows <- list()
baseline_desc_rows <- list()

for (pred in baseline_predictors) {
  for (outcome in usage_outcomes) {
    dat <- feedback_dat[!is.na(feedback_dat[[pred]]) & !is.na(feedback_dat[[outcome]]), c(pred, outcome)]
    if (nrow(dat) < 8) {
      next
    }
    fit <- lm(as.formula(paste0(outcome, " ~ ", pred)), data = dat)
    tab <- extract_lm_terms(
      fit,
      model_name = "baseline_usage_prediction",
      subset_name = "feedback_only",
      outcome_name = outcome,
      moderator_name = pred
    )
    tab <- subset(tab, term == pred)
    baseline_model_rows[[paste(pred, outcome, sep = "_")]] <- tab
  }
}

for (outcome in usage_outcomes) {
  desc <- mean_sd(feedback_dat[[outcome]])
  baseline_desc_rows[[outcome]] <- data.frame(
    outcome = outcome,
    n = desc["n"],
    mean = desc["mean"],
    sd = desc["sd"]
  )
}

if (length(baseline_model_rows) > 0) {
  write.csv(do.call(rbind, baseline_model_rows), file.path(output_dir, "baseline_confidence_predicts_feedback_usage.csv"), row.names = FALSE)
}
write.csv(do.call(rbind, baseline_desc_rows), file.path(output_dir, "feedback_usage_descriptives.csv"), row.names = FALSE)

# 2. PA / NA moderation analyses
attitude_df <- merge(
  wide[, c(
    "full_id", "group_label", "hs_mean_t6", "sm_mean_t6", "cc_mean_t6",
    "self_efficacy_mean_t6", "pa_mean_t6", "na_mean_t6"
  )],
  pm[, c("full_id", "content_exploration_norm_post")],
  by = "full_id",
  all.x = TRUE
)

attitude_df <- merge(
  attitude_df,
  sav[, c("full_id", "RQR_new_post", "questioning_post")],
  by = "full_id",
  all.x = TRUE
)

attitude_df$group_zh <- group_factor_zh(attitude_df$group_label)
attitude_df$group_en <- group_factor_en(attitude_df$group_label)
attitude_df$question_ratio_post <- attitude_df$questioning_post
attitude_df$confidence_post_mean <- rowMeans(
  attitude_df[, c("hs_mean_t6", "sm_mean_t6", "cc_mean_t6")],
  na.rm = TRUE
)
attitude_df$mi_style_post <- (zscore(log1p(attitude_df$RQR_new_post)) + zscore(attitude_df$question_ratio_post)) / 2
attitude_df$confidence_post_z <- rowMeans(
  cbind(
    zscore(attitude_df$hs_mean_t6),
    zscore(attitude_df$sm_mean_t6),
    zscore(attitude_df$cc_mean_t6)
  ),
  na.rm = TRUE
)
attitude_df$style_conf_gap_post <- attitude_df$mi_style_post - attitude_df$confidence_post_z

moderators <- c("pa_mean_t6", "na_mean_t6")
outcomes <- c(
  "hs_mean_t6",
  "sm_mean_t6",
  "cc_mean_t6",
  "content_exploration_norm_post",
  "RQR_new_post",
  "question_ratio_post",
  "mi_style_post",
  "style_conf_gap_post"
)

subsets <- list(
  list(name = "all_groups", filter = rep(TRUE, nrow(attitude_df)), group_var = "group_zh"),
  list(name = "ai_only", filter = attitude_df$group_label %in% c("无反馈组", "实验组"), group_var = "group_zh")
)

interaction_rows <- list()
coefficient_rows <- list()
pred_rows <- list()
band_rows <- list()

for (subset_spec in subsets) {
  subset_name <- subset_spec$name
  group_var <- subset_spec$group_var
  for (moderator in moderators) {
    for (outcome in outcomes) {
      dat <- attitude_df[subset_spec$filter & !is.na(attitude_df[[group_var]]) & !is.na(attitude_df[[moderator]]) & !is.na(attitude_df[[outcome]]), ]
      if (nrow(dat) < 20) {
        next
      }
      dat$group <- dat[[group_var]]
      if (length(unique(dat$group)) < 2) {
        next
      }
      fit <- lm(as.formula(paste0(outcome, " ~ group * ", moderator)), data = dat)
      interaction_rows[[paste(subset_name, moderator, outcome, sep = "_")]] <- extract_anova_term(
        fit = fit,
        term_name = paste0("group:", moderator),
        model_name = "group_by_attitude_moderation",
        subset_name = subset_name,
        outcome_name = outcome,
        moderator_name = moderator
      )
      coefficient_rows[[paste(subset_name, moderator, outcome, sep = "_")]] <- extract_lm_terms(
        fit = fit,
        model_name = "group_by_attitude_moderation",
        subset_name = subset_name,
        outcome_name = outcome,
        moderator_name = moderator
      )
      pred_rows[[paste(subset_name, moderator, outcome, sep = "_")]] <- predict_group_means(
        fit = fit,
        group_levels = levels(droplevels(dat$group)),
        moderator_mean = mean(dat[[moderator]], na.rm = TRUE),
        moderator_sd = sd(dat[[moderator]], na.rm = TRUE),
        moderator_name = moderator,
        subset_name = subset_name,
        outcome_name = outcome
      )

      med <- median(dat[[moderator]], na.rm = TRUE)
      dat$moderator_band <- ifelse(dat[[moderator]] <= med, "low", "high")
      desc <- aggregate(dat[[outcome]], by = list(group = dat$group, moderator_band = dat$moderator_band), FUN = mean_sd)
      band_rows[[paste(subset_name, moderator, outcome, sep = "_")]] <- data.frame(
        subset = subset_name,
        moderator = moderator,
        outcome = outcome,
        group = desc$group,
        moderator_band = desc$moderator_band,
        n = desc$x[, "n"],
        mean = desc$x[, "mean"],
        sd = desc$x[, "sd"],
        median_cut = med
      )
    }
  }
}

if (length(interaction_rows) > 0) {
  write.csv(do.call(rbind, interaction_rows), file.path(output_dir, "attitude_moderation_interactions.csv"), row.names = FALSE)
}
if (length(coefficient_rows) > 0) {
  write.csv(do.call(rbind, coefficient_rows), file.path(output_dir, "attitude_moderation_coefficients.csv"), row.names = FALSE)
}
if (length(pred_rows) > 0) {
  write.csv(do.call(rbind, pred_rows), file.path(output_dir, "attitude_moderation_predicted_means.csv"), row.names = FALSE)
}
if (length(band_rows) > 0) {
  write.csv(do.call(rbind, band_rows), file.path(output_dir, "attitude_moderation_band_descriptives.csv"), row.names = FALSE)
}

notes <- c(
  "# Attitude and Baseline Usage Analyses",
  "",
  "- `PA/NA` are only available at T6 in the cleaned datasets.",
  "- Therefore, PA/NA moderation results are concurrent/post-hoc interactions, not true pre-treatment moderation.",
  "- Baseline confidence -> feedback usage models are run only within the feedback group, because `feedback_count` is only meaningful there.",
  "- `question_ratio_post` here corresponds to `% open questions`.",
  "- `mi_style_post` is the mean of standardized `log1p(RQR_new_post)` and standardized `% open questions`."
)
writeLines(notes, file.path(output_dir, "attitude_and_usage_notes.md"))
