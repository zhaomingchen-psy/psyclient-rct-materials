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
output_dir <- file.path(project_root, "output", "r_storylines")
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

extract_simple_anova <- function(fit, outcome_name) {
  tab <- as.data.frame(summary(fit)[[1]])
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  names(tab) <- make.names(names(tab))
  tab$outcome <- outcome_name
  tab
}

extract_tukey <- function(fit, outcome_name) {
  tab <- as.data.frame(TukeyHSD(fit)$group)
  tab$comparison <- rownames(tab)
  rownames(tab) <- NULL
  names(tab) <- c("diff", "lwr", "upr", "p_adj", "comparison")
  tab$outcome <- outcome_name
  tab
}

run_post_group_models <- function(data, outcome_col, prefix) {
  dat <- data[!is.na(data[[outcome_col]]) & !is.na(data$group), ]
  if (nrow(dat) < 12) {
    return(NULL)
  }
  desc <- aggregate(dat[[outcome_col]], by = list(group = dat$group), FUN = mean_sd)
  desc_out <- data.frame(
    group = desc$group,
    n = desc$x[, "n"],
    mean = desc$x[, "mean"],
    sd = desc$x[, "sd"],
    outcome = outcome_col
  )
  fit <- aov(as.formula(paste0("`", outcome_col, "` ~ group")), data = dat)
  anova_out <- extract_simple_anova(fit, outcome_col)
  tukey_out <- extract_tukey(fit, outcome_col)
  write.csv(desc_out, file.path(output_dir, paste0(prefix, "_descriptives.csv")), row.names = FALSE)
  write.csv(anova_out, file.path(output_dir, paste0(prefix, "_anova.csv")), row.names = FALSE)
  write.csv(tukey_out, file.path(output_dir, paste0(prefix, "_tukey.csv")), row.names = FALSE)
  invisible(list(descriptives = desc_out, anova = anova_out, tukey = tukey_out))
}

extract_lm_terms <- function(fit, keep_terms, analysis_name, outcome_name) {
  tab <- as.data.frame(summary(fit)$coefficients)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  names(tab)[1:4] <- c("estimate", "std_error", "t_value", "p_value")
  tab <- subset(tab, term %in% keep_terms)
  tab$analysis <- analysis_name
  tab$outcome <- outcome_name
  tab
}

sav <- read.spss(
  file.path(project_root, "Psyclient_data.sav"),
  to.data.frame = TRUE,
  use.value.labels = FALSE
)
names(sav) <- trimws(names(sav))
sav$full_id <- normalize_id(sav$ID)
sav$group <- group_factor_sav(sav$Group)

master <- read.csv(file.path(analysis_dir, "participant_master.csv"), check.names = FALSE)
master$full_id <- normalize_id(master$full_id)
master$group <- group_factor_zh(master$group_label)

behavior <- read.csv(file.path(analysis_dir, "behavior_clean.csv"), check.names = FALSE)
behavior$full_id <- normalize_id(behavior$full_id)
behavior$group <- group_factor_zh(behavior$group_label)

practice <- read.csv(file.path(analysis_dir, "practice_student_summary_clean.csv"), check.names = FALSE)
practice$matched_full_id <- normalize_id(practice$matched_full_id)
practice$practice_student_id <- normalize_id(practice$practice_student_id)
practice <- subset(practice, !(practice_student_id %in% c("TA01", "TA02")))
practice$matched_full_id[practice$practice_student_id == "021"] <- "<STUDENT_ID_021>"
practice$matched_group_label[practice$practice_student_id == "021"] <- "实验组"
practice$roster_match_status[practice$practice_student_id == "021"] <- "manual_feedback_resolution"
practice$group <- group_factor_zh(practice$matched_group_label)

story_df <- merge(
  sav,
  behavior[, c("full_id", "content_exploration_norm_post", "content_exploration_norm_pre")],
  by = "full_id",
  all.x = TRUE
)

# MI style composite and confidence composite
story_df$mi_style_pre <- zscore((zscore(log1p(story_df$RQR_new_PRE)) + zscore(story_df$questioning_pre)) / 2)
story_df$mi_style_post <- zscore((zscore(log1p(story_df$RQR_new_post)) + zscore(story_df$questioning_post)) / 2)
story_df$confidence_pre_mean <- rowMeans(story_df[, c("HST1", "SMT1", "CCT1")], na.rm = TRUE)
story_df$confidence_post_mean <- rowMeans(story_df[, c("HST6", "SMT6", "CCT6")], na.rm = TRUE)
story_df$confidence_pre_z <- rowMeans(cbind(zscore(story_df$HST1), zscore(story_df$SMT1), zscore(story_df$CCT1)), na.rm = TRUE)
story_df$confidence_post_z <- rowMeans(cbind(zscore(story_df$HST6), zscore(story_df$SMT6), zscore(story_df$CCT6)), na.rm = TRUE)
story_df$style_conf_gap_pre <- story_df$mi_style_pre - story_df$confidence_pre_z
story_df$style_conf_gap_post <- story_df$mi_style_post - story_df$confidence_post_z
story_df$style_conf_gap_change <- story_df$style_conf_gap_post - story_df$style_conf_gap_pre

run_post_group_models(story_df, "mi_style_post", "mi_style_post")
run_post_group_models(story_df, "confidence_post_mean", "confidence_post_mean")
run_post_group_models(story_df, "style_conf_gap_post", "style_conf_gap_post")
run_post_group_models(story_df, "style_conf_gap_change", "style_conf_gap_change")

# RQR decomposition proxies / estimates
story_df$reflection_proxy_post <- story_df$ReflectionofFeelings_raw_post + story_df$Restatement_raw_post
story_df$reflection_proxy_pre <- story_df$ReflectionofFeelings_raw_pre + story_df$Restatement_raw_pre
story_df$total_questions_est_post <- ifelse(story_df$RQR_new_post > 0, story_df$reflection_proxy_post / story_df$RQR_new_post, NA)
story_df$total_questions_est_pre <- ifelse(story_df$RQR_new_PRE > 0, story_df$reflection_proxy_pre / story_df$RQR_new_PRE, NA)
story_df$open_questions_est_post <- story_df$total_questions_est_post * story_df$questioning_post / 100
story_df$closed_questions_est_post <- story_df$total_questions_est_post - story_df$open_questions_est_post

validation_rows <- data.frame(
  metric = c("RQR derived proxy vs stored RQR", "question proxy vs stored % open"),
  r = c(
    cor(
      ifelse(
        (story_df$InquiringSubjectiveInformation_raw_post + story_df$InquiringObjectiveInformation_raw_post) > 0,
        story_df$reflection_proxy_post / (story_df$InquiringSubjectiveInformation_raw_post + story_df$InquiringObjectiveInformation_raw_post),
        NA
      ),
      story_df$RQR_new_post,
      use = "pairwise.complete.obs"
    ),
    cor(
      ifelse(
        (story_df$InquiringSubjectiveInformation_raw_post + story_df$InquiringObjectiveInformation_raw_post) > 0,
        story_df$InquiringSubjectiveInformation_raw_post / (story_df$InquiringSubjectiveInformation_raw_post + story_df$InquiringObjectiveInformation_raw_post) * 100,
        NA
      ),
      story_df$questioning_post,
      use = "pairwise.complete.obs"
    )
  )
)
write.csv(validation_rows, file.path(output_dir, "mi_decomposition_validation.csv"), row.names = FALSE)

run_post_group_models(story_df, "reflection_proxy_post", "mi_reflection_proxy_post")
run_post_group_models(story_df, "total_questions_est_post", "mi_total_questions_est_post")
run_post_group_models(story_df, "open_questions_est_post", "mi_open_questions_est_post")
run_post_group_models(story_df, "closed_questions_est_post", "mi_closed_questions_est_post")

# MI benchmark pre-post transitions
story_df$rqr_basic_pre <- ifelse(!is.na(story_df$RQR_new_PRE), story_df$RQR_new_PRE >= 1, NA)
story_df$rqr_basic_post <- ifelse(!is.na(story_df$RQR_new_post), story_df$RQR_new_post >= 1, NA)
story_df$rqr_advanced_pre <- ifelse(!is.na(story_df$RQR_new_PRE), story_df$RQR_new_PRE >= 2, NA)
story_df$rqr_advanced_post <- ifelse(!is.na(story_df$RQR_new_post), story_df$RQR_new_post >= 2, NA)
story_df$oq_basic_pre <- ifelse(!is.na(story_df$questioning_pre), story_df$questioning_pre >= 50, NA)
story_df$oq_basic_post <- ifelse(!is.na(story_df$questioning_post), story_df$questioning_post >= 50, NA)
story_df$oq_advanced_pre <- ifelse(!is.na(story_df$questioning_pre), story_df$questioning_pre >= 70, NA)
story_df$oq_advanced_post <- ifelse(!is.na(story_df$questioning_post), story_df$questioning_post >= 70, NA)

benchmarks <- list(
  c("rqr_basic_pre", "rqr_basic_post", "rqr_basic"),
  c("rqr_advanced_pre", "rqr_advanced_post", "rqr_advanced"),
  c("oq_basic_pre", "oq_basic_post", "oq_basic"),
  c("oq_advanced_pre", "oq_advanced_post", "oq_advanced")
)

benchmark_counts <- list()
benchmark_tests <- list()
for (spec in benchmarks) {
  pre_col <- spec[1]
  post_col <- spec[2]
  bench_name <- spec[3]
  for (g in levels(story_df$group)) {
    dat <- subset(story_df, group == g & !is.na(story_df[[pre_col]]) & !is.na(story_df[[post_col]]))
    if (nrow(dat) == 0) {
      next
    }
    tab <- table(dat[[pre_col]], dat[[post_col]])
    count_df <- as.data.frame(tab)
    names(count_df) <- c("pre_met", "post_met", "n")
    count_df$group <- g
    count_df$benchmark <- bench_name
    benchmark_counts[[paste(bench_name, g)]] <- count_df

    test <- suppressWarnings(mcnemar.test(tab))
    benchmark_tests[[paste(bench_name, g)]] <- data.frame(
      benchmark = bench_name,
      group = g,
      n = nrow(dat),
      statistic = unname(test$statistic),
      p_value = test$p.value,
      pre_rate = mean(dat[[pre_col]]),
      post_rate = mean(dat[[post_col]])
    )
  }
}
write.csv(do.call(rbind, benchmark_counts), file.path(output_dir, "mi_benchmark_transition_counts.csv"), row.names = FALSE)
write.csv(do.call(rbind, benchmark_tests), file.path(output_dir, "mi_benchmark_transition_tests.csv"), row.names = FALSE)

# Baseline heterogeneity
story_df$baseline_conf_t1 <- story_df$confidence_pre_mean
story_df$baseline_mi_pre <- story_df$mi_style_pre
story_df$baseline_conf_band <- ifelse(story_df$baseline_conf_t1 <= median(story_df$baseline_conf_t1, na.rm = TRUE), "low", "high")
story_df$baseline_mi_band <- ifelse(story_df$baseline_mi_pre <= median(story_df$baseline_mi_pre, na.rm = TRUE), "low", "high")

hetero_rows <- list()
band_desc_rows <- list()
for (moderator in c("baseline_conf_t1", "baseline_mi_pre")) {
  for (outcome in c("HST6", "SMT6", "CCT6", "mi_style_post", "content_exploration_norm_post")) {
    dat <- story_df[complete.cases(story_df[, c(moderator, outcome, "group")]), c(moderator, outcome, "group")]
    if (nrow(dat) < 18) {
      next
    }
    fit <- lm(as.formula(paste0("`", outcome, "` ~ ", moderator, " * group")), data = dat)
    tab <- as.data.frame(anova(fit))
    tab$term <- rownames(tab)
    rownames(tab) <- NULL
    names(tab) <- make.names(names(tab))
    tab$moderator <- moderator
    tab$outcome <- outcome
    hetero_rows[[paste(moderator, outcome)]] <- tab
  }
}

for (band_var in c("baseline_conf_band", "baseline_mi_band")) {
  for (outcome in c("HST6", "SMT6", "CCT6", "mi_style_post", "content_exploration_norm_post")) {
    dat <- story_df[complete.cases(story_df[, c(band_var, outcome, "group")]), c(band_var, outcome, "group")]
    if (nrow(dat) < 18) {
      next
    }
    desc <- aggregate(dat[[outcome]], by = list(group = dat$group, band = dat[[band_var]]), FUN = mean_sd)
    band_desc <- data.frame(
      moderator_band = band_var,
      group = desc$group,
      band = desc$band,
      outcome = outcome,
      n = desc$x[, "n"],
      mean = desc$x[, "mean"],
      sd = desc$x[, "sd"]
    )
    band_desc_rows[[paste(band_var, outcome)]] <- band_desc
  }
}
write.csv(do.call(rbind, hetero_rows), file.path(output_dir, "heterogeneity_interactions.csv"), row.names = FALSE)
write.csv(do.call(rbind, band_desc_rows), file.path(output_dir, "heterogeneity_band_descriptives.csv"), row.names = FALSE)

# Quadratic copying analyses
feedback_df <- merge(
  subset(practice, group == "实验组"),
  story_df[, c("full_id", "HST6", "SMT6", "CCT6", "mi_style_post", "content_exploration_norm_post")],
  by.x = "matched_full_id",
  by.y = "full_id",
  all.x = TRUE
)
feedback_df$copy_centered <- feedback_df$overall_mean_score_weighted - mean(feedback_df$overall_mean_score_weighted, na.rm = TRUE)
feedback_df$copy_centered_sq <- feedback_df$copy_centered^2

quad_rows <- list()
for (outcome in c("HST6", "SMT6", "CCT6", "mi_style_post", "content_exploration_norm_post")) {
  dat <- feedback_df[complete.cases(feedback_df[, c(outcome, "copy_centered", "copy_centered_sq")]), ]
  if (nrow(dat) < 10) {
    next
  }
  fit <- lm(as.formula(paste0("`", outcome, "` ~ copy_centered + copy_centered_sq")), data = dat)
  quad_rows[[outcome]] <- extract_lm_terms(
    fit,
    keep_terms = c("copy_centered", "copy_centered_sq"),
    analysis_name = "copy_quadratic",
    outcome_name = outcome
  )
}
write.csv(do.call(rbind, quad_rows), file.path(output_dir, "copy_quadratic_models.csv"), row.names = FALSE)

# Early setback analyses
story_df$confidence_t2_mean <- rowMeans(story_df[, c("HST2", "SMT2", "CCT2")], na.rm = TRUE)
story_df$early_drop <- story_df$confidence_t2_mean - story_df$confidence_pre_mean

setback_rows <- list()
for (outcome in c("mi_style_post", "content_exploration_norm_post", "HST6", "SMT6", "CCT6")) {
  dat <- story_df[complete.cases(story_df[, c("early_drop", outcome)]), c("early_drop", outcome)]
  if (nrow(dat) < 15) {
    next
  }
  fit <- lm(as.formula(paste0("`", outcome, "` ~ early_drop")), data = dat)
  setback_rows[[outcome]] <- extract_lm_terms(
    fit,
    keep_terms = c("early_drop"),
    analysis_name = "early_setback",
    outcome_name = outcome
  )
}
write.csv(do.call(rbind, setback_rows), file.path(output_dir, "early_setback_models.csv"), row.names = FALSE)

# Feedback user types
feedback_df$copy_band <- ifelse(feedback_df$overall_mean_score_weighted > median(feedback_df$overall_mean_score_weighted, na.rm = TRUE), "high_copy", "low_copy")
feedback_df$feedback_band <- ifelse(feedback_df$feedback_count > median(feedback_df$feedback_count, na.rm = TRUE), "high_feedback", "low_feedback")
feedback_df$user_type <- paste(feedback_df$feedback_band, feedback_df$copy_band, sep = "__")

user_type_rows <- list()
for (outcome in c("HST6", "SMT6", "CCT6", "mi_style_post", "content_exploration_norm_post")) {
  dat <- feedback_df[complete.cases(feedback_df[, c("user_type", outcome)]), c("user_type", outcome)]
  if (nrow(dat) < 10) {
    next
  }
  desc <- aggregate(dat[[outcome]], by = list(user_type = dat$user_type), FUN = mean_sd)
  user_type_rows[[outcome]] <- data.frame(
    outcome = outcome,
    user_type = desc$user_type,
    n = desc$x[, "n"],
    mean = desc$x[, "mean"],
    sd = desc$x[, "sd"]
  )
}
write.csv(do.call(rbind, user_type_rows), file.path(output_dir, "feedback_user_type_descriptives.csv"), row.names = FALSE)

# MI vs exploration trade-off
trade_rows <- list()
for (g in c(levels(story_df$group), "all")) {
  dat <- story_df[complete.cases(story_df[, c("mi_style_post", "content_exploration_norm_post")]), ]
  if (g != "all") {
    dat <- subset(dat, group == g)
  }
  if (nrow(dat) < 10) {
    next
  }
  test <- cor.test(dat$mi_style_post, dat$content_exploration_norm_post)
  trade_rows[[g]] <- data.frame(
    group = g,
    n = nrow(dat),
    r = unname(test$estimate),
    p_value = test$p.value
  )
}
write.csv(do.call(rbind, trade_rows), file.path(output_dir, "mi_exploration_tradeoff_correlations.csv"), row.names = FALSE)

story_df$style_minus_exploration <- zscore(story_df$mi_style_post) - zscore(story_df$content_exploration_norm_post)
run_post_group_models(story_df, "style_minus_exploration", "style_minus_exploration")

notes <- c(
  "# Storyline Analyses",
  "",
  "- `mi_style_post`, `mi_style_pre`: composite of log-transformed RQR and percent open questions.",
  "- `style_conf_gap`: standardized MI-style performance minus standardized confidence composite.",
  "- RQR decomposition uses proxy reflections from `Reflection of Feelings + Restatement` and estimated total/open/closed question counts based on the stored RQR and percent-open variables. These are exploratory approximations.",
  "- `early_drop` is defined as T2 confidence mean minus T1 confidence mean. More negative values indicate a larger early confidence drop.",
  "- `user_type` in the feedback group is based on median splits of feedback clicks and weighted copy score."
)
writeLines(notes, file.path(output_dir, "storyline_analysis_notes.md"))

cat(file.path(output_dir, "storyline_analysis_notes.md"), "\n")
