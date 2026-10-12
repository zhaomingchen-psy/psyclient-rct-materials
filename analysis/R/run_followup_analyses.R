options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(foreign)
  library(nlme)
})

project_root <- getwd()
analysis_dir <- file.path(project_root, "output", "analysis")
output_dir <- file.path(project_root, "output", "r_followup")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

normalize_id <- function(x) {
  x <- as.character(x)
  x <- sub("\\.0$", "", x)
  x[x %in% c("NA", "NaN", "", "NULL")] <- NA_character_
  x
}

group_factor <- function(x) {
  factor(
    x,
    levels = c("空白组", "无反馈组", "实验组"),
    labels = c("blank", "no_feedback", "feedback")
  )
}

fmt_p <- function(x) {
  ifelse(is.na(x), "", ifelse(x < 0.001, "< .001", sub("^0", "", sprintf("%.3f", x))))
}

extract_lm_contrasts <- function(data, outcome, pre = NULL) {
  dat <- data[!is.na(data[[outcome]]) & !is.na(data$group), ]
  dat$c_ai_vs_blank <- ifelse(dat$group == "blank", -2, 1)
  dat$c_feedback_vs_nofb <- ifelse(dat$group == "feedback", 1, ifelse(dat$group == "no_feedback", -1, 0))
  if (!is.null(pre)) {
    dat <- dat[!is.na(dat[[pre]]), ]
  }
  rhs <- c()
  if (!is.null(pre)) {
    rhs <- c(rhs, pre)
  }
  rhs <- c(rhs, "c_ai_vs_blank", "c_feedback_vs_nofb")
  fit <- lm(as.formula(paste0(outcome, " ~ ", paste(rhs, collapse = " + "))), data = dat)
  coef_table <- as.data.frame(summary(fit)$coefficients)
  coef_table$term <- rownames(coef_table)
  rownames(coef_table) <- NULL
  out <- subset(coef_table, term %in% c("c_ai_vs_blank", "c_feedback_vs_nofb"))
  out$contrast <- c("ai_vs_blank", "feedback_vs_nofb")
  out$outcome <- outcome
  out$n <- nrow(dat)
  out$adjusted <- !is.null(pre)
  names(out)[1:4] <- c("estimate", "std_error", "t_value", "p_value")
  out[, c("outcome", "contrast", "adjusted", "n", "estimate", "std_error", "t_value", "p_value")]
}

extract_anova_lme <- function(model, outcome_name, analysis_name) {
  tab <- as.data.frame(anova(model))
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  names(tab) <- make.names(names(tab))
  tab$outcome <- outcome_name
  tab$analysis <- analysis_name
  tab
}

rct <- read.csv(file.path(analysis_dir, "rct_scale_totals_long.csv"), check.names = FALSE)
pm <- read.csv(file.path(analysis_dir, "participant_master.csv"), check.names = FALSE)
practice <- read.csv(file.path(analysis_dir, "practice_student_summary_clean.csv"), check.names = FALSE)
practice_session <- read.csv(file.path(analysis_dir, "practice_session_level_clean.csv"), check.names = FALSE)
sav <- read.spss(file.path(project_root, "Psyclient_data.sav"), to.data.frame = TRUE, use.value.labels = FALSE)
names(sav) <- trimws(names(sav))

rct$full_id <- normalize_id(rct$full_id)
rct$group <- group_factor(rct$group_label)
rct$front_half <- pmin(rct$timepoint_index - 1, 1)
rct$back_half <- pmax(rct$timepoint_index - 2, 0)

pm$full_id <- normalize_id(pm$full_id)
pm$group <- group_factor(pm$group_label)

practice$matched_full_id <- normalize_id(practice$matched_full_id)
practice$practice_student_id <- normalize_id(practice$practice_student_id)
practice <- subset(practice, !(practice_student_id %in% c("TA01", "TA02")))
practice$matched_full_id[practice$practice_student_id == "021"] <- "<STUDENT_ID_021>"
practice$matched_group_label[practice$practice_student_id == "021"] <- "实验组"
practice$roster_match_status[practice$practice_student_id == "021"] <- "manual_feedback_resolution"
practice <- subset(practice, roster_match_status %in% c("unique_short_id_match", "manual_feedback_resolution"))
practice$group <- group_factor(practice$matched_group_label)

practice_session$matched_full_id <- normalize_id(practice_session$matched_full_id)
practice_session$practice_student_id <- normalize_id(practice_session$practice_student_id)
practice_session <- subset(practice_session, !(practice_student_id %in% c("TA01", "TA02")))
practice_session$matched_full_id[practice_session$practice_student_id == "021"] <- "<STUDENT_ID_021>"
practice_session$matched_group_label[practice_session$practice_student_id == "021"] <- "实验组"
practice_session$roster_match_status[practice_session$practice_student_id == "021"] <- "manual_feedback_resolution"
practice_session <- subset(practice_session, roster_match_status %in% c("unique_short_id_match", "manual_feedback_resolution"))
practice_session$group <- group_factor(practice_session$matched_group_label)
practice_session$date <- as.Date(practice_session$date)

sav$full_id <- normalize_id(sav$ID)
sav <- merge(sav, unique(pm[, c("full_id", "group_label")]), by = "full_id", all.x = TRUE)
sav$group <- group_factor(sav$group_label)
sav$question_ratio_pre <- sav$questioning_pre
sav$question_ratio_post <- sav$questioning_post

# 1. Planned contrasts at T6 / post
contrast_rows <- list()
contrast_specs <- list(
  list(outcome = "HST6", pre = "HST1"),
  list(outcome = "SMT6", pre = "SMT1"),
  list(outcome = "CCT6", pre = "CCT1"),
  list(outcome = "RQR_new_post", pre = "RQR_new_PRE"),
  list(outcome = "question_ratio_post", pre = "question_ratio_pre")
)
for (spec in contrast_specs) {
  contrast_rows[[paste0(spec$outcome, "_raw")]] <- extract_lm_contrasts(sav, spec$outcome, NULL)
  contrast_rows[[paste0(spec$outcome, "_adj")]] <- extract_lm_contrasts(sav, spec$outcome, spec$pre)
}
planned_contrasts <- do.call(rbind, contrast_rows)
planned_contrasts$p_value_holm <- ave(planned_contrasts$p_value, planned_contrasts$adjusted, FUN = function(x) p.adjust(x, method = "holm"))
write.csv(planned_contrasts, file.path(output_dir, "planned_contrasts_t6.csv"), row.names = FALSE)

# 2. Missing data sensitivity
missing_counts <- as.data.frame(with(rct, table(timepoint, group)))
write.csv(missing_counts, file.path(output_dir, "missing_counts_by_group_time.csv"), row.names = FALSE)

t1_base <- subset(
  rct,
  timepoint == "T1",
  select = c("full_id", "group", "hs_total", "sm_total", "cc_total", "self_efficacy_total")
)
t6_presence <- subset(rct, timepoint == "T6", select = c("full_id"))
t6_presence$has_t6 <- 1
t1_missing <- merge(t1_base, t6_presence, by = "full_id", all.x = TRUE)
t1_missing$has_t6[is.na(t1_missing$has_t6)] <- 0

missing_rows <- list()
for (outcome in c("hs_total", "sm_total", "cc_total", "self_efficacy_total")) {
  fit <- lm(as.formula(paste0(outcome, " ~ has_t6 + group")), data = t1_missing)
  tab <- as.data.frame(anova(fit))
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  names(tab) <- make.names(names(tab))
  tab$outcome <- outcome
  missing_rows[[outcome]] <- tab
}
write.csv(do.call(rbind, missing_rows), file.path(output_dir, "baseline_predictors_of_t6_completion.csv"), row.names = FALSE)

complete_ids <- intersect(
  unique(subset(rct, timepoint == "T1" & !is.na(hs_total))$full_id),
  unique(subset(rct, timepoint == "T6" & !is.na(hs_total))$full_id)
)

complete_case_rows <- list()
for (outcome in c("hs_total", "sm_total", "cc_total")) {
  dat <- subset(rct, full_id %in% complete_ids & !is.na(group) & !is.na(rct[[outcome]]))
  fit <- lme(
    as.formula(paste0(outcome, " ~ front_half * group + back_half * group")),
    random = ~1 | full_id,
    data = dat,
    na.action = na.omit,
    method = "REML"
  )
  complete_case_rows[[outcome]] <- extract_anova_lme(fit, outcome, "complete_case_hlm")
}
write.csv(do.call(rbind, complete_case_rows), file.path(output_dir, "complete_case_hlm.csv"), row.names = FALSE)

# 3. Feedback mechanism analyses
feedback_practice <- subset(practice, group == "feedback")
feedback_merged <- merge(
  feedback_practice,
  pm[, c("full_id", "hs_total_t1", "sm_total_t1", "cc_total_t1", "self_efficacy_total_t1", "hs_total_t6", "sm_total_t6", "cc_total_t6", "self_efficacy_total_t6")],
  by.x = "matched_full_id",
  by.y = "full_id",
  all.x = TRUE
)

copy_fit <- lm(overall_mean_score_weighted ~ feedback_count, data = feedback_merged)
copy_tab <- as.data.frame(summary(copy_fit)$coefficients)
copy_tab$term <- rownames(copy_tab)
rownames(copy_tab) <- NULL
names(copy_tab)[1:4] <- c("estimate", "std_error", "t_value", "p_value")
copy_tab$model <- "weighted_copy_from_feedback_count"
write.csv(copy_tab, file.path(output_dir, "feedback_to_copy_model.csv"), row.names = FALSE)

mechanism_rows <- list()
for (outcome in c("hs_total", "sm_total", "cc_total", "self_efficacy_total")) {
  fit <- lm(
    as.formula(
      paste0(outcome, "_t6 ~ ", outcome, "_t1 + feedback_count + overall_mean_score_weighted")
    ),
    data = feedback_merged
  )
  tab <- as.data.frame(summary(fit)$coefficients)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  tab <- subset(tab, term %in% c("feedback_count", "overall_mean_score_weighted"))
  names(tab)[1:4] <- c("estimate", "std_error", "t_value", "p_value")
  tab$outcome <- outcome
  mechanism_rows[[outcome]] <- tab[, c("outcome", "term", "estimate", "std_error", "t_value", "p_value")]
}
write.csv(do.call(rbind, mechanism_rows), file.path(output_dir, "copying_and_feedback_predict_t6.csv"), row.names = FALSE)

# 4. Baseline moderation
moderation_rows <- list()
band_rows <- list()
for (outcome in c("hs_total", "sm_total", "cc_total")) {
  post <- paste0(outcome, "_t6")
  pre <- paste0(outcome, "_t1")
  dat <- pm[, c("group", pre, post)]
  dat <- dat[complete.cases(dat), ]
  fit <- lm(as.formula(paste0(post, " ~ ", pre, " * group")), data = dat)
  tab <- as.data.frame(anova(fit))
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  names(tab) <- make.names(names(tab))
  tab$outcome <- outcome
  moderation_rows[[outcome]] <- tab

  med <- median(dat[[pre]], na.rm = TRUE)
  dat$baseline_band <- ifelse(dat[[pre]] <= med, "low", "high")
  desc <- aggregate(dat[[post]], by = list(group = dat$group, baseline_band = dat$baseline_band), FUN = function(x) c(n = sum(!is.na(x)), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)))
  band_df <- data.frame(
    outcome = outcome,
    group = desc$group,
    baseline_band = desc$baseline_band,
    n = desc$x[, "n"],
    mean = desc$x[, "mean"],
    sd = desc$x[, "sd"],
    median_cut = med
  )
  band_rows[[outcome]] <- band_df
}
write.csv(do.call(rbind, moderation_rows), file.path(output_dir, "baseline_moderation_anova.csv"), row.names = FALSE)
write.csv(do.call(rbind, band_rows), file.path(output_dir, "baseline_moderation_band_descriptives.csv"), row.names = FALSE)

# 5. MI benchmark analyses
sav$rqr_basic <- ifelse(!is.na(sav$RQR_new_post), sav$RQR_new_post >= 1, NA)
sav$oq_basic <- ifelse(!is.na(sav$question_ratio_post), sav$question_ratio_post >= 50, NA)
sav$oq_advanced <- ifelse(!is.na(sav$question_ratio_post), sav$question_ratio_post >= 70, NA)

benchmark_counts <- list()
benchmark_tests <- list()
pairwise_rows <- list()

for (var_name in c("rqr_basic", "oq_basic", "oq_advanced")) {
  dat <- sav[!is.na(sav[[var_name]]) & !is.na(sav$group), c("group", var_name)]
  tab <- table(dat$group, dat[[var_name]])
  count_df <- as.data.frame(tab)
  names(count_df) <- c("group", "benchmark_met", "n")
  count_df$benchmark <- var_name
  benchmark_counts[[var_name]] <- count_df

  chi <- suppressWarnings(chisq.test(tab))
  benchmark_tests[[var_name]] <- data.frame(
    benchmark = var_name,
    test = "chi_square",
    statistic = unname(chi$statistic),
    df = unname(chi$parameter),
    p_value = chi$p.value
  )

  pairs <- list(
    c("feedback", "no_feedback"),
    c("feedback", "blank"),
    c("no_feedback", "blank")
  )
  for (pair in pairs) {
    pair_dat <- subset(dat, group %in% pair)
    pair_tab <- table(pair_dat$group, factor(pair_dat[[var_name]], levels = c(FALSE, TRUE)))
    if (all(dim(pair_tab) == c(2, 2))) {
      ft <- fisher.test(pair_tab)
      pairwise_rows[[paste(var_name, paste(pair, collapse = "_"), sep = "__")]] <- data.frame(
        benchmark = var_name,
        comparison = paste(pair, collapse = "_vs_"),
        odds_ratio = unname(ft$estimate),
        p_value = ft$p.value
      )
    }
  }
}

write.csv(do.call(rbind, benchmark_counts), file.path(output_dir, "mi_benchmark_counts.csv"), row.names = FALSE)
write.csv(do.call(rbind, benchmark_tests), file.path(output_dir, "mi_benchmark_tests.csv"), row.names = FALSE)
write.csv(do.call(rbind, pairwise_rows), file.path(output_dir, "mi_benchmark_pairwise_fisher.csv"), row.names = FALSE)

# 6. Session-level feedback effects in feedback group
psf <- subset(practice_session, group == "feedback")
psf <- psf[order(psf$matched_full_id, psf$date), ]
psf$next_mean_score <- ave(psf$mean_score, psf$matched_full_id, FUN = function(x) c(x[-1], NA))
psf$lag_mean_score <- ave(psf$mean_score, psf$matched_full_id, FUN = function(x) c(NA, x[-length(x)]))
psf$next_total_score <- ave(psf$total_score, psf$matched_full_id, FUN = function(x) c(x[-1], NA))
psf$lag_total_score <- ave(psf$total_score, psf$matched_full_id, FUN = function(x) c(NA, x[-length(x)]))

session_specs <- list(
  list(outcome = "next_mean_score", rhs = c("feedback_count", "mean_score", "lag_mean_score")),
  list(outcome = "next_total_score", rhs = c("feedback_count", "total_score", "lag_total_score"))
)

session_rows <- list()
session_anova <- list()
for (spec in session_specs) {
  dat <- psf[, c("matched_full_id", spec$outcome, spec$rhs)]
  dat <- dat[complete.cases(dat), ]
  names(dat)[1] <- "full_id"
  fit <- lme(
    as.formula(paste0(spec$outcome, " ~ ", paste(spec$rhs, collapse = " + "))),
    random = ~1 | full_id,
    data = dat,
    na.action = na.omit,
    method = "REML"
  )
  tab <- as.data.frame(summary(fit)$tTable)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  names(tab)[1:5] <- c("estimate", "std_error", "df", "t_value", "p_value")
  tab$outcome <- spec$outcome
  session_rows[[spec$outcome]] <- tab
  session_anova[[spec$outcome]] <- extract_anova_lme(fit, spec$outcome, "feedback_session_effect")
}

write.csv(do.call(rbind, session_rows), file.path(output_dir, "feedback_session_effects_coefficients.csv"), row.names = FALSE)
write.csv(do.call(rbind, session_anova), file.path(output_dir, "feedback_session_effects_anova.csv"), row.names = FALSE)

# Markdown summary
lines <- c(
  "# Follow-up Analyses",
  "",
  "## Files",
  "- planned_contrasts_t6.csv: planned contrasts comparing AI groups vs blank and feedback vs no-feedback",
  "- missing_counts_by_group_time.csv, baseline_predictors_of_t6_completion.csv, complete_case_hlm.csv: missing-data sensitivity analyses",
  "- feedback_to_copy_model.csv, copying_and_feedback_predict_t6.csv: feedback mechanism models",
  "- baseline_moderation_anova.csv, baseline_moderation_band_descriptives.csv: baseline moderation analyses",
  "- mi_benchmark_counts.csv, mi_benchmark_tests.csv, mi_benchmark_pairwise_fisher.csv: MI benchmark analyses",
  "- feedback_session_effects_coefficients.csv, feedback_session_effects_anova.csv: session-level feedback effect models"
)
writeLines(lines, con = file.path(output_dir, "followup_output_guide.md"))
