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
output_dir <- file.path(project_root, "output", "r_calibration_groups")
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

compare_independent_rs <- function(r1, n1, r2, n2, label1, label2, outcome_name) {
  if (is.na(r1) || is.na(r2) || n1 <= 3 || n2 <= 3 || abs(r1) >= 1 || abs(r2) >= 1) {
    return(NULL)
  }
  z1 <- atanh(r1)
  z2 <- atanh(r2)
  z_stat <- (z1 - z2) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  p_value <- 2 * pnorm(-abs(z_stat))
  data.frame(
    outcome = outcome_name,
    group1 = label1,
    group2 = label2,
    r1 = r1,
    n1 = n1,
    r2 = r2,
    n2 = n2,
    z = z_stat,
    p_value = p_value
  )
}

pm <- read.csv(file.path(analysis_dir, "participant_master.csv"), check.names = FALSE)
pm$full_id <- normalize_id(pm$full_id)
pm$group <- group_factor_zh(pm$group_label)

practice <- read.csv(file.path(analysis_dir, "practice_student_summary_clean.csv"), check.names = FALSE)
practice$matched_full_id <- normalize_id(practice$matched_full_id)
practice$practice_student_id <- normalize_id(practice$practice_student_id)
practice <- subset(practice, !(practice_student_id %in% c("TA01", "TA02")))
practice$matched_full_id[practice$practice_student_id == "021"] <- "<STUDENT_ID_021>"
practice$matched_group_label[practice$practice_student_id == "021"] <- "实验组"
practice$roster_match_status[practice$practice_student_id == "021"] <- "manual_feedback_resolution"
practice <- subset(practice, roster_match_status %in% c("unique_short_id_match", "manual_feedback_resolution"))

sav <- read.spss(
  file.path(project_root, "Psyclient_data.sav"),
  to.data.frame = TRUE,
  use.value.labels = FALSE
)
names(sav) <- trimws(names(sav))
sav$full_id <- normalize_id(sav$ID)
sav$group <- group_factor_sav(sav$Group)

calibration_df <- merge(
  sav[, c("full_id", "group", "HST6", "SMT6", "CCT6", "RQR_new_post", "questioning_post")],
  practice[, c("matched_full_id", "overall_mean_score_weighted", "total_score", "feedback_count", "session_count")],
  by.x = "full_id",
  by.y = "matched_full_id",
  all.x = TRUE
)

calibration_df$confidence_t6_mean <- rowMeans(
  calibration_df[, c("HST6", "SMT6", "CCT6")],
  na.rm = TRUE
)
calibration_df$objective_mi_post <- (zscore(log1p(calibration_df$RQR_new_post)) + zscore(calibration_df$questioning_post)) / 2
calibration_df$confidence_z <- zscore(calibration_df$confidence_t6_mean)
calibration_df$objective_z <- zscore(calibration_df$objective_mi_post)
calibration_df$discrepancy_signed <- calibration_df$confidence_z - calibration_df$objective_z
calibration_df$discrepancy_abs <- abs(calibration_df$discrepancy_signed)

# 1. T6 self-efficacy and objective performance correlations by group
corr_rows <- list()
for (g in c(levels(calibration_df$group), "all")) {
  if (g == "all") {
    dat <- subset(calibration_df, !is.na(confidence_t6_mean) & !is.na(objective_mi_post))
  } else {
    dat <- subset(calibration_df, group == g & !is.na(confidence_t6_mean) & !is.na(objective_mi_post))
  }
  if (nrow(dat) < 4) {
    next
  }
  test <- cor.test(dat$confidence_t6_mean, dat$objective_mi_post)
  corr_rows[[g]] <- data.frame(
    group = g,
    n = nrow(dat),
    r = unname(test$estimate),
    p_value = test$p.value
  )
}
corr_out <- do.call(rbind, corr_rows)
write.csv(corr_out, file.path(output_dir, "t6_selfeff_objective_correlations.csv"), row.names = FALSE)

group_corrs <- subset(corr_out, group != "all")
pair_rows <- list()
if (nrow(group_corrs) >= 2) {
  for (i in 1:(nrow(group_corrs) - 1)) {
    for (j in (i + 1):nrow(group_corrs)) {
      pair_rows[[paste(i, j, sep = "_")]] <- compare_independent_rs(
        r1 = group_corrs$r[i],
        n1 = group_corrs$n[i],
        r2 = group_corrs$r[j],
        n2 = group_corrs$n[j],
        label1 = as.character(group_corrs$group[i]),
        label2 = as.character(group_corrs$group[j]),
        outcome_name = "confidence_t6_mean_vs_objective_mi_post"
      )
    }
  }
}
if (length(pair_rows) > 0) {
  write.csv(do.call(rbind, pair_rows), file.path(output_dir, "group_correlation_comparisons.csv"), row.names = FALSE)
}

# 2. Discrepancy score comparisons
disc_rows <- list()
for (outcome in c("confidence_t6_mean", "objective_mi_post", "discrepancy_signed", "discrepancy_abs")) {
  dat <- subset(calibration_df, !is.na(group) & !is.na(calibration_df[[outcome]]))
  desc <- aggregate(dat[[outcome]], by = list(group = dat$group), FUN = mean_sd)
  disc_rows[[outcome]] <- data.frame(
    outcome = outcome,
    group = desc$group,
    n = desc$x[, "n"],
    mean = desc$x[, "mean"],
    sd = desc$x[, "sd"]
  )
  fit <- aov(as.formula(paste0(outcome, " ~ group")), data = dat)
  write.csv(extract_simple_anova(fit, outcome), file.path(output_dir, paste0(outcome, "_anova.csv")), row.names = FALSE)
  write.csv(extract_tukey(fit, outcome), file.path(output_dir, paste0(outcome, "_tukey.csv")), row.names = FALSE)
}
write.csv(do.call(rbind, disc_rows), file.path(output_dir, "calibration_group_descriptives.csv"), row.names = FALSE)

# 3. Feedback-group copying and calibration
feedback_df <- subset(
  calibration_df,
  group == "实验组" &
    !is.na(overall_mean_score_weighted) &
    !is.na(confidence_t6_mean) &
    !is.na(objective_mi_post) &
    !is.na(discrepancy_abs)
)

copy_rows <- list()
for (outcome in c("objective_mi_post", "confidence_t6_mean", "discrepancy_abs", "discrepancy_signed")) {
  dat <- subset(feedback_df, !is.na(feedback_df[[outcome]]))
  if (nrow(dat) < 8) {
    next
  }
  fit <- lm(as.formula(paste0(outcome, " ~ overall_mean_score_weighted")), data = dat)
  tab <- as.data.frame(summary(fit)$coefficients)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  tab <- subset(tab, term == "overall_mean_score_weighted")
  names(tab)[1:4] <- c("estimate", "std_error", "t_value", "p_value")
  tab$outcome <- outcome
  tab$r_squared <- summary(fit)$r.squared
  tab$n <- nobs(fit)
  copy_rows[[outcome]] <- tab[, c("outcome", "term", "estimate", "std_error", "t_value", "p_value", "r_squared", "n")]
}
if (length(copy_rows) > 0) {
  write.csv(do.call(rbind, copy_rows), file.path(output_dir, "feedback_copying_calibration_models.csv"), row.names = FALSE)
}

if (nrow(feedback_df) >= 8) {
  med <- median(feedback_df$overall_mean_score_weighted, na.rm = TRUE)
  feedback_df$copy_band <- ifelse(feedback_df$overall_mean_score_weighted <= med, "low_copy", "high_copy")
  band_rows <- list()
  for (outcome in c("objective_mi_post", "confidence_t6_mean", "discrepancy_abs", "discrepancy_signed")) {
    desc <- aggregate(feedback_df[[outcome]], by = list(copy_band = feedback_df$copy_band), FUN = mean_sd)
    band_rows[[outcome]] <- data.frame(
      outcome = outcome,
      copy_band = desc$copy_band,
      n = desc$x[, "n"],
      mean = desc$x[, "mean"],
      sd = desc$x[, "sd"],
      median_cut = med
    )
  }
  write.csv(do.call(rbind, band_rows), file.path(output_dir, "feedback_copy_band_descriptives.csv"), row.names = FALSE)
}

notes <- c(
  "# Calibration Group Analyses",
  "",
  "- `confidence_t6_mean` is the mean of HST6, SMT6, and CCT6.",
  "- `objective_mi_post` is the mean of standardized `log1p(RQR_new_post)` and standardized `% open questions` (`questioning_post`).",
  "- `discrepancy_signed = confidence_z - objective_z`; positive values indicate confidence is higher than objective MI-style performance.",
  "- `discrepancy_abs` is the absolute value of that gap; smaller values indicate better calibration.",
  "- Copying analyses are limited to the feedback group and use `overall_mean_score_weighted`."
)
writeLines(notes, file.path(output_dir, "calibration_group_notes.md"))
