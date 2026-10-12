options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(foreign)
})

root <- getwd()  # project folder that contains output/analysis
output_dir <- file.path(root, "output", "r_feedback_rqr_variability")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

normalize_id <- function(x) {
  x <- as.character(x)
  x <- sub("\\.0$", "", x)
  x[x %in% c("NA", "NaN", "", "NULL")] <- NA_character_
  x
}

extract_lm_term <- function(fit, term, outcome, predictor, subset_name) {
  tab <- as.data.frame(summary(fit)$coefficients)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  names(tab)[1:4] <- c("estimate", "std_error", "t_value", "p_value")
  row <- tab[tab$term == term, ]
  row$outcome <- outcome
  row$predictor <- predictor
  row$subset <- subset_name
  row$r_squared <- summary(fit)$r.squared
  row$n <- nobs(fit)
  row[, c("subset", "outcome", "predictor", "term", "estimate", "std_error", "t_value", "p_value", "r_squared", "n")]
}

safe_cor <- function(x, y, method, outcome, predictor, subset_name) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 4) {
    return(NULL)
  }
  test <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  data.frame(
    subset = subset_name,
    outcome = outcome,
    predictor = predictor,
    method = method,
    n = sum(ok),
    estimate = unname(test$estimate),
    p_value = test$p.value
  )
}

practice <- read.csv(file.path(root, "output", "analysis", "practice_student_summary_clean.csv"), check.names = FALSE)
practice$matched_full_id <- normalize_id(practice$matched_full_id)
practice$practice_student_id <- normalize_id(practice$practice_student_id)
practice <- subset(practice, !(practice_student_id %in% c("TA01", "TA02")))
practice$matched_full_id[practice$practice_student_id == "021"] <- "<STUDENT_ID_021>"
practice$matched_group_label[practice$practice_student_id == "021"] <- "实验组"

sav <- read.spss(file.path(root, "Psyclient_data.sav"), to.data.frame = TRUE, use.value.labels = FALSE)
names(sav) <- trimws(names(sav))
sav$full_id <- normalize_id(sav$ID)
sav$group <- factor(sav$Group, levels = c(3, 2, 1), labels = c("空白组", "无反馈组", "实验组"))

feedback <- merge(
  subset(sav, group == "实验组", select = c(full_id, RQR_new_post, RQR_new_PRE, questioning_post, HST6, SMT6, CCT6)),
  practice[, c(
    "matched_full_id", "feedback_count", "session_count", "total_score", "valid_scores_n",
    "valid_scores_sum", "overall_mean_score_weighted"
  )],
  by.x = "full_id",
  by.y = "matched_full_id",
  all.x = TRUE
)

feedback$feedback_rate <- feedback$feedback_count / feedback$session_count
feedback$copy_rate <- feedback$overall_mean_score_weighted
feedback$log_rqr_post <- log1p(feedback$RQR_new_post)
feedback$confidence_t6_mean <- rowMeans(feedback[, c("HST6", "SMT6", "CCT6")], na.rm = TRUE)

write.csv(feedback[order(feedback$RQR_new_post, na.last = TRUE), ], file.path(output_dir, "feedback_rqr_case_listing.csv"), row.names = FALSE)

dist <- data.frame(
  metric = c("RQR_new_post", "log1p_RQR_new_post"),
  n = c(sum(!is.na(feedback$RQR_new_post)), sum(!is.na(feedback$log_rqr_post))),
  mean = c(mean(feedback$RQR_new_post, na.rm = TRUE), mean(feedback$log_rqr_post, na.rm = TRUE)),
  median = c(median(feedback$RQR_new_post, na.rm = TRUE), median(feedback$log_rqr_post, na.rm = TRUE)),
  sd = c(sd(feedback$RQR_new_post, na.rm = TRUE), sd(feedback$log_rqr_post, na.rm = TRUE)),
  se = c(
    sd(feedback$RQR_new_post, na.rm = TRUE) / sqrt(sum(!is.na(feedback$RQR_new_post))),
    sd(feedback$log_rqr_post, na.rm = TRUE) / sqrt(sum(!is.na(feedback$log_rqr_post)))
  ),
  min = c(min(feedback$RQR_new_post, na.rm = TRUE), min(feedback$log_rqr_post, na.rm = TRUE)),
  max = c(max(feedback$RQR_new_post, na.rm = TRUE), max(feedback$log_rqr_post, na.rm = TRUE))
)
write.csv(dist, file.path(output_dir, "feedback_rqr_distribution.csv"), row.names = FALSE)

predictors <- c("feedback_count", "feedback_rate", "session_count", "total_score", "copy_rate")
outcomes <- c("RQR_new_post", "log_rqr_post", "questioning_post", "confidence_t6_mean")
cor_rows <- list()
lm_rows <- list()

run_subset <- function(data, subset_name) {
  cor_rows_local <- list()
  lm_rows_local <- list()
  for (outcome in outcomes) {
    for (predictor in predictors) {
      if (!all(c(outcome, predictor) %in% names(data))) next
      cor_rows_local[[paste(subset_name, outcome, predictor, "pearson")]] <- safe_cor(data[[outcome]], data[[predictor]], "pearson", outcome, predictor, subset_name)
      cor_rows_local[[paste(subset_name, outcome, predictor, "spearman")]] <- safe_cor(data[[outcome]], data[[predictor]], "spearman", outcome, predictor, subset_name)
      dat <- data[complete.cases(data[, c(outcome, predictor)]), ]
      if (nrow(dat) >= 8) {
        fit <- lm(as.formula(paste(outcome, "~", predictor)), data = dat)
        lm_rows_local[[paste(subset_name, outcome, predictor)]] <- extract_lm_term(fit, predictor, outcome, predictor, subset_name)
      }
    }
  }
  list(cor = cor_rows_local, lm = lm_rows_local)
}

res_all <- run_subset(feedback, "all_feedback")
max_rqr <- max(feedback$RQR_new_post, na.rm = TRUE)
res_no_max <- run_subset(subset(feedback, is.na(RQR_new_post) | RQR_new_post < max_rqr), "without_max_rqr")
res_no_top3 <- run_subset(subset(feedback, is.na(RQR_new_post) | RQR_new_post < sort(feedback$RQR_new_post, decreasing = TRUE, na.last = NA)[3]), "without_top3_rqr")

cor_rows <- c(res_all$cor, res_no_max$cor, res_no_top3$cor)
lm_rows <- c(res_all$lm, res_no_max$lm, res_no_top3$lm)

cor_out <- do.call(rbind, cor_rows[!vapply(cor_rows, is.null, logical(1))])
lm_out <- do.call(rbind, lm_rows[!vapply(lm_rows, is.null, logical(1))])

write.csv(cor_out, file.path(output_dir, "feedback_rqr_usage_correlations.csv"), row.names = FALSE)
write.csv(lm_out, file.path(output_dir, "feedback_rqr_usage_lm_models.csv"), row.names = FALSE)

notes <- c(
  "# Feedback RQR Variability Analysis",
  "",
  "- Analyses are restricted to the feedback group.",
  "- `feedback_rate = feedback_count / session_count`.",
  "- `copy_rate = overall_mean_score_weighted`, the weighted mean copying proportion.",
  "- `log_rqr_post = log1p(RQR_new_post)` is used as a sensitivity check because RQR is highly right-skewed.",
  "- Sensitivity subsets remove either the single highest RQR case or the top three RQR cases."
)
writeLines(notes, file.path(output_dir, "feedback_rqr_variability_notes.md"))
