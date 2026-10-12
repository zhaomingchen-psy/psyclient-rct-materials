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

output_dir <- file.path(project_root, "output", "r_question_progression")
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
    levels = c(3, 2, 1),
    labels = c("空白组", "无反馈组", "实验组")
  )
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

extract_pairwise_t <- function(data, outcome_name) {
  dat <- subset(data, group %in% c("无反馈组", "实验组") & !is.na(data[[outcome_name]]))
  dat$group2 <- factor(dat$group, levels = c("无反馈组", "实验组"))
  fit <- t.test(dat[[outcome_name]] ~ dat$group2)
  m1 <- mean(dat[[outcome_name]][dat$group2 == "无反馈组"], na.rm = TRUE)
  m2 <- mean(dat[[outcome_name]][dat$group2 == "实验组"], na.rm = TRUE)
  s1 <- sd(dat[[outcome_name]][dat$group2 == "无反馈组"], na.rm = TRUE)
  s2 <- sd(dat[[outcome_name]][dat$group2 == "实验组"], na.rm = TRUE)
  n1 <- sum(dat$group2 == "无反馈组" & !is.na(dat[[outcome_name]]))
  n2 <- sum(dat$group2 == "实验组" & !is.na(dat[[outcome_name]]))
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  d <- (m1 - m2) / sp
  data.frame(
    outcome = outcome_name,
    group1 = "无反馈组",
    group2 = "实验组",
    n1 = n1,
    n2 = n2,
    mean1 = m1,
    mean2 = m2,
    t = unname(fit$statistic),
    df = unname(fit$parameter),
    p_value = fit$p.value,
    cohen_d = d
  )
}

sav <- read.spss(
  file.path(project_root, "Psyclient_data.sav"),
  to.data.frame = TRUE,
  use.value.labels = FALSE
)
names(sav) <- trimws(names(sav))
sav$full_id <- normalize_id(sav$ID)
sav$group <- group_factor(sav$Group)

raw_post_cols <- c(
  "AffirmationandReassurance_raw_post",
  "Answer_raw_post",
  "Confrontation_raw_post",
  "InquiringObjectiveInformation_raw_post",
  "InquiringSubjectiveInformation_raw_post",
  "Interpretation_raw_post",
  "InvitetoExploreNewActions_raw_post",
  "InvitetoTakeNewPerspectives_raw_post",
  "MinimalEncouragement_raw_post",
  "Others_raw_post",
  "ReflectionofFeelings_raw_post",
  "Restatement_raw_post",
  "Selfdisclosure_raw_post"
)

for (col in raw_post_cols) {
  sav[[col]][is.na(sav[[col]])] <- 0
}

sav$total_question_count <- sav$InquiringObjectiveInformation_raw_post + sav$InquiringSubjectiveInformation_raw_post
sav$reflection_count <- sav$ReflectionofFeelings_raw_post + sav$Restatement_raw_post
sav$total_technique_count <- rowSums(sav[, raw_post_cols], na.rm = TRUE)
sav$question_density <- ifelse(sav$total_technique_count > 0, sav$total_question_count / sav$total_technique_count, NA)

desc_rows <- list()
pair_rows <- list()
for (outcome in c("total_question_count", "question_density", "reflection_count", "total_technique_count")) {
  dat <- subset(sav, !is.na(group) & !is.na(sav[[outcome]]))
  desc <- aggregate(dat[[outcome]], by = list(group = dat$group), FUN = mean_sd)
  desc_rows[[outcome]] <- data.frame(
    outcome = outcome,
    group = desc$group,
    n = desc$x[, "n"],
    mean = desc$x[, "mean"],
    sd = desc$x[, "sd"]
  )
  fit <- aov(as.formula(paste0(outcome, " ~ group")), data = dat)
  write.csv(extract_simple_anova(fit, outcome), file.path(output_dir, paste0(outcome, "_anova.csv")), row.names = FALSE)
  write.csv(extract_tukey(fit, outcome), file.path(output_dir, paste0(outcome, "_tukey.csv")), row.names = FALSE)
  pair_rows[[outcome]] <- extract_pairwise_t(dat, outcome)
}

write.csv(do.call(rbind, desc_rows), file.path(output_dir, "question_progression_descriptives.csv"), row.names = FALSE)
write.csv(do.call(rbind, pair_rows), file.path(output_dir, "no_feedback_vs_feedback_pairwise.csv"), row.names = FALSE)

notes <- c(
  "# Question-Driven Progression Analyses",
  "",
  "- `total_question_count = InquiringSubjectiveInformation_raw_post + InquiringObjectiveInformation_raw_post`.",
  "- `reflection_count = ReflectionofFeelings_raw_post + Restatement_raw_post`.",
  "- `question_density = total_question_count / total_technique_count`.",
  "- `total_technique_count` is the sum of all 13 technique-level raw post counts in `Psyclient_data.sav`.",
  "- Pairwise tests focus on `无反馈组` vs `实验组` to evaluate the question-driven progression interpretation."
)
writeLines(notes, file.path(output_dir, "question_progression_notes.md"))
