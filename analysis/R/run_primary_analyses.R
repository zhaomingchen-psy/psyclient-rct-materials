options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(foreign)
  library(nlme)
})

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

project_root <- getwd()
output_dir <- file.path(project_root, "output", "r_analysis")
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

extract_fixed_effects <- function(model, analysis_name, outcome_name) {
  table <- as.data.frame(summary(model)$tTable)
  table$term <- rownames(table)
  rownames(table) <- NULL
  names(table) <- c("estimate", "std_error", "df", "t_value", "p_value", "term")
  table$analysis <- analysis_name
  table$outcome <- outcome_name
  table[, c("analysis", "outcome", "term", "estimate", "std_error", "df", "t_value", "p_value")]
}

extract_anova_table <- function(model, analysis_name, outcome_name) {
  table <- as.data.frame(anova(model))
  table$term <- rownames(table)
  rownames(table) <- NULL
  names(table) <- make.names(names(table))
  table$analysis <- analysis_name
  table$outcome <- outcome_name
  table
}

safe_lme <- function(formula, data) {
  tryCatch(
    lme(formula, random = ~1 | full_id, data = data, na.action = na.omit, method = "REML"),
    error = function(e) e
  )
}

safe_lm <- function(formula, data) {
  tryCatch(
    lm(formula, data = data),
    error = function(e) e
  )
}

make_mean_sd_n <- function(x) {
  c(
    n = sum(!is.na(x)),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE)
  )
}

extract_simple_anova <- function(fit, outcome_name) {
  table <- as.data.frame(summary(fit)[[1]])
  table$term <- rownames(table)
  rownames(table) <- NULL
  names(table) <- make.names(names(table))
  table$outcome <- outcome_name
  table
}

extract_tukey <- function(fit, outcome_name) {
  table <- as.data.frame(TukeyHSD(fit)$group)
  table$comparison <- rownames(table)
  rownames(table) <- NULL
  names(table) <- c("diff", "lwr", "upr", "p_adj", "comparison")
  table$outcome <- outcome_name
  table
}

extract_ancova_pairwise <- function(data, post_col, pre_col, outcome_name) {
  refs <- levels(droplevels(data$group))
  rows <- list()

  for (ref in refs) {
    dat <- data
    dat$group_ref <- relevel(dat$group, ref = ref)
    fit <- lm(as.formula(paste0("`", post_col, "` ~ `", pre_col, "` + group_ref")), data = dat)
    coef_table <- as.data.frame(summary(fit)$coefficients)
    coef_table$term <- rownames(coef_table)
    rownames(coef_table) <- NULL
    group_rows <- subset(coef_table, grepl("^group_ref", term))
    if (nrow(group_rows) == 0) {
      next
    }
    group_rows$other_group <- sub("^group_ref", "", group_rows$term)
    group_rows$comparison <- paste(group_rows$other_group, ref, sep = "-")
    group_rows$pair_key <- vapply(
      strsplit(group_rows$comparison, "-", fixed = TRUE),
      function(parts) paste(sort(parts), collapse = "__"),
      character(1)
    )
    group_rows$ref_group <- ref
    rows[[ref]] <- group_rows[, c("comparison", "pair_key", "Estimate", "Std. Error", "t value", "Pr(>|t|)", "ref_group")]
  }

  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0) {
    return(data.frame())
  }

  out <- out[!duplicated(out$pair_key), ]
  names(out) <- c("comparison", "pair_key", "estimate", "std_error", "t_value", "p_value", "ref_group")
  out$outcome <- outcome_name
  out$p_value_holm <- p.adjust(out$p_value, method = "holm")
  out[order(out$p_value), c("outcome", "comparison", "estimate", "std_error", "t_value", "p_value", "p_value_holm", "ref_group")]
}

run_single_post_models <- function(data, post_col, analysis_prefix, outcome_name, pre_col = NULL) {
  post_dat <- subset(data, !is.na(data[[post_col]]) & !is.na(group))
  descriptives <- aggregate(post_dat[[post_col]], by = list(group = post_dat$group), FUN = make_mean_sd_n)
  names(descriptives) <- c("group", "stats")
  descriptives <- data.frame(
    group = descriptives$group,
    n = descriptives$stats[, "n"],
    mean = descriptives$stats[, "mean"],
    sd = descriptives$stats[, "sd"],
    outcome = outcome_name
  )
  write.csv(
    descriptives,
    file.path(output_dir, paste0(analysis_prefix, "_descriptives.csv")),
    row.names = FALSE
  )

  anova_fit <- aov(as.formula(paste0("`", post_col, "` ~ group")), data = post_dat)
  write.csv(
    extract_simple_anova(anova_fit, outcome_name),
    file.path(output_dir, paste0(analysis_prefix, "_anova.csv")),
    row.names = FALSE
  )
  write.csv(
    extract_tukey(anova_fit, outcome_name),
    file.path(output_dir, paste0(analysis_prefix, "_tukey.csv")),
    row.names = FALSE
  )

  if (!is.null(pre_col)) {
    ancova_dat <- subset(data, !is.na(data[[post_col]]) & !is.na(data[[pre_col]]) & !is.na(group))
    ancova_fit <- lm(as.formula(paste0("`", post_col, "` ~ `", pre_col, "` + group")), data = ancova_dat)
    ancova_table <- as.data.frame(anova(ancova_fit))
    ancova_table$term <- rownames(ancova_table)
    rownames(ancova_table) <- NULL
    names(ancova_table) <- make.names(names(ancova_table))
    ancova_table$outcome <- outcome_name
    write.csv(
      ancova_table,
      file.path(output_dir, paste0(analysis_prefix, "_ancova.csv")),
      row.names = FALSE
    )
    ancova_pairwise <- extract_ancova_pairwise(ancova_dat, post_col, pre_col, outcome_name)
    if (nrow(ancova_pairwise) > 0) {
      write.csv(
        ancova_pairwise,
        file.path(output_dir, paste0(analysis_prefix, "_ancova_pairwise.csv")),
        row.names = FALSE
      )
    }
  }
}

source_dir <- file.path(project_root, "output", "analysis")
rct <- read.csv(file.path(source_dir, "rct_scale_totals_long.csv"), stringsAsFactors = FALSE, check.names = FALSE)
behavior <- read.csv(file.path(source_dir, "behavior_clean.csv"), stringsAsFactors = FALSE, check.names = FALSE)
practice <- read.csv(file.path(source_dir, "practice_student_summary_clean.csv"), stringsAsFactors = FALSE, check.names = FALSE)
practice_session <- read.csv(file.path(source_dir, "practice_session_level_clean.csv"), stringsAsFactors = FALSE, check.names = FALSE)
participant_master <- read.csv(file.path(source_dir, "participant_master.csv"), stringsAsFactors = FALSE, check.names = FALSE)
sav <- read.spss(
  file.path(project_root, "Psyclient_data.sav"),
  to.data.frame = TRUE,
  use.value.labels = FALSE
)
names(sav) <- trimws(names(sav))

rct$full_id <- normalize_id(rct$full_id)
rct$group <- group_factor(rct$group_label)
rct$time_c <- rct$timepoint_index - 1
rct$decline_time <- pmin(rct$timepoint_index - 1, 3)
rct$recovery_time <- pmax(rct$timepoint_index - 4, 0)
rct$front_half <- pmin(rct$timepoint_index - 1, 1)
rct$back_half <- pmax(rct$timepoint_index - 2, 0)
rct$timepoint <- factor(rct$timepoint, levels = c("T1", "T2", "T3", "T4", "T5", "T6"))

behavior$full_id <- normalize_id(behavior$full_id)
behavior$group <- group_factor(behavior$group_label)
behavior <- subset(behavior, roster_match_status == "exact_full_id_match")

sav$full_id <- normalize_id(sav$ID)
sav <- merge(
  sav,
  unique(rct[, c("full_id", "group_label")]),
  by = "full_id",
  all.x = TRUE
)
sav$group <- group_factor(sav$group_label)
sav$question_ratio_pre <- sav$questioning_pre
sav$question_ratio_post <- sav$questioning_post

practice$matched_full_id <- normalize_id(practice$matched_full_id)
practice$practice_student_id <- normalize_id(practice$practice_student_id)
practice <- subset(practice, !(practice_student_id %in% c("TA01", "TA02")))
# 021 is ambiguous by short ID, but the record contains many feedback clicks,
# which is only plausible for the feedback arm.
practice$matched_full_id[practice$practice_student_id == "021"] <- "<STUDENT_ID_021>"
practice$matched_group_label[practice$practice_student_id == "021"] <- "实验组"
practice$roster_match_status[practice$practice_student_id == "021"] <- "manual_feedback_resolution"
practice$matched_group <- group_factor(practice$matched_group_label)

practice_session$matched_full_id <- normalize_id(practice_session$matched_full_id)
practice_session$practice_student_id <- normalize_id(practice_session$practice_student_id)
practice_session <- subset(practice_session, !(practice_student_id %in% c("TA01", "TA02")))
practice_session$matched_full_id[practice_session$practice_student_id == "021"] <- "<STUDENT_ID_021>"
practice_session$matched_group_label[practice_session$practice_student_id == "021"] <- "实验组"
practice_session$roster_match_status[practice_session$practice_student_id == "021"] <- "manual_feedback_resolution"
practice_session <- subset(
  practice_session,
  roster_match_status %in% c("unique_short_id_match", "manual_feedback_resolution")
)
practice_session$matched_group <- group_factor(practice_session$matched_group_label)
practice_session$date <- as.Date(practice_session$date)

participant_master$full_id <- normalize_id(participant_master$full_id)
participant_master$group <- group_factor(participant_master$group_label)

practice$ai_emotional_support <- practice$strategy__情感反映 + practice$strategy__认可与安慰
practice$ai_insight <- practice$strategy__解释 + practice$strategy__质疑 + practice$strategy__邀请采取新视角
practice$ai_content_exploration <- (
  practice$strategy__探询主观信息 +
    practice$strategy__探询客观信息 +
    practice$strategy__重述 +
    practice$strategy__最小鼓励
)
practice$ai_action <- practice$strategy__回答 + practice$strategy__邀请探索新行动

practice_session$ai_emotional_support <- (
  practice_session$strategy__情感反映 + practice_session$strategy__认可与安慰
)
practice_session$ai_insight <- (
  practice_session$strategy__解释 + practice_session$strategy__质疑 + practice_session$strategy__邀请采取新视角
)
practice_session$ai_content_exploration <- (
  practice_session$strategy__探询主观信息 +
    practice_session$strategy__探询客观信息 +
    practice_session$strategy__重述 +
    practice_session$strategy__最小鼓励
)
practice_session$ai_action <- (
  practice_session$strategy__回答 + practice_session$strategy__邀请探索新行动
)

write.csv(
  aggregate(
    cbind(hs_total, sm_total, cc_total) ~ group + timepoint + timepoint_index,
    data = rct,
    FUN = make_mean_sd_n
  ),
  file.path(output_dir, "rct_descriptives.csv"),
  row.names = FALSE
)

growth_outcomes <- c("hs_total", "sm_total", "cc_total")
growth_fixed <- list()
growth_anova <- list()
growth_errors <- data.frame()

for (outcome in growth_outcomes) {
  dat <- subset(rct, !is.na(rct[[outcome]]) & !is.na(group) & !is.na(full_id))
  model <- safe_lme(as.formula(paste(outcome, "~ time_c * group")), dat)
  if (inherits(model, "error")) {
    growth_errors <- rbind(
      growth_errors,
      data.frame(analysis = "growth", outcome = outcome, error = conditionMessage(model))
    )
  } else {
    growth_fixed[[outcome]] <- extract_fixed_effects(model, "growth", outcome)
    growth_anova[[outcome]] <- extract_anova_table(model, "growth", outcome)
  }
}

if (length(growth_fixed) > 0) {
  write.csv(do.call(rbind, growth_fixed), file.path(output_dir, "growth_fixed_effects.csv"), row.names = FALSE)
}
if (length(growth_anova) > 0) {
  write.csv(do.call(rbind, growth_anova), file.path(output_dir, "growth_anova.csv"), row.names = FALSE)
}
if (nrow(growth_errors) > 0) {
  write.csv(growth_errors, file.path(output_dir, "growth_errors.csv"), row.names = FALSE)
}

piecewise_fixed <- list()
piecewise_anova <- list()
piecewise_errors <- data.frame()

for (outcome in growth_outcomes) {
  dat <- subset(rct, !is.na(rct[[outcome]]) & !is.na(group) & !is.na(full_id))
  model <- safe_lme(
    as.formula(paste(outcome, "~ decline_time * group + recovery_time * group")),
    dat
  )
  if (inherits(model, "error")) {
    piecewise_errors <- rbind(
      piecewise_errors,
      data.frame(analysis = "growth_piecewise", outcome = outcome, error = conditionMessage(model))
    )
  } else {
    piecewise_fixed[[outcome]] <- extract_fixed_effects(model, "growth_piecewise", outcome)
    piecewise_anova[[outcome]] <- extract_anova_table(model, "growth_piecewise", outcome)
  }
}

if (length(piecewise_fixed) > 0) {
  write.csv(
    do.call(rbind, piecewise_fixed),
    file.path(output_dir, "growth_piecewise_fixed_effects.csv"),
    row.names = FALSE
  )
}
if (length(piecewise_anova) > 0) {
  write.csv(
    do.call(rbind, piecewise_anova),
    file.path(output_dir, "growth_piecewise_anova.csv"),
    row.names = FALSE
  )
}
if (nrow(piecewise_errors) > 0) {
  write.csv(
    piecewise_errors,
    file.path(output_dir, "growth_piecewise_errors.csv"),
    row.names = FALSE
  )
}

front_back_fixed <- list()
front_back_anova <- list()
front_back_errors <- data.frame()

for (outcome in growth_outcomes) {
  dat <- subset(rct, !is.na(rct[[outcome]]) & !is.na(group) & !is.na(full_id))
  model <- safe_lme(
    as.formula(paste(outcome, "~ front_half * group + back_half * group")),
    dat
  )
  if (inherits(model, "error")) {
    front_back_errors <- rbind(
      front_back_errors,
      data.frame(analysis = "growth_front_back", outcome = outcome, error = conditionMessage(model))
    )
  } else {
    front_back_fixed[[outcome]] <- extract_fixed_effects(model, "growth_front_back", outcome)
    front_back_anova[[outcome]] <- extract_anova_table(model, "growth_front_back", outcome)
  }
}

if (length(front_back_fixed) > 0) {
  write.csv(
    do.call(rbind, front_back_fixed),
    file.path(output_dir, "growth_front_back_fixed_effects.csv"),
    row.names = FALSE
  )
}
if (length(front_back_anova) > 0) {
  write.csv(
    do.call(rbind, front_back_anova),
    file.path(output_dir, "growth_front_back_anova.csv"),
    row.names = FALSE
  )
}
if (nrow(front_back_errors) > 0) {
  write.csv(
    front_back_errors,
    file.path(output_dir, "growth_front_back_errors.csv"),
    row.names = FALSE
  )
}

t6 <- subset(rct, timepoint == "T6")
t1 <- subset(rct, timepoint == "T1")[, c("full_id", "hs_total", "sm_total", "cc_total")]
names(t1) <- c("full_id", "hs_total_t1", "sm_total_t1", "cc_total_t1")
t6_for_ancova <- merge(t6, t1, by = "full_id", all.x = TRUE)

t6_anova_rows <- list()
t6_tukey_rows <- list()
t6_ancova_rows <- list()

for (outcome in growth_outcomes) {
  dat <- subset(t6, !is.na(t6[[outcome]]) & !is.na(group))
  anova_fit <- aov(as.formula(paste(outcome, "~ group")), data = dat)
  anova_table <- as.data.frame(summary(anova_fit)[[1]])
  anova_table$term <- rownames(anova_table)
  rownames(anova_table) <- NULL
  names(anova_table) <- make.names(names(anova_table))
  anova_table$outcome <- outcome
  t6_anova_rows[[outcome]] <- anova_table

  tukey_table <- as.data.frame(TukeyHSD(anova_fit)$group)
  tukey_table$comparison <- rownames(tukey_table)
  rownames(tukey_table) <- NULL
  names(tukey_table) <- c("diff", "lwr", "upr", "p_adj", "comparison")
  tukey_table$outcome <- outcome
  t6_tukey_rows[[outcome]] <- tukey_table

  base_col <- paste0(outcome, "_t1")
  ancova_dat <- subset(
    t6_for_ancova,
    !is.na(t6_for_ancova[[outcome]]) & !is.na(t6_for_ancova[[base_col]]) & !is.na(group)
  )
  ancova_fit <- lm(as.formula(paste(outcome, "~", base_col, "+ group")), data = ancova_dat)
  ancova_table <- as.data.frame(anova(ancova_fit))
  ancova_table$term <- rownames(ancova_table)
  rownames(ancova_table) <- NULL
  names(ancova_table) <- make.names(names(ancova_table))
  ancova_table$outcome <- outcome
  t6_ancova_rows[[outcome]] <- ancova_table
}

write.csv(do.call(rbind, t6_anova_rows), file.path(output_dir, "t6_group_anova.csv"), row.names = FALSE)
write.csv(do.call(rbind, t6_tukey_rows), file.path(output_dir, "t6_group_tukey.csv"), row.names = FALSE)
write.csv(do.call(rbind, t6_ancova_rows), file.path(output_dir, "t6_group_ancova.csv"), row.names = FALSE)

make_behavior_long <- function(metric_type) {
  composites <- c("emotional_support", "insight", "content_exploration", "action")
  pieces <- lapply(composites, function(comp) {
    pre_col <- paste0(comp, "_", metric_type, "_pre")
    post_col <- paste0(comp, "_", metric_type, "_post")
    pre_df <- data.frame(
      full_id = behavior$full_id,
      group = behavior$group,
      composite = comp,
      phase = "pre",
      score = behavior[[pre_col]]
    )
    post_df <- data.frame(
      full_id = behavior$full_id,
      group = behavior$group,
      composite = comp,
      phase = "post",
      score = behavior[[post_col]]
    )
    rbind(pre_df, post_df)
  })
  out <- do.call(rbind, pieces)
  out$phase <- factor(out$phase, levels = c("pre", "post"))
  out
}

behavior_raw_long <- make_behavior_long("raw")
behavior_norm_long <- make_behavior_long("norm")

write.csv(
  aggregate(score ~ group + composite + phase, data = behavior_raw_long, FUN = make_mean_sd_n),
  file.path(output_dir, "behavior_raw_descriptives.csv"),
  row.names = FALSE
)
write.csv(
  aggregate(score ~ group + composite + phase, data = behavior_norm_long, FUN = make_mean_sd_n),
  file.path(output_dir, "behavior_norm_descriptives.csv"),
  row.names = FALSE
)

run_behavior_models <- function(metric_type, long_df) {
  composites <- unique(long_df$composite)
  fixed_tables <- list()
  anova_tables <- list()
  errors <- data.frame()

  for (comp in composites) {
    dat <- subset(long_df, composite == comp & !is.na(score) & !is.na(group) & !is.na(full_id))
    model <- safe_lme(score ~ phase * group, dat)
    if (inherits(model, "error")) {
      errors <- rbind(
        errors,
        data.frame(analysis = paste0("behavior_", metric_type), outcome = comp, error = conditionMessage(model))
      )
    } else {
      fixed_tables[[comp]] <- extract_fixed_effects(model, paste0("behavior_", metric_type), comp)
      anova_tables[[comp]] <- extract_anova_table(model, paste0("behavior_", metric_type), comp)
    }
  }

  if (length(fixed_tables) > 0) {
    write.csv(
      do.call(rbind, fixed_tables),
      file.path(output_dir, paste0("behavior_", metric_type, "_fixed_effects.csv")),
      row.names = FALSE
    )
  }
  if (length(anova_tables) > 0) {
    write.csv(
      do.call(rbind, anova_tables),
      file.path(output_dir, paste0("behavior_", metric_type, "_anova.csv")),
      row.names = FALSE
    )
  }
  if (nrow(errors) > 0) {
    write.csv(
      errors,
      file.path(output_dir, paste0("behavior_", metric_type, "_errors.csv")),
      row.names = FALSE
    )
  }
}

run_behavior_models("raw", behavior_raw_long)
run_behavior_models("norm", behavior_norm_long)

run_behavior_post_models <- function(metric_type) {
  composites <- c("emotional_support", "insight", "content_exploration", "action")
  anova_rows <- list()
  tukey_rows <- list()
  ancova_rows <- list()

  for (comp in composites) {
    post_col <- paste0(comp, "_", metric_type, "_post")
    pre_col <- paste0(comp, "_", metric_type, "_pre")
    dat <- subset(behavior, !is.na(behavior[[post_col]]) & !is.na(group))

    post_fit <- aov(as.formula(paste(post_col, "~ group")), data = dat)
    post_anova <- as.data.frame(summary(post_fit)[[1]])
    post_anova$term <- rownames(post_anova)
    rownames(post_anova) <- NULL
    names(post_anova) <- make.names(names(post_anova))
    post_anova$composite <- comp
    anova_rows[[comp]] <- post_anova

    post_tukey <- as.data.frame(TukeyHSD(post_fit)$group)
    post_tukey$comparison <- rownames(post_tukey)
    rownames(post_tukey) <- NULL
    names(post_tukey) <- c("diff", "lwr", "upr", "p_adj", "comparison")
    post_tukey$composite <- comp
    tukey_rows[[comp]] <- post_tukey

    ancova_dat <- subset(
      behavior,
      !is.na(behavior[[post_col]]) & !is.na(behavior[[pre_col]]) & !is.na(group)
    )
    ancova_fit <- lm(as.formula(paste(post_col, "~", pre_col, "+ group")), data = ancova_dat)
    ancova_table <- as.data.frame(anova(ancova_fit))
    ancova_table$term <- rownames(ancova_table)
    rownames(ancova_table) <- NULL
    names(ancova_table) <- make.names(names(ancova_table))
    ancova_table$composite <- comp
    ancova_rows[[comp]] <- ancova_table
  }

  write.csv(
    do.call(rbind, anova_rows),
    file.path(output_dir, paste0("behavior_", metric_type, "_post_anova.csv")),
    row.names = FALSE
  )
  write.csv(
    do.call(rbind, tukey_rows),
    file.path(output_dir, paste0("behavior_", metric_type, "_post_tukey.csv")),
    row.names = FALSE
  )
  write.csv(
    do.call(rbind, ancova_rows),
    file.path(output_dir, paste0("behavior_", metric_type, "_post_ancova.csv")),
    row.names = FALSE
  )
}

run_behavior_post_models("raw")
run_behavior_post_models("norm")

run_single_post_models(
  sav,
  post_col = "question_ratio_post",
  pre_col = "question_ratio_pre",
  analysis_prefix = "question_ratio_post",
  outcome_name = "question_ratio_post"
)

run_single_post_models(
  sav,
  post_col = "RQR_new_post",
  pre_col = "RQR_new_PRE",
  analysis_prefix = "rqr_new_post",
  outcome_name = "RQR_new_post"
)

t6_confidence_specs <- list(
  list(post = "HST6", pre = "HST1", prefix = "t6_hs_itemmean", outcome = "HST6"),
  list(post = "SMT6", pre = "SMT1", prefix = "t6_sm_itemmean", outcome = "SMT6"),
  list(post = "CCT6", pre = "CCT1", prefix = "t6_cc_itemmean", outcome = "CCT6")
)

for (spec in t6_confidence_specs) {
  run_single_post_models(
    sav,
    post_col = spec$post,
    pre_col = spec$pre,
    analysis_prefix = spec$prefix,
    outcome_name = spec$outcome
  )
}

practice_prediction_master <- merge(
  practice,
  participant_master[
    ,
    c(
      "full_id",
      "hs_total_t1",
      "sm_total_t1",
      "cc_total_t1",
      "self_efficacy_total_t1",
      "hs_total_t6",
      "sm_total_t6",
      "cc_total_t6",
      "self_efficacy_total_t6"
    )
  ],
  by.x = "matched_full_id",
  by.y = "full_id",
  all.x = TRUE
)

run_practice_prediction_models <- function(data, sample_name, predictors, outcomes, include_group = FALSE) {
  rows <- list()

  for (outcome_name in names(outcomes)) {
    post_col <- outcomes[[outcome_name]][["post"]]
    pre_col <- outcomes[[outcome_name]][["pre"]]
    for (predictor in predictors) {
      cols <- c(post_col, pre_col, predictor)
      if (include_group) {
        cols <- c(cols, "matched_group")
      }
      dat <- data[, cols]
      dat <- dat[complete.cases(dat), ]
      if (nrow(dat) < 12) {
        next
      }
      rhs <- c(paste0("`", pre_col, "`"), paste0("`", predictor, "`"))
      if (include_group) {
        rhs <- c(rhs, "matched_group")
      }
      formula_text <- paste0("`", post_col, "` ~ ", paste(rhs, collapse = " + "))
      fit <- lm(as.formula(formula_text), data = dat)
      summary_fit <- summary(fit)
      rows[[paste(outcome_name, predictor, sep = "__")]] <- data.frame(
        sample = sample_name,
        outcome = outcome_name,
        post_col = post_col,
        pre_col = pre_col,
        predictor = predictor,
        include_group = include_group,
        n = nrow(dat),
        estimate = coef(summary_fit)[predictor, "Estimate"],
        std_error = coef(summary_fit)[predictor, "Std. Error"],
        t_value = coef(summary_fit)[predictor, "t value"],
        p_value = coef(summary_fit)[predictor, "Pr(>|t|)"],
        r_squared = summary_fit$r.squared,
        adj_r_squared = summary_fit$adj.r.squared
      )
    }
  }

  out <- do.call(rbind, rows)
  if (!is.null(out) && nrow(out) > 0) {
    out$p_value_bh <- p.adjust(out$p_value, method = "BH")
    out <- out[order(out$p_value), ]
  }
  out
}

practice_prediction_outcomes <- list(
  hs_total_t6 = list(post = "hs_total_t6", pre = "hs_total_t1"),
  sm_total_t6 = list(post = "sm_total_t6", pre = "sm_total_t1"),
  cc_total_t6 = list(post = "cc_total_t6", pre = "cc_total_t1"),
  self_efficacy_total_t6 = list(post = "self_efficacy_total_t6", pre = "self_efficacy_total_t1")
)

practice_ai_groups <- subset(practice_prediction_master, matched_group %in% c("no_feedback", "feedback"))
practice_feedback_only <- subset(practice_prediction_master, matched_group == "feedback")

practice_ai_prediction_table <- run_practice_prediction_models(
  practice_ai_groups,
  sample_name = "ai_groups",
  predictors = c("session_count", "feedback_count", "overall_mean_score_weighted"),
  outcomes = practice_prediction_outcomes,
  include_group = TRUE
)
practice_feedback_prediction_table <- run_practice_prediction_models(
  practice_feedback_only,
  sample_name = "feedback_only",
  predictors = c("feedback_count", "overall_mean_score_weighted", "total_score"),
  outcomes = practice_prediction_outcomes,
  include_group = FALSE
)

if (!is.null(practice_ai_prediction_table) && nrow(practice_ai_prediction_table) > 0) {
  write.csv(
    practice_ai_prediction_table,
    file.path(output_dir, "practice_t6_prediction_ai_groups.csv"),
    row.names = FALSE
  )
}
if (!is.null(practice_feedback_prediction_table) && nrow(practice_feedback_prediction_table) > 0) {
  write.csv(
    practice_feedback_prediction_table,
    file.path(output_dir, "practice_t6_prediction_feedback_only.csv"),
    row.names = FALSE
  )
}

practice_session <- practice_session[order(practice_session$matched_full_id, practice_session$date), ]
for (var_name in c("mean_score", "total_score", "feedback_count", "ai_insight")) {
  practice_session[[paste0("next_", var_name)]] <- ave(
    practice_session[[var_name]],
    practice_session$matched_full_id,
    FUN = function(x) c(x[-1], NA)
  )
  practice_session[[paste0("lag_", var_name)]] <- ave(
    practice_session[[var_name]],
    practice_session$matched_full_id,
    FUN = function(x) c(NA, x[-length(x)])
  )
}

run_lagged_session_models <- function(data, specs) {
  rows <- list()
  anova_rows <- list()

  for (spec in specs) {
    cols <- c("matched_full_id", spec$outcome, spec$predictor, spec$baseline)
    dat <- data[, cols]
    dat <- dat[complete.cases(dat), ]
    if (nrow(dat) < 30) {
      next
    }
    names(dat)[1] <- "full_id"
    fit <- lme(
      as.formula(
        paste0("`", spec$outcome, "` ~ `", spec$predictor, "` + `", spec$baseline, "`")
      ),
      random = ~1 | full_id,
      data = dat,
      na.action = na.omit,
      method = "REML"
    )
    t_table <- as.data.frame(summary(fit)$tTable)
    t_table$term <- rownames(t_table)
    rownames(t_table) <- NULL
    pred_row <- subset(t_table, term == spec$predictor)
    rows[[spec$outcome]] <- data.frame(
      outcome = spec$outcome,
      predictor = spec$predictor,
      baseline = spec$baseline,
      n = nrow(dat),
      n_students = length(unique(dat$full_id)),
      estimate = pred_row$Value,
      std_error = pred_row$Std.Error,
      df = pred_row$DF,
      t_value = pred_row$`t-value`,
      p_value = pred_row$`p-value`
    )

    a_table <- as.data.frame(anova(fit))
    a_table$term <- rownames(a_table)
    rownames(a_table) <- NULL
    names(a_table) <- make.names(names(a_table))
    a_table$outcome <- spec$outcome
    anova_rows[[spec$outcome]] <- a_table
  }

  list(
    coefficients = do.call(rbind, rows),
    anova = do.call(rbind, anova_rows)
  )
}

lagged_specs <- list(
  list(outcome = "next_mean_score", predictor = "mean_score", baseline = "lag_mean_score"),
  list(outcome = "next_total_score", predictor = "total_score", baseline = "lag_total_score"),
  list(outcome = "next_feedback_count", predictor = "feedback_count", baseline = "lag_feedback_count"),
  list(outcome = "next_ai_insight", predictor = "ai_insight", baseline = "lag_ai_insight")
)

lagged_results <- run_lagged_session_models(practice_ai_groups <- subset(practice_session, matched_group %in% c("no_feedback", "feedback")), lagged_specs)
if (!is.null(lagged_results$coefficients) && nrow(lagged_results$coefficients) > 0) {
  write.csv(
    lagged_results$coefficients[order(lagged_results$coefficients$p_value), ],
    file.path(output_dir, "practice_next_session_coefficients.csv"),
    row.names = FALSE
  )
}
if (!is.null(lagged_results$anova) && nrow(lagged_results$anova) > 0) {
  write.csv(
    lagged_results$anova,
    file.path(output_dir, "practice_next_session_anova.csv"),
    row.names = FALSE
  )
}

feedback_practice <- subset(
  practice,
  roster_match_status %in% c("unique_short_id_match", "manual_feedback_resolution") &
    matched_group == "feedback"
)

feedback_behavior <- subset(behavior, group == "feedback")

feedback_df <- merge(
  feedback_behavior,
  feedback_practice,
  by.x = "full_id",
  by.y = "matched_full_id",
  all = FALSE
)

feedback_df$ai_emotional_support <- as.numeric(feedback_df$ai_emotional_support)
feedback_df$ai_insight <- as.numeric(feedback_df$ai_insight)
feedback_df$ai_content_exploration <- as.numeric(feedback_df$ai_content_exploration)
feedback_df$ai_action <- as.numeric(feedback_df$ai_action)
feedback_df$ai_total <- (
  feedback_df$ai_emotional_support +
    feedback_df$ai_insight +
    feedback_df$ai_content_exploration +
    feedback_df$ai_action +
    feedback_df$strategy__自我暴露 +
    feedback_df$strategy__其他
)
feedback_df$ai_emotional_support_prop <- ifelse(
  feedback_df$ai_total > 0,
  feedback_df$ai_emotional_support / feedback_df$ai_total,
  NA_real_
)
feedback_df$ai_insight_prop <- ifelse(
  feedback_df$ai_total > 0,
  feedback_df$ai_insight / feedback_df$ai_total,
  NA_real_
)
feedback_df$ai_content_exploration_prop <- ifelse(
  feedback_df$ai_total > 0,
  feedback_df$ai_content_exploration / feedback_df$ai_total,
  NA_real_
)
feedback_df$ai_action_prop <- ifelse(
  feedback_df$ai_total > 0,
  feedback_df$ai_action / feedback_df$ai_total,
  NA_real_
)

write.csv(feedback_df, file.path(output_dir, "feedback_group_prediction_dataset.csv"), row.names = FALSE)

composite_prediction_specs <- list(
  list(outcome = "emotional_support", post = "emotional_support_raw_post", pre = "emotional_support_raw_pre", predictor = "ai_emotional_support", model_spec = "raw_count_to_raw_post"),
  list(outcome = "emotional_support", post = "emotional_support_norm_post", pre = "emotional_support_norm_pre", predictor = "ai_emotional_support", model_spec = "raw_count_to_norm_post"),
  list(outcome = "emotional_support", post = "emotional_support_raw_post", pre = "emotional_support_raw_pre", predictor = "ai_emotional_support_prop", model_spec = "proportion_to_raw_post"),
  list(outcome = "emotional_support", post = "emotional_support_norm_post", pre = "emotional_support_norm_pre", predictor = "ai_emotional_support_prop", model_spec = "proportion_to_norm_post"),
  list(outcome = "insight", post = "insight_raw_post", pre = "insight_raw_pre", predictor = "ai_insight", model_spec = "raw_count_to_raw_post"),
  list(outcome = "insight", post = "insight_norm_post", pre = "insight_norm_pre", predictor = "ai_insight", model_spec = "raw_count_to_norm_post"),
  list(outcome = "insight", post = "insight_raw_post", pre = "insight_raw_pre", predictor = "ai_insight_prop", model_spec = "proportion_to_raw_post"),
  list(outcome = "insight", post = "insight_norm_post", pre = "insight_norm_pre", predictor = "ai_insight_prop", model_spec = "proportion_to_norm_post"),
  list(outcome = "content_exploration", post = "content_exploration_raw_post", pre = "content_exploration_raw_pre", predictor = "ai_content_exploration", model_spec = "raw_count_to_raw_post"),
  list(outcome = "content_exploration", post = "content_exploration_norm_post", pre = "content_exploration_norm_pre", predictor = "ai_content_exploration", model_spec = "raw_count_to_norm_post"),
  list(outcome = "content_exploration", post = "content_exploration_raw_post", pre = "content_exploration_raw_pre", predictor = "ai_content_exploration_prop", model_spec = "proportion_to_raw_post"),
  list(outcome = "content_exploration", post = "content_exploration_norm_post", pre = "content_exploration_norm_pre", predictor = "ai_content_exploration_prop", model_spec = "proportion_to_norm_post"),
  list(outcome = "action", post = "action_raw_post", pre = "action_raw_pre", predictor = "ai_action", model_spec = "raw_count_to_raw_post"),
  list(outcome = "action", post = "action_norm_post", pre = "action_norm_pre", predictor = "ai_action", model_spec = "raw_count_to_norm_post"),
  list(outcome = "action", post = "action_raw_post", pre = "action_raw_pre", predictor = "ai_action_prop", model_spec = "proportion_to_raw_post"),
  list(outcome = "action", post = "action_norm_post", pre = "action_norm_pre", predictor = "ai_action_prop", model_spec = "proportion_to_norm_post")
)

composite_prediction_rows <- list()
composite_prediction_errors <- data.frame()

for (spec in composite_prediction_specs) {
  dat <- feedback_df[, c("full_id", spec$post, spec$pre, spec$predictor)]
  names(dat) <- c("full_id", "post_score", "pre_score", "ai_feedback_score")
  dat <- dat[complete.cases(dat), ]
  model <- safe_lm(post_score ~ pre_score + ai_feedback_score, dat)

  if (inherits(model, "error")) {
    composite_prediction_errors <- rbind(
      composite_prediction_errors,
      data.frame(
        analysis = "feedback_prediction_composite",
        outcome = spec$outcome,
        model_spec = spec$model_spec,
        predictor = spec$predictor,
        error = conditionMessage(model)
      )
    )
  } else {
    summary_model <- summary(model)
    composite_prediction_rows[[paste(spec$outcome, spec$model_spec, sep = "__")]] <- data.frame(
      outcome = spec$outcome,
      model_spec = spec$model_spec,
      predictor = spec$predictor,
      estimate = coef(summary_model)["ai_feedback_score", "Estimate"],
      std_error = coef(summary_model)["ai_feedback_score", "Std. Error"],
      t_value = coef(summary_model)["ai_feedback_score", "t value"],
      p_value = coef(summary_model)["ai_feedback_score", "Pr(>|t|)"],
      n = nrow(dat),
      r_squared = summary_model$r.squared,
      adj_r_squared = summary_model$adj.r.squared
    )
  }
}

composite_prediction_table <- do.call(rbind, composite_prediction_rows)
if (!is.null(composite_prediction_table) && nrow(composite_prediction_table) > 0) {
  composite_prediction_table$p_value_bh <- p.adjust(composite_prediction_table$p_value, method = "BH")
  composite_prediction_table <- composite_prediction_table[order(composite_prediction_table$p_value), ]
  write.csv(
    composite_prediction_table,
    file.path(output_dir, "feedback_prediction_composite_models.csv"),
    row.names = FALSE
  )
}
if (nrow(composite_prediction_errors) > 0) {
  write.csv(
    composite_prediction_errors,
    file.path(output_dir, "feedback_prediction_composite_errors.csv"),
    row.names = FALSE
  )
}

technique_map <- list(
  "strategy__情感反映" = "Reflection of Feelings",
  "strategy__认可与安慰" = "Affirmation and Reassurance",
  "strategy__解释" = "Interpretation",
  "strategy__质疑" = "Confrontation",
  "strategy__邀请采取新视角" = "Invite to Take New Perspectives",
  "strategy__探询主观信息" = "Inquiring Subjective Information",
  "strategy__探询客观信息" = "Inquiring Objective Information",
  "strategy__重述" = "Restatement",
  "strategy__最小鼓励" = "Minimal Encouragement",
  "strategy__回答" = "Answer",
  "strategy__邀请探索新行动" = "Invite to Explore New Actions",
  "strategy__自我暴露" = "Self-disclosure"
)

technique_prediction_rows <- list()

for (strategy_name in names(technique_map)) {
  behavior_base <- technique_map[[strategy_name]]
  pre_raw <- paste0(behavior_base, "_raw_pre")
  post_raw <- paste0(behavior_base, "_raw_post")
  pre_norm <- paste0(behavior_base, "_norm_pre")
  post_norm <- paste0(behavior_base, "_norm_post")
  prop_name <- paste0("prop__", strategy_name)
  feedback_df[[prop_name]] <- ifelse(
    feedback_df$ai_total > 0,
    feedback_df[[strategy_name]] / feedback_df$ai_total,
    NA_real_
  )

  technique_specs <- list(
    list(model_spec = "raw_to_post_raw", post = post_raw, pre = pre_raw, predictor = strategy_name),
    list(model_spec = "proportion_to_post_raw", post = post_raw, pre = pre_raw, predictor = prop_name),
    list(model_spec = "raw_to_post_norm", post = post_norm, pre = pre_norm, predictor = strategy_name),
    list(model_spec = "proportion_to_post_norm", post = post_norm, pre = pre_norm, predictor = prop_name)
  )

  for (spec in technique_specs) {
    dat <- feedback_df[, c("full_id", spec$post, spec$pre, spec$predictor)]
    dat <- dat[complete.cases(dat), ]
    if (nrow(dat) < 10) {
      next
    }
    model <- safe_lm(
      as.formula(
        paste0("`", spec$post, "` ~ `", spec$pre, "` + `", spec$predictor, "`")
      ),
      dat
    )
    if (inherits(model, "error")) {
      next
    }
    summary_model <- summary(model)
    technique_prediction_rows[[paste(strategy_name, spec$model_spec, sep = "__")]] <- data.frame(
      strategy = strategy_name,
      behavior_base = behavior_base,
      model_spec = spec$model_spec,
      predictor = spec$predictor,
      estimate = coef(summary_model)[spec$predictor, "Estimate"],
      std_error = coef(summary_model)[spec$predictor, "Std. Error"],
      t_value = coef(summary_model)[spec$predictor, "t value"],
      p_value = coef(summary_model)[spec$predictor, "Pr(>|t|)"],
      n = nrow(dat),
      r_squared = summary_model$r.squared,
      adj_r_squared = summary_model$adj.r.squared
    )
  }
}

technique_prediction_table <- do.call(rbind, technique_prediction_rows)
if (!is.null(technique_prediction_table) && nrow(technique_prediction_table) > 0) {
  technique_prediction_table$p_value_bh <- p.adjust(technique_prediction_table$p_value, method = "BH")
  technique_prediction_table <- technique_prediction_table[order(technique_prediction_table$p_value), ]
  write.csv(
    technique_prediction_table,
    file.path(output_dir, "feedback_prediction_technique_models.csv"),
    row.names = FALSE
  )
}

feedback_prediction_coefficients <- subset(
  composite_prediction_table,
  model_spec == "raw_count_to_raw_post"
)
feedback_prediction_model_summary <- aggregate(
  cbind(r_squared, adj_r_squared, n) ~ outcome,
  data = feedback_prediction_coefficients,
  FUN = function(x) x[1]
)
write.csv(
  feedback_prediction_coefficients,
  file.path(output_dir, "feedback_prediction_coefficients.csv"),
  row.names = FALSE
)
write.csv(
  feedback_prediction_model_summary,
  file.path(output_dir, "feedback_prediction_model_summary.csv"),
  row.names = FALSE
)

summary_lines <- c(
  "# Primary Analysis Outputs",
  "",
  "## Question 1: HS / SM / CC trajectories",
  "- growth_fixed_effects.csv: lme coefficient tables for hs_total, sm_total, cc_total",
  "- growth_anova.csv: omnibus tests for time, group, and time x group",
  "- growth_piecewise_fixed_effects.csv and growth_piecewise_anova.csv: segmented HLM with decline (T1-T4) and recovery (T4-T6) slopes",
  "- growth_front_back_fixed_effects.csv and growth_front_back_anova.csv: segmented HLM using front half (T1-T2) and back half (T2-T6)",
  "- t6_group_anova.csv, t6_group_tukey.csv, and t6_group_ancova.csv: T6-only group comparisons",
  "- rct_descriptives.csv: group x timepoint descriptive statistics",
  "",
  "## Question 2: behavior pre-post comparison",
  "- behavior_raw_fixed_effects.csv and behavior_raw_anova.csv: pre-post mixed models on raw frequencies",
  "- behavior_norm_fixed_effects.csv and behavior_norm_anova.csv: pre-post mixed models on normalized frequencies",
  "- behavior_raw_post_anova.csv and behavior_raw_post_ancova.csv: post-only raw group comparisons",
  "- behavior_norm_post_anova.csv and behavior_norm_post_ancova.csv: post-only normalized group comparisons",
  "- question_ratio_post_anova.csv and question_ratio_post_ancova.csv: post-only group comparison for the proportion of all techniques that are questions (open + closed)",
  "- rqr_new_post_anova.csv and rqr_new_post_ancova.csv: post-only group comparison for the closed/open question ratio variable from the SAV file",
  "- t6_*_itemmean_*.csv: T6 confidence comparisons for the HS / SM / CC item means, with and without T1 adjustment",
  "- behavior_raw_descriptives.csv and behavior_norm_descriptives.csv: descriptives by group, phase, and composite",
  "",
  "## Question 3: feedback technique predicting post-test behavior",
  "- feedback_group_prediction_dataset.csv: merged feedback-group dataset used for regressions",
  "- feedback_prediction_composite_models.csv: composite-level prediction models using both raw counts and proportions",
  "- feedback_prediction_technique_models.csv: exploratory single-technique prediction models with BH-adjusted p-values",
  "- feedback_prediction_coefficients.csv: primary composite-level raw-count models",
  "- feedback_prediction_model_summary.csv: compact summary of the primary composite-level raw-count models",
  "- practice_t6_prediction_ai_groups.csv: exploratory regressions linking practice engagement to T6 confidence across the two AI-practice groups, adjusting for T1 and group",
  "- practice_t6_prediction_feedback_only.csv: exploratory regressions linking feedback-group usage intensity to T6 confidence, adjusting for T1",
  "- practice_next_session_coefficients.csv and practice_next_session_anova.csv: lagged multilevel models testing whether one practice predicts the next practice",
  "",
  "## Notes",
  "- TA01 and TA02 were excluded from the practice analyses.",
  "- ID 021 was manually resolved to the feedback group because the record contains many feedback clicks and cannot plausibly belong to the blank group.",
  "- Practice rows with missing roster matches were excluded from the feedback prediction models.",
  "- RQR_new_post was read from Psyclient_data.sav after trimming column-name whitespace.",
  "- question_ratio_post was set to questioning_post, treating questioning_post as the percent-open-questions variable stored in the SAV file.",
  "- T6 confidence files use Psyclient_data.sav so that ID <STUDENT_ID_047> follows the SAV version of T6.",
  "- Exploratory technique-level models are nominal tests; use BH-adjusted p-values for stricter interpretation."
)

writeLines(summary_lines, con = file.path(output_dir, "analysis_output_guide.md"))
