options(stringsAsFactors = FALSE)

project_root <- getwd()
analysis_dir <- file.path(project_root, "output", "analysis")
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

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "f"))
}

fmt_mean_sd <- function(mean, sd) {
  paste0(fmt_num(mean), " (", fmt_num(sd), ")")
}

agg_stats <- function(x) {
  c(n = sum(!is.na(x)), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
}

df_to_md <- function(df) {
  header <- paste(names(df), collapse = " | ")
  sep <- paste(rep("---", ncol(df)), collapse = " | ")
  rows <- apply(df, 1, function(row) paste(row, collapse = " | "))
  paste(c(paste0("| ", header, " |"), paste0("| ", sep, " |"), paste0("| ", rows, " |")), collapse = "\n")
}

roster <- read.csv(file.path(analysis_dir, "participant_roster.csv"), stringsAsFactors = FALSE, check.names = FALSE)
rct <- read.csv(file.path(analysis_dir, "rct_scale_totals_long.csv"), stringsAsFactors = FALSE, check.names = FALSE)
behavior <- read.csv(file.path(analysis_dir, "behavior_clean.csv"), stringsAsFactors = FALSE, check.names = FALSE)
practice <- read.csv(file.path(analysis_dir, "practice_student_summary_clean.csv"), stringsAsFactors = FALSE, check.names = FALSE)

roster$full_id <- normalize_id(roster$full_id)
roster$group <- group_factor(roster$group_label)

rct$full_id <- normalize_id(rct$full_id)
rct$group <- group_factor(rct$group_label)
rct$timepoint <- factor(rct$timepoint, levels = c("T1", "T2", "T3", "T4", "T5", "T6"))

behavior$full_id <- normalize_id(behavior$full_id)
behavior$group <- group_factor(behavior$group_label)
behavior <- subset(behavior, roster_match_status == "exact_full_id_match")

practice$matched_full_id <- normalize_id(practice$matched_full_id)
practice$practice_student_id <- normalize_id(practice$practice_student_id)
practice <- subset(practice, !(practice_student_id %in% c("TA01", "TA02")))
practice$matched_full_id[practice$practice_student_id == "021"] <- "<STUDENT_ID_021>"
practice$matched_group_label[practice$practice_student_id == "021"] <- "实验组"
practice$roster_match_status[practice$practice_student_id == "021"] <- "manual_feedback_resolution"
practice$group <- group_factor(practice$matched_group_label)
practice <- subset(practice, roster_match_status %in% c("unique_short_id_match", "manual_feedback_resolution"))

practice$ai_emotional_support <- practice$`strategy__情感反映` + practice$`strategy__认可与安慰`
practice$ai_insight <- practice$`strategy__解释` + practice$`strategy__质疑` + practice$`strategy__邀请采取新视角`
practice$ai_content_exploration <- (
  practice$`strategy__探询主观信息` +
    practice$`strategy__探询客观信息` +
    practice$`strategy__重述` +
    practice$`strategy__最小鼓励`
)
practice$ai_action <- practice$`strategy__回答` + practice$`strategy__邀请探索新行动`

sample_by_group <- aggregate(full_id ~ group, data = roster, FUN = function(x) length(unique(x)))
names(sample_by_group)[2] <- "n_participants"

rct_counts <- aggregate(full_id ~ group + timepoint, data = rct, FUN = function(x) length(unique(x)))
names(rct_counts)[3] <- "n_responded"

rct_desc <- aggregate(
  cbind(hs_mean, sm_mean, cc_mean, hs_total, sm_total, cc_total) ~ group + timepoint,
  data = rct,
  FUN = agg_stats
)
rct_table <- data.frame(
  group = rct_desc$group,
  timepoint = rct_desc$timepoint,
  hs_n = rct_desc$hs_mean[, "n"],
  hs_item_mean_sd = fmt_mean_sd(rct_desc$hs_mean[, "mean"], rct_desc$hs_mean[, "sd"]),
  sm_n = rct_desc$sm_mean[, "n"],
  sm_item_mean_sd = fmt_mean_sd(rct_desc$sm_mean[, "mean"], rct_desc$sm_mean[, "sd"]),
  cc_n = rct_desc$cc_mean[, "n"],
  cc_item_mean_sd = fmt_mean_sd(rct_desc$cc_mean[, "mean"], rct_desc$cc_mean[, "sd"])
)

rct_total_table <- data.frame(
  group = rct_desc$group,
  timepoint = rct_desc$timepoint,
  hs_n = rct_desc$hs_total[, "n"],
  hs_total_mean_sd = fmt_mean_sd(rct_desc$hs_total[, "mean"], rct_desc$hs_total[, "sd"]),
  sm_n = rct_desc$sm_total[, "n"],
  sm_total_mean_sd = fmt_mean_sd(rct_desc$sm_total[, "mean"], rct_desc$sm_total[, "sd"]),
  cc_n = rct_desc$cc_total[, "n"],
  cc_total_mean_sd = fmt_mean_sd(rct_desc$cc_total[, "mean"], rct_desc$cc_total[, "sd"])
)

make_behavior_long <- function(metric_type) {
  comps <- c("emotional_support", "insight", "content_exploration", "action")
  pieces <- lapply(comps, function(comp) {
    pre_col <- paste0(comp, "_", metric_type, "_pre")
    post_col <- paste0(comp, "_", metric_type, "_post")
    rbind(
      data.frame(group = behavior$group, composite = comp, phase = "pre", score = behavior[[pre_col]]),
      data.frame(group = behavior$group, composite = comp, phase = "post", score = behavior[[post_col]])
    )
  })
  do.call(rbind, pieces)
}

behavior_raw_long <- make_behavior_long("raw")
behavior_norm_long <- make_behavior_long("norm")

build_behavior_table <- function(df) {
  agg <- aggregate(score ~ group + composite + phase, data = df, FUN = agg_stats)
  data.frame(
    group = agg$group,
    composite = agg$composite,
    phase = agg$phase,
    n = agg$score[, "n"],
    mean_sd = fmt_mean_sd(agg$score[, "mean"], agg$score[, "sd"])
  )
}

behavior_raw_table <- build_behavior_table(behavior_raw_long)
behavior_norm_table <- build_behavior_table(behavior_norm_long)

practice_desc <- aggregate(
  cbind(
    session_count,
    counselor_count,
    feedback_count,
    total_score,
    overall_mean_score_weighted,
    ai_emotional_support,
    ai_insight,
    ai_content_exploration,
    ai_action
  ) ~ group,
  data = practice,
  FUN = agg_stats
)

practice_table <- data.frame(
  group = practice_desc$group,
  n_students = practice_desc$session_count[, "n"],
  sessions_mean_sd = fmt_mean_sd(practice_desc$session_count[, "mean"], practice_desc$session_count[, "sd"]),
  counselor_turns_mean_sd = fmt_mean_sd(practice_desc$counselor_count[, "mean"], practice_desc$counselor_count[, "sd"]),
  feedback_clicks_mean_sd = fmt_mean_sd(practice_desc$feedback_count[, "mean"], practice_desc$feedback_count[, "sd"]),
  copy_total_mean_sd = fmt_mean_sd(practice_desc$total_score[, "mean"], practice_desc$total_score[, "sd"]),
  copy_weighted_mean_sd = fmt_mean_sd(practice_desc$overall_mean_score_weighted[, "mean"], practice_desc$overall_mean_score_weighted[, "sd"]),
  ai_emotional_support_mean_sd = fmt_mean_sd(practice_desc$ai_emotional_support[, "mean"], practice_desc$ai_emotional_support[, "sd"]),
  ai_insight_mean_sd = fmt_mean_sd(practice_desc$ai_insight[, "mean"], practice_desc$ai_insight[, "sd"]),
  ai_content_exploration_mean_sd = fmt_mean_sd(practice_desc$ai_content_exploration[, "mean"], practice_desc$ai_content_exploration[, "sd"]),
  ai_action_mean_sd = fmt_mean_sd(practice_desc$ai_action[, "mean"], practice_desc$ai_action[, "sd"])
)

write.csv(sample_by_group, file.path(output_dir, "descriptive_sample_by_group.csv"), row.names = FALSE)
write.csv(rct_counts, file.path(output_dir, "descriptive_rct_counts.csv"), row.names = FALSE)
write.csv(rct_table, file.path(output_dir, "descriptive_rct_hs_sm_cc_item_means.csv"), row.names = FALSE)
write.csv(rct_total_table, file.path(output_dir, "descriptive_rct_hs_sm_cc_totals.csv"), row.names = FALSE)
write.csv(behavior_raw_table, file.path(output_dir, "descriptive_behavior_raw.csv"), row.names = FALSE)
write.csv(behavior_norm_table, file.path(output_dir, "descriptive_behavior_norm.csv"), row.names = FALSE)
write.csv(practice_table, file.path(output_dir, "descriptive_practice_by_group.csv"), row.names = FALSE)

md_lines <- c(
  "# Descriptive Summary",
  "",
  "## Sample Size",
  "",
  df_to_md(sample_by_group),
  "",
  "## RCT Response Counts by Group and Time",
  "",
  df_to_md(rct_counts),
  "",
  "## HS / SM / CC Descriptives (Item Means, 1-10 Scale)",
  "",
  df_to_md(rct_table),
  "",
  "## Behavior Norm Pre-Post Descriptives (Proportion of All Techniques)",
  "",
  df_to_md(behavior_norm_table),
  "",
  "## Behavior Raw Pre-Post Descriptives (Reference Only)",
  "",
  df_to_md(behavior_raw_table),
  "",
  "## Practice Descriptives by Matched Group",
  "",
  df_to_md(practice_table),
  "",
  "## Notes",
  "",
  "- HS / SM / CC are now shown as item-level means; the separate total-score table is also saved if you still need it.",
  "- TA01 and TA02 were excluded.",
  "- Practice ID 021 was manually assigned to the feedback group.",
  "- Practice descriptives include only rows with a usable roster match after the TA exclusion and 021 resolution."
)

writeLines(md_lines, con = file.path(output_dir, "descriptive_summary.md"))
