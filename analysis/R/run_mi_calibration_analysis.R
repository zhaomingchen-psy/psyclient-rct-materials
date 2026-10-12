library(foreign)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
  root_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
} else {
  root_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
out_dir <- file.path(root_dir, "output", "r_calibration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fmt_p <- function(p) {
  if (is.na(p)) {
    return(NA_character_)
  }
  if (p < .001) {
    return("< .001")
  }
  sprintf("= %.3f", p)
}

z_score <- function(x) {
  as.numeric(scale(x))
}

make_quartile_factor <- function(x) {
  qs <- quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE, type = 7)
  qs <- unique(qs)
  if (length(qs) < 5) {
    return(factor(rep(NA_character_, length(x))))
  }
  cut(
    x,
    breaks = qs,
    include.lowest = TRUE,
    labels = c("Q1_low", "Q2", "Q3", "Q4_high")
  )
}

run_cor_summary <- function(df, conf_var, perf_var, time_label) {
  keep <- df[complete.cases(df[, c(conf_var, perf_var)]), c(conf_var, perf_var)]
  if (nrow(keep) < 10) {
    return(NULL)
  }
  pearson <- cor.test(keep[[conf_var]], keep[[perf_var]], method = "pearson")
  spearman <- cor.test(keep[[conf_var]], keep[[perf_var]], method = "spearman", exact = FALSE)
  data.frame(
    timepoint = time_label,
    confidence = conf_var,
    performance = perf_var,
    n = nrow(keep),
    pearson_r = unname(pearson$estimate),
    pearson_p = pearson$p.value,
    spearman_rho = unname(spearman$estimate),
    spearman_p = spearman$p.value
  )
}

run_quartile_summary <- function(df, conf_var, perf_var, time_label) {
  keep <- df[complete.cases(df[, c(conf_var, perf_var)]), c(conf_var, perf_var)]
  if (nrow(keep) < 12) {
    return(list(desc = NULL, omnibus = NULL))
  }
  keep$performance_quartile <- make_quartile_factor(keep[[perf_var]])
  keep <- keep[!is.na(keep$performance_quartile), ]
  if (nrow(keep) < 12) {
    return(list(desc = NULL, omnibus = NULL))
  }
  desc <- aggregate(
    keep[[conf_var]],
    by = list(performance_quartile = keep$performance_quartile),
    FUN = function(x) c(n = length(x), mean = mean(x), sd = sd(x))
  )
  desc_out <- data.frame(
    timepoint = time_label,
    confidence = conf_var,
    performance = perf_var,
    performance_quartile = desc$performance_quartile,
    n = desc$x[, "n"],
    mean = desc$x[, "mean"],
    sd = desc$x[, "sd"]
  )

  fit <- aov(keep[[conf_var]] ~ keep$performance_quartile)
  fit_tab <- summary(fit)[[1]]
  omni <- data.frame(
    timepoint = time_label,
    confidence = conf_var,
    performance = perf_var,
    df1 = fit_tab[1, "Df"],
    df2 = fit_tab[2, "Df"],
    F_value = fit_tab[1, "F value"],
    p_value = fit_tab[1, "Pr(>F)"]
  )

  list(desc = desc_out, omnibus = omni)
}

plot_scatter <- function(df, conf_var, perf_var, title, xlab, ylab, out_path) {
  keep <- df[complete.cases(df[, c(conf_var, perf_var)]), c(conf_var, perf_var)]
  if (nrow(keep) < 10) {
    return(invisible(NULL))
  }
  png(out_path, width = 1600, height = 1200, res = 200)
  par(mar = c(4.5, 4.5, 3, 1))
  plot(
    keep[[perf_var]],
    keep[[conf_var]],
    pch = 19,
    col = "#1b9e77",
    xlab = xlab,
    ylab = ylab,
    main = title
  )
  abline(lm(keep[[conf_var]] ~ keep[[perf_var]]), col = "#333333", lwd = 2)
  grid()
  dev.off()
}

sav_path <- file.path(root_dir, "Psyclient_data.sav")
sav_df <- read.spss(
  sav_path,
  to.data.frame = TRUE,
  use.value.labels = FALSE,
  trim.factor.names = TRUE,
  trim_values = TRUE
)

sav_df$mi_perf_post <- z_score((z_score(log1p(sav_df$RQR_new_post)) + z_score(sav_df$questioning_post)) / 2)
sav_df$mi_perf_pre <- z_score((z_score(log1p(sav_df$RQR_new_PRE)) + z_score(sav_df$questioning_pre)) / 2)

confidence_vars_post <- c("HST6", "SMT6", "CCT6", "COSET6")
confidence_vars_pre <- c("HST1", "SMT1", "CCT1", "COSET1")

cor_rows <- list()
quartile_desc <- list()
quartile_omnibus <- list()

for (conf_var in confidence_vars_post) {
  cor_rows[[length(cor_rows) + 1]] <- run_cor_summary(sav_df, conf_var, "mi_perf_post", "post")
  qres <- run_quartile_summary(sav_df, conf_var, "mi_perf_post", "post")
  quartile_desc[[length(quartile_desc) + 1]] <- qres$desc
  quartile_omnibus[[length(quartile_omnibus) + 1]] <- qres$omnibus
}

for (conf_var in confidence_vars_pre) {
  cor_rows[[length(cor_rows) + 1]] <- run_cor_summary(sav_df, conf_var, "mi_perf_pre", "pre")
  qres <- run_quartile_summary(sav_df, conf_var, "mi_perf_pre", "pre")
  quartile_desc[[length(quartile_desc) + 1]] <- qres$desc
  quartile_omnibus[[length(quartile_omnibus) + 1]] <- qres$omnibus
}

cor_out <- do.call(rbind, cor_rows)
quartile_desc_out <- do.call(rbind, quartile_desc)
quartile_omnibus_out <- do.call(rbind, quartile_omnibus)

write.csv(cor_out, file.path(out_dir, "mi_calibration_correlations.csv"), row.names = FALSE)
write.csv(quartile_desc_out, file.path(out_dir, "mi_calibration_quartile_descriptives.csv"), row.names = FALSE)
write.csv(quartile_omnibus_out, file.path(out_dir, "mi_calibration_quartile_anova.csv"), row.names = FALSE)

plot_scatter(
  sav_df,
  conf_var = "HST6",
  perf_var = "mi_perf_post",
  title = "HST6 and MI-style performance",
  xlab = "MI-style performance composite (post)",
  ylab = "HST6",
  out_path = file.path(out_dir, "hst6_mi_performance_scatter.png")
)

plot_scatter(
  sav_df,
  conf_var = "COSET6",
  perf_var = "mi_perf_post",
  title = "COSET6 and MI-style performance",
  xlab = "MI-style performance composite (post)",
  ylab = "COSET6",
  out_path = file.path(out_dir, "coset6_mi_performance_scatter.png")
)

key_post <- subset(cor_out, timepoint == "post")
key_pre <- subset(cor_out, timepoint == "pre")

lines <- c(
  "# MI Calibration Check",
  "",
  "This is not a formal Dunning-Kruger test, because the confidence variables are self-efficacy scores rather than explicit self-estimated performance rankings.",
  "The analysis instead checks whether lower MI-style performance is associated with higher confidence or whether better performers tend to report lower confidence.",
  "",
  "## Post-Test Correlations",
  ""
)

for (i in seq_len(nrow(key_post))) {
  row <- key_post[i, ]
  lines <- c(
    lines,
    sprintf(
      "- %s vs MI performance (post): r = %.3f, p %s; Spearman rho = %.3f, p %s; n = %d",
      row$confidence,
      row$pearson_r,
      fmt_p(row$pearson_p),
      row$spearman_rho,
      fmt_p(row$spearman_p),
      row$n
    )
  )
}

lines <- c(lines, "", "## Pre-Test Correlations", "")
for (i in seq_len(nrow(key_pre))) {
  row <- key_pre[i, ]
  lines <- c(
    lines,
    sprintf(
      "- %s vs MI performance (pre): r = %.3f, p %s; Spearman rho = %.3f, p %s; n = %d",
      row$confidence,
      row$pearson_r,
      fmt_p(row$pearson_p),
      row$spearman_rho,
      fmt_p(row$spearman_p),
      row$n
    )
  )
}

lines <- c(lines, "", "## Interpretation", "")
lines <- c(
  lines,
  "- A classic Dunning-Kruger-like pattern would require lower performers to show clearly higher confidence than higher performers.",
  "- In the current data, the post-test confidence variables are weakly positive, not negative, in relation to the MI-style performance composite.",
  "- The quartile checks therefore should be treated as a descriptive calibration screen rather than evidence of a Dunning-Kruger effect."
)

writeLines(lines, file.path(out_dir, "mi_calibration_notes.md"))

cat(file.path(out_dir, "mi_calibration_correlations.csv"), "\n")
cat(file.path(out_dir, "mi_calibration_quartile_descriptives.csv"), "\n")
cat(file.path(out_dir, "mi_calibration_quartile_anova.csv"), "\n")
cat(file.path(out_dir, "mi_calibration_notes.md"), "\n")
