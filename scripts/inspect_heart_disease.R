args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "inspect_heart_disease.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

raw_dir <- file.path(root, "data", "raw", "heart_disease")
processed_dir <- file.path(root, "data", "processed")
results_dir <- file.path(root, "results")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

base_url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease"
files <- c(
  cleveland = "processed.cleveland.data",
  hungarian = "processed.hungarian.data",
  switzerland = "processed.switzerland.data",
  va = "processed.va.data"
)

cols <- c(
  "age", "sex", "cp", "trestbps", "chol", "fbs", "restecg",
  "thalach", "exang", "oldpeak", "slope", "ca", "thal", "num"
)

download_one <- function(center, file) {
  dest <- file.path(raw_dir, file)
  if (!file.exists(dest)) {
    download.file(file.path(base_url, file), dest, mode = "wb", quiet = FALSE)
  }
  dat <- read.csv(dest, header = FALSE, na.strings = "?", col.names = cols)
  dat$center <- center
  dat
}

hd <- do.call(rbind, Map(download_one, names(files), files))
hd$target_binary <- as.integer(hd$num > 0)

write.csv(hd, file.path(processed_dir, "heart_disease_uci_4center.csv"), row.names = FALSE)

center_summary <- do.call(rbind, lapply(split(hd, hd$center), function(d) {
  data.frame(
    center = d$center[1],
    n = nrow(d),
    positive = sum(d$target_binary == 1, na.rm = TRUE),
    positive_rate = mean(d$target_binary == 1, na.rm = TRUE),
    missing_cells = sum(is.na(d)),
    missing_rate = mean(is.na(d))
  )
}))

missing_by_col <- data.frame(
  variable = names(hd),
  missing = colSums(is.na(hd)),
  missing_rate = colMeans(is.na(hd))
)

target_by_center <- as.data.frame.matrix(table(hd$center, hd$num, useNA = "ifany"))
target_by_center$center <- rownames(target_by_center)
rownames(target_by_center) <- NULL
target_by_center <- target_by_center[, c("center", setdiff(names(target_by_center), "center"))]

continuous_candidates <- c("age", "trestbps", "chol", "thalach", "oldpeak")
continuous_missing_by_center <- do.call(rbind, lapply(split(hd, hd$center), function(d) {
  data.frame(center = d$center[1], t(colSums(is.na(d[, continuous_candidates]))))
}))

continuous_summary <- do.call(rbind, lapply(split(hd, hd$center), function(d) {
  do.call(rbind, lapply(continuous_candidates, function(v) {
    x <- d[[v]]
    data.frame(
      center = d$center[1],
      variable = v,
      n_nonmissing = sum(!is.na(x)),
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      min = min(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  }))
}))

core_cols <- c(
  "age", "sex", "cp", "trestbps", "chol", "fbs", "restecg",
  "thalach", "exang", "oldpeak", "num", "center"
)
full_cols <- c(
  "age", "sex", "cp", "trestbps", "chol", "fbs", "restecg",
  "thalach", "exang", "oldpeak", "slope", "ca", "thal", "num", "center"
)

complete_core_by_center <- as.data.frame(table(hd[complete.cases(hd[, core_cols]), "center"]))
names(complete_core_by_center) <- c("center", "complete_n")
complete_full_by_center <- as.data.frame(table(hd[complete.cases(hd[, full_cols]), "center"]))
names(complete_full_by_center) <- c("center", "complete_n")

write.csv(center_summary, file.path(results_dir, "heart_disease_center_summary.csv"), row.names = FALSE)
write.csv(missing_by_col, file.path(results_dir, "heart_disease_missing_by_col.csv"), row.names = FALSE)
write.csv(target_by_center, file.path(results_dir, "heart_disease_target_by_center.csv"), row.names = FALSE)
write.csv(continuous_missing_by_center,
          file.path(results_dir, "heart_disease_continuous_missing_by_center.csv"),
          row.names = FALSE)
write.csv(continuous_summary,
          file.path(results_dir, "heart_disease_continuous_summary.csv"),
          row.names = FALSE)
write.csv(complete_core_by_center,
          file.path(results_dir, "heart_disease_complete_core_by_center.csv"),
          row.names = FALSE)
write.csv(complete_full_by_center,
          file.path(results_dir, "heart_disease_complete_full_by_center.csv"),
          row.names = FALSE)

cat("\nCenter summary:\n")
print(center_summary)

cat("\nMissing by column:\n")
print(missing_by_col)

cat("\nTarget severity by center (num: 0=no disease, 1-4=disease severity):\n")
print(target_by_center)

cat("\nMissing counts for candidate continuous responses:\n")
print(continuous_missing_by_center)

cat("\nComplete cases with core predictors, dropping high-missing slope/ca/thal:\n")
print(complete_core_by_center)

cat("\nComplete cases with all original predictors:\n")
print(complete_full_by_center)
