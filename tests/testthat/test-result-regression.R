test_that("existing summary CSV files have finite result metrics", {
  result_dir <- test_path("../../results")
  skip_if_not(dir.exists(result_dir))
  files <- list.files(result_dir, pattern = "summary.*\\.csv$", full.names = TRUE)
  expect_gt(length(files), 0)

  for (file in files) {
    tbl <- read.csv(file)
    if ("objective" %in% names(tbl)) {
      expect_true(all(is.finite(tbl$objective)), info = basename(file))
      if ("gap_to_best_observed" %in% names(tbl)) {
        expect_true(all(tbl$gap_to_best_observed >= -1e-10), info = basename(file))
      }
    }
    if ("final_objective" %in% names(tbl)) {
      expect_true(all(is.finite(tbl$final_objective)), info = basename(file))
    }
    if ("final_gap" %in% names(tbl)) {
      expect_true(all(is.finite(tbl$final_gap)), info = basename(file))
    }
  }
})

test_that("large real-data benchmark conclusions match saved summaries", {
  result_dir <- test_path("../../results")
  skip_if_not(file.exists(file.path(result_dir, "nyc_taxi_large_scale_summary_batch5000.csv")))
  skip_if_not(file.exists(file.path(result_dir, "household_power_large_scale_summary.csv")))
  skip_if_not(file.exists(file.path(result_dir, "hd_baseline_aggregate.csv")))

  nyc <- read.csv(file.path(result_dir, "nyc_taxi_large_scale_summary_batch5000.csv"))
  nyc_best <- nyc$method[which.min(nyc$objective)]
  expect_equal(nyc_best, "QR box-dual long")

  power <- read.csv(file.path(result_dir, "household_power_large_scale_summary.csv"))
  power_best <- power$method[which.min(power$objective)]
  expect_equal(power_best, "QR box-dual long")
  fspg_gap <- power$gap_to_best_observed[power$method == "FSPG-smooth"]
  expect_lt(fspg_gap, 0.001)

  hd <- read.csv(file.path(result_dir, "hd_baseline_aggregate.csv"))
  hd_r1r2 <- hd[hd$scenario == "R1 + R2", ]
  for (tau in unique(hd_r1r2$tau)) {
    qr_gap <- hd_r1r2$final_gap[hd_r1r2$tau == tau & hd_r1r2$method == "QR box-dual"]
    fedspd_gap <- hd_r1r2$final_gap[hd_r1r2$tau == tau & hd_r1r2$method == "FedSPD-check"]
    expect_lt(qr_gap, fedspd_gap)
  }
})

test_that("advanced innovation saved results preserve the main conclusions", {
  result_dir <- test_path("../../results")
  summary_file <- file.path(result_dir, "advanced_innovation_summary.csv")
  client_file <- file.path(result_dir, "advanced_innovation_client_loss.csv")
  large_calibration_file <- file.path(result_dir, "advanced_innovation_large_calibration.csv")
  skip_if_not(file.exists(summary_file))
  skip_if_not(file.exists(client_file))
  skip_if_not(file.exists(large_calibration_file))

  advanced <- read.csv(summary_file)
  extreme_09 <- advanced[advanced$dataset == "simulation" &
    advanced$setting == "extreme" & advanced$tau == 0.9, ]
  extreme_gap <- aggregate(target_gap ~ method, extreme_09, mean)
  stale_gap <- extreme_gap$target_gap[extreme_gap$method == "QR box-dual stale"]
  base_gap <- extreme_gap$target_gap[extreme_gap$method == "QR box-dual"]
  expect_lt(stale_gap, base_gap)

  client_loss <- read.csv(client_file)
  hd <- client_loss[client_loss$dataset == "heart_disease", ]
  hd_worst <- aggregate(worst_client_loss ~ tau + method, hd, mean)
  for (tau in unique(hd_worst$tau)) {
    robust <- hd_worst$worst_client_loss[
      hd_worst$tau == tau & hd_worst$method == "QR box-dual robust"
    ]
    base <- hd_worst$worst_client_loss[
      hd_worst$tau == tau & hd_worst$method == "QR box-dual"
    ]
    expect_lt(robust, base)
  }

  large_cal <- read.csv(large_calibration_file)
  for (dataset in unique(large_cal$dataset)) {
    qr_long <- large_cal[large_cal$dataset == dataset &
      large_cal$method == "QR box-dual long", ]
    raw <- qr_long$global_coverage_error[qr_long$calibration == "raw"]
    global <- qr_long$global_coverage_error[qr_long$calibration == "global_intercept"]
    client <- qr_long$mean_client_coverage_error[qr_long$calibration == "client_offset"]
    raw_client <- qr_long$mean_client_coverage_error[qr_long$calibration == "raw"]
    expect_lte(global, raw)
    expect_lt(client, raw_client)
  }
})

test_that("large-client scale stress winners come from QR box-dual variants", {
  result_dir <- test_path("../../results")
  summary_file <- file.path(result_dir, "scale_stress_summary.csv")
  client_file <- file.path(result_dir, "scale_stress_client_loss.csv")
  skip_if_not(file.exists(summary_file))
  skip_if_not(file.exists(client_file))

  summary_tbl <- read.csv(summary_file)
  gap <- aggregate(target_gap ~ client_count + heterogeneity + tau + method,
                   summary_tbl, mean)
  groups <- interaction(gap$client_count, gap$heterogeneity, gap$tau, drop = TRUE)
  winners <- vapply(split(seq_len(nrow(gap)), groups), function(idx) {
    gap$method[idx[which.min(gap$target_gap[idx])]]
  }, character(1))
  expect_true(all(grepl("^QR box-dual", winners)))

  client_tbl <- read.csv(client_file)
  fair <- aggregate(worst_client_loss ~ client_count + heterogeneity + tau + method,
                    client_tbl, mean)
  fair_groups <- interaction(fair$client_count, fair$heterogeneity, fair$tau, drop = TRUE)
  fair_winners <- vapply(split(seq_len(nrow(fair)), fair_groups), function(idx) {
    fair$method[idx[which.min(fair$worst_client_loss[idx])]]
  }, character(1))
  expect_gte(sum(fair_winners == "QR box-dual adaptive"), 4)
})
