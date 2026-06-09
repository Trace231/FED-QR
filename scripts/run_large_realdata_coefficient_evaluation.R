args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_large_realdata_coefficient_evaluation.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "scripts", "load_rfedqr.R"))
make_experiment_dirs(root)

evaluate_coefficients <- function(dataset, X, y, client_indices, tau, coef_path) {
  coef_tbl <- read.csv(coef_path, check.names = FALSE)
  methods <- setdiff(names(coef_tbl), "term")

  client_loss <- list()
  calibration <- list()
  k <- 1
  for (method in methods) {
    beta <- coef_tbl[[method]][match(colnames(X), coef_tbl$term)]
    if (any(!is.finite(beta))) {
      warning("Skipping method with incompatible coefficients: ", method)
      next
    }
    loss <- client_loss_summary(X, y, beta, client_indices, tau = tau)
    loss$dataset <- dataset
    loss$method <- method
    client_loss[[k]] <- loss[, c(
      "dataset", "method", "global_mean_loss", "client_mean_loss",
      "worst_client_loss", "client_loss_sd", "client_q90_loss", "min_client_loss"
    )]

    raw <- calibration_summary(X, y, beta, tau = tau, client_indices = client_indices)
    raw$calibration <- "raw"
    global_cal <- calibrate_quantile_intercept(X, y, beta, tau = tau, mode = "global")
    global <- calibration_summary(X, y, global_cal$beta, tau = tau, client_indices = client_indices)
    global$calibration <- "global_intercept"
    client_cal <- calibrate_quantile_intercept(
      X, y, beta,
      tau = tau,
      client_indices = client_indices,
      mode = "client_offset"
    )
    client <- calibration_summary(
      X, y, client_cal$beta,
      tau = tau,
      client_indices = client_indices,
      offsets = client_cal$offsets
    )
    client$calibration <- "client_offset"
    cal <- rbind(raw, global, client)
    cal$dataset <- dataset
    cal$method <- method
    calibration[[k]] <- cal[, c(
      "dataset", "method", "calibration",
      "global_coverage", "global_coverage_error",
      "mean_client_coverage_error", "worst_client_coverage_error",
      "coverage_sd", "min_client_coverage", "max_client_coverage"
    )]
    k <- k + 1
  }

  list(
    client_loss = do.call(rbind, client_loss),
    calibration = do.call(rbind, calibration)
  )
}

outputs <- list()

nyc_data <- file.path(root, "data", "processed", "nyc_taxi_yellow_qr_2024q1_sample.csv")
nyc_coef <- file.path(root, "results", "nyc_taxi_large_scale_coefficients.csv")
if (file.exists(nyc_data) && file.exists(nyc_coef)) {
  message("Evaluating NYC Taxi saved coefficients...")
  taxi <- read.csv(nyc_data)
  taxi <- taxi[complete.cases(taxi), ]
  taxi$client_id <- as.integer(taxi$client_id)
  for (v in c("log_trip_distance", "duration_min", "passenger_count", "congestion_surcharge", "airport_fee")) {
    taxi[[paste0(v, "_z")]] <- as.numeric(scale(taxi[[v]]))
  }
  taxi$duration_log_z <- as.numeric(scale(log1p(taxi$duration_min)))
  taxi$hour_sin <- sin(2 * pi * taxi$pickup_hour / 24)
  taxi$hour_cos <- cos(2 * pi * taxi$pickup_hour / 24)
  taxi$dow_sin <- sin(2 * pi * taxi$pickup_dow / 7)
  taxi$dow_cos <- cos(2 * pi * taxi$pickup_dow / 7)
  taxi$pickup_borough <- factor(taxi$pickup_borough)
  taxi$vendor_id <- factor(taxi$vendor_id)
  taxi$rate_code <- factor(taxi$rate_code)
  taxi$payment_type <- factor(taxi$payment_type)
  formula_nyc <- log_total_amount ~
    log_trip_distance_z + duration_log_z + passenger_count_z +
    congestion_surcharge_z + airport_fee_z +
    hour_sin + hour_cos + dow_sin + dow_cos +
    pickup_borough + vendor_id + rate_code + payment_type
  X <- model.matrix(formula_nyc, data = taxi)
  y <- taxi$log_total_amount
  client_raw <- split(seq_along(y), taxi$client_id)
  client_indices <- client_raw[lengths(client_raw) >= 100]
  client_indices <- client_indices[order(as.integer(names(client_indices)))]
  keep_idx <- sort(unlist(client_indices, use.names = FALSE))
  X <- X[keep_idx, , drop = FALSE]
  y <- y[keep_idx]
  client_indices <- split(seq_along(y), taxi$client_id[keep_idx])
  client_indices <- client_indices[order(as.integer(names(client_indices)))]
  outputs[["nyc"]] <- evaluate_coefficients("nyc_taxi", X, y, client_indices, tau = 0.9, nyc_coef)
}

power_data <- file.path(root, "data", "processed", "household_power_consumption.csv")
power_coef <- file.path(root, "results", "household_power_large_scale_coefficients.csv")
if (file.exists(power_data) && file.exists(power_coef)) {
  message("Evaluating Household Power saved coefficients...")
  power <- read.csv(power_data, na.strings = c("?", ""), stringsAsFactors = FALSE)
  names(power) <- make.names(names(power))
  numeric_cols <- c(
    "Global_active_power", "Global_reactive_power", "Voltage",
    "Global_intensity", "Sub_metering_1", "Sub_metering_2", "Sub_metering_3"
  )
  for (v in numeric_cols) {
    power[[v]] <- as.numeric(power[[v]])
  }
  power <- power[complete.cases(power[, c("Date", "Time", numeric_cols)]), ]
  date <- as.Date(power$Date, format = "%d/%m/%Y")
  hour <- as.integer(substr(power$Time, 1, 2))
  minute <- as.integer(substr(power$Time, 4, 5))
  dow <- as.POSIXlt(date)$wday
  month <- as.integer(format(date, "%m"))
  year_month <- format(date, "%Y-%m")
  power$log_active_power <- log1p(power$Global_active_power)
  power$reactive_z <- as.numeric(scale(power$Global_reactive_power))
  power$voltage_z <- as.numeric(scale(power$Voltage))
  power$sub1_z <- as.numeric(scale(power$Sub_metering_1))
  power$sub2_z <- as.numeric(scale(power$Sub_metering_2))
  power$sub3_z <- as.numeric(scale(power$Sub_metering_3))
  power$hour_sin <- sin(2 * pi * (hour + minute / 60) / 24)
  power$hour_cos <- cos(2 * pi * (hour + minute / 60) / 24)
  power$dow_sin <- sin(2 * pi * dow / 7)
  power$dow_cos <- cos(2 * pi * dow / 7)
  power$month <- factor(month)
  power$year_month <- factor(year_month)
  formula_power <- log_active_power ~
    reactive_z + voltage_z + sub1_z + sub2_z + sub3_z +
    hour_sin + hour_cos + dow_sin + dow_cos + month
  X <- model.matrix(formula_power, data = power)
  y <- power$log_active_power
  client_raw <- split(seq_along(y), droplevels(power$year_month), drop = TRUE)
  client_indices <- client_raw[lengths(client_raw) >= 1000]
  client_indices <- client_indices[order(names(client_indices))]
  keep_idx <- sort(unlist(client_indices, use.names = FALSE))
  X <- X[keep_idx, , drop = FALSE]
  y <- y[keep_idx]
  client_indices <- split(seq_along(y), droplevels(power$year_month[keep_idx]), drop = TRUE)
  client_indices <- client_indices[order(names(client_indices))]
  outputs[["household_power"]] <- evaluate_coefficients("household_power", X, y, client_indices, tau = 0.9, power_coef)
}

if (length(outputs) == 0) {
  stop("No large-data coefficient files were available for evaluation.")
}

client_loss_tbl <- do.call(rbind, lapply(outputs, `[[`, "client_loss"))
calibration_tbl <- do.call(rbind, lapply(outputs, `[[`, "calibration"))

write.csv(client_loss_tbl,
          file.path(root, "results", "advanced_innovation_large_client_loss.csv"),
          row.names = FALSE)
write.csv(calibration_tbl,
          file.path(root, "results", "advanced_innovation_large_calibration.csv"),
          row.names = FALSE)

cat("\nLarge-data coefficient client loss:\n")
print(client_loss_tbl[order(client_loss_tbl$dataset, client_loss_tbl$worst_client_loss), ])
cat("\nLarge-data coefficient calibration:\n")
print(calibration_tbl[order(calibration_tbl$dataset, calibration_tbl$method, calibration_tbl$calibration), ])
