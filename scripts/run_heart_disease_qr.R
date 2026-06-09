args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_heart_disease_qr.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "R", "prox.R"))
source(file.path(root, "R", "objective.R"))
source(file.path(root, "R", "simulate.R"))
source(file.path(root, "R", "qr_pdhg.R"))
source(file.path(root, "R", "fed_qr_spd.R"))

dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)

data_path <- file.path(root, "data", "processed", "heart_disease_uci_4center.csv")
if (!file.exists(data_path)) {
  stop("Heart Disease data not found. Run scripts/inspect_heart_disease.R first.")
}

hd <- read.csv(data_path)

vars <- c(
  "center", "thalach", "age", "sex", "cp", "trestbps", "chol",
  "fbs", "restecg", "exang", "oldpeak", "num"
)
hd_model <- hd[, vars]
hd_model <- hd_model[complete.cases(hd_model), ]

factor_vars <- c("sex", "cp", "fbs", "restecg", "exang", "num")
for (v in factor_vars) {
  hd_model[[v]] <- factor(hd_model[[v]])
}

numeric_vars <- c("age", "trestbps", "chol", "oldpeak")
for (v in numeric_vars) {
  hd_model[[paste0(v, "_z")]] <- as.numeric(scale(hd_model[[v]]))
}
hd_model$thalach_z <- as.numeric(scale(hd_model$thalach))

formula_hd <- thalach_z ~ age_z + sex + cp + trestbps_z + chol_z +
  fbs + restecg + exang + oldpeak_z + num

X <- model.matrix(formula_hd, data = hd_model)
y <- hd_model$thalach_z
center <- hd_model$center
client_indices <- split(seq_along(y), center)

write.csv(hd_model,
          file.path(root, "data", "processed", "heart_disease_qr_model_data.csv"),
          row.names = FALSE)

response_summary <- do.call(rbind, lapply(split(hd_model, hd_model$center), function(d) {
  data.frame(
    center = d$center[1],
    n = nrow(d),
    thalach_mean = mean(d$thalach),
    thalach_sd = sd(d$thalach),
    thalach_q10 = as.numeric(quantile(d$thalach, 0.10)),
    thalach_q50 = as.numeric(quantile(d$thalach, 0.50)),
    thalach_q90 = as.numeric(quantile(d$thalach, 0.90))
  )
}))
write.csv(response_summary,
          file.path(root, "results", "heart_disease_thalach_by_center.csv"),
          row.names = FALSE)

fit_tau <- function(tau) {
  central_pdhg <- qr_pdhg(
    X, y,
    tau = tau,
    max_iter = 3000,
    step_rule = "box",
    trace_every = 50,
    seed = 20260526
  )

  fed_full <- fed_qr_spd(
    X, y,
    client_indices = client_indices,
    tau = tau,
    rounds = 1200,
    clients_per_round = length(client_indices),
    batch_size = max(lengths(client_indices)),
    step_rule = "box",
    trace_every = 20,
    seed = 20260526
  )

  fed_r1 <- fed_qr_spd(
    X, y,
    client_indices = client_indices,
    tau = tau,
    rounds = 1200,
    clients_per_round = 2,
    batch_size = max(lengths(client_indices)),
    step_rule = "box",
    trace_every = 20,
    seed = 20260526
  )

  out <- data.frame(
    tau = tau,
    method = c("central_qr_pdhg", "fed_full_clients", "fed_r1_clients"),
    objective = c(central_pdhg$objective, fed_full$objective, fed_r1$objective),
    eta = c(central_pdhg$eta, fed_full$eta, fed_r1$eta),
    sigma = c(central_pdhg$sigma, fed_full$sigma, fed_r1$sigma)
  )

  if (requireNamespace("quantreg", quietly = TRUE)) {
    rq_fit <- quantreg::rq(formula_hd, data = hd_model, tau = tau)
    rq_beta <- as.numeric(coef(rq_fit))
    rq_obj <- qr_objective(X, y, rq_beta, tau = tau)
    out <- rbind(out, data.frame(
      tau = tau,
      method = "quantreg_rq",
      objective = rq_obj,
      eta = NA_real_,
      sigma = NA_real_
    ))

    coef_cmp <- data.frame(
      tau = tau,
      term = colnames(X),
      quantreg = rq_beta,
      central_qr_pdhg = central_pdhg$beta,
      fed_full_clients = fed_full$beta,
      fed_r1_clients = fed_r1$beta
    )
    coef_cmp$central_abs_diff <- abs(coef_cmp$quantreg - coef_cmp$central_qr_pdhg)
    coef_cmp$fed_full_abs_diff <- abs(coef_cmp$quantreg - coef_cmp$fed_full_clients)
    coef_cmp$fed_r1_abs_diff <- abs(coef_cmp$quantreg - coef_cmp$fed_r1_clients)
  } else {
    coef_cmp <- data.frame()
  }

  local_summary <- do.call(rbind, lapply(names(client_indices), function(center_name) {
    idx <- client_indices[[center_name]]
    fit <- qr_pdhg(
      X[idx, , drop = FALSE],
      y[idx],
      tau = tau,
      max_iter = 2500,
      step_rule = "box",
      trace_every = 50,
      seed = 20260526
    )
    data.frame(
      tau = tau,
      center = center_name,
      n = length(idx),
      local_objective = fit$objective,
      global_objective_of_local_beta = qr_objective(X, y, fit$beta, tau = tau)
    )
  }))

  trace <- rbind(
    data.frame(
      tau = tau,
      method = "central_qr_pdhg",
      iter = central_pdhg$trace$iter,
      objective = central_pdhg$trace$objective
    ),
    data.frame(
      tau = tau,
      method = "fed_full_clients",
      iter = fed_full$trace$round,
      objective = fed_full$trace$objective
    ),
    data.frame(
      tau = tau,
      method = "fed_r1_clients",
      iter = fed_r1$trace$round,
      objective = fed_r1$trace$objective
    )
  )

  list(summary = out, coef = coef_cmp, local = local_summary, trace = trace)
}

taus <- c(0.5, 0.75, 0.9)
runs <- lapply(taus, fit_tau)

summary_tbl <- do.call(rbind, lapply(runs, `[[`, "summary"))
coef_tbl <- do.call(rbind, lapply(runs, `[[`, "coef"))
local_tbl <- do.call(rbind, lapply(runs, `[[`, "local"))
trace_tbl <- do.call(rbind, lapply(runs, `[[`, "trace"))

write.csv(summary_tbl,
          file.path(root, "results", "heart_disease_qr_summary.csv"),
          row.names = FALSE)
write.csv(coef_tbl,
          file.path(root, "results", "heart_disease_qr_coef_comparison.csv"),
          row.names = FALSE)
write.csv(local_tbl,
          file.path(root, "results", "heart_disease_local_qr_summary.csv"),
          row.names = FALSE)
write.csv(trace_tbl,
          file.path(root, "results", "heart_disease_qr_trace.csv"),
          row.names = FALSE)

cat("\nHeart Disease response summary:\n")
print(response_summary)

cat("\nHeart Disease QR summary:\n")
print(summary_tbl)

if (nrow(coef_tbl) > 0) {
  cat("\nMax coefficient absolute differences vs quantreg:\n")
  print(aggregate(
    cbind(central_abs_diff, fed_full_abs_diff, fed_r1_abs_diff) ~ tau,
    data = coef_tbl,
    FUN = max
  ))
}

cat("\nLocal QR summary:\n")
print(local_tbl)
