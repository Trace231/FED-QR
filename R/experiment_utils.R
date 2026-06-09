#' Create standard experiment output directories
#'
#' @param root Project root.
#' @return Invisibly returns the created directory paths.
#' @export
make_experiment_dirs <- function(root = ".") {
  dirs <- file.path(root, c("results", "figures"))
  for (dir in dirs) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(dirs)
}

#' Save standard experiment outputs
#'
#' @param name Output file prefix.
#' @param result Result returned by [run_fedqr_methods()].
#' @param root Project root.
#' @return Invisibly returns paths written.
#' @export
save_experiment_outputs <- function(name, result, root = ".") {
  make_experiment_dirs(root)
  paths <- c(
    summary = file.path(root, "results", paste0(name, "_summary.csv")),
    trace = file.path(root, "results", paste0(name, "_trace.csv")),
    coefficients = file.path(root, "results", paste0(name, "_coefficients.csv"))
  )
  write.csv(result$summary, paths[["summary"]], row.names = FALSE)
  write.csv(result$trace, paths[["trace"]], row.names = FALSE)
  write.csv(result$coefficients, paths[["coefficients"]], row.names = FALSE)
  invisible(paths)
}

#' Validate a method summary table
#'
#' @param table Data frame with `method` and `objective`.
#' @return Invisibly returns `TRUE`; errors on invalid tables.
#' @export
validate_result_table <- function(table) {
  required <- c("method", "objective")
  missing <- setdiff(required, names(table))
  if (length(missing) > 0) {
    stop("Result table is missing required columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (any(!is.finite(table$objective))) {
    stop("Result table contains non-finite objective values.", call. = FALSE)
  }
  if ("gap_to_best_observed" %in% names(table) &&
      any(table$gap_to_best_observed < -sqrt(.Machine$double.eps), na.rm = TRUE)) {
    stop("Result table contains negative best-observed gaps.", call. = FALSE)
  }
  invisible(TRUE)
}

fedqr_default_method_cols <- function(methods) {
  palette <- c(
    "QR box-dual" = "#2A6FBB",
    "QR box-dual long" = "#1B9E77",
    "QR box-dual stale" = "#0072B2",
    "QR box-dual robust" = "#CC79A7",
    "QR box-dual stale+robust" = "#009E73",
    "QR box-dual adaptive" = "#E7298A",
    "QR box-dual adaptive+VR" = "#66A61E",
    "FSPG-smooth" = "#984EA3",
    "FedSubGrad" = "#4DAF4A",
    "FedSPD-check" = "#D95F02",
    "FedSPD-smooth" = "#E6AB02",
    "FedQR-ADMM" = "#A6761D"
  )
  out <- palette[methods]
  missing <- is.na(out)
  if (any(missing)) {
    out[missing] <- grDevices::rainbow(sum(missing))
  }
  out
}

#' Plot objective-gap convergence
#'
#' @param trace Trace table from [run_fedqr_methods()].
#' @param path PNG path.
#' @param title Plot title.
#' @param method_cols Optional named color vector.
#' @return Invisibly returns `path`.
#' @export
plot_convergence <- function(trace, path, title = "Federated QR convergence",
                             method_cols = NULL) {
  stopifnot(all(c("round", "method", "gap_to_best_observed") %in% names(trace)))
  methods <- unique(trace$method)
  if (is.null(method_cols)) {
    method_cols <- fedqr_default_method_cols(methods)
  }
  grDevices::png(path, width = 1500, height = 1000, res = 180)
  op <- graphics::par(mar = c(5, 5, 4, 1))
  on.exit({
    graphics::par(op)
    grDevices::dev.off()
  }, add = TRUE)
  y <- pmax(trace$gap_to_best_observed, 0)
  graphics::plot(
    NA,
    xlim = range(trace$round),
    ylim = c(0, stats::quantile(y, 0.98, na.rm = TRUE)),
    xlab = "Communication round",
    ylab = "Objective gap to best observed",
    main = title
  )
  for (method in methods) {
    d <- trace[trace$method == method, , drop = FALSE]
    graphics::lines(d$round, pmax(d$gap_to_best_observed, 0),
                    col = method_cols[[method]], lwd = 2)
  }
  graphics::legend("topright", legend = methods, col = method_cols[methods],
                   lwd = 2, bty = "n", cex = 0.8)
  invisible(path)
}

#' Plot final objective gaps
#'
#' @param summary Summary table from [run_fedqr_methods()].
#' @param path PNG path.
#' @param title Plot title.
#' @param method_cols Optional named color vector.
#' @return Invisibly returns `path`.
#' @export
plot_final_gap <- function(summary, path, title = "Final objective",
                           method_cols = NULL) {
  validate_result_table(summary)
  methods <- summary$method
  if (is.null(method_cols)) {
    method_cols <- fedqr_default_method_cols(methods)
  }
  grDevices::png(path, width = 1300, height = 900, res = 180)
  op <- graphics::par(mar = c(9, 5, 4, 1))
  on.exit({
    graphics::par(op)
    grDevices::dev.off()
  }, add = TRUE)
  vals <- pmax(summary$gap_to_best_observed, 1e-8)
  names(vals) <- methods
  graphics::barplot(
    log10(vals),
    col = method_cols[methods],
    las = 2,
    cex.names = 0.75,
    ylab = "log10(gap to best observed)",
    main = title
  )
  invisible(path)
}
