if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
  library(rfedqr)
  test_check("rfedqr")
} else {
  message("testthat is not installed; skipping testthat tests.")
}
