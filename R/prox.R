soft_threshold <- function(z, gamma) {
  sign(z) * pmax(abs(z) - gamma, 0)
}

prox_mcp <- function(z, step_size, lambda_vec, mcp_gamma = 3) {
  az <- abs(z)
  out <- numeric(length(z))
  zero <- az <= step_size * lambda_vec
  middle <- az > step_size * lambda_vec & az <= mcp_gamma * lambda_vec
  large <- az > mcp_gamma * lambda_vec

  denom <- pmax(1 - step_size / mcp_gamma, .Machine$double.eps)
  out[zero] <- 0
  out[middle] <- sign(z[middle]) *
    (az[middle] - step_size * lambda_vec[middle]) / denom
  out[large] <- z[large]
  out
}

prox_scad <- function(z, step_size, lambda_vec, scad_a = 3.7) {
  az <- abs(z)
  out <- numeric(length(z))

  region1 <- az <= lambda_vec * (1 + step_size)
  region2 <- az > lambda_vec * (1 + step_size) & az <= scad_a * lambda_vec
  region3 <- az > scad_a * lambda_vec

  out[region1] <- soft_threshold(z[region1], step_size * lambda_vec[region1])
  denom <- pmax(scad_a - 1 - step_size, .Machine$double.eps)
  out[region2] <- ((scad_a - 1) * z[region2] -
    sign(z[region2]) * scad_a * step_size * lambda_vec[region2]) / denom
  out[region3] <- z[region3]
  out
}

prox_penalty <- function(z, gamma, penalty = c("none", "l1", "mcp", "scad"),
                         penalty_factor = NULL, step_size = NULL,
                         lambda_value = NULL, mcp_gamma = 3, scad_a = 3.7) {
  penalty <- match.arg(penalty)

  if (is.null(penalty_factor)) {
    penalty_factor <- rep(1, length(z))
  }

  if (penalty == "none" || gamma <= 0) {
    return(z)
  }

  if (penalty == "l1") {
    out <- soft_threshold(z, gamma * penalty_factor)
    out[penalty_factor == 0] <- z[penalty_factor == 0]
    return(out)
  }

  if (penalty %in% c("mcp", "scad")) {
    if (is.null(step_size)) {
      step_size <- 1
    }
    if (is.null(lambda_value)) {
      lambda_value <- gamma / step_size
    }
    lambda_vec <- lambda_value * penalty_factor
    if (penalty == "mcp") {
      out <- prox_mcp(z, step_size, lambda_vec, mcp_gamma = mcp_gamma)
    } else {
      out <- prox_scad(z, step_size, lambda_vec, scad_a = scad_a)
    }
    out[penalty_factor == 0] <- z[penalty_factor == 0]
    return(out)
  }

  stop("Unsupported penalty: ", penalty)
}
