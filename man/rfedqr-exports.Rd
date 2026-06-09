\name{rfedqr-exports}
\alias{rfedqr-exports}
\alias{check_loss}
\alias{qr_objective}
\alias{client_qr_objective}
\alias{client_loss_summary}
\alias{quantile_coverage}
\alias{calibrate_quantile_intercept}
\alias{calibration_summary}
\alias{qr_pdhg}
\alias{fed_qr_spd}
\alias{fedspd_dp_qr}
\alias{qr_box_fed_pdhg}
\alias{qr_admm}
\alias{fed_qr_admm}
\alias{fed_subgrad_qr}
\alias{fed_smooth_qr}
\alias{fedqr_methods}
\alias{fit_fedqr}
\alias{run_fedqr_methods}
\alias{make_experiment_dirs}
\alias{save_experiment_outputs}
\alias{plot_convergence}
\alias{plot_final_gap}
\alias{validate_result_table}
\alias{make_qr_sim}
\alias{iid_partition}
\alias{dirichlet_partition}
\title{rfedqr exported functions}
\description{
Public interfaces for federated penalized quantile regression experiments.
}
\details{
Core solvers include centralized QR-PDHG, QR box-dual federated PDHG,
FedSPD-DP-style QR, federated subgradient/smoothing baselines, and ADMM
baselines. The wrapper functions \code{fit_fedqr()} and
\code{run_fedqr_methods()} provide a normalized interface for experiment
scripts.
}
\keyword{models}
