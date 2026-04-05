#' Simulated microbiome dataset for scaPCoA
#'
#' A simulated dataset for demonstrating and testing scaPCoA methods.
#'
#' @format A list with 4 elements:
#' \describe{
#'   \item{Z}{Data frame of 3 covariate variables (Z1, Z2, Z3) with 100 observations}
#'   \item{Y}{Integer vector of binary outcome (0/1) with 100 observations}
#'   \item{X}{Numeric matrix of microbiome features (100 x 196)}
#'   \item{features_signals}{Character vector of 15 signal feature names}
#' }
#' @source Simulated using the zinLDA model
#'   (\url{https://github.com/rebeccadeek/zinLDA}).
"data_microbe"
