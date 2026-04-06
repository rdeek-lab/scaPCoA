#' Calculate Feature Importance from PCoA Eigenvectors
#'
#' Computes feature importance scores by calculating the covariance between
#' original features and principal coordinate axes, weighted by eigenvalues.
#'
#' @param eigen_vectors A matrix of eigenvectors from PCoA decomposition.
#' @param eigen_values A numeric vector of corresponding eigenvalues.
#' @param X_train A numeric matrix of training data (rows = samples,
#'   columns = features).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{feature}{Feature name from column names of \code{X_train}.}
#'     \item{sum}{Absolute sum of standardized covariances across all axes.}
#'     \item{sum_weighted}{Absolute eigenvalue-weighted sum of covariances.}
#'     \item{rank_sum}{Rank based on \code{sum} (1 = most important).}
#'     \item{rank_sum_weighted}{Rank based on \code{sum_weighted} (1 = most important).}
#'   }
#'   Additional columns correspond to the standardized covariance with each
#'   principal coordinate axis.
#'
#' @details
#' Standardization of covariances follows the approach described in
#' Legendre & Legendre (1998). Each covariance is scaled by the inverse
#' square root of the eigenvalue divided by \code{n - 1}.
#'
#' @references
#' Legendre, P. & Legendre, L. (1998). \emph{Numerical Ecology}. Elsevier.
#'
#' @importFrom stats cov
#' @export
calc_feature_importance <- function(eigen_vectors, eigen_values, X_train) {

  # Standardize eigenvectors
  eigen_vectors_standard <- scale(eigen_vectors)

  # Compute covariance of variables with all axes
  S <- cov(X_train, eigen_vectors_standard)

  # Weighted sum (eigenvalue-weighted)
  sum_weighted <- as.vector(t(eigen_values %*% t(S[, 1:ncol(S)])))

  # Standardize covariance (see Legendre & Legendre 1998)
  U <- S %*% diag((eigen_values / (nrow(X_train) - 1))^(-0.5))
  colnames(U) <- colnames(eigen_vectors)

  # Convert to data frame (base R, no tibble/dplyr)
  U_df <- as.data.frame(U)

  # Add feature column from row names
  U_df$feature <- rownames(U_df)

  # Add sum and sum_weighted columns
  U_df$sum          <- abs(rowSums(U))
  U_df$sum_weighted <- abs(sum_weighted)

  # Add rank columns
  U_df$rank_sum          <- rank(-U_df$sum,          ties.method = "min")
  U_df$rank_sum_weighted <- rank(-U_df$sum_weighted, ties.method = "min")

  # Reorder columns: feature first, then everything else
  col_order <- c("feature", colnames(U), "sum", "sum_weighted",
                 "rank_sum", "rank_sum_weighted")
  U_df <- U_df[, col_order]

  # Reset row names
  rownames(U_df) <- NULL

  return(U_df)
}
