#' @title eigen_dist
#' @description Computes PCoA eigendecomposition using a dissimilarity matrix
#'   calculated via \code{\link[vegan]{vegdist}}. Supports all distance methods
#'   available in \code{vegdist}, with Bray-Curtis as the default.
#' @param X n x p numeric matrix of feature data (e.g., microbiome abundance table)
#' @param pos Logical. Default is \code{TRUE}. If \code{TRUE}, only positive eigenvalues (and their
#'   corresponding eigenvectors) are returned. If \code{FALSE}, all eigenvalues
#'   and eigenvectors are returned.
#' @param method Character. Dissimilarity method passed to
#'   \code{\link[vegan]{vegdist}}. Default is \code{"bray"} (Bray-Curtis).
#'   All supported methods: \code{"manhattan"}, \code{"euclidean"},
#'   \code{"canberra"}, \code{"bray"}, \code{"kulczynski"}, \code{"gower"},
#'   \code{"morisita"}, \code{"horn"}, \code{"mountford"}, \code{"jaccard"},
#'   \code{"raup"}, \code{"binomial"}, \code{"chao"}, \code{"altGower"},
#'   \code{"cao"}, \code{"mahalanobis"}, \code{"clark"}, \code{"chisq"},
#'   \code{"chord"}, \code{"hellinger"}, \code{"aitchison"},
#'   \code{"robust.aitchison"}.
#' @param binary Logical. If \code{TRUE}, data are converted to
#'   presence/absence before computing dissimilarity. Default is \code{FALSE}.
#' @param na.rm Logical. If \code{TRUE}, missing values are handled gracefully.
#'   Default is \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link[vegan]{vegdist}}
#'   (e.g., \code{pseudocount} for Aitchison distance).
#' @return A list containing:
#'   \item{eigenvectors}{Matrix of eigenvectors (columns correspond to PCoA axes)}
#'   \item{eigenvalues}{Numeric vector of eigenvalues}
#' @importFrom vegan vegdist
#' @export
#' @examples
#' # Default Bray-Curtis
#' eigen_dist(X, pos = TRUE)
#'
#' # Jaccard distance
#' eigen_dist(X, pos = TRUE, method = "jaccard")
#'
#' # Aitchison distance with pseudocount for zero handling
#' eigen_dist(X, pos = TRUE, method = "aitchison", pseudocount = 0.5)
#'
#' # Binary Bray-Curtis (presence/absence)
#' eigen_dist(X, pos = TRUE, method = "bray", binary = TRUE)

eigen_dist <- function(X, pos = TRUE, method = "bray", binary = FALSE, na.rm = FALSE, ...) {

  # Validate method
  VALID_METHODS <- c("manhattan", "euclidean", "canberra", "bray",
                     "kulczynski", "gower", "morisita", "horn",
                     "mountford", "jaccard", "raup", "binomial", "chao",
                     "altGower", "cao", "mahalanobis", "clark", "chisq",
                     "chord", "hellinger", "aitchison", "robust.aitchison")

  if (!method %in% VALID_METHODS) {
    stop(paste("Invalid method '", method, "'. Must be one of: ",
               paste(VALID_METHODS, collapse = ", "), sep = ""))
  }

  # Step 1: Calculate dissimilarity using specified method
  D <- vegdist(X, method = method, binary = binary, na.rm = na.rm, ...)
  D <- as.matrix(D)

  # Step 2: Center and square D
  S <- -0.5 * (D^2)

  # Step 3: Double center S
  N <- nrow(S)
  one_vector <- matrix(1, nrow = N, ncol = 1)
  k <- one_vector / N
  I <- diag(N)
  A_normalized <- (I - k %*% t(one_vector)) %*% S %*% (I - one_vector %*% t(k))

  # Step 4: Eigen Decomposition
  eigen_decomp <- eigen(A_normalized)

  if (pos == TRUE) {
    # This epsilon is referencing the epsilon selected by acPCoA method. To be consistent when adapting their function,
    # we set the same epsilon.
    # https://github.com/YuWang28/acPCoA/blob/master/R/acPCoA.R
    epsilon <- sqrt(.Machine$double.eps) * 10000
    positive_indices <- which(eigen_decomp$values > epsilon)

    if (length(positive_indices) == 0) {
      stop("No positive eigenvalues found.")
    }

    # Step 5: Obtain reduced dimension X_hat with dimension of positive eigenvalues
    eigenvalues  <- eigen_decomp$values[positive_indices]
    eigenvectors <- eigen_decomp$vectors[, positive_indices]
    print(paste("Number of features in X_hat:", dim(eigenvectors)[2]))

    return(list(eigenvectors = eigenvectors,
                eigenvalues  = eigenvalues))

  } else {

    # Return all eigenvectors and eigenvalues
    eigenvalues  <- eigen_decomp$values
    eigenvectors <- eigen_decomp$vectors
    print(paste("Number of features in X_hat:", dim(eigenvectors)[2]))

    return(list(eigenvectors = eigenvectors,
                eigenvalues  = eigenvalues))
  }
}

#' @title eigen_bray (deprecated)
#' @description Deprecated wrapper for \code{\link{eigen_dist}} with
#'   \code{method = "bray"}. Please use \code{eigen_dist()} directly.
#' @param X n x p numeric matrix of feature data
#' @param pos Logical. If \code{TRUE}, only positive eigenvalues returned.
#' @param binary Logical. Default is \code{FALSE}.
#' @param na.rm Logical. Default is \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link{eigen_dist}}.
#' @export
eigen_bray <- function(X, pos, binary = FALSE, na.rm = FALSE, ...) {
  .Deprecated("eigen_dist", msg = "'eigen_bray' is deprecated. Use 'eigen_dist(method = \"bray\")' instead.")
  eigen_dist(X, pos, method = "bray", binary = binary, na.rm = na.rm, ...)
}
