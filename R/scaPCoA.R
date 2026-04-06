
#' PLS Maximization Function
#'
#' Computes the matrix-vector product for partial least squares eigenvalue
#' problem. Called internally by \code{\link{scapcoa}}.
#'
#' @param v A numeric vector.
#' @param args A list containing \code{X} (feature matrix) and \code{Y}
#'   (outcome vector).
#'
#' @return A numeric vector resulting from \code{X'YY'X * v}.
#'
#' @keywords internal
plsMax = function(v,args){
  X = args$X
  Y = args$Y

  out = tcrossprod(crossprod(X,Y))%*%matrix(v,ncol=1)
  # Replace NA values with 0
  out[is.na(out)] <- 0

  return(out)
}

#' scapcoa Maximization Function
#'
#' Computes the matrix-vector product for the penalized decomposition
#' eigenvalue problem. Called internally by \code{\link{scapcoa}}.
#'
#' @param v A numeric vector.
#' @param args A list containing \code{X} (feature matrix),
#'   \code{X.penalize} (penalization matrix), \code{Y} (outcome),
#'   \code{Z} (covariates), and \code{lambda} (tuning parameter).
#'
#' @return A numeric vector resulting from the penalized objective function
#'
#' @keywords internal
scapcoaMax = function(v,args){
  # maximization function for scapcoa (called within function scapcoa()).
  # output: v which maximizes objective function
  X = args$X
  X.penalize = args$X.penalize
  Y = args$Y
  Z = args$Z
  lambda = args$lambda # tuning parameter lambda

  if (is.null(Z)) stop("Error: Z is NULL")

  XTY=crossprod(X,Y)
  K=Z%*%t(Z)


  out = XTY%*%t(XTY)%*%v - lambda*t(X.penalize)%*%K%*%X.penalize%*%v # based on output from calAv function in acPCA R package (Lin et al., 2016)

  # Replace NA values with 0
  out[is.na(out)] <- 0


  return(out)
}

#' Penalized Decomposition for Confounding Removal (scapcoa)
#'
#' Performs penalized eigen decomposition that maximizes association with the
#' outcome Y while penalizing association with covariates Z.
#'
#' @param X A numeric matrix of features (rows = samples, columns = features).
#' @param X.penalize A numeric matrix used for the penalization term. Must have
#'   the same number of columns as \code{X}. Can be \code{NULL} if
#'   \code{lambda = 0}.
#' @param Z Z numeric matrix of covariates to adjust for.
#' @param Y A numeric vector or matrix of outcome values.
#' @param lambda Numeric. Penalization tuning parameter. When \code{lambda = 0},
#'   no penalization is applied (equivalent to PLS).
#' @param nPC Integer. Number of principal components to compute.
#' @param centerX Logical. Whether to center X. Default is \code{FALSE}.
#' @param centerZ Logical. Whether to center Z. Default is \code{FALSE}.
#' @param centerY Logical. Whether to center Y. Default is \code{FALSE}.
#' @param scaleX Logical. Whether to scale X. Default is \code{FALSE}.
#' @param scaleZ Logical. Whether to scale Z. Default is \code{FALSE}.
#' @param scaleY Logical. Whether to scale Y. Default is \code{FALSE}.
#'
#' @return A list with:
#'   \describe{
#'     \item{values}{A numeric vector of eigenvalues.}
#'     \item{vectors}{A matrix of eigenvectors.}
#'   }
#'   Returns \code{NULL} if the penalized matrix is not symmetric.
#'
#'
#' @importFrom RSpectra eigs_sym
#' @keywords internal
scapcoa_decomp <- function(X, X.penalize = NULL, Z, Y, lambda, nPC, centerX = FALSE, centerZ = FALSE, centerY = FALSE, scaleX = FALSE, scaleZ = FALSE, scaleY = FALSE) {
  # Center and scale input matrices
  Z <- as.matrix(Z)
  X <- as.matrix(X)
  X <- scale(X, center = centerX, scale = scaleX)


  if (!is.null(X.penalize) & lambda> 0) { # check for dimension
    if (ncol(X) != ncol(X.penalize)) {
      stop("X and X.penalize should have the same number of columns")
    }


    X.penalize <- as.matrix(X.penalize)
    X.penalize <- scale(X.penalize, center = centerX, scale = scaleX)


    args <- list(X = X, X.penalize = X.penalize, Y = Y, Z = Z, lambda = lambda)


    X = args$X
    X.penalize = args$X.penalize
    Y = args$Y
    Z = args$Z
    Z <- as.matrix(Z)
    lambda= args$lambda

    K <- Z %*% t(Z)  # Kernel matrix for penalization
    XTY <- crossprod(X, Y)  # Compute X'Y

    penalized_matrix <- XTY %*% t(XTY) - lambda* t(X.penalize) %*% K %*% X.penalize
    penalized_matrix[is.na(penalized_matrix)] <- 0  # Replace NAs with 0

    is_symmetric <- function(mat) all(abs(mat - t(mat)) < 1e-10)
    #print(paste("Penalized matrix is symmetric:",is_symmetric(penalized_matrix)))

    if (!is_symmetric(penalized_matrix)) {
      message(paste("Skipping lambda =", lambda, "because penalized_matrix is not symmetric"))
      return(NULL)
    }

    # Handle eigenvalue computation
    if (nPC == ncol(X)) {
      print("Full spectrum computation using eigen()")
      # Full spectrum computation using eigen()

      eigen_out <- eigen(penalized_matrix, symmetric = TRUE)

      eig_out <- list(values = eigen_out$values, vectors = eigen_out$vectors)
    } else {
      # Use eigs_sym for partial spectrum computation
      eig_out <- eigs_sym(
        scapcoaMax,
        k = nPC,
        which = "LA",
        n = ncol(X),
        args = args,
        opts = list(tol = 1e-10, maxitr = 1e6)
      )
      #niter = eig_out$niter
      #nops = eig_out$nops
      if(all(eig_out$values == 0)){
        message("Warning: All eigenvalues are zero, returning NA")
      }
    }
  } else {
    # No penalty term (maximize X'YY'X)
    args <- list(X = X, Y = Y)

    if (nPC == ncol(X)) {
      # Full spectrum computation using eigen()
      X = args$X
      Y = args$Y

      XTY <- crossprod(X, Y)  # Compute X'Y
      pls_matrix <- XTY %*% t(XTY)  # Compute full matrix

      eigen_out <- eigen(pls_matrix, symmetric = TRUE)
      eig_out <- list(values = eigen_out$values, vectors = eigen_out$vectors)
    }  else {
      # Use eigs_sym for partial spectrum computation
      eig_out <- eigs_sym(
        plsMax,
        k = nPC,
        which = "LA", #largest eigenvalues (positive)
        n = ncol(X),
        args = args,
        opts = list(tol = 1e-10, maxitr = 1e6)
      )
    }
  }
  return(eig_out)
}


#' Tune Lambda for scapcoa
#'
#' Selects the optimal penalization parameter lambda by maximizing the
#' eigenvalue-weighted difference between partial correlations with the
#' outcome Y and partial correlations with covariates Z.
#'
#' @param X.orig A numeric matrix of original features.
#' @param X.max A numeric matrix used for the maximization term.
#' @param X.penalize A numeric matrix used for the penalization term.
#' @param lambdas A numeric vector of candidate lambda values to evaluate.
#' @param Z A numeric matrix of covariates.
#' @param Y A numeric vector or matrix of outcome values.
#' @param nPC Integer. Number of principal components to compute.
#' @param centerX Logical. Whether to center X. Default is \code{FALSE}.
#' @param scaleX Logical. Whether to scale X. Default is \code{FALSE}.
#' @param centerY Logical. Whether to center Y. Default is \code{FALSE}.
#' @param scaleY Logical. Whether to scale Y. Default is \code{FALSE}.
#' @param centerZ Logical. Whether to center Z. Default is \code{FALSE}.
#' @param scaleZ Logical. Whether to scale Z. Default is \code{FALSE}.
#' @param n.cores Integer. Number of cores for parallel computation.
#'   Default is 1.
#' @param plot Logical. Currently unused. Default is \code{FALSE}.
#' @param thresh Numeric. Currently unused. Default is 0.001.
#'
#' @return A list with:
#'   \describe{
#'     \item{weighted_sum_abs_diff}{Numeric vector of tuning criterion values
#'       for each valid lambda.}
#'     \item{lambdas}{The candidate lambda values evaluated.}
#'     \item{lambda_tune}{The optimal lambda that maximizes the tuning
#'       criterion. \code{NA} if no valid lambda was found.}
#'   }
#'
#' @details
#' For each candidate lambda, the function runs \code{\link{scapcoa}},
#' computes in-sample PC scores, and evaluates partial correlations via
#' \code{\link{partial_cor}}. The tuning criterion is the mean across
#' covariates of the eigenvalue-weighted difference:
#' \code{|pcor(Y, PC)| - |pcor(Z, PC)|}.
#'
#' @importFrom parallel mclapply
#' @export
scapcoa.tune = function(X.orig, X.max, X.penalize, lambdas, Z, Y, nPC, centerX = F, scaleX = F, centerY = F, scaleY = F, centerZ = F, scaleZ = F, n.cores = 1, plot = FALSE, thresh = 0.001){

  sum_abs_diff = numeric(length = length(lambdas))
  names_Z = colnames(Z)

  #lambdas = sort(lambdas)  # sort in ascending order

  best_lambda = best_weighted_sum_abs_diff = numeric(length = length(lambdas))
  all_weighted_sum_abs_diff = c()
  valid_lambdas = c()  # Track only lambdas that contribute to all_weighted_sum_abs_diff

  # Iterate through all lambdas without splitting into groups
  for (lambda_index in 1:length(lambdas)) {
    lambda_i = lambdas[lambda_index]
    print(paste("lambda=", lambda_i))

    weighted_sum_abs_diff = unlist(mclapply(lambda_i, FUN = function(lambda_i) {
      tryCatch({

        # Penalized decomposition using residuals for a given lambda:
        temp_out = scapcoa_decomp(X = X.max, X.penalize = X.penalize, Z = Z, Y = Y, lambda= lambda_i, nPC = nPC,
                            centerX = centerX, scaleX = scaleX, centerZ = centerZ, scaleZ = scaleZ,
                            centerY = centerY, scaleY = scaleY)

        if (is.null(temp_out)){
          return(NA)
        }
        # In-sample scores for lambda = l:
        temp_scores = X.orig %*% temp_out$vectors

        # Partial correlations:
        temp_pcor = partial_cor(scores = temp_scores, Z = Z, Y = Y)$partial$estimates

        # Weighted sum of difference of absolute values of partial correlations
        # weighted_sum_abs_diff_pcor = sapply(1:ncol(A), FUN = function(a) {
        #   (temp_out$values / sum(temp_out$values)) %*% (abs(temp_pcor["Y", ]) - abs(temp_pcor[names_Z[a], ]))
        # })

        weighted_sum_abs_diff_pcor = sapply(1:ncol(Z), FUN = function(z) {
          if (sum(temp_out$values) == 0) {
            return(0)  # Avoid division by zero, return NA
          } else {
            return((temp_out$values / sum(temp_out$values)) %*%
                     (abs(temp_pcor["Y", ]) - abs(temp_pcor[names_Z[z], ])))
          }
        })

        return(mean(weighted_sum_abs_diff_pcor))
      }, error = function(e){
        message(paste("Error for lambda =", lambda_i, ":", e$message))
        return(NA)  # Return NA if there's an error
      })

    },mc.cores = n.cores))

    # Store results for each lambda, only if values are valid
    if (!all(is.na(weighted_sum_abs_diff))) {
      all_weighted_sum_abs_diff = c(all_weighted_sum_abs_diff, weighted_sum_abs_diff)
      valid_lambdas = c(valid_lambdas, lambda_i)
    }
  }
  df <- data.frame(all_weighted_sum_abs_diff,valid_lambdas)
  # Check if df is empty before proceeding
  if (nrow(df) > 0) {
    max_index <- which.max(df$all_weighted_sum_abs_diff)
    lambda_max <- df$valid_lambdas[max_index]
  } else {
    lambda_max <- NA
    message("No valid lambda tuning values were found.")
  }

  return(list(weighted_sum_abs_diff = all_weighted_sum_abs_diff,
              lambdas = lambdas[1:length(all_weighted_sum_abs_diff)],
              lambda_tune = lambda_max))
}


#' Check Matrix Symmetry
#'
#' @param mat A numeric matrix.
#'
#' @return Logical. \code{TRUE} if the matrix is symmetric.
#'
#' @keywords internal
is_symmetric <- function(mat) {
  all.equal(mat, t(mat)) == TRUE
}




#' Supervised and Covariate-Adjusted Principal Coordinates Analysis (scaPCoA)
#'
#' Main function for performing scaPCoA or scaPCoA-R (with residualization).
#' Finds principal coordinates that maximize association with the outcome Y
#' while minimizing association with confounding covariates Z.
#'
#' @param obj A list containing:
#'   \describe{
#'     \item{Z}{A data frame or matrix of covariates.}
#'     \item{Y}{A numeric vector of outcome values.}
#'     \item{X}{A numeric matrix of features (e.g., microbiome counts).}
#'   }
#' @param method Character. Dissimilarity method passed to
#'   \code{\link[vegan]{vegdist}}. Default is \code{"bray"} (Bray-Curtis).
#'   All supported methods: \code{"manhattan"}, \code{"euclidean"},
#'   \code{"canberra"}, \code{"bray"}, \code{"kulczynski"}, \code{"gower"},
#'   \code{"morisita"}, \code{"horn"}, \code{"mountford"}, \code{"jaccard"},
#'   \code{"raup"}, \code{"binomial"}, \code{"chao"}, \code{"altGower"},
#'   \code{"cao"}, \code{"mahalanobis"}, \code{"clark"}, \code{"chisq"},
#'   \code{"chord"}, \code{"hellinger"}, \code{"aitchison"},
#'   \code{"robust.aitchison"}.
#' @param lambda_ls A numeric vector of candidate lambda values for tuning.
#' @param nPC Integer. Number of principal components to compute. Default is 3.
#' @param residualization Logical. If \code{TRUE}, use the scaPCoA-R variant
#'   with residualized matrices. Default is \code{FALSE}.
#'
#' @return A list with:
#'   \describe{
#'     \item{PC.scores}{A matrix of principal coordinate scores.}
#'     \item{v_hat}{A matrix of estimated eigenvectors.}
#'     \item{best.lambda}{The optimal lambda selected by tuning.}
#'     \item{eigen_values}{A numeric vector of eigenvalues.}
#'     \item{scores.partial_cor}{A matrix of partial correlations between
#'       PC scores and Y/Z.}
#'   }
#' @importFrom compositions clr
#' @export
scapcoa <- function(obj,method = "bray",
                    lambda_ls = seq(0, 10, 0.5),
                    nPC=3, residualization = FALSE){
    results = list()

    Z = obj$Z
    Y = obj$Y
    X = obj$X

    Y <- as.matrix(Y)
    Z <- as.matrix(Z)
    X <- as.matrix(X)

    n = nrow(X)

    if(length(unique(Y)) > 2){
      y_cont = TRUE
    }else{
      y_cont = FALSE
    }

    eigen_output <-  eigen_dist(X, method = method)
    eigen_vectors <- eigen_output$eigenvectors
    eigen_values = eigen_output$eigenvalues

    X_hat <- eigen_vectors %*% diag(sqrt(eigen_values))
    X_hat_train = X_hat
    X.orig = X_hat

    m = ncol(X_hat)

    if(residualization == TRUE){
      print("Run scaPCoA-R")
      resid.dat = get_resid(X.orig,Y, Z)
      X.max = resid.dat$X.star
      X.penalize = resid.dat$X.tilde
    }else{
      print("Run scaPCoA")
      X.max = X_hat
      X.penalize =X_hat
    }


    lambda.tune = scapcoa.tune(X.orig = X.orig,
                                  X.max = X.max,
                                  X.penalize = X.penalize,
                                  lambdas = lambda_ls,
                                  Z = Z,
                                  Y = Y,
                                  nPC = m, centerX = T, scaleX = T)

      best.lambda = lambda.tune$lambda_tune
      print(paste('Best Lambda from tuning:',best.lambda))

      scapcoa.out = scapcoa_decomp(X = X.max,
                              X.penalize = X.penalize,
                              Z = Z,
                              Y = Y,
                              lambda = best.lambda,
                              nPC = m, centerX = T, scaleX = T)

      v_hat = scapcoa.out$vectors
      PC.scores = X.orig%*%v_hat
      eigen_values <- scapcoa.out$values

      # Partial Correlation
      # Look at correlations between the first few PC scores and the nuisance variables (Z1, Z2, Y)
      cor.scores = partial_cor(PC.scores, Z, Y)

      scores.partial_cor = cor.scores$partial$estimates

      # Loop for y and PCs
      for (j in 1:m) {  # Assuming 3 principal components (PC1, PC2, PC3)
        results[[paste0("pcor_y_pc", j)]] <- abs(scores.partial_cor[ncol(Z) + 1, j])
      }

      # Loop for all Z and PCs
      for (i in 1:ncol(Z)) {  # Loop through Z1, Z2, ..., Zl
        for (j in 1:m) {  # Loop through PC1, PC2, PC3
          results[[paste0("pcor_z", i, "_pc", j)]] <- abs(scores.partial_cor[i, j])
        }
      }

      return(list(PC.scores = PC.scores, v_hat = v_hat, best.lambda = best.lambda, eigen_values = eigen_values, scores.partial_cor= scores.partial_cor))
  }



