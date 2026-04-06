#' Check That a Correlation Value Is Within Valid Bounds
#'
#' @param cor A numeric correlation value.
#'
#' @return Invisible NULL. Throws an error if \code{cor} is outside the range.
#'
#' @keywords internal
check_correlation <- function(cor) {
  if (cor < -1 | cor > 1) {
    stop("Error: Correlation is out of bounds (-1 to 1). Check your data.")
  }
}


#' Compute Partial Correlations Between PC Scores, Outcome, and Covariates
#'
#' Estimates partial correlations between each principal coordinate score and
#' the outcome Y conditional on covariates A, and between each score and each
#' covariate conditional on Y. Automatically selects the appropriate correlation
#' measure based on variable types (continuous, binary, categorical).
#'
#' @param scores A numeric matrix of principal coordinate scores (rows = samples,
#'   columns = principal coordinates).
#' @param Z A numeric matrix of covariates. Column names are used to label
#'   results.
#' @param Y A numeric vector of outcome values (binary or continuous).
#'
#' @return A list containing:
#'   \describe{
#'     \item{partial}{A list with:
#'       \describe{
#'         \item{estimates}{A matrix of partial correlations. Rows correspond to
#'           each covariate and Y. Columns correspond to each PC.}
#'         \item{statistics}{A matrix of t-statistics (currently not populated).}
#'         \item{df}{Degrees of freedom for partial correlation tests
#'           (\code{n - 2 - ncol(Z)}).}
#'       }
#'     }
#'   }
#'
#' @details
#' The function constructs a mixed correlation matrix for each PC score using
#' the appropriate measure for each variable pair:
#' \describe{
#'   \item{Continuous vs Continuous}{Pearson correlation}
#'   \item{Binary vs Continuous}{Point-biserial correlation via
#'     \code{\link[ltm]{biserial.cor}}}
#'   \item{Binary vs Binary}{Matthews Correlation Coefficient via
#'     \code{\link[mltools]{mcc}}}
#'   \item{Categorical vs Continuous/Binary}{Eta squared via
#'     \code{\link[rstatix]{anova_test}}}
#' }
#'
#' Partial correlations are derived by inverting the correlation matrix. If the
#' matrix is nearly singular, the Moore-Penrose pseudoinverse is used via
#' \code{\link[MASS]{ginv}}.
#'
#' @importFrom MASS ginv
#' @importFrom mltools mcc
#' @importFrom ltm biserial.cor
#' @importFrom rstatix anova_test
#' @importFrom stats cor sd
#' @export
partial_cor = function(scores, Z, Y){
  if(is.character(Y)){
    print("Change Y from character to numeric")
  }
  # input: scores:matrix of PC scores (e.g., Xv1, Xv2, Xv3)
  # estimate partial (+ marginal) correlations between each PC score and Y|A and A|Y
  # output: list of partial and marginal correlations between each PC score and A and Y
  scores = as.matrix(scores)
  Z = as.matrix(Z)
  colnames(scores)  = paste0("PC",1:ncol(scores))
  n = nrow(scores)

  pcor.Xv.j.Y.given.Z = vector(mode = "numeric",length = ncol(scores))
  pcor.Xv.j.Z.given.Y = rep(list(vector(mode = "numeric", length = ncol(scores))), ncol(Z))
  names(pcor.Xv.j.Z.given.Y) = colnames(Z)

  t.statistics = matrix(nrow=(ncol(Z) + 1),ncol=ncol(scores))
  rownames(t.statistics) = c(paste0(colnames(Z)),"Y")
  colnames(t.statistics) = paste0("PC",1:ncol(scores))

  for (j in 1:ncol(scores)) {

    omega.j <- matrix(0, ncol(Z) + 2, ncol(Z) + 2)  # Initialize omega.j
    colnames(omega.j) <- rownames(omega.j) <- c(colnames(Z), "Y", paste0("Xv", j))

    # A and Y
    for (z in 1:ncol(Z)) {
      # Check if Z[, z] has zero standard deviation
      if (is.numeric(Z[, z]) && sd(Z[, z]) == 0) {
        stop(paste("Error: Column", z, "in Z has no standard deviation (all values are the same)."))
      }

      is_Y_binary <- length(unique(Y)) == 2
      is_Z_binary <- length(unique(Z[, z])) == 2
      is_Z_categorical <- is.factor(Z[, z]) || is.character(Z[, z]) || (!is_Z_binary && length(unique(Z[, z])) <= 10)


      if (is_Z_binary && is_Y_binary) {
        #print("Binary Zi and Binary Y: MCC")
        Y_vector = as.vector(Y)
        Y_vector <- as.numeric(as.character(Y_vector))

        cor_i = mcc(Z[,z], Y_vector)
      }else if(is_Z_binary && !is_Y_binary){
        #print("Z is Binary, Y is Continuous: Point-Biserial Correlation")
        cor_i = biserial.cor(Y, Z[, z],use = "all.obs", level = 2)

      }else if (!is_Z_binary && is_Y_binary) {
        if (is_Z_categorical) {
          # Categorical Z, Binary Y: Eta Squared
          #print("Categorical Z, Binary Y: Eta Squared")
          df_tmp <- data.frame(score = as.vector(Y), group = as.factor(Z[, z]))

          test <- anova_test(data = df_tmp, dv = score, between = group)
          cor_i <-  test$ges[1]
        } else {
          #print("Continuous Z, Binary Y: Point-Biserial")
          cor_i <- biserial.cor(Z[, z], Y, use = "all.obs", level = 2)
        }

      } else if (is_Z_categorical && !is_Y_binary) {
        # Categorical A, Continuous Y: Eta Squared
        df_tmp <- data.frame(score = as.vector(Y), group = as.factor(Z[, z]))
        test <- anova_test(data = df_tmp, dv = score, between = group)

        cor_i <- test$ges[1]

      }else {
        # Both Ai and Y are continuous: Use Pearson Correlation
        #print("Continuous Ai and Continuous Y: Pearson Correlation")
        cor_i = cor(Z[, z], Y, method = "pearson")
      }
      check_correlation(cor_i)
      omega.j[z, ncol(Z) + 1] <- omega.j[ncol(Z) + 1, z] <- cor_i
    }


    # Compute correlation between scores and Y
    if (length(table(Y)) == 2) {
      #print("Y is Binary: Point-Biserial Correlation with scores")
      cor_i = biserial.cor( scores[, j],Y, use = "all.obs", level = 2)
      check_correlation(cor_i)
      omega.j[ncol(Z) + 1, ncol(Z) + 2] <- omega.j[ncol(Z) + 2, ncol(Z) + 1] <- cor_i
    } else {
      #print("Y is Continuous: Pearson Correlation with scores")
      cor_i = cor(Y, scores[, j], method = "pearson")
      check_correlation(cor_i)
      omega.j[ncol(Z) + 1, ncol(Z) + 2] <- omega.j[ncol(Z) + 2, ncol(Z) + 1] <- cor_i
    }

    # Compute correlation between scores and A
    for (z in 1:ncol(Z)) {
      if (length(table(Z[, z])) == 2) {
        #print("Zi is Binary: Point-Biserial Correlation with scores")
        cor_i =biserial.cor(scores[, j], Z[, z], use = "all.obs", level = 2)
        check_correlation(cor_i)
        omega.j[z, ncol(Z) + 2] <- omega.j[ncol(Z) + 2, z] <- cor_i
      } else if(length(table(Z[, z])) > 2 && length(table(Z[, z])) < 10 ){
        # Zi is Categorical with >2 levels: use eta squared (0 to 1)
        #print("Z is categorical with > 2 levels. Use ANOVA eta score.")
        df_tmp <- data.frame(score = scores[, j], group = Z[, z])

        result <- tryCatch({
          test <- anova_test(data = df_tmp, dv = score, between = group)
          eta_sq <- test$ges[1]  # Generalized eta squared
          check_correlation(eta_sq)
          eta_sq
        }, error = function(e) {
          warning(paste("ANOVA failed for Z column", z))
          NA
        })

        omega.j[z, ncol(Z) + 2] <- omega.j[ncol(Z) + 2, z] <- result

      }else {
        #print("Ai is Continuous: Pearson Correlation with scores")
        cor_i = cor(Z[, z], scores[, j], method = "pearson")
        check_correlation(cor_i)

        omega.j[z, ncol(Z) + 2] <- omega.j[ncol(Z) + 2, z] <- cor_i
      }
    }

    # Fill in the correlations among the A variables
    for (z1 in 1:ncol(Z)) {
      for (z2 in z1:ncol(Z)) {  # We only fill in the upper triangle (or lower triangle since it's symmetric)
        if (z1 != z2) {  # Skip diagonal as it's already filled with variance
          if (length(table(Z[, z1])) == 2 && length(table(Z[, z2])) == 2) {
            #print("Both Z1 and Z2 are Binary: Matthews Correlation Coefficient (MCC)")
            omega.j[z1, z2] <- omega.j[z2, z1] <- mcc(Z[, z1], Z[, z2])
          } else if (length(table(Z[, z1])) == 2 && length(table(Z[, z2])) > 2) {
            #print("Z1 is Binary, Z2 is Continuous: Point-Biserial Correlation")
            cor_i = biserial.cor(Z[, z2], Z[, z1], use = "all.obs", level = 2)
            check_correlation(cor_i)
            omega.j[z1, z2] <- omega.j[z2, z1] <- cor_i
          } else if (length(table(Z[, z1])) > 2 && length(table(Z[, z2])) == 2) {
            #print("Z1 is Continuous, Z2 is Binary: Point-Biserial Correlation")
            cor_i = biserial.cor(Z[, z1], Z[, z2], use = "all.obs", level = 2)
            check_correlation(cor_i)
            omega.j[z1, z2] <- omega.j[z2, z1] <- cor_i
          } else {
            #print("Both Z1 and Z2 are Continuous: Pearson Correlation")
            cor_i = cor(Z[, z1], Z[, z2], method = "pearson")
            check_correlation(cor_i)
            omega.j[z1, z2] <- omega.j[z2, z1] <- cor_i
          }
        }
      }
    }

    diag(omega.j) = 1

    diag(omega.j) <- diag(omega.j) + 1e-10

    # Calculate the inverse of omega.j. (If not singular, a pseudoinverse is calculated)
    det_omega <- det(omega.j)
    if (abs(det_omega) < 1e-10) {
      print(paste("Determinant =", det_omega, "Matrix is nearly singular, performing Moore-Penrose Pseudoinverse"))
      #P.j <- ginv(omega.j)
      omega.j_reg <- omega.j + diag(1e-6, nrow(omega.j))
      P.j <- ginv(omega.j_reg)
      dimnames(P.j) <- dimnames(omega.j)
    } else {
      P.j <- solve(omega.j)
    }

    # Compute partial correlation between scores and Y given A
    pcor.Xv.j.Y.given.Z[j] <- -P.j["Y", ncol(Z) + 2] / sqrt(P.j["Y", "Y"] * P.j[ncol(Z) + 2, ncol(Z) + 2])

    # Compute partial correlations between scores and each A given Y
    for (z in 1:ncol(Z)) {
      pcor.Xv.j.Z.given.Y[[colnames(Z)[z]]][j] <- -P.j[colnames(Z)[z], ncol(Z) + 2] / sqrt(P.j[colnames(Z)[z], colnames(Z)[z]] * P.j[ncol(Z) + 2, ncol(Z) + 2])
    }
  }

  partial.correlations = rbind(do.call("rbind",pcor.Xv.j.Z.given.Y),
                               Y = pcor.Xv.j.Y.given.Z)

  colnames(partial.correlations) = paste0("PC",1:ncol(scores))

  out.partial = list(estimates = partial.correlations,
                     statistics = t.statistics)

  return(list(partial = out.partial))

}
