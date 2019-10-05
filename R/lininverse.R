#' Linear Inverse Problem with Noise
#'
#' @param y  An (Nx1) vector representing the dependent variable where N is the number of observations
#' @param X  An (NxK) matrix representing a set of independent variables where K is number of regressors
#' @param v  An optional argument representing a support space for error terms:
#'           (a) if missing then \strong{v} is a (5x1) vector of equally spaced points in \code{[a,b]} interval;
#'           (b) if a scalar (e.g. \strong{H}) then \strong{v} is a (Hx1) vector of equally spaced points in \code{[a,b]} interval;
#'           (c) can be a user-supplied vector; (d) can be a user-supplied matrix.
#'           Please note that in case (a) and (b) the \code{[a,b]} interval is centered around zero, and \strong{a} and \strong{b} are calculated
#'           using the empirical three-sigma rule Pukelsheim (1994).
#' @param nu Optional: A weight parameter representing the trade-off between prediction and precision
#' @param p0 Optional: Prior probabilities associated with the regression coefficients
#' @param w0 Optional: Prior probabilities associated with the error terms
#' @param optim_method Optional: same as the "method" argument for the "optim" function in "stats"
#'
#' @return This function returns a list which has the following elements.
#' \itemize{
#'   \item lambda - Estimated Lagrange Multipliers.
#'   \item hess - Hessian matrix associated with the Lagrange Multipliers.
#'   \item p - Estimated probabilities associated with the regressions coefficients.
#'   \item p_e - Error bounds based on the weaker version of the Fano's inequality
#'   \item w - Estimated probabilities associated with the error terms.
#'   \item e - Estimated Residuals.
#'   \item Sp - The (signal) information of the whole system.
#'   \item H_p_w - Value of the joint entropies of \code{p} and \code{w} at the final iteration.
#'   \item ER - Entropy Ratio Statistic.
#'   \item Pseudo_R2 - Pseudo R-squared.
#'   \item conv - convergence (same as in the optim function).
#' }
#' @export
#'
#' @importFrom stats optim sd qchisq
#'
#' @examples
#' set.seed(123)
#' y <- runif(100)
#' X <- matrix(runif(200), nrow = 100, ncol = 2)
#' lin_inv(y, X)
lin_inv <- function(y, X, v, nu, p0, w0, optim_method = "BFGS") {

  dimX <- dim(X)
  J <- dimX[1]
  N <- dimX[2]
  if (missing(v)) {
    M <- 5
    sd_y <- sd(y)
    v <- seq(from = -3 * sd_y, to = 3 * sd_y, length.out = M)
  } else {
    if (is.vector(v)) {
      len_v <- length(v)
      if (len_v == 1) {
        M <- v
        sd_y <- sd(y)
        v <- seq(from = -3 * sd_y, to = 3 * sd_y, length.out = M)
      } else if (len_v > 1) {
        M <- len_v
      }
    } else if (is.matrix(v) && dim(v)[1] == J) {
      M <- dim(v)[2]
    }
  }
  if (missing(nu)) nu <- 0.5
  if (missing(p0)) {
    p0 <- rep(1 / N, N)
    p0_is_uniform <- TRUE
  } else {
    p0_is_uniform <- FALSE
  }
  if (missing(w0)) w0 <- matrix(1 / M, nrow = J, ncol = M)

  lambda0 <- rep(0, J)

  gce_optim <- optim(lambda0, lin_inv_obj, lin_inv_grad, y = y, X = X,
                     v = v, nu = nu, p0 = p0, w0 = w0, N = N, M = M,
                     J = J, method = optim_method, hessian = T)

  lambda <- gce_optim$par

  p <- rep(0, N)
  for (n in 1:N) {
    p[n] <- p0[n] * exp(sum(lambda * X[, n]) / (1 - nu))
  }
  Omega <- sum(p)
  p <- p / Omega

  w   <- matrix(0, J, M)
  e   <- rep(0, J)
  Psi <- rep(0, J)
  for (j in 1:J) {
    if (is.vector(v)) {
      w[j, ] <- w0[j, ] * exp(lambda[j] * v)
      Psi[j] <- sum(w[j, ])
      w[j, ] <- w[j, ] / Psi[j]
      e[j]   <- sum(v * w[j, ])
    } else if (is.matrix(v) && dim(v)[1] == J) {
      w[j, ] <- w0[j, ] * exp(lambda[j] * v[j, ])
      Psi[j] <- sum(w[j, ])
      w[j, ] <- w[j, ] / Psi[j]
      e[j]   <- sum(v[j, ] * w[j, ])
    }
  }

  # Information measures
  if (p0_is_uniform == TRUE) {
    Sp <- -sum(p * log(p)) / log(N)
  } else {
    Sp <-  sum(p * log(p)) / sum(p0 * log(p0))
  }

  # INFERENCE AND DISGNOSTICS
  # Error bounds for weaker version of the Fano's inequality
  p_e <- Sp - 1 / log(N)
  # Entropy ratio (assuming Ho[beta = 0])
  ER <- 2 * log(N) * (1 - Sp)
  # Pseudo R-squared
  R2 <- 1 - Sp
  # All outputs
  info_estim_all <- list("lambda" = lambda, "hess" = gce_optim$hessian,
                         "p" = p, "p_e" = p_e, "w" = w, "e" = e, "Sp" = Sp,
                         "H_p_w" = - gce_optim$value, "ER" = ER, "Pseudo_R2" = R2,
                         "conv" = gce_optim$convergence)
  return(info_estim_all)

}

lin_inv_obj <- function(lambda, y, X, v, nu, p0, w0, N, J, M) {

  p <- rep(0, N)
  for (n in 1:N) {
    p[n] <- p0[n] * exp(sum(lambda * X[, n]) / (1 - nu))
  }
  Omega <- sum(p)

  w   <- matrix(0, J, M)
  Psi <- rep(0, J)
  for (j in 1:J) {
    if (is.vector(v)) {
      w[j, ] <- w0[j, ] * exp(lambda[j] * v)
    } else if (is.matrix(v) && dim(v)[1] == J) {
      w[j, ] <- w0[j, ] * exp(lambda[j] * v[j, ])
    }
    Psi[j] <- sum(w[j, ])
  }

  l <- -sum(lambda * y) + (1 - nu) * log(Omega) + nu * sum(log(Psi))
  return(l)

}

lin_inv_grad <- function(lambda, y, X, v, nu, p0, w0, N, J, M) {

  p <- rep(0, N)
  for (n in 1:N) {
    p[n] <- p0[n] * exp(sum(lambda * X[, n]) / (1 - nu))
  }
  Omega <- sum(p)
  p <- p / Omega

  w   <- matrix(0, J, M)
  e   <- rep(0, J)
  Psi <- rep(0, J)
  for (j in 1:J) {
    if (is.vector(v)) {
      w[j, ] <- w0[j, ] * exp(lambda[j] * v)
      Psi[j] <- sum(w[j, ])
      w[j, ] <- w[j, ] / Psi[j]
      e[j]   <- sum(v * w[j, ])
    } else if (is.matrix(v) && dim(v)[1] == J) {
      w[j, ] <- w0[j, ] * exp(lambda[j] * v[j, ])
      Psi[j] <- sum(w[j, ])
      w[j, ] <- w[j, ] / Psi[j]
      e[j]   <- sum(v[j, ] * w[j, ])
    }
  }

  lambda_grad <- - y + X %*% matrix(p, ncol = 1) + e
  return(lambda_grad)

}

