#' Linear Regression Model
#'
#' @param y  An (Nx1) vector representing the dependent variable where N is the number of observations
#' @param X  An (NxK) matrix representing a set of independent variables where K is number of regressors
#' @param Z  An (KxM) matrix representing support spaces  the for regression coefficients where M is the dimension of the support spaces
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
#' @return A list
#' \itemize{
#'   \item lambda - Estimated Lagrange Multipliers.
#'   \item beta - Regression coefficients.
#'   \item var_beta - Variance-covariance matrix of the regression coefficients.
#'   \item p - Estimated probabilities associated with the regressions coefficients.
#'   \item w - Estimated probabilities associated with the error terms.
#'   \item e - Estimated Residuals.
#'   \item Sp - The (signal) information of the whole system.
#'   \item Sp_k - The (signal) information associated with the choice k-th regression coefficient.
#'   \item H_p_w - Value of the joint entropies of \code{p} and \code{w} at the final iteration.
#'   \item dH - Delta-H from the Entropy Concentration Theorem.
#'   \item ER - Entropy Ratio Statistic.
#'   \item Pseudo_R2 - Pseudo R-squared.
#'   \item conv - Hessian matrix associated with the estimated Lagrange multipliers.
#' }
#' @export
#'
#' @examples
#' set.seed(123)
#' y <- runif(100)
#' X <- matrix(runif(200), nrow = 100, ncol = 2)
#' Z <- matrix(c(-1, 0, 1, -1, 0, 1), nrow = 2, byrow = TRUE)
#' gce_lin(y, X, Z)
#' gce_lin(y, X, Z, 3)
#' gce_lin(y, X, Z, 3, 0.6)
gce_lin <- function(y, X, Z, v, nu, p0, w0, optim_method = "BFGS") {

  dimX <- dim(X)
  N <- dimX[1]
  K <- dimX[2]
  M <- dim(Z)[2]
  if (missing(v)) {
    dimV <- 5
    sd_y <- sd(y)
    v <- seq(from = -3 * sd_y, to = 3 * sd_y, length.out = dimV)
  } else {
    if (is.vector(v)) {
      len_v <- length(v)
      if (len_v == 1) {
        dimV <- v
        sd_y <- sd(y)
        v <- seq(from = -3 * sd_y, to = 3 * sd_y, length.out = dimV)
      } else if (len_v > 1) {
        dimV <- len_v
      }
    } else if (is.matrix(v) && dim(v)[1] == N) {
      dimV <- dim(v)[2]
    }
  }
  if (missing(nu)) nu <- 0.5
  if (missing(p0)) p0 <- matrix(1 / M, nrow = K, ncol = M)
  if (missing(w0)) w0 <- matrix(1 / dimV, nrow = N, ncol = dimV)

  lambda0 <- rep(0, N)

  gce_optim <- optim(lambda0, gce_lin_obj, gce_lin_grad, y = y, X = X, Z = Z,
                     v = v, nu = nu, p0 = p0, w0 = w0, N = N, K = K, M = M,
                     J = dimV, method = optim_method, hessian = T)

  lambda_hat <- gce_optim$par

  p <- matrix(0, K, M)
  Omega <- rep(0, K)
  for (k in 1:K) {
    temp <- sum(lambda_hat * X[, k])
    for (m in 1:M) {
      p[k, m] <- p0[k, m] * exp(Z[k, m] * temp / (1 - nu))
    }
    Omega[k] <- sum(p[k, ])
    p[k, ]   <- p[k, ] / Omega[k]
  }
  beta_hat <- matrix(apply(Z * p, 1, sum), ncol = 1)

  Psi <- rep(0, N)
  if (is.vector(v)) {
    J <- length(v)
    w <- matrix(0, N, J)
    for (n in 1:N) {
      for (j in 1:J) {
        w[n, j] <- w0[n, j] * exp(v[j] * lambda_hat[n] / nu)
      }
      Psi[n]  <- sum(w[n, ])
      w[n, ] <- w[n, ] / Psi[n]
    }
    e <- w %*% matrix(v, ncol = 1)
  } else if (is.matrix(v) && dim(v)[1] == N) {
    J <- dim(v)[2]
    w <- matrix(0, N, J)
    for (n in 1:N) {
      for (j in 1:J) {
        w[n, j] <- w0[n, j] * exp(v[n, j] * lambda_hat[n] / nu)
      }
      Psi[n]  <- sum(w[n, ])
      w[n, ] <- w[n, ] / Psi[n]
    }
    e <- apply(w * v, 1, sum)
  }

  # Normalized entropy
  Sp <- sum(p * log(p)) / sum(p0 * log(p0))
  S_pk <- rep(0, K)
  for (k in 1:K) {
    S_pk[k] <- sum(p[k, ] * log(p[k, ])) / sum(p0[k, ] * log(p0[k, ]))
  }

  # INFERENCE AND DISGNOSTICS
  # Asymptotic variance
  sigma2_beta <- sum(lambda_hat * lambda_hat) / N
  sigma2_e <- rep(0, N)
  if (is.vector(v)) {
    for (n in 1:N) {
      sigma2_e[n] <- sum((v * v) * w[n, ]) - (sum(v * w[n, ])) ^ 2
    }
  } else if (is.matrix(v) && dim(v)[1] == N) {
    # !!! IMPORTANT this part of the code will be adjusted further
    for (n in 1:N) {
      sigma2_e[n] <- sum((v[n, ] * v[n, ]) * w[n, ]) - (sum(v[n, ] * w[n, ])) ^ 2
    }
  }
  w2_beta <- (sum(1 / sigma2_e) / N) ^ 2
  var_beta <- (sigma2_beta / w2_beta) * solve(t(X) %*% X)
  # Entropy ratio (assuming Ho[beta = 0])
  ER <- rep(0, K)
  p_temp <- rep(1 / M, M)
  H_R <- -sum(p_temp * log(p_temp))
  for (k in 1:K) {
    H_U <- -sum(p[k, ] * log(p[k, ]))
    ER[k] <- 2 * H_R - 2 * H_U
  }
  # Pseudo R-squared
  R2 <- 1 - Sp
  # Interval associated with the Entropy Concentration Theorem
  CC <- K * (M - 1) + N * (J - 2)
  dH <- qchisq(c(0.9, 0.95, 0.99), df = CC) / (2 * N)
  # All outputs
  info_estim_all <- list("lambda" = lambda_hat, "beta" = beta_hat, "var_beta" = var_beta,
                         "p" = p, "w" = w, "e" = e, "Sp" = Sp, "S_pk" = S_pk,
                         "H_p_w" = gce_optim$H, "dH" = dH, "ER" = ER, "Pseudo_R2" = R2,
                         "conv" = gce_optim$convergence)
  return(info_estim_all)

}

gce_lin_obj <- function(lambda, y, X, Z, v, nu, p0, w0, N, K, ...) {

  Omega <- rep(0, K)
  for (k in 1:K) {
    Omega[k] <- sum(p0[k, ] * exp(Z[k, ] * sum(lambda * X[, k]) / (1 - nu)))
  }

  Psi <- rep(0, N)
  if (is.vector(v)) {
    for (n in 1:N) {
      Psi[n] <- sum(w0[n, ] * exp(lambda[n] * v / nu))
    }
  } else if (is.matrix(v) && dim(v)[1]==N) {
    for (n in 1:N) {
      Psi[n] <- sum(w0[n, ] * exp(lambda[n] * v[n, ] / nu))
    }
  }

  l <- -sum(lambda * y) + (1 - nu) * sum(log(Omega)) + nu * sum(log(Psi))
  return(l)

}

gce_lin_grad <- function(lambda, y, X, Z, v, nu, p0, w0, N, K, M, J) {

  p <- matrix(0, K, M)
  Omega <- rep(0, K)
  for (k in 1:K) {
    temp <- sum(lambda * X[, k])
    for (m in 1:M) {
      p[k, m] <- p0[k, m] * exp(Z[k, m] * temp / (1 - nu))
    }
    Omega[k] <- sum(p[k, ])
    p[k, ]   <- p[k, ] / Omega[k]
  }
  beta <- matrix(apply(Z * p, 1, sum), ncol = 1)

  w   <- matrix(0, N, J)
  Psi <- rep(0, N)
  if (is.vector(v)) {
    for (n in 1:N) {
      for (j in 1:J) {
        w[n, j] <- w0[n, j] * exp(v[j] * lambda[n] / nu)
      }
      Psi[n]  <- sum(w[n, ])
      w[n, ] <- w[n, ] / Psi[n]
    }
    e <- w %*% matrix(v, ncol = 1)
  } else if (is.matrix(v) && dim(v)[1]==N) {
    for (n in 1:N) {
      for (j in 1:J) {
        w[n, j] <- w0[n, j] * exp(v[n, j] * lambda[n] / nu)
      }
      Psi[n]  <- sum(w[n, ])
      w[n, ] <- w[n, ] / Psi[n]
    }
    e <- apply(w * v, 1, sum)
  }

  lambda_grad <- - y + X %*% beta + e
  return(lambda_grad)

}


