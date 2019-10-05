#' Matrix Balancing
#'
#' This function
#'
#' @param y A vector representing the sum of all the cells in respective rows of the unknown P matrix
#' @param x A vector representing the sum of all the cells in respective columns of the unknown P matrix
#' @param dimV An optional argument (scalar) representing dimensions of the support space for the error terms, default value is 5.
#' @param nu An optional argument (scalar in \code{[0,1]}) representing the trade-off between prediction and precision
#' @param p0 Optional: Prior probabilities associated with the P matrix
#' @param w0 Optional: Prior probabilities associated with the error terms
#' @param optim_method Optional: same as the "method" argument for the "optim" function in "stats"
#'
#' @return This function returns a list which has the following elements.
#' \itemize{
#'   \item lambda - Estimated Lagrange Multipliers.
#'   \item hess - Hessian matrix associated with the Lagrange Multipliers.
#'   \item p - Estimated transition matrix.
#'   \item w - Estimated probabilities associated with the error terms.
#'   \item e - Estimated Residuals.
#'   \item Sp - The (signal) information of the whole system.
#'   \item S_p_j - Reflects the (signal) information in each moment.
#'   \item p_e_j - Error bounds based on the weaker version of the Fano's inequality
#'   \item H_p_w - Value of the joint entropies of \code{p} and \code{w} at the final iteration.
#'   \item ER - Entropy Ratio Statistic.
#'   \item Pseudo_R2 - Pseudo R-squared.
#'   \item conv - convergence (same as in the optim function).
#' }
#' @export
#'
#' @examples
#' y <- c(0.3, 0.2, 0.5)
#' x <- c(0.1, 0.25, 0.3, 0.2, 0.15)
#' gce_matrix(y, x)
#' gce_matrix(y, x, 3)
#' gce_matrix(y, x, 3, 0.6)
gce_matrix <- function(y, x, dimV, nu, p0, w0, optim_method = "BFGS") {
  N <- length(y)
  J <- length(x)
  if (missing(dimV)) {
    M <- 5
  } else {
    M <- dimV
  }
  v <- matrix(seq(from = -1, to = 1, length.out = M), ncol = 1)
  if (missing(nu)) nu <- 0.5
  if (missing(p0)) {
    p0 <- matrix(1 / N, N, J)
    p0_is_uniform <- TRUE
  } else {
    p0_is_uniform <- FALSE
  }
  if (missing(w0)) w0 <- matrix(1 / M, N, M)

  lambda0 <- rep(0, N)
  gce_optim <- optim(lambda0, gce_matrix_obj, gce_matrix_grad,
                     y = y, x = x, v = v, nu = nu, p0 = p0,
                     w0 = w0, N = N, J = J, M = M,
                     method = optim_method, hessian = TRUE)
  lambda <- gce_optim$par

  p     <- matrix(0, N, J)
  Omega <- rep(0, J)
  for (j in 1:J) {
    for (i in 1:N) {
      p[i, j] <- p0[i, j] * exp(lambda[i] * x[j] / (1 - nu))
    }
    Omega[j] <- sum(p[, j])
    p[, j]   <- p[, j] / Omega[j]
  }

  w   <- matrix(0, N, M)
  e   <- rep(0, N)
  Psi <- rep(0, N)
  for (i in 1:N) {
    for (m in 1:M) {
      w[i, m] <- w0[i, m] * exp(lambda[i] * v[m] / nu)
    }
    Psi[i] <- sum(w[i, ])
    w[i, ] <- w[i, ] / Psi[i]
    e[i]   <- sum(v * w[i, ])
  }

  # Information measures, Golan (1988)
  if (p0_is_uniform == TRUE) {
    Sp <- -sum(p * log(p)) / (J * log(N))
  } else {
    Sp <-  sum(p * log(p)) / sum(p0 * log(p0))
  }

  S_p_j <- rep(0, J)
  for (j in 1:J) {
    if (p0_is_uniform == TRUE) {
      S_p_j[j] <- -sum(p[, j] * log(p[, j])) / log(N)
    } else {
      S_p_j[j] <-  sum(p[, j] * log(p[, j])) / sum(p0[, j] * log(p0[, j]))
    }
  }

  # Error bounds based on the Fano's [weaker] inequality
  p_e_j <- rep(0, J)
  for (j in 1:J) {
    p_e_j[j] <- S_p_j[j] - log(N)
  }

  # Entropy ratio Statistic
  if (p0_is_uniform == TRUE) {
    ER <- 2 * J * log(N) * (1 - Sp)
  } else {
    R <- 2 * (-sum(p0[, j] * log(p0[, j]))) * (1 - Sp)
  }
  # Pseudo R-squared
  R2 <- 1 - Sp

  info_estim_all <- list("lambda" = lambda, "p" = p, "w" = w, "e" = e, "Sp" = Sp,
                         "S_p_j" = S_p_j, "p_e_j" = p_e_j, "H_p_w" = - gce_optim$value,
                         "ER" = ER, "Pseudo_R2" = R2, "conv" = gce_optim$convergence,
                         "hess" = gce_optim$hessian)

  return(info_estim_all)

}

gce_matrix_obj <- function(lambda, y, x, v, nu, p0, w0, N, J, M) {

  p     <- matrix(0, N, J)
  Omega <- rep(0, J)
  for (j in 1:J) {
    for (i in 1:N) {
      p[i, j] <- p0[i, j] * exp(lambda[i] * x[j] / (1 - nu))
    }
    Omega[j] <- sum(p[, j])
  }

  w   <- matrix(0, N, M)
  Psi <- rep(0, N)
  for (i in 1:N) {
    for (m in 1:M) {
      w[i, m] <- w0[i, m] * exp(lambda[i] * v[m] / nu)
    }
    Psi[i] <- sum(w[i, ])
  }

  like  <- - sum(lambda * y) + (1 - nu) * sum(log(Omega)) + nu * sum(log(Psi))
  return(like)

}

gce_matrix_grad <- function(lambda, y, x, v, nu, p0, w0, N, J, M) {

  p     <- matrix(0, N, J)
  Omega <- rep(0, J)
  for (j in 1:J) {
    for (i in 1:N) {
      p[i, j] <- p0[i, j] * exp(lambda[i] * x[j] / (1 - nu))
    }
    Omega[j] <- sum(p[, j])
    p[, j]   <- p[, j] / Omega[j]
  }

  w   <- matrix(0, N, M)
  Psi <- rep(0, N)
  for (i in 1:N) {
    for (m in 1:M) {
      w[i, m] <- w0[i, m] * exp(lambda[i] * v[m] / nu)
    }
    Psi[i] <- sum(w[i, ])
    w[i, ] <- w[i, ] / Psi[i]
  }

  dlambda <- - y + p %*% x + w %*% v
  return(dlambda)

}


