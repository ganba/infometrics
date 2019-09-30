#' Matrix Balancing
#'
#' @param y A vector representing the sum of all the cells in respective rows of the unknown P matrix ma
#' @param x A vector representing the sum of all the cells in respective columns of the unknown P matrix ma
#' @param dimV An optional argument (scalar) representing dimensions of the support space for the error terms, default value is 5.
#' @param nu An optional argument (scalar) representing the trade-off between prediction and precision
#' @param p0 Optional: Prior probabilities associated with the betas
#' @param w0 Optional: Prior probabilities associated with the error terms
#' @param optim_method Optional: same as "method" argument for "optim" in "stats"
#'
#' @return A list which includes estimated Largange Multipliers (LMs), Hessian matrix associated with LMs
#' @export
#'
#' @examples
#' \dontrun{
#' gce_matrix(y, x)
#' gce_matrix(y, x, 3)
#' gce_matrix(y, x, 3, 0.6)
#' }
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
  if (missing(p0)) p0 <- matrix(1 / N, N, J)
  if (missing(w0)) w0 <- matrix(1 / M, N, M)

  lambda0 <- rep(0, N)
  gce_optim <- optim(lambda0, gce_matrix_obj, gce_matrix_grad,
                     y = y, x = x, v = v, nu = nu, p0 = p0,
                     w0 = w0, N = N, J = J, M = M,
                     method = optim_method)
  lambda_hat <- gce_optim$par
  gce_all <- list("par" = lambda_hat, "convergence" = gce_optim$convergence,
                  "H" = - gce_optim$value, "y" = y, "x" = x, "v" = v, "nu" = nu,
                  "p0"= p0, "w0" = w0, "type" = "gce")
  return(gce_all)

}

#' Summary Statistics for the Matrix Balancing Problem
#'
#' @param gme_mult_results Output from function gce_matrix()
#'
#' @return Estimated parameters of the model
#' @export
#'
#' @examples
#' \dontrun{
#' gce_matrix(y, x)
#' gce_matrix(y, x, 3)
#' gce_matrix(y, x, 3, 0.6)
#' }
summary_matrix <- function(gme_mult_results) {
  y  <- gme_mult_results$y
  x  <- gme_mult_results$x
  v  <- gme_mult_results$v
  nu <- gme_mult_results$nu
  p0 <- gme_mult_results$p0
  w0 <- gme_mult_results$w0
  lambda <- gme_mult_results$par

  N <- length(y)
  J <- length(x)
  M <- length(v)

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

  y_hat <- p %*% x + w %*% v

  dp_dx <- matrix(0, N, J)
  for (i in 1:N) {
    for (j in 1:J) {
      dp_dx[i, j] <- p[i, j] * (lambda[i] - sum(lambda * p[, j]))
    }
  }
  ave_dp_dx <- apply(dp_dx, 2, mean)

  est_sum_all <- list("y_hat" = y_hat,"p_hat" = p,
                      "w_hat" = w, "marg_eff" = ave_dp_dx)
  return(est_sum_all)
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


