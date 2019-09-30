#' Multinominal Choice Model
#'
#' @param Y Dependent Variable: an (NxJ) matrix where N represents the number of individuals and J is the number of choices
#' @param X Independent Variable(s): (NxK) where K represents number of independent variables
#' @param dimV An optional argument (scalar) representing dimensions of the support space for the error terms, default value is 5.
#' @param nu An optional argument (scalar) representing the trade-off between prediction and precision
#' @param p0 Optional: Prior probabilities associated with the unknown P matrix
#' @param w0 Optional: Prior probabilities associated with the error terms
#' @param optim_method Optional: same as the "method" argument for the "optim" function in "stats"
#'
#' @return A list which includes estimated Largange Multipliers (LMs), Hessian matrix associated with LMs
#'
#' @export
#'
#' @examples
#' Y <- c(0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1)
#' Y <- matrix(Y, nrow = 4, byrow = TRUE)
#' X <- c(1, 3, 1, 2, 1, 5, 1, 4)
#' X <- matrix(X, nrow = 4, byrow = TRUE)
#' gce_mult(Y, X)
#' gce_mult(Y, X, 3)
#' gce_mult(Y, X, 3, 0.6)
gce_mult <- function(Y, X, dimV, nu, p0, w0, optim_method = "BFGS") {
  dimX <- dim(X)
  N <- dimX[1]
  K <- dimX[2]
  J <- dim(Y)[2]
  if (missing(dimV)) {
    M <- 5
  } else {
    M <- dimV
  }
  v <- matrix(seq(from = -1, to = 1, length.out = M), nrow = 1)
  if (missing(nu)) nu <- 0.5
  if (missing(p0)) p0 <- matrix(1 / J, N, J)
  if (missing(w0)) w0 <- array(1 / M, c(N, J, M))

  lambda0 <- rep(0, K * (J - 1))
  gce_optim <- optim(lambda0, gce_mult_obj, gce_mult_grad,
                     Y = Y, X = X, v = v, nu = nu, p0 = p0,
                     w0 = w0, N = N, K = K, J = J, M = M,
                     method = optim_method)
  lambda_hat <- matrix(gce_optim$par, K, J - 1)
  lambda_hat <- cbind(rep(0, K), lambda_hat)
  gce_all <- list("par" = lambda_hat, "convergence" = gce_optim$convergence,
                  "H" = - gce_optim$value, "Y" = Y, "X" = X, "v" = v, "nu" = nu,
                  "p0"= p0, "w0" = w0, "type" = "gce")
  return(gce_all)

}

#' Title
#'
#' @param gme_mult_results b
#'
#' @return z
#' @export
#'
#' @examples
#' \dontrun{
#' gce_mult(Y, X)
#' gce_mult(Y, X, 3)
#' gce_mult(Y, X, 3, 0.6)
#' }
summary_mult <- function(gme_mult_results) {
  Y  <- gme_mult_results$Y
  X  <- gme_mult_results$X
  v  <- gme_mult_results$v
  nu <- gme_mult_results$nu
  p0 <- gme_mult_results$p0
  w0 <- gme_mult_results$w0
  beta <- -gme_mult_results$par
  lambda <- gme_mult_results$par

  dimX <- dim(X)
  N <- dimX[1]
  K <- dimX[2]
  J <- dim(Y)[2]
  M <- length(v)

  temp  <- X %*% lambda
  p     <- p0 * exp(temp / (1 - nu))
  Omega <- apply(p, 1, sum)
  for (i in 1:N) {
    p[i, ] <- p[i, ] / Omega[i]
  }

  w   <- array(0, c(N, J, M))
  Psi <- matrix(0, N, J)
  e   <- matrix(0, N, J)
  for (i in 1:N) {
    for (j in 1:J) {
      w[i, j, ] <- exp(temp[i, j] * v / nu)
      Psi[i, j] <- sum(w[i, j, ])
      w[i, j, ] <- w[i, j, ] / Psi[i, j]
      e[i, j]   <- sum(v * w[i, j, ])
    }
  }

  y_hat <- p + e

  CM <- matrix(0, J, J)
  for (n in 1:N) {
    t1 <- match(1, Y[n, ])
    t2 <- which.max(y_hat[n, ])
    CM[t1, t2] <- CM[t1, t2] + 1
  }

  dp_dx <- array(0, c(N, J, K))
  for (i in 1:N) {
    for (j in 1:J) {
      for (k in 1:K) {
        dp_dx[i, j, k] <- p[i, j] * (lambda[k, j] - sum(lambda[k, ] * p[i, ]))
      }
    }
  }
  ave_dp_dx <- apply(dp_dx, c(2, 3), mean)

  est_sum_all <- list("y_hat" = y_hat,"p_hat" = p, "w_hat" = w, "e_hat" = e,
                      "CM" = CM, "beta" = beta, "marg_eff" = ave_dp_dx)
  return(est_sum_all)
}

.gce_mult_obj <- function(lambda_vector, Y, X, v, nu, p0, w0, N, K, J, M) {

  lambda <- matrix(lambda_vector, K, (J - 1))
  lambda <- cbind(rep(0, K), lambda)

  temp  <- X %*% lambda
  p <- p0 * exp(temp / (1 - nu))
  Omega <- apply(p, 1, sum)

  w <- array(0, c(N, J, M))
  Psi <- matrix(0, N, J)
  for (i in 1:N) {
    for (j in 1:J) {
      w[i, j, ] <- w0[i, j, ] * exp(temp[i, j] * v / nu)
      Psi[i, j] <- sum(w[i, j, ])
    }
  }

  like  <- - sum(lambda * (t(X) %*% Y)) + (1 - nu) * sum(log(Omega)) + nu * sum(log(Psi))
  return(like)

}

.gce_mult_grad <- function(lambda_vector, Y, X, v, nu, p0, w0, N, K, J, M) {

  lambda <- matrix(lambda_vector, K, (J - 1))
  lambda <- cbind(rep(0, K), lambda)

  temp  <- X %*% lambda
  p <- p0 * exp(temp / (1 - nu))
  Omega <- apply(p, 1, sum)
  for (i in 1:N) {
    p[i, ] <- p[i, ] / Omega[i]
  }

  w   <- array(0, c(N, J, M))
  Psi <- matrix(0, N, J)
  e   <- matrix(0, N, J)
  for (i in 1:N) {
    for (j in 1:J) {
      w[i, j, ] <- w0[i, j, ] * exp(temp[i, j] * v / nu)
      Psi[i, j] <- sum(w[i, j, ])
      w[i, j, ] <- w[i, j, ] / Psi[i, j]
      e[i, j]   <- sum(as.vector(v) * w[i, j, ])
    }
  }

  dlambda <- rep(0, K * (J -1) )
  index_lambda <- 1
  for (j in 2:J) {
    for (k in 1:K) {
      dlambda[index_lambda] <- - sum(Y[, j] * X[, k]) +
        sum(p[, j] * X[, k]) +
        sum(e[, j] * X[, k])
      index_lambda <- index_lambda + 1
    }
  }
  return(dlambda)

}

