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
#' @return This function returns a list which has the following elements.
#' \itemize{
#'   \item lambda - Estimated Lagrange Multipliers.
#'   \item hess - Hessian matrix associated with the Lagrange Multipliers.
#'   \item p - Estimated transition matrix.
#'   \item w - Estimated probabilities associated with the error terms.
#'   \item e - Estimated Residuals.
#'   \item marg_eff - Marginal Effects of independent variables on the individual choices.
#'   \item Sp - The (signal) information of the whole system.
#'   \item S_p_j - Reflects the (signal) information in each moment.
#'   \item p_e_j - Error bounds based on the weaker version of the Fano's inequality
#'   \item H_p_w - Value of the joint entropies of \code{p} and \code{w} at the final iteration.
#'   \item ER - Entropy Ratio Statistic.
#'   \item Pseudo_R2 - Pseudo R-squared.
#'   \item CM - Confusion matrix.
#'   \item conv - convergence (same as in the optim function).
#' }
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
  v <- seq(from = -1, to = 1, length.out = M)
  if (missing(nu)) nu <- 0.5
  if (missing(p0)) {
    p0 <- matrix(1 / J, N, J)
    p0_is_uniform <- TRUE
  } else {
    p0_is_uniform <- FALSE
  }
  if (missing(w0)) w0 <- array(1 / M, c(N, J, M))

  lambda0 <- rep(0, K * (J - 1))
  gce_optim <- optim(lambda0, gce_mult_obj, gce_mult_grad,
                     Y = Y, X = X, v = v, nu = nu, p0 = p0,
                     w0 = w0, N = N, K = K, J = J, M = M,
                     method = optim_method, hessian = TRUE)

  lambda <- matrix(gce_optim$par, K, J - 1)
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
      e[i, j]   <- sum(v * w[i, j, ])
    }
  }

  # Marginal effects
  dp_dx <- array(0, c(N, J, K))
  for (i in 1:N) {
    for (j in 1:J) {
      for (k in 1:K) {
        dp_dx[i, j, k] <- (p[i, j] / (1 - nu)) * (lambda[k, j] - sum(lambda[k, ] * p[i, ]))
      }
    }
  }
  ave_dp_dx <- apply(dp_dx, c(2, 3), mean)

  # Confusion matrix
  y_hat <- p + e
  CM    <- matrix(0, J, J)
  for (n in 1:N) {
    t1 <- match(1, Y[n, ])
    t2 <- which.max(y_hat[n, ])
    CM[t1, t2] <- CM[t1, t2] + 1
  }

  # Information measures, Golan (1988)
  if (p0_is_uniform == TRUE) {
    Sp <- -sum(p * log(p)) / (N * log(J))
  } else {
    Sp <-  sum(p * log(p)) / sum(p0 * log(p0))
  }

  S_p_i <- rep(0, N)
  for (i in 1:N) {
    if (p0_is_uniform == TRUE) {
      S_p_i[i] <- -sum(p[i, ] * log(p[i, ])) / log(J)
    } else {
      S_p_i[i] <-  sum(p[i, ] * log(p[i, ])) / sum(p0[i, ] * log(p0[i, ]))
    }
  }

  # Error bounds based on the Fano's [weaker] inequality
  p_e_i <- rep(0, N)
  for (i in 1:N) {
    p_e_i[i] <- S_p_i[i] - log(J)
  }

  # Entropy ratio Statistic
  if (p0_is_uniform == TRUE) {
    ER <- 2 * N * log(J) * (1 - Sp)
  } else {
    ER <- 2 * (-sum(p0[i, ] * log(p0[i, ]))) * (1 - Sp)
  }

  # Pseudo R-squared
  R2 <- 1 - Sp

  info_estim_all <- list("lambda" = lambda, "hess" = gce_optim$hessian, "p" = p, "w" = w, "e" = e,
                         "marg_eff" = ave_dp_dx,"Sp" = Sp, "S_p_i" = S_p_i, "p_e_i" = p_e_i,
                         "H_p_w" = - gce_optim$value, "ER" = ER, "Pseudo_R2" = R2,
                         "CM" = CM, "conv" = gce_optim$convergence)

  return(info_estim_all)

}

gce_mult_obj <- function(lambda_vector, Y, X, v, nu, p0, w0, N, K, J, M) {

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

gce_mult_grad <- function(lambda_vector, Y, X, v, nu, p0, w0, N, K, J, M) {

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
      e[i, j]   <- sum(v * w[i, j, ])
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

