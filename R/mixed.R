#' Mixed Discrete Choice Model
#'
#' @param Y Dependent Variable: an (NxJ) matrix where N represents the number of individuals and J is the number of choices
#' @param X Independent Variable(s): (NxJxK) where J represents the number of choices, K represents the number of independent variables
#' @param dimV An optional argument (scalar) representing dimensions of the support space for the error terms, default value is 5.
#' @param optim_method Optional: same as the "method" argument for the "optim" function in "stats"
#'
#' @return A list which includes estimated Largange Multipliers (LMs), Hessian matrix associated with LMs
#'
#' @export
#'
#' @examples
gme_mixed <- function(Y, X1, X2, dimS, dimV, optim_method = "BFGS") {
  dY <- dim(Y)
  N  <- dY[1]
  J  <- dY[2]
  K1 <- dim(X1)[2]
  K2 <- dim(X2)[3]

  if (missing(dimS)) {
    M <- 3
  } else {
    M <- dimS
  }
  s <- matrix(seq(from = -1, to = 1, length.out = M), nrow = 1)

  if (missing(dimV)) {
    H <- 5
  } else {
    H <- dimV
  }
  v <- matrix(seq(from = -1, to = 1, length.out = H), nrow = 1)

  param0 <- rep(0, (K1 * (J - 1) + K2 + N))
  gce_optim <- optim(param0, gme_mixed_obj, gme_mixed_grad,
                     Y = Y, X1 = X1, X2 = X2, s = s, v = v,
                     N = N, J = J, K1 = K1, K2 = K2, M = M,
                     H = H, method = optim_method)
  lambda_hat <- matrix(gce_optim$par, K, J - 1)
  lambda_hat <- cbind(rep(0, K), lambda_hat)
  gce_all <- list("par" = lambda_hat, "convergence" = gce_optim$convergence,
                  "H" = - gce_optim$value, "Y" = Y, "X" = X, "v" = v, "nu" = nu,
                  "p0"= p0, "w0" = w0, "type" = "gce")
  return(gce_all)

}

summary_mixed <- function(gme_mult_results) {
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

.gme_mixed_obj <- function(param, Y, X1, X2, s, v, N, J, K1, K2, M, H) {

  lambda <- param[1:K1 * (J - 1)]
  lambda <- matrix(lambda, K1, (J - 1))
  lambda <- cbind(rep(0, K1), lambda)
  alpha  <- param((K1 * (J - 1) + 1):(K1 * (J - 1) + K2))
  rho    <- param[-1:(K1 * (J - 1) + K2)]

  temp1  <- X1 %*% lambda
  temp2  <- matrix(0, N, J)
  for (i in 1:N) {
    for (j in 1:J) {
      temp2[i, j] <- sum(alpha * X[i, j, ])
    }
  }

  phi    <- array(0, c(N, J, M))
  Omega  <- matrix(0, N, J)
  for (i in 1:N) {
    for (j in 1:J) {
      for (m in 1:M) {
        phi[i, j, m] <- exp(s[m] * (- temp1[i, j] - temp2[i, j] - rho[i]))
      }
      Omega[i, j] <- sum(phi[i, j, ])
    }
  }

  w   <- array(0, c(N, J, H))
  Psi <- matrix(0, N, J)
  for (i in 1:N) {
    for (j in 1:J) {
      for (h in 1:H) {
        w[i, j, h] <- exp(v[h] * (- temp1[i, j] - temp2[i, j]))
      }
      Psi[i, j] <- sum(w[i, j, ])
    }
  }

  temp3 <- rep(0, K2)
  for (k in 1:K2) {
    temp3[k] <- sum(Y * X2[, , k])
  }
  like  <- sum(lambda * (t(X1) %*% Y)) + sum(alpha * temp3) + sum(log(Omega)) + sum(log(Psi))
  return(like)

}

.gme_mixed_grad <- function(param, Y, X1, X2, s, v, N, J, K1, K2, M, H) {

  lambda <- param[1:K1 * (J - 1)]
  lambda <- matrix(lambda, K1, (J - 1))
  lambda <- cbind(rep(0, K1), lambda)
  alpha  <- param((K1 * (J - 1) + 1):(K1 * (J - 1) + K2))
  rho    <- param[-1:(K1 * (J - 1) + K2)]

  temp1  <- X1 %*% lambda
  temp2  <- matrix(0, N, J)
  for (i in 1:N) {
    for (j in 1:J) {
      temp2[i, j] <- sum(alpha * X2[i, j, ])
    }
  }

  phi    <- array(0, c(N, J, M))
  p      <- matrix(0, N, J)
  Omega  <- matrix(0, N, J)
  for (i in 1:N) {
    for (j in 1:J) {
      for (m in 1:M) {
        phi[i, j, m] <- exp(s[m] * (- temp1[i, j] - temp2[i, j] - rho[i]))
      }
      Omega[i, j] <- sum(phi[i, j, ])
      phi[i, j, ] <- phi[i, j, ] / Omega[i, j]
      p[i, j]     <- sum(s * phi[i, j, ])
    }
  }

  w   <- array(0, c(N, J, H))
  e   <- matrix(0, N, J)
  Psi <- matrix(0, N, J)
  for (i in 1:N) {
    for (j in 1:J) {
      for (h in 1:H) {
        w[i, j, h] <- exp(v[h] * (- temp1[i, j] - temp2[i, j]))
      }
      Psi[i, j] <- sum(w[i, j, ])
      w[i, j, ] <- w[i, j, ] / Psi[i, j]
      e[i, j]   <- sum(v * w[i, j, ])
    }
  }

  dparam <- rep(0, (K1 * (J -1) + K2 + N))
  index_param <- 1
  for (j in 2:J) {
    for (k in 1:K1) {
      dparam[index_param] <- sum(Y[, j] * X1[, k]) -
        sum(p[, j] * X1[, k]) -
        sum(e[, j] * X1[, k])
      index_param <- index_param + 1
    }
  }

  for (k in 1:K2) {
    dparam[index_param] <- sum(Y * X2[, , k]) - sum(p * X2[, , k]) - sum(e * X2[, , k])
    index_param <- index_param + 1
  }

  for (i in 1:N) {
    dparam[index_param] <- 1 - sum(p[i, ])
    index_param <- index_param + 1
  }
  return(dparam)

}

