#' Mixed Discrete Choice Model
#'
#' @param Y  Dependent Variable: an (NxJ) matrix where N represents the number of individuals and J is the number of choices
#' @param X1 Alternative-invariant regressors: (NxK1) matrix where K1 represents the number of regressors
#' @param X2 Alternative-variant regressors: (NxJxK2) matrix where K2 represents the number of regressors
#' @param dimS An optional argument (scalar) representing dimensions of the support space for the phi probabilities, default value is 3
#' @param dimV An optional argument (scalar) representing dimensions of the support space for the error terms, default value is 5
#' @param optim_method Optional: same as the "method" argument for the "optim" function in "stats"
#'
#' @return A list which has the following elements:
#' \itemize{
#'   \item lambda - Lagrange multipliers associated with the alternative-invariant regressors.
#'   \item alpha - Lagrange multipliers associated with the alternative-variant regressors.
#'   \item rho - Lagrange multipliers associated with the normalization constraint.
#'   \item hess - Hessian matrix associated with the estimated Lagrange multipliers.
#'   \item H_phi_w - Value of the joint entropies of \code{phi} and \code{w} at the final iteration.
#'   \item marg_eff_x1 - Marginal effects associated with the alternative-invariant regressors.
#'   \item marg_eff_x2 - Marginal effects associated with the alternative-variant regressors.
#'   \item S1 - The (signal) information of the whole system.
#'   \item S2 - The (signal) information associated with the choice \code{j} of the individual \code{i}.
#'   \item Pseudo_R2 - Pseudo R-squared.
#'   \item ER - Entropy Ratio Statistic.
#'   \item conv - Convergence (details in the help file of the optim function).
#' }
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
  gme_optim <- optim(param0, gme_mixed_obj, gme_mixed_grad,
                     Y = Y, X1 = X1, X2 = X2, s = s, v = v,
                     N = N, J = J, K1 = K1, K2 = K2, M = M,
                     H = H, method = optim_method, hessian = T)
  param <- gme_optim$par
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

  dphi_dx1 <- array(0, c(N, J, M, K1))
  for (i in 1:N) {
    for (j in 1:J) {
      for (m in 1:M) {
        for (k in 1:K1) {
          dphi_dx1[i, j, m, k] <- phi[i, j, m] * lambda[k, j] * (s[m] - sum(s * phi[i, j, ]))
        }
      }
    }
  }

  dp_dx1 <- array(0, c(N, J, K1))
  for (i in 1:N) {
    for (j in 1:J) {
      for (k in 1:K1) {
        dp_dx1[i, j, k] <- sum(s * dphi_dx1[i, j, , k])
      }
    }
  }
  me_x1 <- apply(dp_dx, c(2, 3), sum)

  dphi_dx2 <- array(0, c(N, J, M, K2))
  for (i in 1:N) {
    for (j in 1:J) {
      for (m in 1:M) {
        for (k in 1:K2) {
          dphi_dx2[i, j, m, k] <- phi[i, j, m] * alpha[k] * (s[m] - sum(s * phi[i, j, ]))
        }
      }
    }
  }

  dp_dx2 <- array(0, c(N, J, K2))
  for (i in 1:N) {
    for (j in 1:J) {
      for (k in 1:K2) {
        dp_dx2[i, j, k] <- sum(s * dphi_dx2[i, j, , k])
      }
    }
  }
  me_x2 <- apply(dp_dx, c(3), sum)

  H_phi <-  phi * log(phi)
  S1    <- -sum(H_phi) / (N * J * log(M))
  S2    <- -apply(H_phi, 3, sum) / log(M)
  R2    <- 1 - S1
  ER    <- 2 * N * J * log(M) * R2

  gme_all <- list("lambda" = lambda, "alpha" = alpha, "rho" = rho,
                  "hess" = gme_optim$hessian, "H_phi_w" = gme_optim$value,
                  "marg_eff_x1" = me_x1, "marg_eff_x2" = me_x2,
                  "S1" = S1, "S2" = S2, "Pseudo_R2" = R2,
                  "Entropy_ratio" = ER, "convergence" = gme_optim$convergence)

  return(gme_all)

}

gme_mixed_obj <- function(param, Y, X1, X2, s, v, N, J, K1, K2, M, H) {

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

gme_mixed_grad <- function(param, Y, X1, X2, s, v, N, J, K1, K2, M, H) {

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

