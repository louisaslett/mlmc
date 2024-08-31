# This code is derived and adapted from the original GPL-2 Matlab version by
# Mike Giles.  See https://people.maths.ox.ac.uk/~gilesm/mlmc/

sig_dW <- function(x, dW, h) {
  dW[2,] <- -0.5*dW[1,] + sqrt(0.75)*dW[2,]

  rbind(sqrt(pmax(0,x[2,]))*x[1,]*dW[1,],
        exp(-5*h)*0.25*sqrt(pmax(0,x[2,]))*dW[2,]);
}

mu <- function(x, h) {
  rbind(0.05*x[1,],
        ((1-exp(-5*h))/h)*(0.04-x[2,]))
}

#' Financial options using an Euler-Maruyama discretisation
#'
#' Financial options based on scalar geometric Brownian motion and Heston models, similar to Mike Giles' original 2008 Operations Research paper, Giles (2008), using an Euler-Maruyama discretisation
#'
#' This function is based on GPL-2 'Matlab' code by Mike Giles.
#'
#' @param l
#'        the level to be simulated.
#' @param N
#'        the number of samples to be computed.
#' @param option
#'        the option type, between 1 and 5.
#'        The options are:
#'        \describe{
#'          \item{1 = European call;}{}
#'          \item{2 = Asian call;}{}
#'          \item{3 = lookback call;}{}
#'          \item{4 = digital call;}{}
#'          \item{5 = Heston model.}{}
#'        }
#'
#' @return A named list containing: \describe{
#'           \item{\code{sums}}{is a vector of length six \eqn{\left(\sum Y_i, \sum Y_i^2, \sum Y_i^3, \sum Y_i^4, \sum X_i, \sum X_i^2\right)} where \eqn{Y_i} are iid simulations with expectation \eqn{E[P_0]} when \eqn{l=0} and expectation \eqn{E[P_l-P_{l-1}]} when \eqn{l>0}, and \eqn{X_i} are iid simulations with expectation \eqn{E[P_l]}.
#'                              Note that only the first two components of this are used by the main \code{\link[=mlmc]{mlmc()}} driver, the full vector is used by \code{\link[=mlmc.test]{mlmc.test()}} for convergence tests etc;}
#'           \item{\code{cost}}{is a scalar with the total cost of the paths simulated, computed as \eqn{N \times 4^l} for level \eqn{l}.}
#'         }
#'
#' @author Louis Aslett <louis.aslett@durham.ac.uk>
#' @author Mike Giles <Mike.Giles@maths.ox.ac.uk>
#' @author Tigran Nagapetyan <nagapetyan@stats.ox.ac.uk>
#'
#' @references
#' Giles, M.B. (2008) 'Multilevel Monte Carlo Path Simulation', \emph{Operations Research}, 56(3), pp. 607–617. Available at: \doi{10.1287/opre.1070.0496}.
#'
#' @examples
#' \donttest{
#' # These are similar to the MLMC tests for the original
#' # 2008 Operations Research paper, using an Euler-Maruyama
#' # discretisation with 4^l timesteps on level l.
#' #
#' # The differences are:
#' # -- the plots do not have the extrapolation results
#' # -- two plots are log_2 rather than log_4
#' # -- the new MLMC driver is a little different
#' # -- switch to X_0=100 instead of X_0=1
#' #
#' # Note the following takes quite a while to run, for a toy example see after
#' # this block.
#'
#' N0   <- 1000 # initial samples on coarse levels
#' Lmin <- 2 # minimum refinement level
#' Lmax <- 6 # maximum refinement level
#'
#' test.res <- list()
#' for(option in 1:5) {
#'   if(option == 1) {
#'     cat("\n ---- Computing European call ---- \n")
#'     N      <- 1000000 # samples for convergence tests
#'     L      <- 5 # levels for convergence tests
#'     Eps    <- c(0.005, 0.01, 0.02, 0.05, 0.1)
#'   } else if(option == 2) {
#'     cat("\n ---- Computing Asian call ---- \n")
#'     N      <- 1000000 # samples for convergence tests
#'     L      <- 5 # levels for convergence tests
#'     Eps    <- c(0.005, 0.01, 0.02, 0.05, 0.1)
#'   } else if(option == 3) {
#'     cat("\n ---- Computing lookback call ---- \n")
#'     N      <- 1000000 # samples for convergence tests
#'     L      <- 5 # levels for convergence tests
#'     Eps    <- c(0.01, 0.02, 0.05, 0.1, 0.2)
#'   } else if(option == 4) {
#'     cat("\n ---- Computing digital call ---- \n")
#'     N      <- 4000000 # samples for convergence tests
#'     L      <- 5 # levels for convergence tests
#'     Eps    <- c(0.02, 0.05, 0.1, 0.2, 0.5)
#'   } else if(option == 5) {
#'     cat("\n ---- Computing Heston model ---- \n")
#'     N      <- 2000000 # samples for convergence tests
#'     L      <- 5 # levels for convergence tests
#'     Eps    <- c(0.005, 0.01, 0.02, 0.05, 0.1)
#'   }
#'
#'   test.res[[option]] <- mlmc.test(opre_l, N, L, N0, Eps, Lmin, Lmax, option = option)
#'
#'   # print exact analytic value, based on S0=K
#'   T   <- 1
#'   r   <- 0.05
#'   sig <- 0.2
#'   K   <- 100
#'
#'   k   <- 0.5*sig^2/r;
#'   d1  <- (r+0.5*sig^2)*T / (sig*sqrt(T))
#'   d2  <- (r-0.5*sig^2)*T / (sig*sqrt(T))
#'
#'   if(option == 1) {
#'     val <- K*( pnorm(d1) - exp(-r*T)*pnorm(d2) )
#'   } else if(option == 2) {
#'     val <- NA
#'   } else if(option == 3) {
#'     val <- K*( pnorm(d1) - pnorm(-d1)*k - exp(-r*T)*(pnorm(d2) - pnorm(d2)*k) )
#'   } else if(option == 4) {
#'     val <- K*exp(-r*T)*pnorm(d2)
#'   } else if(option == 5) {
#'     val <- NA
#'   }
#'
#'   if(is.na(val)) {
#'     cat(sprintf("\n Exact value unknown, MLMC value: %f \n", test.res[[option]]$P[1]))
#'   } else {
#'     cat(sprintf("\n Exact value: %f, MLMC value: %f \n", val, test.res[[option]]$P[1]))
#'   }
#'
#'   # plot results
#'   plot(test.res[[option]])
#' }
#' }
#'
#' # The level sampler can be called directly to retrieve the relevant level sums:
#' opre_l(l = 7, N = 10, option = 1)
#'
#' @importFrom stats rnorm
#' @export
opre_l <- function(l, N, option) {
  M <- 4

  T   <- 1
  r   <- 0.05
  sig <- 0.2
  K   <- 100

  nf <- M^l
  nc <- nf/M

  hf <- T/nf
  hc <- T/nc

  sums <- rep(0, 6)

  for(N1 in seq(1, N, by=10000)) {
    N2 <- min(10000, N-N1+1)

    #
    # GBM model
    #
    if(option<5) {
      X0 <- K

      Xf <- rep(X0, N2)
      Xc <- Xf

      Af <- 0.5*hf*Xf
      Ac <- 0.5*hc*Xc

      Mf <- Xf
      Mc <- Xc

      if(l==0) {
        dWf <- sqrt(hf)*rnorm(N2)
        Xf  <- Xf + r*Xf*hf + sig*Xf*dWf
        Af <- Af + 0.5*hf*Xf
        Mf <- pmin(Mf,Xf)
      } else {
        for (n in 1:nc){
          dWc <- rep(0, N2)
          for (m in 1:M){
            dWf <- sqrt(hf)*rnorm(N2)
            dWc <- dWc + dWf
            Xf  <- Xf + r*Xf*hf + sig*Xf*dWf
            Af  <- Af + hf*Xf
            Mf  <- pmin(Mf,Xf)
          }
          Xc <- Xc + r*Xc*hc + sig*Xc*dWc
          Ac <- Ac + hc*Xc
          Mc <- pmin(Mc,Xc)
        }
        Af <- Af - 0.5*hf*Xf
        Ac <- Ac - 0.5*hc*Xc
      }

      if(option==1) {
        Pf <- pmax(0,Xf-K)
        Pc <- pmax(0,Xc-K)
      } else if(option==2) {
        Pf <- pmax(0,Af-K)
        Pc <- pmax(0,Ac-K)
      } else if(option==3) {
        beta <- 0.5826 # special factor for offset correction
        Pf <- Xf - Mf*(1-beta*sig*sqrt(hf))
        Pc <- Xc - Mc*(1-beta*sig*sqrt(hc))
      } else if(option==4) {
        Pf <- K * 0.5*(sign(Xf-K)+1)
        Pc <- K * 0.5*(sign(Xc-K)+1)
      }
      #
      # Heston model
      #
    } else {
      Xf <- matrix(c(K, 0.04), nrow=2, ncol=N2)
      Xc <- Xf

      if(l==0) {
        dWf <- sqrt(hf)*matrix(rnorm(2*N2), nrow=2, ncol=N2)
        Xf  <- Xf + mu(Xf,hf)*hf + sig_dW(Xf,dWf,hf)
      } else {
        for(n in 1:nc) {
          dWc <- matrix(0, nrow=2, ncol=N2)
          for(m in 1:M) {
            dWf <- sqrt(hf)*matrix(rnorm(2*N2), nrow=2, ncol=N2)
            dWc <- dWc + dWf
            Xf  <- Xf + mu(Xf,hf)*hf + sig_dW(Xf,dWf,hf)
          }
          Xc <- Xc + mu(Xc,hc)*hc + sig_dW(Xc,dWc,hc)
        }
      }

      Pf <- pmax(0,Xf[1,]-K)
      Pc <- pmax(0,Xc[1,]-K)
    }

    Pf <- exp(-r*T)*Pf
    Pc <- exp(-r*T)*Pc

    if(l==0) {
      Pc <- 0
    }

    sums[1] <- sums[1] + sum(Pf-Pc)
    sums[2] <- sums[2] + sum((Pf-Pc)^2)
    sums[3] <- sums[3] + sum((Pf-Pc)^3)
    sums[4] <- sums[4] + sum((Pf-Pc)^4)
    sums[5] <- sums[5] + sum(Pf)
    sums[6] <- sums[6] + sum(Pf^2)
  }

  cost <- N*nf # cost defined as number of fine timesteps

  list(sums = sums, cost = cost)
}
