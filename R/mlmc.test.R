# This code is derived and adapted from the original GPL-2 Matlab version by
# Mike Giles.  See https://people.maths.ox.ac.uk/~gilesm/mlmc/
# There is one key change, whereby alpha, beta, and gamma can be specified if
# desired. If not specified then the behaviour is the same as the original
# Matlab code.

#' Multi-level Monte Carlo estimation test suite
#'
#' Computes a suite of diagnostic values for an MLMC estimation problem.
#'
#' See one of the example level sampler functions (e.g. \code{\link[=opre_l]{opre_l()}}) for example usage.
#'
#' This function is based on GPL-2 'Matlab' code by Mike Giles.
#'
#' @param mlmc_l
#'        a user supplied function which provides the estimate for level \eqn{l}.
#'        It must take at least two arguments, the first is the level number to be simulated and the second the number of paths.
#'        Additional arguments can be taken if desired: all additional \code{...} arguments to this function are forwarded to the user defined \code{mlmc_l} function.
#'
#'        The user supplied function should return a named list containing one element named \code{sums} and second named \code{cost}, where:
#'        \describe{
#'          \item{\code{sums}}{is a vector of length six \eqn{\left(\sum Y_i, \sum Y_i^2, \sum Y_i^3, \sum Y_i^4, \sum X_i, \sum X_i^2\right)} where \eqn{Y_i} are iid simulations with expectation \eqn{E[P_0]} when \eqn{l=0} and expectation \eqn{E[P_l-P_{l-1}]} when \eqn{l>0}, and \eqn{X_i} are iid simulations with expectation \eqn{E[P_l]}.
#'                             Note that this differs from the main \code{\link[=mlmc]{mlmc()}} driver, which only requires the first two of these elements in order to calculate the estimate.
#'                             The remaining elements are required by \code{mlmc.test()} since they are used for convergence tests, kurtosis, and telescoping sum checks.}
#'          \item{\code{cost}}{is a scalar with the total cost of the paths simulated.
#'                             For example, in the financial options samplers included in this package, this is calculated as \eqn{NM^l}, where \eqn{N} is the number of paths requested in the call to the user function \code{mlmc_l}, \eqn{M} is the refinement cost factor (\eqn{M=2} for \code{\link[=mcqmc06_l]{mcqmc06_l()}} and \eqn{M=4} for \code{\link[=opre_l]{opre_l()}}), and \eqn{l} is the level being sampled.}
#'        }
#'
#'        See the function (and source code of) \code{\link[=opre_l]{opre_l()}} and \code{\link[=mcqmc06_l]{mcqmc06_l()}} in this package for an example of user supplied level samplers.
#' @param N
#'        number of samples to use in convergence tests, kurtosis, telescoping sum check.
#' @param L
#'        number of levels to use in convergence tests, kurtosis, telescoping sum check.
#' @param N0
#'        initial number of samples which are used for the first 3 levels and for any subsequent levels which are automatically added in the complexity tests.
#'        Must be \eqn{> 0}.
#' @param eps.v
#'        a vector of one or more target accuracies for the complexity tests.
#'        Must all be \eqn{> 0}.
#' @param Lmin
#'        the minimum level of refinement for complexity tests.
#'        Must be \eqn{\ge 2}.
#' @param Lmax
#'        the maximum level of refinement for complexity tests.
#'        Must be \eqn{\ge} \code{Lmin}.
#' @param alpha
#'        the weak error, \eqn{O(2^{-\alpha l})}.
#'        Must be \eqn{> 0} if specified.
#'        If \code{NA} then \code{alpha} will be estimated.
#' @param beta
#'        the variance, \eqn{O(2^{-\beta l})}.
#'        Must be \eqn{> 0} if specified.
#'        If \code{NA} then \code{beta} will be estimated.
#' @param gamma
#'        the sample cost, \eqn{O(2^{\gamma l})}.
#'        Must be \eqn{> 0} if specified.
#'        If \code{NA} then \code{gamma} will be estimated.
#' @param parallel
#'        if an integer is supplied, R will fork \code{parallel} parallel processes.
#'        This is done for the convergence tests section by splitting the \code{N} samples as evenly as possible across cores when sampling each level.
#'        This is also done for the MLMC complexity tests by passing the \code{parallel} argument on to the \code{\link[=mlmc]{mlmc()}} driver when targeting each accuracy level in \code{eps}.
#' @param silent
#'        set to TRUE to supress running output (identical output can still be printed by printing the return result)
#' @param ...
#'        additional arguments which are passed on when the user supplied \code{mlmc_l} function is called
#'
#' @return An \code{mlmc.test} object which contains all the computed diagnostic values.
#'         This object can be printed or plotted (see \code{\link{plot.mlmc.test}}).
#'
#' @author Louis Aslett <louis.aslett@durham.ac.uk>
#' @author Mike Giles <Mike.Giles@maths.ox.ac.uk>
#' @author Tigran Nagapetyan <nagapetyan@stats.ox.ac.uk>
#'
#' @examples
#' \donttest{
#' # Example calls with realistic arguments
#' # Financial options using an Euler-Maruyama discretisation
#' tst <- mlmc.test(opre_l, N = 2000000,
#'                  L = 5, N0 = 1000,
#'                  eps.v = c(0.005, 0.01, 0.02, 0.05, 0.1),
#'                  Lmin = 2, Lmax = 6,
#'                  option = 1)
#' tst
#' plot(tst)
#'
#' # Financial options using a Milstein discretisation
#' tst <- mlmc.test(mcqmc06_l, N = 20000,
#'                  L = 8, N0 = 200,
#'                  eps.v = c(0.005, 0.01, 0.02, 0.05, 0.1),
#'                  Lmin = 2, Lmax = 10,
#'                  option = 1)
#' tst
#' plot(tst)
#' }
#'
#' # Toy versions for CRAN tests
#' tst <- mlmc.test(opre_l, N = 10000,
#'                  L = 5, N0 = 1000,
#'                  eps.v = c(0.025, 0.1),
#'                  Lmin = 2, Lmax = 6,
#'                  option = 1)
#'
#' tst <- mlmc.test(mcqmc06_l, N = 10000,
#'                  L = 8, N0 = 1000,
#'                  eps.v = c(0.025, 0.1),
#'                  Lmin = 2, Lmax = 10,
#'                  option = 1)
#'
#' @importFrom stats lm
#' @export
mlmc.test <- function(mlmc_l, N, L, N0, eps.v, Lmin, Lmax, alpha = NA, beta = NA, gamma = NA, parallel = NA, silent = FALSE, ...) {
  if(silent)
    cat <- function(...) { }

  # check parameters are acceptable
  if(N<=0) {
    stop("N must be > 0.")
  }
  if(L<0) {
    stop("L must be >= 0.")
  }
  if(Lmin<2) {
    stop("Lmin must be >= 2.")
  }
  if(Lmax<Lmin) {
    stop("must have Lmax >= Lmin.")
  }
  if(N0<=0 || any(eps.v<=0)) {
    stop("N0 and eps.v must be greater than zero.")
  }
  if(!is.na(alpha) && alpha<=0) {
    stop("if specified, alpha must be greater than zero.  Set alpha to NA to automatically estimate.")
  }
  if(!is.na(beta) && beta<=0) {
    stop("if specified, beta must be greater than zero.  Set beta to NA to automatically estimate.")
  }
  if(!is.na(gamma) && gamma<=0) {
    stop("if specified, gamma must be greater than zero.  Set gamma to NA to automatically estimate.")
  }
  if(!is.na(parallel) && parallel<=0) {
    stop("if specified, parallel must be greater than zero.")
  }

  # check user supplied mlmc function on fastest level for just 1 simulation is
  # returning the correct shape of list(sums, cost)
  res.tst <- mlmc_l(0, 1, ...)
  if(!is.list(res.tst)) {
    stop("user suppled mlmc sampler must return a list.")
  } else if(!identical(sort(names(res.tst)), c("cost", "sums"))) {
    stop("user suppled mlmc sampler must return a list with elements named sums and cost.")
  } else if(!is.numeric(res.tst$sums) || length(res.tst$sums) != 6) {
    stop("sums returned by user suppled mlmc sampler must be length 6.")
  } else if(!is.numeric(res.tst$cost) || length(res.tst$cost) != 1) {
    stop("cost returned by user suppled mlmc sampler must be length 1.")
  }

  res <- within(list(), {
    N <- N
    L <- L
    eps.v <- eps.v

    cat("\n")
    cat("**********************************************************\n")
    cat("*** Convergence tests, kurtosis, telescoping sum check ***\n")
    cat("*** using N =", N, "samples    ", paste0(rep(' ', 21 - nchar(as.character(N))), collapse = ''), "      ***\n")
    cat("**********************************************************\n")
    cat("\n l   ave(Pf-Pc)    ave(Pf)   var(Pf-Pc)    var(Pf)")
    cat("    kurtosis      check       cost\n-------------------------")
    cat("--------------------------------------------------------------\n")

    del1 <- del2 <- var1 <- var2 <- kur1 <- chk1 <- cost <- c()

    for(l in 0:L) {
      sums <- 0
      cst  <- 0
      if(is.na(parallel)) {
        res  <- mlmc_l(l, N, ...)
        sums <- sums + res$sums/N
        cst  <- cst  + res$cost/N
      } else {
        dN <- c(rep(N %/% parallel + 1, N %% parallel), rep(N %/% parallel, parallel - N %% parallel)) # Spread simulation of N samples over parallel slots
        res <- mcmapply(function(l, dN, ...) {
          x <- mlmc_l(l, dN, ...)
          c(x$sums, x$cost)
        }, l = rep(l, parallel), dN = dN, ..., mc.preschedule = FALSE, mc.cores = parallel)
        sums <- sums + rowSums(res[-nrow(res),]/N)
        cst  <- cst  +     sum(res[ nrow(res),]/N)
      }
      if(l==0) {
        kurt <- 0.0
      } else {
        kurt <- (sums[4] -
                   4*sums[3]*sums[1] +
                   6*sums[2]*sums[1]^2 -
                   3*sums[1]*sums[1]^3 ) /
          (sums[2]-sums[1]^2)^2
      }

      cost <- c(cost, cst)
      del1 <- c(del1, sums[1])
      del2 <- c(del2, sums[5])
      var1 <- c(var1, sums[2]-sums[1]^2)
      var1 <- pmax(var1, 1e-10) # fix for cases with var=0
      var2 <- c(var2, sums[6]-sums[5]^2)
      var2 <- pmax(var2, 1e-10) # fix for cases with var=0
      kur1 <- c(kur1, kurt)

      if(l==0) {
        check <- 0
      } else {
        check <- abs(del1[l+1] + del2[l] - del2[l+1]) /
          (3.0*(sqrt(var1[l+1]) + sqrt(var2[l]) + sqrt(var2[l+1]) )/sqrt(N))
      }
      chk1 <- c(chk1, check)

      cat(sprintf("%2d   %8.4e  %8.4e  %8.4e  %8.4e  %8.4e  %8.4e  %8.4e \n",
                  l, del1[l+1], del2[l+1], var1[l+1], var2[l+1], kur1[l+1], chk1[l+1], cost[l+1]))
    }

    # print out a warning if kurtosis or consistency check looks bad

    if( kur1[length(kur1)] > 100.0 ) {
      cat(sprintf("\n WARNING: kurtosis on finest level = %f \n", kur1[length(kur1)]))
      cat(" indicates MLMC correction dominated by a few rare paths; \n")
      cat(" for information on the connection to variance of sample variances,\n")
      cat(" see http://mathworld.wolfram.com/SampleVarianceDistribution.html\n\n")
    }

    if( max(chk1) > 1.0 ) {
      cat(sprintf("\n WARNING: maximum consistency error = %f \n", max(chk1)))
      cat(" indicates identity E[Pf-Pc] = E[Pf] - E[Pc] not satisfied; \n")
      cat(" to be more certain, re-run mlmc.test with larger N \n\n")
    }

    # use linear regression to estimate alpha, beta and gamma
    # if they are not specified

    L1 <- 2
    L2 <- L+1
    msg <- NULL
    if(is.na(alpha)) {
      alpha <- max(0.5, -lm(y~x, data.frame(y=log2(abs(del1[L1:L2])), x=L1:L2))$coefficients["x"])
      msg <- c(msg, "alpha")
    }
    if(is.na(beta)) {
      beta <- max(0.5, -lm(y~x, data.frame(y=log2(abs(var1[L1:L2])), x=L1:L2))$coefficients["x"])
      msg <- c(msg, "beta")
    }
    if(is.na(gamma)) {
      gamma <- max(0.5, lm(y~x, data.frame(y=log2(abs(cost[L1:L2])), x=L1:L2))$coefficients["x"])
      msg <- c(msg, "gamma")
    }
    if(is.null(msg)) {
      msg <- "All MLMC parameters specified by user"
    } else {
      if(length(msg) == 3) {
        msg <- "Linear regression estimates of MLMC parameters"
      } else {
        msg <- paste0(paste0(msg, collapse = ", "),
                      " estimated by linear regression; ",
                      paste0(setdiff(c("alpha", "beta", "gamma"), msg), collapse = ", "),
                      " user supplied")
      }
    }

    cat("\n****", paste0(rep("*", nchar(msg)), collapse = ""), "****\n", sep = "")
    cat("*** ", msg, " ***\n", sep = "")
    cat("****", paste0(rep("*", nchar(msg)), collapse = ""), "****\n", sep = "")
    cat(sprintf("\n alpha = %f  (exponent for MLMC weak convergence)\n", alpha))
    cat(sprintf(" beta  = %f  (exponent for MLMC variance) \n", beta))
    cat(sprintf(" gamma = %f  (exponent for MLMC cost) \n", gamma))

    # second, mlmc complexity tests

    cat("\n")
    cat("***************************** \n")
    cat("*** MLMC complexity tests *** \n")
    cat("***************************** \n\n")
    cat("  eps      value    mlmc_cost   std_cost  savings     N_l \n")
    cat("----------------------------------------------------------- \n")

    alpha <- max(alpha, 0.5)
    beta  <- max(beta, 0.5)
    theta <- 0.25

    P <- mlmc_cost <- std_cost <- c()
    Nl <- list()
    for(i in 1:length(eps.v)) {
      mlmcres <- mlmc(Lmin, Lmax, N0, eps.v[i], mlmc_l, alpha, beta, gamma, parallel, ...)
      P <- c(P, mlmcres$P)
      Nl[[i]] <- mlmcres$Nl
      mlmc_cost <- c(mlmc_cost, sum(mlmcres$Nl*mlmcres$Cl))
      std_cost <- c(std_cost, var2[min(length(var2), length(mlmcres$Nl))]*mlmcres$Cl[length(mlmcres$Cl)] / ((1.0-theta)*eps.v[i]^2))

      cat(sprintf("%.4f  %.4e  %.3e  %.3e  %7.2f ",
                  eps.v[i], P[i], mlmc_cost[i], std_cost[i], std_cost[i]/mlmc_cost[i]))
      cat(sprintf("%9d", Nl[[i]]))
      cat("\n")
    }

    cat("\n")
    rm(mlmcres, res, i, theta, L1, L2, check, kurt, l, sums, cst)
  })
  class(res) <- "mlmc.test"
  invisible(res)
}
