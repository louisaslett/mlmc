# This code is derived and adapted from the original GPL-2 Matlab version by
# Mike Giles.  See https://people.maths.ox.ac.uk/~gilesm/mlmc/

#' Multi-level Monte Carlo estimation
#'
#' This function is the Multi-level Monte Carlo driver which will sample from the levels of user specified function.
#'
#' The Multilevel Monte Carlo Method method originated in the works Giles (2008) and Heinrich (1998).
#'
#' Consider a sequence \eqn{P_0, P_1, \ldots}, which approximates \eqn{P_L} with increasing accuracy, but also increasing cost, we have the simple identity
#' \deqn{E[P_L] = E[P_0] + \sum_{l=1}^L E[P_l-P_{l-1}],}
#' and therefore we can use the following unbiased estimator for \eqn{E[P_L]},
#' \deqn{N_0^{-1} \sum_{n=1}^{N_0} P_0^{(0,n)} + \sum_{l=1}^L \left\{ N_l^{-1} \sum_{n=1}^{N_l} \left(P_l^{(l,n)} - P_{l-1}^{(l,n)}\right) \right\}}
#' where \eqn{N_l} samples are produced at level \eqn{l}.
#' The inclusion of the level \eqn{l} in the superscript \eqn{(l,n)} indicates that the samples used at each level of correction are independent.
#'
#' Set \eqn{C_0}, and \eqn{V_0} to be the cost and variance of one sample of \eqn{P_0}, and \eqn{C_l, V_l} to be the cost and variance of one sample of \eqn{P_l - P_{l-1}}, then the overall cost and variance of the multilevel estimator is \eqn{\sum_{l=0}^L N_l C_l} and \eqn{\sum_{l=0}^L N_l^{-1} V_l}, respectively.
#'
#' The idea behind the method, is that provided that the product \eqn{V_l C_l} decreases with \eqn{l}, i.e. the cost increases with level slower than the variance decreases, then one can achieve significant computational savings, which can be formalised as in Theorem 1 of Giles (2015).
#'
#' For further information on multilevel Monte Carlo methods, see the webpage \url{https://people.maths.ox.ac.uk/gilesm/mlmc_community.html} which lists the research groups working in the area, and their main publications.
#'
#' This function is based on GPL-2 'Matlab' code by Mike Giles.
#'
#' @author Louis Aslett <louis.aslett@durham.ac.uk>
#' @author Mike Giles <Mike.Giles@maths.ox.ac.uk>
#' @author Tigran Nagapetyan <nagapetyan@stats.ox.ac.uk>
#'
#' @references
#' Giles, M.B. (2008) 'Multilevel Monte Carlo Path Simulation', \emph{Operations Research}, 56(3), pp. 607–617. Available at: \doi{10.1287/opre.1070.0496}.
#'
#' Giles, M.B. (2015) 'Multilevel Monte Carlo methods', \emph{Acta Numerica}, 24, pp. 259–328. Available at: \doi{10.1017/S096249291500001X}.
#'
#' Heinrich, S. (1998) 'Monte Carlo Complexity of Global Solution of Integral Equations', \emph{Journal of Complexity}, 14(2), pp. 151–175. Available at: \doi{10.1006/jcom.1998.0471}.
#'
#' @param Lmin
#'        the minimum level of refinement.  Must be \eqn{\ge 2}.
#' @param Lmax
#'        the maximum level of refinement.  Must be \eqn{\ge} \code{Lmin}.
#' @param N0
#'        initial number of samples which are used for the first 3 levels and for any subsequent levels which are automatically added.
#'        Must be \eqn{> 0}.
#' @param eps
#'        the target accuracy of the estimate (root mean square error).
#'        Must be \eqn{> 0}.
#' @param mlmc_l
#'        a user supplied function which provides the estimate for level \eqn{l}.
#'        It must take at least two arguments, the first is the level number to be simulated and the second the number of paths.
#'        Additional arguments can be taken if desired: all additional \code{...} arguments to this function are forwarded to the user defined \code{mlmc_l} function.
#'
#'        The user supplied function should return a named list containing one element named \code{sums} and second named \code{cost}, where:
#'        \describe{
#'          \item{\code{sums}}{is a vector of length at least two.
#'                             The first two elements should be \eqn{\left(\sum Y_i, \sum Y_i^2\right)} where \eqn{Y_i} are iid simulations with expectation \eqn{E[P_0]} when \eqn{l=0} and expectation \eqn{E[P_l-P_{l-1}]} when \eqn{l>0}.
#'                             Note that typically the user supplied level sampler will actually return a vector of length six, also enabling use of the \code{\link[=mlmc.test]{mlmc.test()}} function to perform convergence tests, kurtosis, and telescoping sum checks.
#'                             See \code{\link[=mlmc.test]{mlmc.test()}} for the definition of these remaining four elements.}
#'          \item{\code{cost}}{is a scalar with the total cost of the paths simulated.
#'                             For example, in the financial options samplers included in this package, this is calculated as \eqn{NM^l}, where \eqn{N} is the number of paths requested in the call to the user function \code{mlmc_l}, \eqn{M} is the refinement cost factor (\eqn{M=2} for \code{\link[=mcqmc06_l]{mcqmc06_l()}} and \eqn{M=4} for \code{\link[=opre_l]{opre_l()}}), and \eqn{l} is the level being sampled.}
#'        }
#'
#'        See the function (and source code of) \code{\link[=opre_l]{opre_l()}} and \code{\link[=mcqmc06_l]{mcqmc06_l()}} in this package for an example of user supplied level samplers.
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
#'        if an integer is supplied, R will fork \code{parallel} parallel processes and spread the simulations required at each level as evenly as possible across all cores.
#' @param ...
#'        additional arguments which are passed on when the user supplied \code{mlmc_l} function is called.
#'
#' @return A named list containing: \describe{
#'   \item{\code{P}}{The MLMC estimate;}
#'   \item{\code{Nl}}{A vector of the number of samples performed on each level;}
#'   \item{\code{Cl}}{Per sample cost at each level.}
#' }
#'
#' @examples
#' mlmc(2, 6, 1000, 0.01, opre_l, option = 1)
#'
#' mlmc(2, 10, 1000, 0.01, mcqmc06_l, option = 1)
#'
#' @importFrom parallel mclapply
#' @importFrom stats lm aggregate
#' @export
mlmc <- function(Lmin, Lmax, N0, eps, mlmc_l, alpha = NA, beta = NA, gamma = NA, parallel = NA, ...) {
  # check parameters are acceptable
  if(Lmin<2) {
    stop("Lmin must be >= 2.")
  }
  if(Lmax<Lmin) {
    stop("must have Lmax >= Lmin.")
  }
  if(N0<=0 || eps<=0){
    stop("N0 and eps must be greater than zero.")
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

  # check user supplied mlmc function on fastest level for just 1 simulation is
  # returning the correct shape of list(sums, cost)
  res.tst <- mlmc_l(0, 1, ...)
  if(!is.list(res.tst)) {
    stop("user suppled mlmc sampler must return a list.")
  } else if(!identical(sort(names(res.tst)), c("cost", "sums"))) {
    stop("user suppled mlmc sampler must return a list with elements named sums and cost.")
  } else if(!is.numeric(res.tst$sums) || length(res.tst$sums) < 2) {
    stop("sums returned by user suppled mlmc sampler must be at least length 2.")
  } else if(!is.numeric(res.tst$cost) || length(res.tst$cost) != 1) {
    stop("cost returned by user suppled mlmc sampler must be length 1.")
  }

  # initialise the MLMC run
  est.alpha <- is.na(alpha)
  alpha <- ifelse(is.na(alpha), 0, alpha)
  est.beta  <- is.na(beta)
  beta <- ifelse(is.na(beta), 0, beta)
  est.gamma  <- is.na(gamma)
  gamma <- ifelse(is.na(gamma), 0, gamma)

  theta <- 0.25

  L <- Lmin

  Nl <- rep(0, L+1)
  suml <- matrix(0, nrow=2, ncol=L+1)
  costl <- rep(0, L+1)
  dNl <- rep(N0, L+1)

  while(sum(dNl) > 0) {
    # update sample sums from each level
    if(is.na(parallel)) {
      for(l in 0:L) {
        if(dNl[l+1] > 0) {
          res         <- mlmc_l(l, dNl[l+1], ...)
          sums        <- res$sums
          cost        <- res$cost
          Nl[l+1]     <- Nl[l+1] + dNl[l+1]
          suml[1,l+1] <- suml[1,l+1] + sums[1]
          suml[2,l+1] <- suml[2,l+1] + sums[2]
          costl[l+1]  <- costl[l+1] + cost
        }
      }
    } else if(is.numeric(parallel)) {
      # Figure out how to distribute across cores in the most even way possible
      # by giving the same number of sims on each level to each core
      cores.work <- list()
      for(i in 1:parallel) {
        core.work <- data.frame()
        for(l in 0:L) {
          core.dNl <- c(rep(dNl[l+1] %/% parallel + 1, dNl[l+1] %% parallel), rep(dNl[l+1] %/% parallel, parallel - dNl[l+1] %% parallel))[i] # Spread simulation of dNl samples over parallel slots and extract for this core
          if(core.dNl > 0)
            core.work <- rbind(core.work,
                               data.frame(l = l, dNl = core.dNl))
        }
        cores.work[[i]] <- core.work
      }
      rm(core.work)

      # Verify the above gives the correct overall number of sims
      par.vars <- data.frame(l=0:L, dNl=dNl)[!(dNl==0),]
      row.names(par.vars) <- NULL
      if(!identical(par.vars, aggregate(dNl ~ l, data = do.call("rbind", cores.work), sum)))
        stop("incorrect parallelism split, please file a bug report with example, this should never happen! https://github.com/louisaslett/mlmc/issues")

      # Execute parallel simulations
      res <- mclapply(cores.work, function(core.work, ...) {
        if(nrow(core.work) == 0) return(data.frame())
        res <- data.frame(l = core.work$l, sums1 = NA, sums2 = NA, cost = NA)
        for(i in 1:nrow(core.work)) {
          x <- mlmc_l(core.work[i,"l"], core.work[i,"dNl"], ...)
          res[i,2:4] <- c(x$sums[1], x$sums[2], x$cost)
        }
        res
      }, ..., mc.preschedule = FALSE, mc.cores = parallel)
      # Gather results and combine sums/costs from different cores
      res <- do.call("rbind", res)
      res <- aggregate(cbind(sums1, sums2, cost) ~ l, data = res, sum)
      sums <- unname(t(res[,-c(1,4),drop = FALSE]))
      cost <- unname(res[,4])

      # Update running totals
      Nl <- Nl+dNl
      suml[,res$l+1] <- suml[,res$l+1,drop = FALSE] + sums
      costl[res$l+1] <- costl[res$l+1] + cost
    }

    # compute absolute average, variance and cost
    ml <- abs(    suml[1,]/Nl)
    Vl <- pmax(0, suml[2,]/Nl - ml^2)
    Cl <- costl/Nl

    # fix to cope with possible zero values for ml and Vl
    # (can happen in some applications when there are few samples)
    for(l in 3:(L+1)) {
      ml[l] <- max(ml[l], 0.5*ml[l-1]/2^alpha)
      Vl[l] <- max(Vl[l], 0.5*Vl[l-1]/2^beta)
    }

    # use linear regression to estimate alpha, beta, gamma if not given
    if(est.alpha) {
      alpha <- max(0.5, -lm(y~x, data.frame(y=log2(ml[-1]), x=1:L))$coefficients["x"])
    }
    if(est.beta) {
      beta <- max(0.5, -lm(y~x, data.frame(y=log2(Vl[-1]), x=1:L))$coefficients["x"])
    }
    if(est.gamma) {
      gamma <- max(0.5, lm(y~x, data.frame(y=log2(Cl[-1]), x=1:L))$coefficients["x"])
    }
    # set optimal number of additional samples
    Ns  <- ceiling(sqrt(Vl/Cl) * sum(sqrt(Vl*Cl)) / ((1-theta)*eps^2))
    dNl <- pmax(0, Ns-Nl)
    # if (almost) converged, estimate remaining error and decide
    # whether a new level is required
    if( !any(dNl > 0.01*Nl) ) {
      # 23/3/18 this change copes with cases with erratic ml values
      range <- 0:min(2, L-1);
      rem   <- max(ml[L+1-range] / 2^(range*alpha)) / (2^alpha - 1)
      #rem <- ml[L+1] / (2^alpha - 1)

      if(rem > sqrt(theta)*eps) {
        if(L==Lmax) {
          warning("failed to achieve weak convergence.")
        } else {
          L       <- L+1
          Vl[L+1] <- Vl[L] / 2^beta
          Cl[L+1] <- Cl[L] * 2^gamma
          Nl[L+1] <- 0
          suml <- cbind(suml, 0)
          costl[L+1] <- 0

          Ns  <- ceiling(sqrt(Vl/Cl) * sum(sqrt(Vl*Cl)) / ((1-theta)*eps^2))
          dNl <- pmax(0, Ns-Nl)
        }
      }
    }
  }

  # finally, evaluate multilevel estimator
  P <- sum(suml[1,]/Nl)
  list(P = P, Nl = Nl, Cl = Cl)
}
