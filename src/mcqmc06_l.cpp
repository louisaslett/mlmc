// This code is derived and adapted from the original GPL-2 C++ version by
// Mike Giles.  See https://people.maths.ox.ac.uk/~gilesm/mlmc/

#include <Rcpp.h>
using namespace Rcpp;

#define NORMCDF(a) (0.5*erfc(-(a)/sqrt(2.0)))

//' Financial options using a Milstein discretisation
//'
//' Financial options based on scalar geometric Brownian motion, similar to Mike Giles' MCQMC06 paper, Giles (2008), using a Milstein discretisation.
//'
//' This function is based on GPL-2 C++ code by Mike Giles.
//'
//' @param l
//'        the level to be simulated.
//' @param N
//'        the number of samples to be computed.
//' @param option
//'        the option type, between 1 and 5.
//'        The options are:
//'        \describe{
//'          \item{1 = European call;}{}
//'          \item{2 = Asian call;}{}
//'          \item{3 = lookback call;}{}
//'          \item{4 = digital call;}{}
//'          \item{5 = barrier call.}{}
//'        }
//'
//' @return A named list containing: \describe{
//'           \item{\code{sums}}{is a vector of length six \eqn{\left(\sum Y_i, \sum Y_i^2, \sum Y_i^3, \sum Y_i^4, \sum X_i, \sum X_i^2\right)} where \eqn{Y_i} are iid simulations with expectation \eqn{E[P_0]} when \eqn{l=0} and expectation \eqn{E[P_l-P_{l-1}]} when \eqn{l>0}, and \eqn{X_i} are iid simulations with expectation \eqn{E[P_l]}.
//'                              Note that only the first two components of this are used by the main \code{\link[=mlmc]{mlmc()}} driver, the full vector is used by \code{\link[=mlmc.test]{mlmc.test()}} for convergence tests etc;}
//'           \item{\code{cost}}{is a scalar with the total cost of the paths simulated, computed as \eqn{N \times 2^l} for level \eqn{l}.}
//'         }
//'
//' @author Louis Aslett <louis.aslett@durham.ac.uk>
//' @author Mike Giles <Mike.Giles@maths.ox.ac.uk>
//'
//' @references
//' Giles, M. (2008) 'Improved Multilevel Monte Carlo Convergence using the Milstein Scheme', in A. Keller, S. Heinrich, and H. Niederreiter (eds) \emph{Monte Carlo and Quasi-Monte Carlo Methods 2006}. Berlin, Heidelberg: Springer, pp. 343–358. Available at: \doi{10.1007/978-3-540-74496-2_20}.
//'
//' @examples
//' \donttest{
//' # These are similar to the MLMC tests for the MCQMC06 paper
//' # using a Milstein discretisation with 2^l timesteps on level l
//' #
//' # The figures are slightly different due to:
//' # -- change in MSE split
//' # -- change in cost calculation
//' # -- different random number generation
//' # -- switch to S_0=100
//' #
//' # Note the following takes quite a while to run, for a toy example see after
//' # this block.
//'
//' N0   <- 200 # initial samples on coarse levels
//' Lmin <- 2 # minimum refinement level
//' Lmax <- 10 # maximum refinement level
//'
//' test.res <- list()
//' for(option in 1:5) {
//'   if(option == 1) {
//'     cat("\n ---- Computing European call ---- \n")
//'     N      <- 20000 # samples for convergence tests
//'     L      <- 8 # levels for convergence tests
//'     Eps    <- c(0.005, 0.01, 0.02, 0.05, 0.1)
//'   } else if(option == 2) {
//'     cat("\n ---- Computing Asian call ---- \n")
//'     N      <- 20000 # samples for convergence tests
//'     L      <- 8 # levels for convergence tests
//'     Eps    <- c(0.005, 0.01, 0.02, 0.05, 0.1)
//'   } else if(option == 3) {
//'     cat("\n ---- Computing lookback call ---- \n")
//'     N      <- 20000 # samples for convergence tests
//'     L      <- 10 # levels for convergence tests
//'     Eps    <- c(0.005, 0.01, 0.02, 0.05, 0.1)
//'   } else if(option == 4) {
//'     cat("\n ---- Computing digital call ---- \n")
//'     N      <- 200000 # samples for convergence tests
//'     L      <- 8 # levels for convergence tests
//'     Eps    <- c(0.01, 0.02, 0.05, 0.1, 0.2)
//'   } else if(option == 5) {
//'     cat("\n ---- Computing barrier call ---- \n")
//'     N      <- 200000 # samples for convergence tests
//'     L      <- 8 # levels for convergence tests
//'     Eps    <- c(0.005, 0.01, 0.02, 0.05, 0.1)
//'   }
//'
//'   test.res[[option]] <- mlmc.test(mcqmc06_l, N, L, N0, Eps, Lmin, Lmax, option = option)
//'
//'   # print exact analytic value, based on S0=K
//'   T   <- 1
//'   r   <- 0.05
//'   sig <- 0.2
//'   K   <- 100
//'   B   <- 0.85*K
//'
//'   k   <- 0.5*sig^2/r;
//'   d1  <- (r+0.5*sig^2)*T / (sig*sqrt(T))
//'   d2  <- (r-0.5*sig^2)*T / (sig*sqrt(T))
//'   d3  <- (2*log(B/K) + (r+0.5*sig^2)*T) / (sig*sqrt(T))
//'   d4  <- (2*log(B/K) + (r-0.5*sig^2)*T) / (sig*sqrt(T))
//'
//'   if(option == 1) {
//'     val <- K*( pnorm(d1) - exp(-r*T)*pnorm(d2) )
//'   } else if(option == 2) {
//'     val <- NA
//'   } else if(option == 3) {
//'     val <- K*( pnorm(d1) - pnorm(-d1)*k - exp(-r*T)*(pnorm(d2) - pnorm(d2)*k) )
//'   } else if(option == 4) {
//'     val <- K*exp(-r*T)*pnorm(d2)
//'   } else if(option == 5) {
//'     val <- K*(                             pnorm(d1) - exp(-r*T)*pnorm(d2) -
//'               ((K/B)^(1-1/k))*((B^2)/(K^2)*pnorm(d3) - exp(-r*T)*pnorm(d4)) )
//'   }
//'
//'   if(is.na(val)) {
//'     cat(sprintf("\n Exact value unknown, MLMC value: %f \n", test.res[[option]]$P[1]))
//'   } else {
//'     cat(sprintf("\n Exact value: %f, MLMC value: %f \n", val, test.res[[option]]$P[1]))
//'   }
//'
//'   # plot results
//'   plot(test.res[[option]])
//' }
//' }
//'
//' # The level sampler can be called directly to retrieve the relevant level sums:
//' mcqmc06_l(l = 7, N = 10, option = 1)
//'
//' @export
// [[Rcpp::export]]
List mcqmc06_l(int l, int N, int option) {
  RNGScope scope;

  NumericVector sums(6);
  NumericVector cost(1);

  if(l < 0) {
    stop("l must be >= 0\n");
  }
  if(N < 1) {
    stop("N must be > 0\n");
  }
  if(option < 1 || option > 5) {
    stop("option must be between 1 and 5 inclusive\n");
  }

  int   nf, nc;
  double T, r, sig, B, hf, hc, X0, Xf, Xc, Af, Ac, Mf, Mc, Bf, Bc,
  Xf0 = 0.0, Xc0 = 0.0, Xc1, vf, vc, dWc, ddW, Pf, Pc, dP, K;

  double dWf[2], dIf[2], Lf[2];

  // ull   v1[3], v2[3];       // needed for RNG
  // float x1, x2 = nanf("");  // needed for Normal RNG

  // initialise seeds

  // for (int m=0; m<3; m++) {
  //   v1[m] = CPU_mrg32k3a_v1[m];
  //   v2[m] = CPU_mrg32k3a_v2[m];
  // }

  K   = 100.0;
  T   = 1.0;
  r   = 0.05;
  sig = 0.2;
  B   = 0.85*K;

  nf = 1<<l;
  nc = nf/2;

  hf = T / ((double) nf);
  hc = T / ((double) nc);

  for (int np = 0; np<N; np++) {
    X0 = K;

    Xf = X0;
    Xc = Xf;

    Af  = 0.5*hf*Xf;
    Ac  = 0.5*hc*Xc;

    Mf  = Xf;
    Mc  = Xc;

    Bf  = 1.0;
    Bc  = 1.0;

    if (l==0) {
      // CPU_mrg32k3a_next_normal(v1,v2,x1,x2);
      dWf[0] = sqrt(hf)*R::rnorm(0,1);

      // CPU_mrg32k3a_next_uniform(v1,v2,x1);
      Lf[0] = log(R::runif(0,1));

      // CPU_mrg32k3a_next_normal(v1,v2,x1,x2);
      dIf[0] = sqrt(hf/12.0)*hf*R::rnorm(0,1);

      Xf0 = Xf;
      Xf  = Xf + r*Xf*hf + sig*Xf*dWf[0]
      + 0.5*sig*sig*Xf*(dWf[0]*dWf[0]-hf);
      vf  = sig*Xf0;
      Af  = Af + 0.5*hf*Xf + vf*dIf[0];
      Mf  = fmin(Mf,
                  0.5*(Xf0+Xf-sqrt((Xf-Xf0)*(Xf-Xf0)-2.0*hf*vf*vf*Lf[0])));
      Bf  = Bf*(1.0-exp(-2.0*fmax(0.0,(Xf0-B)*(Xf-B)/(hf*vf*vf))));
    }

    else {
      for (int n=0; n<nc; n++) {
        // CPU_mrg32k3a_next_normal(v1,v2,x1,x2);
        dWf[0] = sqrt(hf)*R::rnorm(0,1);
        // CPU_mrg32k3a_next_normal(v1,v2,x1,x2);
        dWf[1] = sqrt(hf)*R::rnorm(0,1);

        // CPU_mrg32k3a_next_uniform(v1,v2,x1);
        Lf[0] = log(R::runif(0,1));
        // CPU_mrg32k3a_next_uniform(v1,v2,x1);
        Lf[1] = log(R::runif(0,1));

        // CPU_mrg32k3a_next_normal(v1,v2,x1,x2);
        dIf[0] = sqrt(hf/12.0)*hf*R::rnorm(0,1);
        // CPU_mrg32k3a_next_normal(v1,v2,x1,x2);
        dIf[1] = sqrt(hf/12.0)*hf*R::rnorm(0,1);

        for (int m=0; m<2; m++) {
          Xf0 = Xf;
          Xf  = Xf + r*Xf*hf + sig*Xf*dWf[m]
          + 0.5*sig*sig*Xf*(dWf[m]*dWf[m]-hf);
          vf  = sig*Xf0;
          Af  = Af + hf*Xf + vf*dIf[m];
          Mf  = fmin(Mf,
                      0.5*(Xf0+Xf-sqrt((Xf-Xf0)*(Xf-Xf0)-2.0*hf*vf*vf*Lf[m])));
          Bf  = Bf*(1.0-exp(-2.0f*fmax(0.0,(Xf0-B)*(Xf-B)/(hf*vf*vf))));
        }

        dWc = dWf[0] + dWf[1];
        ddW = dWf[0] - dWf[1];

        Xc0 = Xc;
        Xc  = Xc + r*Xc*hc + sig*Xc*dWc + 0.5*sig*sig*Xc*(dWc*dWc-hc);

        vc  = sig*Xc0;
        Ac  = Ac + hc*Xc + vc*(dIf[0]+dIf[1] + 0.25*hc*ddW);
        Xc1 = 0.5*(Xc0 + Xc + vc*ddW);
        Mc  = fmin(Mc,
                    0.5*(Xc0+Xc1-sqrt((Xc1-Xc0)*(Xc1-Xc0)-2.0*hf*vc*vc*Lf[0])));
        Mc  = fmin(Mc,
                    0.5*(Xc1+Xc -sqrt((Xc -Xc1)*(Xc -Xc1)-2.0*hf*vc*vc*Lf[1])));
        Bc  = Bc *(1.0-exp(-2.0*fmax(0.0,(Xc0-B)*(Xc1-B)/(hf*vc*vc))));
        Bc  = Bc *(1.0-exp(-2.0*fmax(0.0,(Xc1-B)*(Xc -B)/(hf*vc*vc))));
      }
      Af = Af - 0.5*hf*Xf;
      Ac = Ac - 0.5*hc*Xc;
    }

    if (option==1) {
      Pf  = fmax(0.0,Xf-K);
      Pc  = fmax(0.0,Xc-K);
    }
    else if (option==2) {
      Pf  = fmax(0.0,Af-K);
      Pc  = fmax(0.0,Ac-K);
    }
    else if (option==3) {
      Pf  = Xf - Mf;
      Pc  = Xc - Mc;
    }
    else if (option==4) {
      Pf  = K*NORMCDF((Xf0+r*Xf0*hf-K)/(sig*Xf0*sqrt(hf)));
      if (l==0)
        Pc = Pf;
      else
        Pc = K
        * NORMCDF((Xc0+r*Xc0*hc+sig*Xc0*dWf[0]-K)/(sig*Xc0*sqrt(hf)));
    }
    else if (option==5) {
      Pf  = Bf*fmax(0.0f,Xf-K);
      Pc  = Bc*fmax(0.0f,Xc-K);
    } else {
      // Should be impossible to reach here, but adding to hint to compiler
      // these variables are never uninitialised
      Pf = 0.0; Pc = 0.0;
      stop("option must be between 1 and 5 inclusive\n");
    }

    dP  = exp(-r*T)*(Pf-Pc);
    Pf  = exp(-r*T)*Pf;

    if (l==0) dP = Pf;

    //sums[0] += nf;     // add number of timesteps as cost
    sums[0] += dP;
    sums[1] += dP*dP;
    sums[2] += dP*dP*dP;
    sums[3] += dP*dP*dP*dP;
    sums[4] += Pf;
    sums[5] += Pf*Pf;

    cost[0] += nf; // add number of timesteps as cost
  }

  return(List::create(Named("sums") = sums,
                      Named("cost") = cost));
}
