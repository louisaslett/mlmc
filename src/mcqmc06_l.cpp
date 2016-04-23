// This code is derived and adapted from the original GPL-2 C++ version by
// Mike Giles.  See http://people.maths.ox.ac.uk/~gilesm/mlmc/

#include <Rcpp.h>
using namespace Rcpp;

#define NORMCDF(a) (0.5*erfc(-(a)/sqrt(2.0)))

//' Financial options using a Milstein discretisation
//'
//' Financial options based on scalar geometric Brownian motion, similar to
//' Mike Giles' MCQMC06 paper, using a Milstein discretisation
//'
//' This function is based on GPL-2 C++ code by Mike Giles.
//'
//' @param l the level to be simulated.
//' @param N the number of samples to be computed.
//' @param option the option type, between 1 and 5.  The options are: \describe{
//'   \item{1 = European call;}{}
//'   \item{2 = Asian call;}{}
//'   \item{3 = lookback call;}{}
//'   \item{4 = digital call;}{}
//'   \item{5 = barrier call.}{}
//' }
//'
//' @author Louis Aslett <aslett@stats.ox.ac.uk>
//' @author Mike Giles <Mike.Giles@maths.ox.ac.uk>
//'
//' @examples
//' \dontrun{
//' # These are similar to the MLMC tests for the MCQMC06 paper
//' # using a Milstein discretisation with 2^l timesteps on level l
//' #
//' # The figures are slightly different due to:
//' # -- change in MSE split
//' # -- change in cost calculation
//' # -- different random number generation
//' # -- switch to S_0=100
//'
//' M    <- 2 # refinement cost factor
//' N0   <- 200 # initial samples on coarse levels
//' Lmin <- 2 # minimum refinement level
//' Lmax <- 10 # maximum refinement level
//'
//' test.res <- list()
//' for(option in 1:5) {
//'   if(option==1) {
//'     cat("\n ---- Computing European call ---- \n")
//'     N      <- 20000 # samples for convergence tests
//'     L      <- 8 # levels for convergence tests
//'     Eps    <- c(0.005, 0.01, 0.02, 0.05, 0.1)
//'   } else if(option==2) {
//'     cat("\n ---- Computing Asian call ---- \n")
//'     N      <- 20000 # samples for convergence tests
//'     L      <- 8 # levels for convergence tests
//'     Eps    <- c(0.005, 0.01, 0.02, 0.05, 0.1)
//'   } else if(option==3) {
//'     cat("\n ---- Computing lookback call ---- \n")
//'     N      <- 20000 # samples for convergence tests
//'     L      <- 10 # levels for convergence tests
//'     Eps    <- c(0.005, 0.01, 0.02, 0.05, 0.1)
//'   } else if(option==4) {
//'     cat("\n ---- Computing digital call ---- \n")
//'     N      <- 200000 # samples for convergence tests
//'     L      <- 8 # levels for convergence tests
//'     Eps    <- c(0.01, 0.02, 0.05, 0.1, 0.2)
//'   } else if(option==5) {
//'     cat("\n ---- Computing barrier call ---- \n")
//'     N      <- 200000 # samples for convergence tests
//'     L      <- 8 # levels for convergence tests
//'     Eps    <- c(0.005, 0.01, 0.02, 0.05, 0.1)
//'   }
//'
//'   test.res[[option]] <- mlmc.test(mcqmc06_l, M, N, L, N0, Eps, Lmin, Lmax, option=option)
//'
//'   # plot results
//'   plot(test.res[[option]])
//' }
//' }
//' @export
// [[Rcpp::export]]
NumericVector mcqmc06_l(int l, int N, int option) {
  RNGScope scope;

  int   nf, nc;
  double T, r, sig, B, hf, hc, X0, Xf, Xc, Af, Ac, Mf, Mc, Bf, Bc,
  Xf0, Xc0, Xc1, vf, vc, dWc, ddW, Pf, Pc, dP, K;

  double dWf[2], dIf[2], Lf[2];

  // ull   v1[3], v2[3];       // needed for RNG
  // float x1, x2 = nanf("");  // needed for Normal RNG

  NumericVector sums(6);

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
  }
  return(sums);
}
