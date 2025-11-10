# Multi-level Monte Carlo estimation

This function is the Multi-level Monte Carlo driver which will sample
from the levels of user specified function.

## Usage

``` r
mlmc(
  Lmin,
  Lmax,
  N0,
  eps,
  mlmc_l,
  alpha = NA,
  beta = NA,
  gamma = NA,
  parallel = NA,
  ...
)
```

## Arguments

- Lmin:

  the minimum level of refinement. Must be \\\ge 2\\.

- Lmax:

  the maximum level of refinement. Must be \\\ge\\ `Lmin`.

- N0:

  initial number of samples which are used for the first 3 levels and
  for any subsequent levels which are automatically added. Must be \\\>
  0\\.

- eps:

  the target accuracy of the estimate (root mean square error). Must be
  \\\> 0\\.

- mlmc_l:

  a user supplied function which provides the estimate for level \\l\\.
  It must take at least two arguments, the first is the level number to
  be simulated and the second the number of paths. Additional arguments
  can be taken if desired: all additional `...` arguments to this
  function are forwarded to the user defined `mlmc_l` function.

  The user supplied function should return a named list containing one
  element named `sums` and second named `cost`, where:

  `sums`

  :   is a vector of length at least two. The first two elements should
      be \\\left(\sum Y_i, \sum Y_i^2\right)\\ where \\Y_i\\ are iid
      simulations with expectation \\E\[P_0\]\\ when \\l=0\\ and
      expectation \\E\[P_l-P\_{l-1}\]\\ when \\l\>0\\. Note that
      typically the user supplied level sampler will actually return a
      vector of length six, also enabling use of the
      [`mlmc.test()`](https://mlmc.louisaslett.com/reference/mlmc.test.md)
      function to perform convergence tests, kurtosis, and telescoping
      sum checks. See
      [`mlmc.test()`](https://mlmc.louisaslett.com/reference/mlmc.test.md)
      for the definition of these remaining four elements.

  `cost`

  :   is a scalar with the total cost of the paths simulated. For
      example, in the financial options samplers included in this
      package, this is calculated as \\NM^l\\, where \\N\\ is the number
      of paths requested in the call to the user function `mlmc_l`,
      \\M\\ is the refinement cost factor (\\M=2\\ for
      [`mcqmc06_l()`](https://mlmc.louisaslett.com/reference/mcqmc06_l.md)
      and \\M=4\\ for
      [`opre_l()`](https://mlmc.louisaslett.com/reference/opre_l.md)),
      and \\l\\ is the level being sampled.

  See the function (and source code of)
  [`opre_l()`](https://mlmc.louisaslett.com/reference/opre_l.md) and
  [`mcqmc06_l()`](https://mlmc.louisaslett.com/reference/mcqmc06_l.md)
  in this package for an example of user supplied level samplers.

- alpha:

  the weak error, \\O(2^{-\alpha l})\\. Must be \\\> 0\\ if specified.
  If `NA` then `alpha` will be estimated.

- beta:

  the variance, \\O(2^{-\beta l})\\. Must be \\\> 0\\ if specified. If
  `NA` then `beta` will be estimated.

- gamma:

  the sample cost, \\O(2^{\gamma l})\\. Must be \\\> 0\\ if specified.
  If `NA` then `gamma` will be estimated.

- parallel:

  if an integer is supplied, R will fork `parallel` parallel processes
  and spread the simulations required at each level as evenly as
  possible across all cores.

- ...:

  additional arguments which are passed on when the user supplied
  `mlmc_l` function is called.

## Value

A named list containing:

- `P`:

  The MLMC estimate;

- `Nl`:

  A vector of the number of samples performed on each level;

- `Cl`:

  Per sample cost at each level.

## Details

The Multilevel Monte Carlo Method method originated in the works Giles
(2008) and Heinrich (1998).

Consider a sequence \\P_0, P_1, \ldots\\, which approximates \\P_L\\
with increasing accuracy, but also increasing cost, we have the simple
identity \$\$E\[P_L\] = E\[P_0\] + \sum\_{l=1}^L E\[P_l-P\_{l-1}\],\$\$
and therefore we can use the following unbiased estimator for
\\E\[P_L\]\\, \$\$N_0^{-1} \sum\_{n=1}^{N_0} P_0^{(0,n)} + \sum\_{l=1}^L
\left\\ N_l^{-1} \sum\_{n=1}^{N_l} \left(P_l^{(l,n)} -
P\_{l-1}^{(l,n)}\right) \right\\\$\$ where \\N_l\\ samples are produced
at level \\l\\. The inclusion of the level \\l\\ in the superscript
\\(l,n)\\ indicates that the samples used at each level of correction
are independent.

Set \\C_0\\, and \\V_0\\ to be the cost and variance of one sample of
\\P_0\\, and \\C_l, V_l\\ to be the cost and variance of one sample of
\\P_l - P\_{l-1}\\, then the overall cost and variance of the multilevel
estimator is \\\sum\_{l=0}^L N_l C_l\\ and \\\sum\_{l=0}^L N_l^{-1}
V_l\\, respectively.

The idea behind the method, is that provided that the product \\V_l
C_l\\ decreases with \\l\\, i.e. the cost increases with level slower
than the variance decreases, then one can achieve significant
computational savings, which can be formalised as in Theorem 1 of Giles
(2015).

For further information on multilevel Monte Carlo methods, see the
webpage <https://people.maths.ox.ac.uk/gilesm/mlmc_community.html> which
lists the research groups working in the area, and their main
publications.

This function is based on GPL-2 'Matlab' code by Mike Giles.

## References

Giles, M.B. (2008) 'Multilevel Monte Carlo Path Simulation', *Operations
Research*, 56(3), pp. 607–617. Available at:
[doi:10.1287/opre.1070.0496](https://doi.org/10.1287/opre.1070.0496) .

Giles, M.B. (2015) 'Multilevel Monte Carlo methods', *Acta Numerica*,
24, pp. 259–328. Available at:
[doi:10.1017/S096249291500001X](https://doi.org/10.1017/S096249291500001X)
.

Heinrich, S. (1998) 'Monte Carlo Complexity of Global Solution of
Integral Equations', *Journal of Complexity*, 14(2), pp. 151–175.
Available at:
[doi:10.1006/jcom.1998.0471](https://doi.org/10.1006/jcom.1998.0471) .

## Author

Louis Aslett \<louis.aslett@durham.ac.uk\>

Mike Giles \<Mike.Giles@maths.ox.ac.uk\>

Tigran Nagapetyan \<nagapetyan@stats.ox.ac.uk\>

## Examples

``` r
mlmc(2, 6, 1000, 0.01, opre_l, option = 1)
#> $P
#> [1] 10.43702
#> 
#> $Nl
#> [1] 4285886  356016  107236   22251
#> 
#> $Cl
#> [1]  1  4 16 64
#> 

mlmc(2, 10, 1000, 0.01, mcqmc06_l, option = 1)
#> $P
#> [1] 10.46893
#> 
#> $Nl
#> [1] 2966949   60337   25167    8103    2922    1056     382     144
#> 
#> $Cl
#> [1]   1   2   4   8  16  32  64 128
#> 
```
