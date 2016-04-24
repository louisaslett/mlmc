# mlmc R package
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![license](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![metacran version](http://www.r-pkg.org/badges/version/mlmc)](http://cran.r-project.org/web/packages/mlmc/index.html)
[![metacran downloads](http://cranlogs.r-pkg.org/badges/mlmc?color=brightgreen)](http://cran.r-project.org/web/packages/mlmc/index.html)
[CRAN check result](http://cran.r-project.org/web/checks/check_results_mlmc.html)

An implementation of Multi-level Monte Carlo for R.  This package builds on the original GPL-2 Matlab and C++ implementations by Mike Giles (see <http://people.maths.ox.ac.uk/~gilesm/mlmc/>) to provide a full MLMC driver and example level samplers.  Multi-core parallel sampling of levels is provided built-in.

## Contact

Please feel free to:

* submit suggestions and bug-reports at: <https://github.com/louisaslett/mlmc/issues>
* compose an e-mail to: <aslett@stats.ox.ac.uk>, <nagapetyan@stats.ox.ac.uk> or <vollmer@stats.ox.ac.uk>

## Install

You can install the latest release directly from
[CRAN](http://cran.r-project.org/web/packages/mlmc/index.html).

```r
install.packages("mlmc")
```

## Install development version (not recommended)

Installing directly from [GitHub](https://github.com) is not supported by the
`install.packages` command. You could use the
[devtools](http://cran.r-project.org/web/packages/devtools/index.html) package
to install the development version if desired.

```r
install.packages("devtools")
library("devtools")
install_github("louisaslett/mlmc")
```

Under releases, the tree/commit from which CRAN releases were made are recorded,
so historic source can be downloaded from there.

## Acknowledgements

Louis Aslett is supported by the i-like programme grant (EPSRC grant reference number EP/K014463/1 <http://www.i-like.org.uk>).  Tigran Nagapetyan and Sebastian Vollmer are supported by EPSRC Grant EP/N000188/1.

## Citation

If you use this software, please cite:

Aslett, L. J. M., Nagapetyan, T. and Vollmer, S. J. (2016), *mlmc: Tools for Multilevel Monte Carlo*.  R package.

Thank-you.
