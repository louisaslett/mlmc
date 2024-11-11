# mlmc 2.1.1

* Bug fix in parallel processing for main driver and `mlmc.test` (thanks to Qian Xin, University of Bristol, for bug report).
* At the same time, improve the method of splitting simulations in parallel for the main `mlmc` driver, so that work is more evenly distributed to keep all cores busy.

# mlmc 2.1.0

* Add parameter value checks in `mlmc.test`.
* Allow user to specify `alpha`, `beta`, and `gamma` to `mlmc.test`, rather than forcing estimation by linear regression.
  Note this is a departure from the original Matlab code, but if they are left unspecified then the same results as under Matlab are reproduced.
* Improve specificity of some argument documentation in `mlmc.test`.

# mlmc 2.0.2

* Package was removed from CRAN because I didn't notice my old Oxford email address wasn't forwarding any longer.
  In order to comply with CRAN changes, the C++ routines are now registered and maintainer info updated to my Durham email.
* The Matlab driver code by Mike Giles has been quite substantially updated, so this major version bump in the R package addresses updating this code to match the new driver API.
  None of these sub-bullets are bug fixes, merely changing to match the new best-practice for the MLMC driver designed by Mike Giles.
  In particular:
    * User level sampling functions must now also return the total cost of all samples simulated at that level.
      Therefore user level sampler functions must return a list with a `sums` and `cost` element.
    * The `gamma` argument is no longer required, since it is not used in automatic cost computation, and can be estimated as for `alpha` and `beta`.
    * `mlmc.test()` no longer takes `M`, a level refinement factor, since this was only used to calculate the cost as `N*M^l`.
      Per above comment, the user now defines cost completely via the return from the level sampler function.
    * Along these lines, `mlmc.test()` now uses the user returned cost in all places: previously CPU time was measured as cost in the convergence tests section, whilst the MLMC complexity tests previously forced costs to be `N*M^l`.
* Some (very) old bugs were squashed in the Euler-Maruyama discretisation level sampler, `opre_l()` which affected lookback call and Heston model options.
* I managed to get hold of a Matlab license, so have now confirmed that the examples in the docs return (within sampling variability) the same results for both Euler-Maruyama and Milstein discretisation example level sampler functions.
* There is now a hex sticker!
  It is hopefully fairly self explanatory: many fast simulations are done at low levels (lots of dice, with the hare running at the bottom of the stairs); fewer simulations are done at higher levels (fewer dice as you go up each step, with a tortoise and fewest dice on top step)!
* There is now a documentation website at <https://mlmc.louisaslett.com/>

# mlmc 1.0.0

* Initial CRAN submission.
