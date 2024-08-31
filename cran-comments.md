## Resubmission

This is a resubmission in light of the kind CRAN review by Beni Altmann.
In this version I have:

* Added acronym explanation and single quoted 'C++' in DESCRIPTION, as requested.
* Added return value information to `opre_l`, `mcqml06_l`, and `plot.mlmc.test`. Thank you very much for identifying this (CRAN review is so important and works!), I thought I'd been so careful and am horrified these slipped through.
* Changed all \dontrun to \donttest since the examples are indeed merely long-running, not unable to be run.


## R CMD check results

0 errors | 0 warnings | 1 note

* The note arises because the package was archived due to email to the maintainer being undeliverable. I have updating my maintainer email address, at the same time as making a substantial number of updates, bug fixes and improvements to the package.
