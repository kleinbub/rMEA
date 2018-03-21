## Resubmission
This is a resubmission. In this version I have:

* un-capitalized "Motion Energy" in my description
* added a reference for the method in the 'Description' field in the DESCRIPTION file
* Included a small example for almost each of the functions included in the package. _In this regard please note that motion energy analysis requires high density data and results in relatively large objects. The most efficient way to provide examples is to start from the raw data import and run the analyses instead of saving a completely analyzed .RData of several megabytes. This leads to some example during a bit longer than 10s, these have been put in \donttest{} sections._
* I removed from the package any possible attempt to write in the user's home filespace by default. At least to the best of my knowledge.
* small bug fixes.

Thank you


## Previous resubmissions
* provided correct image links for the README documents

## Test environments
* local OS X install, R 3.4.3
* ubuntu 14.04 (on travis-ci), R 3.4.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
