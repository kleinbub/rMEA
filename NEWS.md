---
title: "rMEA NEWS"
pagetitle: rMEA
---
### v1.2.1
#### [2021/03/31] minor bugfixes
  
  * lagSec=0 and lagSec=1 are now supported
  * documentation update

### v1.2.0
#### [2020/07/29] Yearly update to CRAN

  * the package is now compatible with R version 4.0
  * New function: shuffle_segment(), offering within subject shuffling. This function is included for replication of older studies. shuffle() is still more conservative, and suggested approach to pseudosynchrony.
  * New function CCFartefacts() to identify sequences with extremely high correlation, which may require inspection of the videos.
  * New function: MEAreplace() allows to automatically delete (set to NA) or replace windows of MEA data, based on a data frame created by hand, or with functions like CCFartefacts().
  
  
##### Minor changes

  * various improvements in the documentation
  * all function now dynamically use ccfResNames.
  * writeMEA now reports all ccfResNames.
  * ccfRes now stores as well the 'bestLag' values, and the start and end of the synchronization windows (thanks to anonymous reviewer 2 for the suggestion).
  * MEAheatmap received a new parameter 'mirror' that offers finer control on the color scale.
  * shuffle() now allows to specify size = "max", to use all possible combinations.
  * automatic detection of 'na.rm' arguments in MEAscale FUNs was buggy. So now there is a new argument 'removeNA' which applies na.omit on the data before submitting it to the scaling FUN.
  * improved the MEAmap use of dots (...), which are now recycled and iterating correctly (when not a function)
  

### v1.1.0.9002
  * fixed a bug preventing heatmaps to be plotted with missing data
  * improved MEAdistplot() to allow grouping by "id", "session", and "group" 

### v1.1.0.9001
  * improved MEAdistplot() to support single session groups
  * fixed a bug in MEAdistplot for plotting single groups

### v1.1.0 
#### [2019/03/22] Yearly update to CRAN

  * included LICENSE file
  * improved fisher's transform performance
  * updated rangeRescale function
  * fixed minor bug in heatmap scaling
  * added heatmap "rescale" parameter to highlight small trends in data
  * improved axis label in various plots
  * plot.MEA can now use "duration" instead of "to"
  * various minor bug fixes
  * added plot.MEAlist


### v1.0.0.9012
* Implemented "[.MEAlist" method to extract from MEAlist while retaining metadata.
* included a function to extract CCF values (getCCF())
* shuffle now retains the names of the original signals in its _uid_

### v1.0.0.9011

* fixed a bug in plot.MEA when interpolation results in NA vals
* updated documentation in readMEA
* all functions accepting _MEAlist_ objects now should accept a list instead and try to parse as _MEAlist_ internally
* new getter functions for MEA and MEAlist attributes


### v1.0.0.9010
* more informative error messages
* fixed smaller bugs

### v1.0.0 Initial release
* Everything should work
