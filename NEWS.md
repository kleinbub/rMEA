# rMEA

### v1.1.0 
[2019/03/22] Yearly update to CRAN 

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
