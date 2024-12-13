#' Moving average smoothing for motion energy data
#'
#' This function applies a moving average filter, based on SAS "proc expand" procedure. The moving average is applied independently on each subject's motion energy data.
#' NA values are set to 0.
#'
#' @param mea an object of class \code{MEA} or a list of \code{MEA} objects (see function \code{\link{readMEA}})
#' @param moving.average.win numeric. The size of the filter window, in seconds or fractions of seconds.
#'
#' @return The filtered object(s)
#' @examples
#' ## read the first 4 minutes of the normal sample
#' ##   (intake interviews of patients that carried on therapy)
#' path_normal <- system.file("extdata/normal", package = "rMEA")
#' mea_raw <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
#'                      s1Name = "Patient", s2Name = "Therapist",
#'                      idOrder = c("id","session"), idSep="_", skip=1, nrow = 6000)
#'
#' ## filter with moving average
#' mea_filter = MEAsmooth(mea_raw)
#'
#' ## assign groups names
#' mea_raw <- setGroup(mea_raw, "raw")
#' mea_filter <- setGroup(mea_filter, "filtered")
#'
#' ## Compute ccf
#' mea_raw <- MEAccf(mea_raw, lagSec = 5, winSec = 60, incSec = 30, r2Z = TRUE, ABS = FALSE)
#' mea_filter <- MEAccf(mea_filter, lagSec = 5, winSec = 60, incSec = 30, r2Z = TRUE, ABS = FALSE)
#'
#' ## Compare the effect of filtering on ccf
#' MEAlagplot(mea_filter, contrast = mea_raw)
#'
#' @export
#'
#'
#
MEAsmooth = function(mea, moving.average.win = 0.5){
  cat("\r\nMoving average smoothing:\r\n")
  res = MEAmap(mea, function(x){
    na1 = which(is.na(x[[1]]))
    na2 = which(is.na(x[[2]]))
    x[[1]] = sasFiltB(x[[1]], winSec = moving.average.win, sampRate = attr(mea,"sampRate"))
    x[[2]] = sasFiltB(x[[2]], winSec = moving.average.win, sampRate = attr(mea,"sampRate"))
    x[[1]][na1] = NA
    x[[2]][na2] = NA
    x
  }, label="smooth" )
}

sasFiltB = function(a, sampRate, winSec=0.5){
  n = winSec*sampRate
  a = unlist(lapply(seq_along(a), function(t){
    i1 = t-(n-1)/2; if(i1<1) i1=1;
    i2 = t+(n-1)/2;
    res = sum(a[i1:i2])/n
    res[is.na(res)]=0
    return(res)
  }))
}



#' Scaling (and centering) of motion energy time-series
#'
#' @param mea an object of class \code{MEA} or a list of \code{MEA} objects (see function \code{\link{readMEA}})
#' @param scale either a numeric value or a function to be applied to each motion
#' energy time-series to calculate a scaling factor. Default is standard deviation.
#' @param ... further arguments passed to \code{scale} if it is a function.
#' @param center either a logical value or a numeric vector of length 2 specifying separate centering values for s1 and s2.
#' @param removeNA logical. If \code{scale} is a function, defines whether NAs be removed prior to calculating the scaling factor.
#' @details If \code{scale} is a function, it is found by a call to  \code{\link[base]{match.fun}} and typically is either a function
#' or a symbol (e.g., a backquoted name) or a character string specifying a function
#' to be searched for from the environment of the call to apply. Note that the chosen function must return a single numeric value.
#'
#' \code{center} is directly passed to \code{\link[base]{scale}}. If \code{center} is \code{TRUE} then centering
#' is done by subtracting the means (omitting NAs) from the motion energy time-series. If \code{center} is a numeric vector,
#' the first value will be subtracted from s1 and the second from s2.
#' Note: the s1 and s2 signals are scaled independently.
#' @return returns the same \code{MEA} or \code{MEAlist} object, with all motion energy data rescaled
# #' @importFrom stats sd
#' @examples
#' ## read the first 4 minutes of the normal sample
#' ##   (intake interviews of patients that carried on therapy)
#' path_normal <- system.file("extdata/normal", package = "rMEA")
#' mea_raw <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
#'                      s1Name = "Patient", s2Name = "Therapist",
#'                      idOrder = c("id","session"), idSep="_", skip=1, nrow = 6000)
#'
#' ## rescale by factor 0.7
#' mea_scaled = MEAscale(mea_raw, scale = 0.7)
#'
#' ## rescale with standard deviation
#' mea_scaled = MEAscale(mea_raw, scale = "sd", removeNA = TRUE)
#'
#' ## assign groups names
#' mea_raw <- setGroup(mea_raw, "raw")
#' mea_scaled <- setGroup(mea_scaled, "scaled")
#'
#' ## Compute ccf
#' mea_raw <- MEAccf(mea_raw, lagSec = 5, winSec = 60, incSec = 30, r2Z = TRUE, ABS = FALSE)
#' mea_scaled <- MEAccf(mea_scaled, lagSec = 5, winSec = 60, incSec = 30, r2Z = TRUE, ABS = FALSE)
#'
#' ## Compare the effect of scaling on ccf
#' MEAlagplot(mea_scaled, contrast = mea_raw)
#'
#' @export
#'

MEAscale = function(mea, scale="sd", ..., center=F, removeNA=T){
  cat("\r\nRescaling data:\r\n")
  if(is.numeric(scale)){
    scaleVAL = scale
    fil = paste0("rescaled (",scaleVAL,")")
    scaleFUN = function(x){scaleVAL}
    dots = list()
  } else {
    scaleFUN = match.fun(scale)
    fil = paste0("rescaled (", as.character(substitute(scale)),")" )
    if (length(fil)>1) fil = "rescaled (custom)"
    dots = list(...)
  }
  res = MEAmap(mea, function(x){
    if(removeNA) {
      scaleX1 = stats::na.omit(x[[1]]); scaleX2 = stats::na.omit(x[[2]])
    } else {
      scaleX1 = x[[1]]; scaleX2 = x[[2]]
    }
    x[[1]]=base::scale(x[[1]], scale = do.call(scaleFUN, c(list(scaleX1),dots)), center = center)
    x[[2]]=base::scale(x[[2]], scale = do.call(scaleFUN, c(list(scaleX2),dots)), center = center)
    x
  },
  label=fil)
}

#' Replace outliers with given values
#'
#' Sometimes motion energy analysis generates excessively high peaks resulting from video artifacts or other anomalies in the video source.
#'
# This function allows to substitute the values greater or less than a specific threshold. The default threshold is 10 times the standard deviation of the time-series.
# **NOTE: As of version 1.2.0.0 this function is deprecated and may be removed in future versions.
# Please use \code{\link{MEAthreshold}} to identify values and \code{link{MEAreplace}} to perform the substitutions**
#'
#' @param mea an object of class \code{MEA} or a list of \code{MEA} objects (see function \code{\link{readMEA}})
#' @param threshold a numeric value, or a function returning the threshold value to consider data as outliers.
#' @param direction a text string. One of "\code{greater}" or "\code{less}": can be abbreviated.
#' @param replace a numeric, NULL, or NA value to use as substitution.
#'
#' @return The same \code{mea} object with all extreme values substituted.
#' @examples
#' ## read the first 4 minutes of the normal sample
#' ##   (intake interviews of patients that carried on therapy)
#' path_normal <- system.file("extdata/normal", package = "rMEA")
#' mea_raw <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
#'                      s1Name = "Patient", s2Name = "Therapist",
#'                      idOrder = c("id","session"), idSep="_", skip=1, nrow = 6000)
#'
#' ## Remove extreme values, higher than 10 times the standard deviation
#' mea_clean = MEAoutlier(mea_raw, threshold=function(x){sd(x)*10}, direction = "greater")
#' @export
MEAoutlier = function(mea, threshold=function(x){stats::sd(x)*10}, direction=c("greater", "less") ,replace=NA){
  # message("As of version 1.2.0.0 this function is deprecated and may be removed in future versions.
  #         Please use MEAthreshold() to identify values and MEAreplace() to perform the substitutions")
  cat0("\r\nSuppressing outliers to ",as.character(replace),":\r\n")
  direction = match.arg(direction)
  res = MEAmap(mea, function(x){
    if(is.function(threshold)){
      threshold = match.fun(threshold)
      threshold = threshold(c(x[[1]],x[[2]]))
    }
    if (length(threshold)>1) stop("Threshold must be a single value, not ",length(threshold), call.=F)
    if(direction == "greater"){
    x[[1]][x[[1]]>threshold] =replace
    x[[2]][x[[2]]>threshold] =replace
    } else if (direction =="less") {
      x[[1]][x[[1]]<threshold] =replace
      x[[2]][x[[2]]<threshold] =replace
    } else stop("'direction' argument must be either 'greater' or 'less'.",call.=F)
    x
  },
  label= paste("outliers to",as.character(replace))
  )
}



#' Apply a function to a single or a list of MEA objects
#'
#' MEApply is a wrapper to \code{\link[base]{do.call}} that allows to apply a function on the motion energy data of one or multiple \code{MEA} objects. Complex constructions are possible, see details.
#'
#' @param mea an object of class \code{MEA} or a list of \code{MEA} objects (see function \code{\link{readMEA}})
#' @param FUN function to apply, found via \code{\link[base]{match.fun}}.
#' @param label a character vector to update the 'filter' attribute of \code{mea}.
#' @param ... further arguments passed to \code{FUN}. If a function is provided,
#' it will be run on each MEA object and then passed as an argument to \code{FUN}.
#'
#' @details \code{FUN} will be applied on the motion energy time-series of \code{MEA} objects, which is stored as a data frame with 2 columns, respectively for s1 and s2.
#' @return an object of the same class of the provided \code{'mea'} object, with the transformed motion energy data
#' @export

MEAmap = function(mea, FUN, label="", ...){
  UseMethod("MEAmap",mea)
}
#' @export
MEAmap.default = function(mea,FUN,label,...){
  if(is.list(mea)){
    mea = MEAlist(mea)
    MEAmap(mea,FUN,label,...)
  } else {
    stop("This function accepts only MEA objects (individual or a list of them). Please use readMEA() to import files")
  }
}

#' @export
MEAmap.MEAlist = function(mea,FUN,label,...){
  if(any(!sapply(mea, is.MEA)))
    stop("Some elements in the provided MEAlist object are not of class MEA. Please try reimporting your data using readMEA()")
  res = Map(function (x,i){
    dots = list(...)
    #iterate and recycle items in the dots (beacuse MEAmap is not really a map() call at the MEA object level)
    other = which(!sapply(dots, is.function))
    for(k in seq_along(other)){
      if (length(dots[[other[k]]] ) ==length(mea)){
        dots[[other[k]]] = dots[[other[k]]][i]
      } else if (length(other[[k]])==1){
        dots[[other[k]]] = dots[[other[k]]]
      } else stop ("... must contain either functions, or object of length 1 or of length == length(mea")
    }

    prog(i,length(mea))
    do.call(MEAmap, c(list(x, FUN, label),dots))
    # MEAmap(x,FUN,label,...)
  },mea,seq_along(mea))
  # attributes(res) = attributes(mea)
  res = MEAlist(res)
  # attr(res,"filter") = paste(attr(mea,"filter"),"-->", label)
  return(res)
}

#' @export
MEAmap.MEA = function(mea,FUN,label,...){

  # #debug
  # FUN = function(x){
  #     x[[1]] = sasFiltB(x[[1]], winSec = moving.average.win, sampRate = attr(mea,"sampRate"))
  #     x[[2]] = sasFiltB(x[[2]], winSec = moving.average.win, sampRate = attr(mea,"sampRate"))
  #     x
  #   }
  # # label="supahsmooth"

  # FUN = function(x ,mySD){x[[1]] = x[[1]]/2+mySD }
  # dots = list(mySD =sd)
  # label="asd"

  FUN = match.fun(FUN)
  dots = list(...)
  if(length(dots)>0){
    for(d in 1:length(dots)){
      if(is.function(dots[[d]])){
        myF = match.fun(dots[[d]])
        dots[[d]] = myF(mea$MEA)
      }
    }
  }
  res = mea
  res$MEA = do.call(FUN, c(list(res$MEA),dots))
  attr(res,"filter") = paste(attr(mea,"filter"),"-->", label)
  return(res)
}
