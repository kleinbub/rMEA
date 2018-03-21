#' Moving average smoothing for motion energy data
#'
#' This function applies a moving average filter, based on SAS "proc expand" procedure. The moving average is applied independently on each subject's motion energy data.
#' NA values are set to 0.
#'
#' @param mea an object of class \code{MEA} or \code{MEAlist} (see function \code{\link{readMEA}})
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
    x[[1]] = sasFiltB(x[[1]], winSec = moving.average.win, sampRate = attr(mea,"sampRate"))
    x[[2]] = sasFiltB(x[[2]], winSec = moving.average.win, sampRate = attr(mea,"sampRate"))
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
#' @param mea an object of class \code{MEA} or \code{MEAlist} (see function \code{\link{readMEA}})
#' @param scaleFUN the function to be applied to each motion energy time-series to calculate a scaling factor. Default is standard deviation.
#' @param ... further arguments passed to \code{FUN}
#' @param center either a logical value or a numeric vector of length 2 specifying separate centering values for s1 and s2.
#' @details \code{FUN} is found by a call to  \code{\link[base]{match.fun}} and typically is either a function
#' or a symbol (e.g., a backquoted name) or a character string specifying a function
#' to be searched for from the environment of the call to apply.
#' If a \code{na.rm} argument is present in \code{FUN}
#' it is automatically set to TRUE.
#'
#' \code{center} is directly passed to \code{\link[base]{scale}}. If \code{center} is \code{TRUE} then centering
#' is done by subtracting the means (omitting NAs) from the motion energy time-series. If a \code{center} is a numeric vector,
#' the first value will be subtracted from s1 and the second from s2.
#' Note: the s1 and s2 signals are scaled independently.
#' @return returns the same \code{MEA} or \code{MEAlist} object, with all motion energy data rescaled
#' @importFrom stats sd
#' @examples
#' ## read the first 4 minutes of the normal sample
#' ##   (intake interviews of patients that carried on therapy)
#' path_normal <- system.file("extdata/normal", package = "rMEA")
#' mea_raw <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
#'                      s1Name = "Patient", s2Name = "Therapist",
#'                      idOrder = c("id","session"), idSep="_", skip=1, nrow = 6000)
#'
#' ## filter with moving average
#' mea_scaled = MEAscale(mea_raw, scaleFUN = sd)
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

MEAscale = function(mea, scaleFUN=sd, ..., center=F){
  cat("\r\nRescaling data:\r\n")
  fil = paste0("rescaled (", as.character(substitute(scaleFUN)),")" )
  scaleFUN = match.fun(scaleFUN)
  dots = list(...)
  if("na.rm" %in%  methods::formalArgs(scaleFUN) && !"na.rm"%in% names(dots)){ #if the function allows na.rm, set it to TRUE
    dots[["na.rm"]] = T
  }
  res = MEAmap(mea, function(x){
    x[[1]]=scale(x[[1]], scale = do.call(scaleFUN, c(list(x[[1]]),dots)), center = center)
    x[[2]]=scale(x[[2]], scale = do.call(scaleFUN, c(list(x[[2]]),dots)), center = center)
    x
  },
  label=fil)
}

#' Replace outliers with given values
#'
#' Sometimes motion energy analysis generates excessively high peaks resulting from video artifacts or other anomalies in the video source.
#'
#' This function allows to substitute the values greater or less than a specific threshold. The default threshold is 10 times the standard deviation of the time-series.
#'
#' @param mea an object of class \code{MEA} or \code{MEAlist} (see function \code{\link{readMEA}}).
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



#' Apply a function to a MEA or MEAlist object
#'
#' MEApply is a wrapper to \code{\link[base]{do.call}} that allows to apply a function on the motion energy data of one or multiple \code{MEA} objects. Complex constructions are possible, see details.
#'
#' @param mea an object of class \code{MEA} or \code{MEAlist} (see function \code{\link{readMEA}})
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
MEAmap.MEAlist = function(mea,FUN,label,...){
  if(any(!sapply(mea, is.MEA)))
    stop("Some elements in the provided MEAlist object are not of class MEA. Please try reimporting your data using readMEA()")
  res = Map(function (x,i){
    prog(i,length(mea))
    MEAmap(x,FUN,label,...)
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
