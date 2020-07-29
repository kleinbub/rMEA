

#' Substitute values from MEA data
#' @description This function allows to substitute MEA data from a list of time epochs. This is useful to mark and remove artefacts,
#' or to substitute extreme values.
#'
#' @param mea an object of class \code{MEA} or a list of \code{MEA} objects (see function \code{\link{readMEA}})
#' @param epochs a data.frame containing a list of epochs that must be changed (see Details)
#' @param replacement the value used to mark artefacts. Use 'NA' to remove artefacts and '0' to apply thresholds.
#' Other values are allowed but should not be used without a good reason.
#' @param filterLabel can be used to update the \code{filter} attribute, to keep track of the data transformations.
#'
#' @details the \code{artefacts} data.frame must contain a "start" and "end" columns, with the boundaries of the epochs that must be
#' marked as artefacts. The start and end values can be either integers (denoting seconds), or string values in the format hh:mm:ss, or mm:ss.
#' Furthermore the data.frame must contain a \code{uid} column containing strings in the format "group_id_session", OR three columns \code{group, id, session}
#' presenting the information separately. Data identifiers must match those of the \code{mea} object.
#'
#' The data.frame can be either hand crafted, for instance by importing a csv file (see \code{\link[utils]{read.table}}), or
#' generated through the packages own artefact detection tools such as \code{\link{CCFartefacts}}
#' @return returns the same \code{MEA} or \code{MEAlist} object, with all artefactual data substituted.
#' @export
#'
#' @examples
#' #' ## read the first 4 minutes of the normal sample
#' ## (intake interviews of patients that carried on therapy)
#' path_normal <- system.file("extdata/normal", package = "rMEA")
#' mea_normal <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
#'                       s1Name = "Patient", s2Name = "Therapist",
#'                       idOrder = c("id","session"), idSep="_", skip=1, nrow = 6000)
#' mea_normal <- MEAccf(mea_normal, lagSec = 5, winSec = 30, incSec = 10, ABS = FALSE)
#'
#' ## find potential artefacts
#' artefacts = CCFartefacts(mea_normal, threshold = 0.2, delta=1)
#'
#' ##replace values
#' mea_replaced <- MEAreplace(mea_normal, epochs = artefacts,
#'                            replace = NA, filterLabel = "artefacts deletion")
#'
#' #visualize results on first case
#' plot(mea_normal$all_200_1)
#' plot(mea_replaced$all_200_1)
#'
#'
MEAreplace = function(mea, epochs, replacement, filterLabel = "replaced"){
  # path_normal <- system.file("extdata/normal", package = "rMEA")
  # mea_normal <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
  #                       s1Name = "Patient", s2Name = "Therapist",
  #                       idOrder = c("id","session"), idSep="_", skip=1, nrow = 6000)
  # mea_normal <- MEAccf(mea_normal, lagSec = 5, winSec = 30, incSec = 10, ABS = FALSE)
  # mea_normal <- setGroup(mea_normal, "normal")
  # epochs = CCFartefacts(mea_normal, threshold = 0.2, delta=1)
  # mea = mea_normal
  # x = mea_normal[[1]]
  # replacement=NA


  #validate epochs
  if(!is.data.frame(epochs)) stop("epochs must be a dataframe")
  if(!"uid" %in%colnames(epochs)){
    if(all(c("group","id","session") %in% colnames(epochs))){
      epochs$uid = apply(epochs, 1, function(x)paste(x["group"],x["id"],x["session"],sep="_"))
    } else stop("epochs must contain either a 'uid' column, or all of 'group', 'id', and 'session' columns")
  }
  uids = sapply(mea, uid )
  srs = sapply(mea, sampRate )
  res = MEAmap(mea, function(x,y,sr){
    ex = epochs[epochs$uid == y,]
    if(nrow(ex)>0){ #there were matches
      ex$st_s = timeMaster(ex$start,0,"sec")*sr
      ex$en_s = timeMaster(ex$end,0,"sec")*sr
      for(i in 1:nrow(ex)){
        x[[1]][ex$st_s[i]:ex$en_s[i]] =replacement
        x[[2]][ex$st_s[i]:ex$en_s[i]] =replacement
      }

    } else {
      #return x without changes
    }
    x
  }, uids, srs,
  label= paste(as.character(replacement))
  )

}


#' Detection of potential artefacts in CCF results.
#'
#' @description High synchronization values for extended time durations may be the effect of artefacts in the MEA data. For instance subject 2 movement may have
#' been captured by subject 1's camera (or ROI) as well, or some ambiental characteristic (e.g. light) is changing for both cameras (or ROIs). This function identifies
#' those moments allowing to inspect the original videos with temporal precision.
#' **Please note that is impossible to discriminate real high synchronization phenomena from artefacts without inspecting the
#' original videos.**
#'
#' @param mea an object of class \code{MEA} or a list of \code{MEA} objects (see function \code{\link{readMEA}})
#' @param threshold A single numeric value specifying the absolute correlation value above (and below) which a window must be considered artefactual.
#' @param delta Integer. The maximum numbers of consecutive CCF windows below threshold to be allowed in an artefactual streak without determining its end. A value of 1 is
#' default and recommended, but it can be increased if the reports become too long (e.g. with very noisy source material, or when working with very small ccf windows),  to achieve
#' the desired level of reporting detail.
#' @param duration Integer. Minimum duration of the artefacts to be reported. Note that artefacts smaller than the CCF windows increments cannot be detected.
#' @details The function only considers lag_zero correlations as MEA artefacts are expected to be non-lagged phenomena.
#' @return a data.frame object with all potential artefact epochs
#' @export
#'
#' @examples
#' ## read the first 4 minutes of the normal sample
#' ## (intake interviews of patients that carried on therapy)
#' path_normal <- system.file("extdata/normal", package = "rMEA")
#' mea_normal <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
#'                       s1Name = "Patient", s2Name = "Therapist",
#'                       idOrder = c("id","session"), idSep="_", skip=1, nrow = 6000)
#' mea_normal <- MEAccf(mea_normal, lagSec = 5, winSec = 30, incSec = 10, ABS = FALSE)
#'
#' ## find potential artefacts with various granularity of reporting
#' print(CCFartefacts(mea_normal, threshold = 0.2, delta=1))
#' print(CCFartefacts(mea_normal, threshold = 0.2, delta=5))
#' print(CCFartefacts(mea_normal, threshold = 0.2, delta=5, duration=60))
#'
CCFartefacts = function(mea, threshold, delta=1, duration=attributes(mea)$ccf$inc){
  if(!hasCCF(mea)) stop("CCFartefacts requires MEA objects with valid CCF informations. Refer to function MEAccf()", call.=F)
  final = data.frame()
  for(i in 1:length(mea)) {
    x = mea[[i]]
    #where lag0 > threshold?
    dud = which(x$ccf$lag0 > threshold)
    #if we find results...
    if(length(dud)>0){
      res = lag = c(max(delta+1,999),diff(dud)) #save the distance between matches (first value must always be huge to become START)
      res[lag>delta & c(lag[c(2:length(lag))],0) > delta ] = "startstop" #isolated high windows
      res[lag>delta & c(lag[c(2:length(lag))],0) <= delta ] = "start" #start of lenghty high sync segments
      res[lag<=delta & c(lag[c(2:length(lag))],0) > delta ] = "stop"  #end of segments
      if(res[length(res)] == "start" ) res[length(res)] = "startstop" #fix last line
      if(!res[length(res)] %in% c("start","stop","startstop") ) res[length(res)] = "stop" #fix last line
      # if(length(res)==1)  res = "startstop"
      a = data.frame(dud,lag,res)
      a= a[a$res %in% c("start","stop","startstop"),]
      if(nrow(a)==1) a$res[1] = "startstop"

      #reshape the results
      RES = data.frame(
        uid = attr(x, "uid"),
        start=c(a$dud[a$res=="start"], a$dud[a$res=="startstop"]),
        end=c(a$dud[which(a$res=="start")+1], a$dud[a$res=="startstop"])

      )
      RES = RES[order(RES$start),] #sort
      RES$start = timeMaster((RES$start-1)*attr(x,"ccf")$inc,out = "h") #dud = 1 --> 00:00
      RES$end = timeMaster((RES$end)*attr(x,"ccf")$inc,out = "h")
      RES$duration = timeMaster(timeMaster(RES$end,out="s")-timeMaster(RES$start,out="s"), out="s")
      final = rbind(final, RES)

    }

  }
  final = final[final$duration>=duration,]

}

# #' Detection of extreme values in MEA data
# #'
# #' Sometimes motion energy analysis generates excessively high peaks resulting from video artifacts or other anomalies in the video source.
# #'
# #' This function allows to identify the values greater or less than a specific threshold.
# #' The default threshold is 10 times the standard deviation of the time-series.
# #'
# #' The identified values can later be substituted with a fixed value (eg. NA or 0) using the \code{\link{MEAreplace}} function.
# #'
# #' @param mea an object of class \code{MEA} or a list of \code{MEA} objects (see function \code{\link{readMEA}})
# #' @param threshold a numeric value, or a function returning the threshold value to consider data as outliers.
# #' @param direction a text string. One of "\code{greater}" or "\code{less}": can be abbreviated.
# #'
# #' @return The same \code{mea} object with all extreme values substituted.
# #' @examples
# #' ## read the first 4 minutes of the normal sample
# #' ##   (intake interviews of patients that carried on therapy)
# #' path_normal <- system.file("extdata/normal", package = "rMEA")
# #' mea_raw <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
# #'                      s1Name = "Patient", s2Name = "Therapist",
# #'                      idOrder = c("id","session"), idSep="_", skip=1, nrow = 6000)
# #'
# #' ## Remove extreme values, higher than 10 times the standard deviation
# #' mea_clean = MEAoutlier(mea_raw, threshold=function(x){sd(x)*10}, direction = "greater")
# #' @export
# MEAoutlier = function(mea, threshold=function(x){stats::sd(x)*10}, direction=c("greater", "less") ,replace=NA){
#   cat0("\r\nSuppressing outliers to ",as.character(replace),":\r\n")
#   direction = match.arg(direction)
#   res = MEAmap(mea, function(x){
#     if(is.function(threshold)){
#       threshold = match.fun(threshold)
#       threshold = threshold(c(x[[1]],x[[2]]))
#     }
#     if (length(threshold)>1) stop("Threshold must be a single value, not ",length(threshold), call.=F)
#     if(direction == "greater"){
#       x[[1]][x[[1]]>threshold] =replace
#       x[[2]][x[[2]]>threshold] =replace
#     } else if (direction =="less") {
#       x[[1]][x[[1]]<threshold] =replace
#       x[[2]][x[[2]]<threshold] =replace
#     } else stop("'direction' argument must be either 'greater' or 'less'.",call.=F)
#     x
#   },
#   label= paste("outliers to",as.character(replace))
#   )
# }


