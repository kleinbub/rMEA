# id = id,
# session = session,
# group = group,
# sampRate = sampRate,
# filter = filter,
# ccf = "",
# s1Name = s1Name,
# s2Name = s2Name,
# uid = paste(group,id,session,sep="_")


#' Get MEA attributes
#'
#' @param mea an object of class \code{MEA} or a list of \code{MEA} objects (see function \code{\link{readMEA}})
#' @return A string or a vector of strings containing the metadata.
#' @details if a well formatted list of MEA objects is provided, the function returns a vector
#' of results for id, session, group and uid. sampRate, s1Name, and s2Name return always a single
#' value, as they are not allowed to be mixed.
#' @export
id <- function(mea) {
  UseMethod("id", mea)
}
#' @export
id.MEA <- function(mea){
  attr(mea, "id")
}
#' @export
id.default <- function(mea){
  if (is.list(mea)) mea = MEAlist(mea)
  sapply(mea, attr, "id")
}

#' @rdname id
#' @export
group <- function(mea) {
  UseMethod("group", mea)
}
#' @export
group.MEA <- function(mea){
  attr(mea, "group")
}
#' @export
group.default <- function(mea){
  if (is.list(mea)) mea = MEAlist(mea)
  sapply(mea, attr, "group")
}

#' @rdname id
#' @export
session <- function(mea) {
  UseMethod("session", mea)
}
#' @export
session.MEA <- function(mea){
  attr(mea, "session")
}
#' @export
session.default <- function(mea){
  if (is.list(mea)) mea = MEAlist(mea)
  sapply(mea, attr, "session")
}

#' @rdname id
#' @export
sampRate <- function(mea) {
  UseMethod("sampRate", mea)
}
#' @export
sampRate.MEA <- function(mea){
  attr(mea, "sampRate")
}
#' @export
sampRate.default <- function(mea){
  if (is.list(mea)) mea = MEAlist(mea)
  attr(mea, "sampRate")
}

#' @rdname id
#' @export
s1Name <- function(mea) {
  UseMethod("s1Name", mea)
}
#' @export
s1Name.MEA <- function(mea){
  attr(mea, "s1Name")
}
#' @export
s1Name.default <- function(mea){
  if (is.list(mea)) mea = MEAlist(mea)
  attr(mea, "s1Name")
}

#' @rdname id
#' @export
s2Name <- function(mea) {
  UseMethod("s2Name", mea)
}
#' @export
s2Name.MEA <- function(mea){
  attr(mea, "s2Name")
}
#' @export
s2Name.default <- function(mea){
  if (is.list(mea)) mea = MEAlist(mea)
  attr(mea, "s2Name")
}

#' @rdname id
#' @export
uid <- function(mea) {
  UseMethod("uid", mea)
}
#' @export
uid.MEA <- function(mea){
  attr(mea, "uid")
}
#' @export
uid.default <- function(mea){
  if (is.list(mea)) mea = MEAlist(mea)
  sapply(mea, attr, "uid")
}


#' Extract ccf values from MEA objects
#'
#' @param mea an object of class \code{MEA} or a list of \code{MEA} objects (see function \code{\link{readMEA}})
#' @param type A character vector defining which ccf must be extracted.
#' Either "matrix", one of the ccfRes indexes identified with \code{\link{ccfResNames}}
#' or the name of one lag value which can be identified with \code{\link{lagNames}}
#'
#' @return If \code{type="matrix"}, the whole ccf matrix is returned. Otherwise a vector containing the ccf
#' time-series for the selected lag, or aggregated values is returned.
#' If \code{mea} is a list, the return value is a list of the individual ccf of each MEA object.
#' @export

getCCF <- function (mea, type) {
  UseMethod("getCCF", mea)
}
#' @export
getCCF.MEA <- function (mea, type) {
  if (!hasCCF(mea)) stop ("No ccf computation found, please refer to MEAccf() function.")
  if (type %in% lagNames(mea)) {
    return (mea$ccf[[type]])
  } else if (type %in% names(mea$ccfRes)) {
    return (mea$ccfRes[[type]])
  } else if (type == "matrix") {
    return (mea$ccf)
  } else stop ("'type' must be either \"matrix\", a lag label extracted with lagNames(), or one of:\r\n\"",paste0(ccfResNames(mea),collapse = "\", \""),"\"", call.=F)
}
#' @export
getCCF.default <- function (mea, type) {
  if (is.list(mea)) mea = MEAlist(mea)
  mea <- MEAlist(mea)
  sapply(mea, getCCF, type)
}


#' Extract the lag names of a ccf analysis in MEA objects
#'
#' @param mea an object of class \code{MEA} or a list of \code{MEA} objects (see function \code{\link{readMEA}})
#'
#' @return a vector containing the labels of the lag values
#' @export
lagNames <- function (mea) {
  UseMethod("lagNames", mea)
}
#' @export
lagNames.MEA <- function (mea) {
  if (!hasCCF(mea)) stop("No ccf computation found, please refer to MEAccf() function.")
  names(mea$ccf)
}
#' @export
lagNames.default <- function (mea){
  if (is.list(mea)) mea = MEAlist(mea)
  mea <- MEAlist(mea)
  names(mea[[1]]$ccf)
}

#' Extract the names of the ccf analysis summaries in a MEA objects
#'
#' @param mea an object of class \code{MEA} or a list of \code{MEA} objects (see function \code{\link{readMEA}})
#'
#' @return a vector containing the labels of the ccfRes indexes
#' @export
ccfResNames <- function (mea) {
  UseMethod("ccfResNames", mea)
}
#' @export
ccfResNames.MEA <- function (mea) {
  if (!hasCCF(mea)) stop("No ccf computation found, please refer to MEAccf() function.")
  names(mea$ccfRes)
}
#' @export
ccfResNames.default <- function (mea){
  if (is.list(mea)) mea = MEAlist(mea)
  mea <- MEAlist(mea)
  names(mea[[1]]$ccfRes)
}


