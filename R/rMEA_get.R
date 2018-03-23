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
#' @param mea a single or a list of MEA objects
#' @return A string or a vector of strings containing the metadata.
#' @details if a well formatted list of MEA objects is provided, the function returns a vector of results.
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
  sapply(mea, attr, "sampRate")
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
  sapply(mea, attr, "s1Name")
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
  sapply(mea, attr, "s2Name")
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
