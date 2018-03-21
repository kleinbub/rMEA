#' \code{MEA} class constructor
#'
#' The preferred way to create an object of class \code{MEA} is through the function \code{\link{readMEA}}.
#'
#' @param dataframe a data frame with 2 columns containing MEA data respectively for subject 1 (s1) and subject 2 (s2).
#' @param sampRate integer. The sampling rate of the MEA data. Normally derived from the framerate of the analyzed video sequence (frames per second; fps)
#' @param filter a string describing the pre-processing that has been applied on the raw data.
#' @param id a string representing a unique identifier of the dyad that the MEA data belong to.
#' @param session an integer representing the session (or experiment, interaction, etc); if each dyad is measured only once, the default value is 1.
#' @param group a string naming the group the dyad belongs to, such as diagnostic group, clinic, etc.
#' @param s1Name a string naming subject 1
#' @param s2Name a string naming subject 2
#' @param x object to be tested.
#' @details It is advised to \strong{not} create the \code{MEA} object manually but to always use the function \code{\link{readMEA}} instead.
#'
#' @return
#' A list containing three objects:
#'
#'     MEA: the dataframe containing the motion energy data
#'
#'     ccf: the matrix of lagged cross-correlations between s1 and s2 (if \code{\link[rMEA]{MEAccf}} was run)
#'
#'     ccfRes: some useful row marginals
#' @export

MEA = function(dataframe, sampRate, filter = "raw",id,
               session, group, s1Name, s2Name
               ){
  x = list("MEA"=dataframe,
           "ccf" = NULL,
           "ccfRes" = NULL
  )
  attributes(x) = c(attributes(x),list(
    id = id,
    session = session,
    group = group,
    sampRate = sampRate,
    filter = filter,
    ccf = "",
    s1Name = s1Name,
    s2Name = s2Name,
    uid = paste(group,id,session,sep="_")
  ))
  class(x) = c("MEA",class(x))
  return(x)
}

#' Well formatted list of \code{MEA} objects
#'
#' This constructor function checks if all the supplied \code{MEA} objects share the same sampling rate, pre-processing, and metadata,
#' and returns an object with additional attributes summarizing the contained \code{MEA} objects.
#'
#' @param listOfMea a list containing \code{MEA} objects
#' @param x object to be tested.
#' @return an object of class \code{MEAlist}
#' @export
MEAlist = function(listOfMea){
  if(any(!sapply(listOfMea, is.MEA) ) ) stop("Only a list of MEA objects can be used to create a MEAlist", call.=F)
  class(listOfMea) = "MEAlist"
  attributes(listOfMea) = c(attributes(listOfMea), list(
    nId = length(unique(sapply(listOfMea, attr, "id"))),
    n = length(listOfMea),
    groups = unique(sapply(listOfMea, attr, "group")),
    #tutte le linee qua sotto controllano che tutti gli oggetti MEA abbiano un valore uguale per certi attributi, e lo salvano nelle info dell'oggetto MEALIST
    sampRate = ifelse(length(unique(sapply(listOfMea, attr, "sampRate")))==1, unique(sapply(listOfMea, attr, "sampRate")), stop("Cannot construct a MEAlist with different sampling rates") ),
    filter = ifelse(length(unique(sapply(listOfMea, attr, "filter")))==1, unique(sapply(listOfMea, attr, "filter")), stop("Cannot construct a MEAlist with different filtering procedures") ),
    s1Name = ifelse(length(unique(sapply(listOfMea, attr, "s1Name")))==1, unique(sapply(listOfMea, attr, "s1Name")), stop("Cannot construct a MEAlist with multiple s1Name labels") ),
    s2Name = ifelse(length(unique(sapply(listOfMea, attr, "s2Name")))==1, unique(sapply(listOfMea, attr, "s2Name")), stop("Cannot construct a MEAlist with multiple s2Name labels") ),
    ccf = if(all(  sapply(listOfMea, function(x) !is.null(x$ccf)) )) {# se tutti hanno CCF
      list(
        filter = ifelse(length(unique(sapply(listOfMea, function(x) attr(x, "ccf")$filter )))==1, unique(sapply(listOfMea, function(x)  attr(x, "ccf")$filter)), stop("Cannot construct a MEAlist with CCF computed with different settings") ),
        lag    = ifelse(length(unique(sapply(listOfMea, function(x) attr(x, "ccf")$lag    )))==1, unique(sapply(listOfMea, function(x)  attr(x, "ccf")$lag   )), stop("Cannot construct a MEAlist with CCF computed with different settings") ),
        win    = ifelse(length(unique(sapply(listOfMea, function(x) attr(x, "ccf")$win    )))==1, unique(sapply(listOfMea, function(x)  attr(x, "ccf")$win   )), stop("Cannot construct a MEAlist with CCF computed with different settings") ),
        inc    = ifelse(length(unique(sapply(listOfMea, function(x) attr(x, "ccf")$inc    )))==1, unique(sapply(listOfMea, function(x)  attr(x, "ccf")$inc   )), stop("Cannot construct a MEAlist with CCF computed with different settings") )
      )
    } else if(all(  sapply(listOfMea, function(x) is.null(x$ccf)) )){
      ""
    } else stop("Cannot construct a MEAlist when only a part of the MEA objects have ccf analyses")

  ))
  names(listOfMea) = sapply(listOfMea, attr, "uid")
  return(listOfMea)
}

#' @rdname MEA
#' @return \code{is.MEA} returns TRUE if and only if its argument is of class \code{MEA}
#' @export
is.MEA = function(x) inherits(x,"MEA") && length(x)

hasCCF = function(x){
  if(is.MEA(x))
   !is.null(x$ccf)
  else if(is.list(x))
    all(sapply(x, hasCCF))
  else stop()
}

#' @rdname MEAlist
#' @return \code{is.MEAlist} returns TRUE if and only if its argument is of class \code{MEAlist}
#' @export
is.MEAlist = function(x) inherits(x,"MEAlist") && length(x)

#' @export
c.MEAlist = function(...){
  dots = list(...)
  # dots= list(group1,group2)
  dots = lapply(dots, function(x){attributes(x)=NULL;x})
  MEAlist(do.call("c",dots))
}


#' Sets the group of MEA objects
#' @param mea a single or a list of \code{MEA} objects (see function \code{\link{readMEA}})
#' @param group a text string specifying a group name
#' @return an object of the same type of \code{'mea'}, with the group attributes set to \code{group}.
#' @examples ## read a sample
#' path_normal <- system.file("extdata/normal", package = "rMEA")
#' mea_normal <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
#'                      s1Name = "Patient", s2Name = "Therapist",
#'                      idOrder = c("id","session"), idSep = "_",  skip = 1)
#' mea_normal <- setGroup(mea_normal, "normal")
#' @export
setGroup = function(mea,group){
  UseMethod("setGroup", mea)
}
#' @export
setGroup.MEA = function(mea,group){
  attr(mea,"uid") = paste(group, attr(mea,"id"), attr(mea, "session"),sep ="_")
  attr(mea,"group") = group
  mea
}

#' @export
setGroup.default = function(mea,group){
  #NB: questa funzione è così perché funzioni anche su liste di MEA normali
  if(is.list(mea)){
    newNames = vector("character",length(mea))
    for(i in 1:length(mea)){
      if(is.MEA(mea[[i]])) {
        mea[[i]] = setGroup(mea[[i]],group)
        newNames[i] = attr(mea[[i]],"uid")
      }
      else stop("An object of a class different from MEA was found in the list")
    }
  } else stop("unrecognized format")
  MEAlist(mea)
}



#' @export
summary.MEAlist = function(object, ...){
  filters = sapply(object, attr, "filter")
  if (length(unique(filters))>1) stop("Different processing pipeline found in MEAlist:\r\n",unique(filters),call. = F)
  s1_perc = sapply(object, function(x){length(x$MEA[[1]][x$MEA[[1]]>0])/nrow(x$MEA)})#
  s2_perc = sapply(object, function(x){length(x$MEA[[2]][x$MEA[[2]]>0])/nrow(x$MEA)})#

  noCCF = any(sapply(object, function(x){is.null(x$ccf)}))
  if(noCCF){ ## output solo per i dati importati
    Q = data.frame(
      "dyad"      = factor(sapply(object, function(mea){attr(mea, "id")})),
      "session"   = factor(sapply(object, function(mea){attr(mea, "session")})),
      "group"     = factor(sapply(object, function(mea){attr(mea, "group")})),
      "s1_perc"  =  round(s1_perc*100,1),
      "s2_perc"  =  round(s2_perc*100,1)
    )
    names(Q)[names(Q) == 's1_perc'] = paste0(attr(object, "s1Name"),"_%")
    names(Q)[names(Q) == 's2_perc'] = paste0(attr(object, "s2Name"),"_%")

    cat("\r\nMEA analysis results:\r\n")
    print(Q)
    cat("\r\nData processing: ",attr(object, "filter"))

  } else { ## output if you have calculated CCF
    pace  = sapply(object, function(mea){mean(unlist(mea$ccf[, 1:(attr(mea, "ccf")$lag*attr(mea, "sampRate"))]),na.rm=T)} )
    zero  = sapply(object, function(mea){mean(unlist(mea$ccf[, (attr(mea, "ccf")$lag*attr(mea, "sampRate")+1)]),na.rm=T)} )
    lead  = sapply(object, function(mea){mean(unlist(mea$ccf[, (attr(mea, "ccf")$lag*attr(mea, "sampRate")+2):(attr(mea, "ccf")$lag*attr(mea, "sampRate")*2+1) ]),na.rm=T)} )
    pace0 = sapply(object, function(mea){mean(unlist(mea$ccf[, 1:(attr(mea, "ccf")$lag*attr(mea, "sampRate") +1)]),na.rm=T)} )
    lead0 = sapply(object, function(mea){mean(unlist(mea$ccf[, (attr(mea, "ccf")$lag*attr(mea, "sampRate")   +1):(attr(mea, "ccf")$lag*attr(mea, "sampRate")*2+1) ]),na.rm=T)} )
    grandAver = unlist(lapply(object, function(mea) { mean(unlist(mea$ccf),na.rm = T)}))


    #  save/report/plot
    Q = data.frame(
      "dyad"      = factor(sapply(object, function(mea){attr(mea, "id")})),
      "session"   = factor(sapply(object, function(mea){attr(mea, "session")})),
      "group"     = factor(sapply(object, function(mea){attr(mea, "group")})),
      "s1_perc"   =  round(s1_perc*100,1),
      "s2_perc"   =  round(s2_perc*100,1),
      "all_lags"  = round(grandAver,4), #the average across all lags
      "s1_lead"   = round(lead,4), #the average of the lags > 0
      "s2_lead"    = round(pace,4), #the average of the lags < 0
      "lag_zero"   = round(zero,4) #the sync values at lag = 0
      #,"pacing_0"  = round(pace0,4) #the average of the lags <= 0
      #,"leading_0" = round(lead0,4) #the average of the lags >= 0
    )
    names(Q)[names(Q) == 's1_perc'] = paste0(attr(object, "s1Name"),"_%")
    names(Q)[names(Q) == 's2_perc'] = paste0(attr(object, "s2Name"),"_%")
    names(Q)[names(Q) == 's1_lead'] = paste0(attr(object, "s1Name"),"_lead")
    names(Q)[names(Q) == 's2_lead'] = paste0(attr(object, "s2Name"),"_lead")

    cat("\r\nMEA analysis results:\r\n")
    print(Q)
    cat("\r\nData processing: ",attr(object, "filter"))
    cat0("\r\nCCF settings:\r\nWindow = ",attr(object, "ccf")$win, " s | Increments = ",attr(object, "ccf")$inc," s | Lag = ",attr(object, "ccf")$lag, " s.")
    cat0("\r\n",attr(object, "ccf")$filter)

  }
  invisible(Q)

}

