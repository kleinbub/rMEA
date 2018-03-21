#######################################################################
########## BOOTSTRAP ANALYSIS
#
# ToDo:
# metti i nomi della lista shuffled = uid1 + uid2
#

#' Shuffle MEA data
#'
#' This function recombines the s1 and s2 motion energy time-series between all \code{MEA} objects in the supplied list.
#' It is typically used to compare genuine synchrony of real data with pseudosynchrony of shuffled (recombined) data.
#'
#' @param mea an object of class \code{MEAlist} (see function \code{\link{readMEA}}).
#' @param size integer. The number of combinations to be returned.
#'
#' @details The shuffling process first creates all possible combinations between s1 and s2 of all \code{MEA} objects in the supplied list,
#' then removes the original pairings, and finally extracts the desired numbers of dyads without replacement.
#'
#' Note: all the ccf data, if present, are discarded from the shuffled objects and have to be calculated again using \code{\link{MEAccf}}
#'
#' @return an object of class \code{MEAlist} containing randomly combined dyads.
#' @examples
#' ## read the first 4 minutes of the normal sample
#' ##   (intake interviews of patients that carried on therapy)
#' path_normal <- system.file("extdata/normal", package = "rMEA")
#' mea_normal <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
#'                      s1Name = "Patient", s2Name = "Therapist",
#'                      idOrder = c("id","session"), idSep="_", skip=1, nrow = 6000)
#' mea_normal <- setGroup(mea_normal, "normal")
#'
#'## Create a shuffled sample
#'mea_rand = shuffle(mea_normal, 50)
#'
#'summary(mea_rand)
#' @export
#'
shuffle = function(mea, size) {
  #if(any(!sapply(mea, is.MEA)))
  if(!is.MEAlist(mea))
    stop("This function accepts only MEAlist objects. Please use readMEA() to import files")

  lrfn = mea
  cat("\r\nShuffling dyads:\r\n")


  #ll è una lista cont tutti i segnali, prima tutti gli s1, poi tutti gli s2
  ll = c(
    lapply(lrfn, function(x) { x$MEA[,1] }),
    lapply(lrfn, function(x) { x$MEA[,2] })
  )
  names(ll)[1:(length(ll)/2)] = paste(names(ll)[1:(length(ll)/2)],substr(attr(mea,"s1Name"),1,3) ,sep="_")
  names(ll)[(length(ll)/2+1):length(ll)] = paste(names(ll)[(length(ll)/2+1):length(ll)],substr(attr(mea,"s2Name"),1,3),sep="_")
  percs = sapply(mea,function(x){attr(x,"s1_%_movement")})
  percs = c(percs,sapply(mea,function(x){attr(x,"s2_%_movement")}) )
  # calcola le possibili combinazioni ed estendi la lista con le pseudo-diadi
  combo = utils::combn(1:length(ll),2)
  #remove real combination
  comborem = unlist(lapply(seq_along(lrfn), function(i){which(combo[1,] == i & combo[2,] == i+length(lrfn))}))
  combo = combo[,-comborem]
  if (length(combo)/2 <=size) {
    com = data.frame(t(combo))
  } else {
    zamp = sample(1:(length(combo)/2), size=size )
    com = data.frame(t(combo[,zamp]))
  }
  n_boot = nrow(com)

  res = lapply(1:n_boot, function(i){
    prog(i,n_boot)
    res = stats::na.omit(unequalCbind(ll[[com[i,1]]], ll[[com[i,2]]]))
    colnames(res) = c(attr(mea,"s1Name"),attr(mea,"s2Name"))
    MEA(res, sampRate=attr(mea,"sampRate"), id=lead0(i,4),
        session="1", group="random", s1Name="s1Random", s2Name="s2Random",
        filter = attr(mea,"filter")
        )

  })
  cat('\r\n',n_boot,"/",length(combo)/2,"possible combinations were randomly selected")
  MEAlist(res)

}
