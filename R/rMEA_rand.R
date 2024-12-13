#######################################################################
########## BOOTSTRAP ANALYSIS
#
# ToDo:
# metti i nomi della lista shuffled = uid1 + uid2
#

#' Shuffle MEA data (between subjects)
#'
#' This function recombines the s1 and s2 motion energy time-series between all \code{MEA} objects in the supplied list.
#' It is typically used to compare genuine synchrony of real data with pseudosynchrony of shuffled (recombined) data.
#'
#' @param mea a list of \code{MEA} objects (see function \code{\link{readMEA}}).
#' @param size either "max" or an integer specifying the number of combinations to be returned.
#' @param keepRoles Boolean. If TRUE the resulting random dyad will preserve the roles, i.e. they will all have a new s1 sampled among all
#' s1s and a new s2 sampled among all s2s. If FALSE (default), the role will be disregared.
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
shuffle = function(mea, size="max", keepRoles=FALSE) {

  lrfn = MEAlist(mea)
  cat("\r\nShuffling dyads:\r\n")


  #ll è una lista cont tutti i segnali, prima tutti gli s1, poi tutti gli s2
  ll = c(
    lapply(lrfn, function(x) { x$MEA[,1] }),
    lapply(lrfn, function(x) { x$MEA[,2] })
  )
  #build names containing the original info
  names(ll)[1:(length(ll)/2)] = paste(names(ll)[1:(length(ll)/2)],substr(attr(mea,"s1Name"),1,3) ,sep="_")
  names(ll)[(length(ll)/2+1):length(ll)] = paste(names(ll)[(length(ll)/2+1):length(ll)],substr(attr(mea,"s2Name"),1,3),sep="_")

  # calcola le possibili combinazioni ed estendi la lista con le pseudo-diadi
  combo = utils::combn(1:length(ll),2)

  #remove real combination
  comborem = unlist(lapply(seq_along(lrfn), function(i){which(combo[1,] == i & combo[2,] == i+length(lrfn))}))
  combo = combo[,-comborem]

  #per request of thekryz on github
  if(keepRoles){
    comborem = c(
                  which(combo[1,] <= length(lrfn) & combo[2,] <= length(lrfn)), #dove s1 <-> s1
                  which(combo[1,] > length(lrfn) & combo[2,] > length(lrfn))
                )
    combo = combo[,-comborem]
  }
  if (size=="max" || length(combo)/2 <=size) {
    com = data.frame(t(combo))
  } else {
    zamp = sample(1:(length(combo)/2), size=size )
    com = data.frame(t(combo[,zamp]))
  }
  n_boot = nrow(com)

  res = lapply(1:n_boot, function(i){
    prog(i,n_boot)
    res = stats::na.omit(unequalCbind(ll[[com[i,1]]], ll[[com[i,2]]], keep=FALSE))
    colnames(res) = c(attr(mea,"s1Name"),attr(mea,"s2Name"))
    MEA(res, sampRate=attr(mea,"sampRate"), id=lead0(i,4),
        session="1", group="random", s1Name="s1Random", s2Name="s2Random",
        filter = attr(mea,"filter"), uid = paste(names(ll)[com[i,1]],names(ll)[com[i,2]],sep="_|_")
        )

  })
  cat('\r\n',n_boot,"/",length(combo)/2,"possible combinations were randomly selected")

  #debug comborem using demo data
  # x =   MEAlist(res)
  # x = names(x)
  # x = gsub('normal_20','',x,fixed = T)
  # x = gsub('_1_P','P',x,fixed = T)
  # x = gsub('_1_T','T',x,fixed = T)
  # x = gsub('_|_',' | ',x,fixed = T)
  # x = gsub('The','_S2',x,fixed = T)
  # x = gsub('Pat','_S1',x,fixed = T)
  # data.frame(x)

  MEAlist(res)

}

#' Shuffle MEA data (within subjects)
#'
#' This function generates fakes dyads to be used for pseudosynchrony calculations following the Ramseyer & Tschacher (2010)
#' within-subject segment shuffling approach.
#' Between subjects shuffling \code{\link{shuffle}} is probably more conservative, and suggested for most cases.
#' This function is provided for replicability of older studies, and can be useful to quickly assess pseudosynchrony in single sessions,
#' or very small samples.
#'
#' @param mea a list of \code{MEA} objects (see function \code{\link{readMEA}}).
#' @param n_each the number of random dyads to be generated from each real dyad.
#' @param segSec the width (in seconds) of the shuffling segments.
#'
#' @details For each \code{MEA} object, the shuffling procedure first divides s1 and s2 MEA data in segments of size \code{segSec},
#' then shuffles them within subject (so that the new segments of s1, are the old segments of s1 in a new order). This is repeated
#' for \code{n_each} times, before getting to the next \code{MEA} object
#'
#' Note: all the ccf data, if present, are discarded from the shuffled objects and have to be calculated again using \code{\link{MEAccf}}
#'
#' @return an object of class \code{MEAlist} containing \code{n_each * length(mea)} random dyads.
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
#'mea_rand = shuffle_segments(mea_normal, n_each=10, segSec=30)
#'
#'summary(mea_rand)
#' @export
#'
shuffle_segments = function(mea, n_each, segSec){
  message("new version. please carefully check data! DUDE!")

  # Capture sampling rate and filter attribute from the original MEA
  sr = sampRate(mea)
  filter_attr = paste(attr(mea, "filter"), "--> within shuffle" ) # Capture the filter attribute
  RES = list()  # Initialize result list for flattened shuffled MEA objects

  # Loop through each MEA object in `mea`
  for(j in seq_along(mea)) {
    # Extract original MEA data, ID, and session for current object
    x <- mea[[j]]$MEA
    id <- attr(mea[[j]], "id")  # Retrieve the ID attribute
    session <- attr(mea[[j]], "session")  # Retrieve the session attribute

    # Perform `n_each` shuffles for the current MEA object
    for(i in 1:n_each) {
      # Define a unique identifier for each shuffle using ID and session
      shuffle_id <- paste0("shuffle_", id, "_session_", session, "_", sprintf("%04d", i))

      # Calculate segment size and number of segments
      seg = (sr * segSec)
      nseg = ceiling((nrow(x)-seg+1)/seg)
      startx = (seq_len(nseg) - 1) * seg + 1
      endx = seq_len(nseg) * seg

      # Shuffle segment indices independently for s1 and s2
      ran1 = sample(1:nseg)
      ran2 = sample(1:nseg)

      # Initialize shuffling index vectors
      val1 = val2 = numeric(nrow(x))
      for(k in 1:nseg) {
        val1[startx[k]:endx[k]] = startx[ran1[k]]:endx[ran1[k]]
        val2[startx[k]:endx[k]] = startx[ran2[k]]:endx[ran2[k]]
      }
      # pad=runif(1,0,nrow(x))
      # plot(val1,t="l",xlim=c(100,500)+pad,ylim=c(0,1))
      # abline(v=startx,col=2,lwd=3)
      # abline(v=endx,col=4  ,lwd=3)
      # lines(val1,lwd=2 )

      # Apply shuffled indices and trim to segment length
      x_shuffled = x[1:nrow(x), ]
      x_shuffled[[1]] = x[[1]][val1]
      x_shuffled[[2]] = x[[2]][val2]

      # Create a new MEA object with the unique name and filter attribute
      shuffled_mea = MEA(x_shuffled, sampRate = sr, id = shuffle_id, session = session, group = "shuffled",
                         s1Name = attr(mea, "s1Name"), s2Name = attr(mea, "s2Name"), filter = filter_attr)

      # Add the shuffled MEA object directly to the flat RES list with unique name
      RES[[shuffle_id]] <- shuffled_mea
    }
  }

  # Apply consistent filter attribute to all MEA objects in RES
  RES = lapply(RES, function(x) {
    attr(x, "filter") <- filter_attr
    if ("ccf" %in% names(x)) x["ccf"] = list(NULL)
    if ("ccfRes" %in% names(x)) x["ccfRes"] = list(NULL)
    return(x)
  })

  cat("Final shuffled list created with", length(RES), "unique shuffled MEA objects.\n")
  return(RES)
}
