# PLOTS for MEA

dotsList = function(parList, ...){
  dots = list(...)
  for(i in names(parList)){
    if (!i %in% names(dots)) dots[[i]] = parList[[i]]
  }
  dots
}

# rangeRescale <- function(x, rangeMin, rangeMax){
#   (rangeMax-rangeMin) * (
#     (x - min(x, na.rm=T))  / (max(x, na.rm=T) - min(x, na.rm=T) )
#   ) + rangeMin
# }


rangeRescale <- function(x, rangeMin, rangeMax, xmin =min(x, na.rm=T), xmax = max(x, na.rm=T), pres.signs=FALSE){
  #rangeMin e rangeMax indicano il minimo e il massimo della nuova scala
  #
  #se xmin e xmax mancano, vengono usati il minimo e il massimo del campione
  if(any(x>xmax, na.rm=T) || any(x<xmin, na.rm=T)) stop ("Value found outside xmax and ymin boundaries. xmax and xmin should be equal or larger than the data range.")
  if(pres.signs){
    # if((!missing(xmin) && !missing(xmax)) || (xmin!=-max(abs(x)) || xmax != max(abs(x)) ) )
    #   stop("Either x")
    if( rangeMin != -rangeMax )
      stop("with pres.sign = TRUE, rangeMin and rangeMax should be opposites (e.g. -1 and 1")
    mightyMax = max(abs(xmax),abs(xmin)) #così centra qualsiasi range?
    xmin = -mightyMax
    xmax = mightyMax
  }
  (rangeMax-rangeMin) * (
    (x - xmin)  / (xmax - xmin )
  ) + rangeMin
}


#' Transform color
#'
#' @param col a color to begin with in hex format
#' @param luminosity numeric. negative numbers darken the color, positive lighten it. Eg: a value of -2 make the color two times darker.
#' @param alpha numeric from 0 to 1. the value of opacity of the resulting color
#'
#' @return a color string
#'
colTrans <- function(col, luminosity=NA, alpha=NA){
  if(!is.na(luminosity)){
    col <- grDevices::col2rgb(col)

    if(luminosity<0)
      col <- col/abs(luminosity)
    else if(luminosity>0)
      col <-col*luminosity
    col <- grDevices::rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
  }
  if(!is.na(alpha)){
    col <- grDevices::col2rgb(col)
    alpha = 255*alpha
    col <- grDevices::rgb(col[1],col[2],col[3],alpha,maxColorValue = 255)
  }
  col
}



#' Plots the average cross-correlation at different lags
#'
#' Provides a graphical representation of the comparison between two lists of \code{MEA} objects.
#' The X-axis represents the lag values over which cross-correlation was calculated (in seconds), the Y-axis represents the averaged strength of the cross-correlation.
#' Typically, the  is useful for a visual inspection of the strength of synchrony from real dyads in relation to synchrony expected by coincidence (pseudosynchrony).
#'
#'
#' @param mea a list of \code{MEA} objects (see function \code{\link{MEAlist}}).
#' @param contrast either FALSE or a list of \code{MEA} objects to be used as a contrast
#' @param by.group logical. Should the different groups of \code{mea} be plotted separately?
#' @param sub.line on which margin line should the 'social presence' subtitle be printed, starting at 0 counting outwards.
#' @param ... further arguments and \code{\link[graphics]{par}} parameters passed to \code{\link[base]{plot}}
#' @param mea.alpha numeric from 0 to 1. The value of opacity of individual lines for the main MEA data. If set to zero, drawing is suppressed to improve performance.
#' @param contrast.alpha numeric from 0 to 1. The value of opacity of individual lines for contrast data. If set to zero, drawing is suppressed to improve performance.
#'
#' @details A typical application of \code{MEAlagplot} is to represent the difference between real dyads and random dyads obtained through a \code{\link{shuffle}} procedure.
#' It may also be used to see the difference among various filtering procedures or different regions of interest (e.g. head-synchrony versus body-synchrony, female vs. male dyads, etc).
#'
#' Percentages indicate the relative amount of synchrony where the values are higher than the contrast sample.
#' @examples
#' \donttest{
#' ## This example is excluded from test as it takes more than 10s to run
#' ## read the first 4 minutes of the normal sample
#' ##   (intake interviews of patients that carried on therapy)
#' path_normal <- system.file("extdata/normal", package = "rMEA")
#' mea_normal <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
#'                      s1Name = "Patient", s2Name = "Therapist",
#'                      idOrder = c("id","session"), idSep="_", skip=1, nrow = 6000)
#' mea_normal <- setGroup(mea_normal, "normal")
#'
#' ## read the dropout sample (intake interviews of patients that dropped out)
#' path_dropout <- system.file("extdata/dropout", package = "rMEA")
#' mea_dropout <- readMEA(path_dropout, sampRate = 25, s1Col = 1, s2Col = 2,
#'                      s1Name = "Patient", s2Name = "Therapist",
#'                      idOrder = c("id","session"), idSep="_", skip=1, nrow = 6000)
#' mea_dropout <- setGroup(mea_dropout, "dropout")
#'
#' ## Combine into a single object
#' mea_all = c(mea_normal, mea_dropout)
#'
#'## Create a shuffled sample
#'mea_rand = shuffle(mea_all, 20)
#'
#'## Compute ccf
#'mea_all = MEAccf(mea_all, lagSec = 5, winSec = 60, incSec = 30, r2Z = TRUE, ABS = TRUE)
#'mea_rand = MEAccf(mea_rand, lagSec = 5, winSec = 60, incSec = 30, r2Z = TRUE, ABS = TRUE)
#'
#'## Visualize the effects:
#'
#'MEAlagplot(mea_all, contrast = mea_rand, by.group = TRUE)
#'MEAlagplot(mea_all, contrast = mea_rand, by.group = FALSE, col=c(2,4))
#'}
#' @export
#'
MEAlagplot = function(mea, contrast=F, by.group=T, sub.line=0.5,
                      mea.lines = TRUE, mea.alpha = 0.8,
                      contrast.lines = TRUE, contrast.alpha=0.5,
                      ...){
  #1 social presence
  # mea = mea3
  # contrast = meaR
  # colz = c(1,2)

  mea = MEAlist(mea)
  if(!hasCCF(mea)) stop("Density plot requires MEA objects with valid CCF informations. Refer to function MEAccf()", call.=F)

  # if(!is.list(mea) && any(!sapply(mea,is.MEA))) stop("mea must be a list of MEA objects")

  if(is.logical(contrast) && contrast == F){
    contrast = F
  } else if(is.list(contrast) && all(sapply(contrast,is.MEA))){
    boot = MEAlist(contrast)
    if(!hasCCF(boot)) stop("'contrast' object had no valid CCF informations. Refer to function MEAccf()", call.=F)
    contrast = T
    if(length(attr(boot,"groups"))>1 ){
      boot = setGroup(boot, "Contrast")
      warning("Contrast data was collapsed to a single group.",call.=F)
    }
  } else stop("contrast must either be FALSE or a list of MEA object")

  if(by.group)
    mea = MEAlist(mea)
  else
    mea = setGroup(mea,"real")

  groups = attr(mea,"groups")
  lagSec = attr(mea,"ccf")$lag
  colz = mycolz(length(groups))
  glist = list()
  for(g in groups){
    glist[[g]] = lapply(mea,function(x){if(attr(x,"group")==g) x else NULL }) #set to NULL the meas of other groups
    glist[[g]] = Filter(Negate(is.null), glist[[g]])#remove null values
  }

  mean_lags = lapply(glist, function(gmea){ lapply(gmea,function(iFile){apply(iFile$ccf,2,mean,na.rm=T)})})
  if(contrast)bmean_lags = lapply(boot,function(iFile){apply(iFile$ccf,2,mean,na.rm=T)})

  # mean_win = lapply(mea,function(iFile){apply(iFile$ccf,1,mean,na.rm=T)})
  sampRate = attr(mea,"sampRate")
  ran = (-attr(mea , "ccf")$lag*sampRate): (attr(mea , "ccf")$lag*sampRate)

  #calculate overall means
  laggers = lapply(mean_lags, function(gmean_lags){ apply(matrix(unlist(gmean_lags),nrow = length(gmean_lags), byrow = T),2,mean)})
  if(contrast)blaggers = apply(matrix(unlist(bmean_lags),nrow = length(bmean_lags), byrow = T),2,mean)

  if(contrast) { myYlim = c(min(c(unlist(mean_lags,recursive = T),unlist(blaggers,recursive = T)),na.rm=T),
                            max(c(unlist(mean_lags,recursive = T),unlist(blaggers,recursive = T), 0.22),na.rm=T) )
  } else myYlim = c(min(unlist(mean_lags,recursive = T),na.rm=T),
                    max(c(unlist(mean_lags,recursive = T), 0.22),na.rm=T) )

  defPar = list(type="n",ylim=myYlim, xlim=c(min(ran), max(ran)), col=c(colz, "gray40" ), xaxt="n",bty='n',     #neu standard y-achse
                main="Cross correlation of real and contrast dyads", ylab=attr(mea,"ccf")$filter ,
                xlab=paste(attr(mea,"s2Name"),"leading    <<< ------   simultaneous  ------ >>>    ",attr(mea,"s1Name"),"leading\nLag (seconds)"))

  resPar = dotsList(defPar,...)
  for(k in c("col")){if(length(resPar[[k]])<3) resPar[[k]] = rep(resPar[[k]], length.out=length(groups)+1) }

  # print(resPar)
  do.call(base::plot, c(list(x = mean_lags[[1]][[1]]),resPar))

  graphics::axis(side = 1,at = (-lagSec:lagSec)*sampRate, labels = -lagSec:lagSec )

  #1. plot bootstrap values
  if(contrast && contrast.alpha != 0) for(i in 1:length(bmean_lags)){
    graphics::lines(ran,bmean_lags[[i]],col=colTrans(resPar$col[length(groups)+1],1.1,contrast.alpha ),lwd=0.3)  # original: gray70 fuer bootstrapped
    #       graphics::lines(ran,bmean_lags[[i]],col="white",lwd=0.2)   # fuer presence: white fuer boot, so dass nicht sichtbar
  }
  #2 plot values for each group
  if(mea.alpha != 0) for(g in 1:length(groups)){
    for(i in 1:length(mean_lags[[g]])){
      graphics::lines(ran,mean_lags[[g]][[i]],col=colTrans(resPar$col[g],1.1,mea.alpha),lwd=0.8)
      #       graphics::lines(ran,mean_lags[[i]],col="blue",lwd=0.2)
    }
  }

  #3. plot the bold average lines
  if(contrast) graphics::lines(ran,blaggers, col=resPar$col[length(groups)+1], lwd=6)
  for(g in 1:length(groups)){
    graphics::lines(ran,laggers[[g]],col=resPar$col[g], lwd=6)
  }

  leglab = if(contrast)c(groups,attr(boot,"groups")) else groups

  graphics::legend("topright", legend=leglab, col=resPar$col,lty=1,lwd = 2 ,bty='n')

  ###SOCIAL PRESENCE
  if(contrast){
    socPres = list()
    for(g in 1:length(groups)){
      socPres[[groups[g]]] = paste0(groups[g],": ", round(sum(laggers[[g]]>blaggers)/length(ran),2)*100,"%")
    }
    graphics::mtext(text = paste0("% above ",attr(boot,"groups"),": ", do.call(paste, c(socPres,list(sep=" | ")))) ,side=3,line = sub.line)

  }
}



#' Distribution of cross-correlations
#'
#' Plots the distribution of the average cross-correlations in a list of \code{MEA} objects.
#'
#' @param mea a well formatted list of \code{MEA} objects (see function \code{\link{MEAlist}}).
#' @param contrast either FALSE or a list of \code{MEA} objects to be used as a contrast
#' @param by Either "none", "group", "id", or "session". Defines the grouping to be used.
#' @param by.group deprecated argument. Use by="group" instead.
#' @param sub.line on which margin line should the Effect Size subtitle be printed, starting at 0 counting outwards.
#' @param ... further graphical parameters passed to  \code{\link[base]{plot}}
#'
#' @details If \code{contrast} is defined, then a normalized difference (Cohen's \emph{d}) between the means of each group and the contrast is provided.
#' Otherwise, if the \code{mea} object has 3 or less groups, Cohen's \emph{d} will be calculated on the group differences.
#' @examples
#' \donttest{
#' ## This example is excluded from test as it may take more than 10s to run
#' ## read the first 4 minutes of the normal sample
#' ##   (intake interviews of patients that carried on therapy)
#' path_normal <- system.file("extdata/normal", package = "rMEA")
#' mea_normal <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
#'                      s1Name = "Patient", s2Name = "Therapist",
#'                      idOrder = c("id","session"), idSep="_", skip=1, nrow = 6000)
#' mea_normal <- setGroup(mea_normal, "normal")
#'
#' ## read the dropout sample (intake interviews of patients that dropped out)
#' path_dropout <- system.file("extdata/dropout", package = "rMEA")
#' mea_dropout <- readMEA(path_dropout, sampRate = 25, s1Col = 1, s2Col = 2,
#'                      s1Name = "Patient", s2Name = "Therapist",
#'                      idOrder = c("id","session"), idSep="_", skip=1, nrow = 6000)
#' mea_dropout <- setGroup(mea_dropout, "dropout")
#'
#'## Combine into a single object
#'mea_all = c(mea_normal, mea_dropout)
#'
#'## Create a shuffled sample
#'mea_rand = shuffle(mea_all, 20)
#'
#'## Compute ccf
#'mea_all = MEAccf(mea_all, lagSec = 5, winSec = 60, incSec = 30, r2Z = TRUE, ABS = TRUE)
#'mea_rand = MEAccf(mea_rand, lagSec = 5, winSec = 60, incSec = 30, r2Z = TRUE, ABS = TRUE)
#'
#'## Visualize the effects:
#'
#'MEAdistplot(mea_all, contrast = mea_rand, by.group = TRUE)
#'}
#' @export
#'
MEAdistplot = function(mea, contrast=FALSE, by=c("none","group","id","session"), by.group=FALSE, sub.line = 0.5, ...) {
  #2 density
  # mea = d
  # contrast = r
  # colz = c(1,2)
  by = match.arg(by,c("none","group","id","session"))
  if(missing(by) && by.group){message("The 'by.group' argument is deprecated and will be removed from future versions. Please use by = 'group' instead, or refer to the help page."); by = "group"}
  mea = MEAlist(mea)
  if(!hasCCF(mea)) stop("Density plot requires MEA objects with valid CCF informations. Refer to function MEAccf()", call.=F)


  if(is.logical(contrast) && contrast == F){
    contrast = F
  } else if(is.list(contrast) && all(sapply(contrast,is.MEA))){
    boot = MEAlist(contrast)
    if(!hasCCF(boot)) stop("'contrast' object had no valid CCF informations. Refer to function MEAccf()", call.=F)
    if(length(boot) < 2) stop("'contrast' object must contain at least two MEA objects");
    contrast = T
    if(length(attr(boot,"groups"))>1 ){
      boot = setGroup(boot, "Contrast")
      warning("Contrast with multiple groups is not supported. Contrast data was collapsed to a single group.",call.=F)
    }
  } else stop("contrast must either be FALSE or a MEAlist object")


  if(by=="group"){
    mea = MEAlist(mea)
    groups = attr(mea,"groups")

  } else if(by == "id"){
    groups = unique(id(mea))

  } else if(by == "session") {
    groups = sort(unique(unlist(lapply(mea, attr, "session"))))

  } else {
    by = "group"
    mea = setGroup(mea,"Original")
    groups = attr(mea,"groups")
  }


  lagSec = attr(mea,"ccf")$lag
  colz = mycolz(length(groups))
  # colz = grey.colors(length(groups))
  glist = list()
  for(g in groups){
    glist[[g]] = lapply(mea,function(x){if(attr(x,by)==g) x else NULL }) #set to NULL the meas of other groups
    glist[[g]] = Filter(Negate(is.null), glist[[g]])#remove null values
  }

  #generate data
  grandAver = lapply(glist, function(iMea) sapply(iMea,function(x)x$ccfRes$grandAver ))
  dreal = lapply(grandAver, function(x){if(length(x)==1) list("x"=x, "y"=NA) else stats::density(x)})
  cohs = list()
  if(contrast){
    bgrandAver = sapply(boot,function(x)x$ccfRes$grandAver )
    dboot = stats::density(bgrandAver)
    for(g in 1:length(groups)){
      cohs[[groups[g]]] = paste0(groups[g]," d = ", round(cohens_d(grandAver[[g]],bgrandAver),2 ))
    }
  } else {
    dboot = list(x=NULL,y=NULL)
    if(length(groups)==2)
      cohs[[1]] = paste0("Cohen's d = ", round(cohens_d(grandAver[[1]],grandAver[[2]]),2 ))
    else if(length(groups)==3){
      shortGroups = substr(groups, 1, 4)
      cohs[[1]] = paste0(shortGroups[1],"-",shortGroups[2]," d = ", round(cohens_d(grandAver[[1]],grandAver[[2]]),2 ))
      cohs[[2]] = paste0(shortGroups[1],"-",shortGroups[2]," d = ", round(cohens_d(grandAver[[2]],grandAver[[3]]),2 ))
      cohs[[3]] = paste0(shortGroups[1],"-",shortGroups[2]," d = ", round(cohens_d(grandAver[[1]],grandAver[[3]]),2 ))
    }
  }
  myYlim = c(0,max(c(  dboot$y, unlist(lapply(dreal,function(x)x$y))  ),na.rm = T) )
  myXlim = c(min(c(  dboot$x, unlist(lapply(dreal,function(x)x$x))  )),max(c(  dboot$x, unlist(lapply(dreal,function(x)x$x))  )))

  defPar = list(
    x = dboot,type="n",ylim=myYlim, xlim = myXlim,col=grDevices::rgb(1,0,0,0.2),  bty='n', ylab="Density",
    main="Density of MEA average cross-correlations",xlab=paste(attr(mea,"ccf")$filter,"grand average")
  )
  resPar = dotsList(defPar,...)
  do.call(base::plot, resPar)
  graphics::lines(dboot, col ="gray40", lwd=4,lty=3)
  for(g in 1:length(groups)){
    if(length(dreal[[g]]$x) == 1)
      graphics::abline(v=dreal[[g]]$x,col=colz[g], lwd=2)
    else
      graphics::lines(dreal[[g]],col=colz[g], lwd=2)
  }

  leglab = if(contrast)c(groups,attr(boot,"groups")) else groups
  legcol = if(contrast)c(colz,"gray40") else colz
  legLTY = c(rep(1,length(groups)),3)
  legLWD = c(rep(2,length(groups)),4)
  graphics::legend("topright", legend=leglab, col=legcol,lty=legLTY,lwd = legLWD ,bty='n')

  if(length(cohs)>0){
    if(contrast) vsLab = paste0("ES vs ",attr(boot,"groups"),": ") else vsLab = ""
    vsLab = paste0(vsLab, do.call(paste, c(cohs,list(sep=" | "))))
    if(nchar(vsLab)>100){
      vsLab = ""
      # vsLab = gsub("(.{100})", "\\1\\\n", vsLab)
    }
    graphics::mtext(text = vsLab ,side=3,line = sub.line)
  }

}
#
# MEAdistplot = function(mea, contrast=F, by=c("none","group","id","session"), by.group=FALSE, ...) {
#   #2 density
#   # mea = d
#   # contrast = r
#   # colz = c(1,2)
#   by = match.arg(by,c("none","group","id","session"))
#   if(missing(by) && by.group){print("ASD"); by = "group"}
#   mea = MEAlist(mea)
#   if(!hasCCF(mea)) stop("Density plot requires MEA objects with valid CCF informations. Refer to function MEAccf()", call.=F)
#
#
#   if(is.logical(contrast) && contrast == F){
#     contrast = F
#   } else if(is.list(contrast) && all(sapply(contrast,is.MEA))){
#     boot = MEAlist(contrast)
#     if(!hasCCF(boot)) stop("'contrast' object had no valid CCF informations. Refer to function MEAccf()", call.=F)
#     if(length(boot) < 2) stop("'contrast' object must contain at least two MEA objects");
#     contrast = T
#     if(length(attr(boot,"groups"))>1 ){
#       boot = setGroup(boot, "Contrast")
#       warning("Contrast with multiple groups is not supported. Contrast data was collapsed to a single group.",call.=F)
#     }
#   } else stop("contrast must either be FALSE or a MEAlist object")
#
#
#   if(by=="group"){
#     mea = MEAlist(mea)
#     groups = attr(mea,"groups")
#
#   } else if(by == "id"){
#     groups = unique(id(mea))
#
#   } else if(by == "session") {
#     groups = sort(unique(unlist(lapply(mea, attr, "session"))))
#
#   } else {
#     by = "group"
#     mea = setGroup(mea,"Original")
#     groups = attr(mea,"groups")
#   }
#
#
#   lagSec = attr(mea,"ccf")$lag
#   colz = mycolz(length(groups))
#   # colz = grey.colors(length(groups))
#   glist = list()
#   for(g in groups){
#     glist[[g]] = lapply(mea,function(x){if(attr(x,by)==g) x else NULL }) #set to NULL the meas of other groups
#     glist[[g]] = Filter(Negate(is.null), glist[[g]])#remove null values
#   }
#
#   #generate data
#   grandAver = lapply(glist, function(iMea) sapply(iMea,function(x)x$ccfRes$grandAver ))
#   dreal = lapply(grandAver, function(x){if(length(x)==1) list("x"=x, "y"=NA) else stats::density(x)})
#   cohs = list()
#   if(contrast){
#     bgrandAver = sapply(boot,function(x)x$ccfRes$grandAver )
#     dboot = stats::density(bgrandAver)
#     for(g in 1:length(groups)){
#       cohs[[groups[g]]] = paste0(groups[g]," d = ", round(cohens_d(grandAver[[g]],bgrandAver),2 ))
#     }
#   } else {
#     dboot = list(x=NULL,y=NULL)
#     if(length(groups)==2)
#       cohs[[1]] = paste0("Cohen's d = ", round(cohens_d(grandAver[[1]],grandAver[[2]]),2 ))
#     else if(length(groups)==3){
#       shortGroups = substr(groups, 1, 4)
#       cohs[[1]] = paste0(shortGroups[1],"-",shortGroups[2]," d = ", round(cohens_d(grandAver[[1]],grandAver[[2]]),2 ))
#       cohs[[2]] = paste0(shortGroups[1],"-",shortGroups[2]," d = ", round(cohens_d(grandAver[[2]],grandAver[[3]]),2 ))
#       cohs[[3]] = paste0(shortGroups[1],"-",shortGroups[2]," d = ", round(cohens_d(grandAver[[1]],grandAver[[3]]),2 ))
#     }
#   }
#   myYlim = c(0,max(c(  dboot$y, unlist(lapply(dreal,function(x)x$y))  ),na.rm = T) )
#   myXlim = c(min(c(  dboot$x, unlist(lapply(dreal,function(x)x$x))  )),max(c(  dboot$x, unlist(lapply(dreal,function(x)x$x))  )))
#
#   defPar = list(
#     x = dboot,type="n",ylim=myYlim, xlim = myXlim,col=grDevices::rgb(1,0,0,0.2),  bty='n', ylab="Density",
#     main="Density of MEA average cross-correlations",xlab=paste(attr(mea,"ccf")$filter,"grand average")
#   )
#   resPar = dotsList(defPar,...)
#   do.call(base::plot, resPar)
#   graphics::lines(dboot, col ="gray40", lwd=4,lty=3)
#   for(g in 1:length(groups)){
#     if(length(dreal[[g]]$x) == 1)
#       graphics::abline(v=dreal[[g]]$x,col=colz[g], lwd=2)
#     else
#       graphics::lines(dreal[[g]],col=colz[g], lwd=2)
#   }
#
#   leglab = if(contrast)c(groups,attr(boot,"groups")) else groups
#   legcol = if(contrast)c(colz,"gray40") else colz
#   legLTY = c(rep(1,length(groups)),3)
#   legLWD = c(rep(2,length(groups)),4)
#   graphics::legend("topright", legend=leglab, col=legcol,lty=legLTY,lwd = legLWD ,bty='n')
#
#   if(length(cohs)>0){
#     if(contrast) vsLab = paste0("ES vs ",attr(boot,"groups"),": ") else vsLab = ""
#     vsLab = paste0(vsLab, do.call(paste, c(cohs,list(sep=" | "))))
#     if(nchar(vsLab)>100){
#       vsLab = ""
#       # vsLab = gsub("(.{100})", "\\1\\\n", vsLab)
#     }
#     graphics::mtext(text = vsLab ,3,line = -5.5)
#   }
#
# }
#



#' Plots the initial, middle and ending part of a MEA object
#'
#' This is typically useful to check if the motion energy time-series are good.
#' The middle section is chosen randomly among possible middle sections.
#'
#' @param mea   an object of class \code{MEA} (see function \code{\link{readMEA}}).
#' @param width integer. The number of seconds to be plotted for each panel
#' @param ...   further arguments passed to \code{plot}
#'
#' @details Motion energy time-series should always be visually inspected for possible artifacts. Periodic peaks or drops in time-series are indicators of e.g. key-frames or duplicated video-frames.
#' For further information regarding the program MEA, please refer to the documentation available at \code{http://www.psync.ch}.
#' @examples
#' ## read a single file
#' path_normal <- system.file("extdata/normal/200_01.txt", package = "rMEA")
#' mea_normal <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
#'                      s1Name = "Patient", s2Name = "Therapist", skip=1,
#'                      idOrder = c("id","session"), idSep="_")
#' ## Visual inspection of the data
#' diagnosticPlot(mea_normal[[1]])
#'
#' @export
#'
diagnosticPlot = function(mea,width=60,...){
  #### debug
  # mea = meaRaw[[1]]
  # width = 60
  #######
  if(!is.MEA(mea)) stop("This function only accepts individual MEA objects")
  sampRate = attr(mea,"sampRate")
  width = width*sampRate

  #find last nonzero value
  zeroes = apply(mea$MEA,1,sum)
  i=length(zeroes)
  k=zeroes[i]
  while(k == 0){
    i = i-1
    k = zeroes[i]
  }
  end = min(i+10*sampRate-1, length(zeroes)) #aggiungi max 10 secondi di zero alla fine
  if(end<width*3+10)
    width = trunc(end/3)

  # get intervals
  seg1 = c(1,(width)) #first  min
  seg3 = c((end-(width)),end) #last  min
  offset = stats::rnorm(1, 0, (end-2*width)/5 )
  b_center = width+mean(c(width,end-width))   #random value toward the center
  if(offset<width*1.5 || offset>end-width*1.5) offset = 0
  b_center = trunc(b_center+offset)
  seg2 = c((b_center - trunc(width/2) ),( b_center + trunc(width/2) ))
  all = c(seg1[1]:seg1[2],rep(NA,20),seg2[1]:seg2[2],rep(NA,20),seg3[1]:seg3[2])

  s1 = mea$MEA[all,1]
  s2 = mea$MEA[all,2]

  ##rescale values
  s1 = rangeRescale(s1,0,10)
  s2 = rangeRescale(s2,0,10)



  #https://coolors.co/ef3e36-17bebb-2e282a-edb88b-fad8d6
  myYlim= c(-0.2,max(s1, s2, na.rm = T))
  colz = c("#2E282A","#EF3E36")
  origPar = graphics::par(c("xpd","mar","font"))
  graphics::par(xpd=TRUE,mar=c(3,0.5,2.5,0.5))
  base::plot(s1,type = "l",col=colz[1] ,axes = F, ann=FALSE, lwd = 2,ylim = myYlim,...)
  graphics::lines(s2,col=colz[2],lwd = 2)
  graphics::segments(width+10, myYlim[1]-0.1, width +10, myYlim[2])
  graphics::segments(width*2+32, myYlim[1]-0.1, width*2+32, myYlim[2])
  graphics::segments(0,-0.1,length(all),-0.1)
  graphics::text(x = c(width/2,width*1.5,width*2.5)+20*(0:2)+1,-0.35, labels = c(
    paste(timeMaster(round(0),out="min"), "to", timeMaster(round(seg1[2]/sampRate),out="min")),
    paste(timeMaster(round(seg2[1]/sampRate),out="min"), "to", timeMaster(round(seg2[2]/sampRate),out="min")),
    paste(timeMaster(round(seg3[1]/sampRate),out="min"), "to", timeMaster(round(seg3[2]/sampRate),out="min"))
    ),cex=1 )
  graphics::title(paste0("Group: ",attr(mea,"group"),", Id: ",attr(mea,"id"),", Session: ",attr(mea,"session") ))
  graphics::par(font=2)
  graphics::legend(length(all)/2,-0.3,xjust = 0.5,legend=c(attr(mea,"s1Name"),attr(mea,"s2Name")),lwd=2,col = colz, bty = "n",horiz=T,cex =1.3)
  graphics::par(origPar)
}



#' Plots an object of class \code{MEA}
#'
#'
#'
#' @param x an object of class \code{MEA} (see function \code{\link{readMEA}}).
#' @param from either an integer or a string in the format hh:mm:ss or mm:ss representing the starting second.
#' @param to if \code{duration} is not specified, either an integer or a string in the format hh:mm:ss or mm:ss representing the ending second.
#' @param duration if \code{to} is not specified, either an integer or a string in the format hh:mm:ss or mm:ss representing the amount of seconds to be plotted.
#' @param ccf either FALSE or a string representing the type of ccf to be overlayed.
#' One of "all_lags"  "s1_lead"   "s2_lead"   "lag_zero"  "s1_lead_0" "s2_lead_0" "bestLag" "grandAver" "winTimes".
#' @param rescale logical. Should the motion energy time-series be rescaled?
 #' @param ... further arguments passed to \code{\link[base]{plot}}
#'
#' @details Note: if more of than 10s of trailing zeroes are found at the end of both s1 and s2 signals they are truncated.
#' @examples
#' ## read a single file
#' path_normal <- system.file("extdata/normal/200_01.txt", package = "rMEA")
#' mea_normal <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
#'                      s1Name = "Patient", s2Name = "Therapist", skip=1,
#'                      idOrder = c("id","session"), idSep="_")
#' mea_normal <- MEAccf(mea_normal, lagSec = 5, winSec = 30, incSec = 10, ABS = FALSE)
#' ## Visual inspection of the data
#' plot(mea_normal[[1]], from = 60, to = "2:00")
#' plot(mea_normal[[1]], from = 0, duration = "5:00")
#'
#' #' ## Visualize CCF inspection of the data
#' plot(mea_normal[[1]], from = 0, duration = "2:00", ccf = "lag_zero", rescale=TRUE)
#'
#' @export
#'
plot.MEA = function(x, from = 0, to = NULL, duration = NULL, ccf = F, rescale = F, ... ){
  #### debug
  # mea = mea3$all_38508_8
  # width = 60
  #######
  mea = x
  sampRate = attr(mea,"sampRate")
  from = timeMaster(from, out="sec") * sampRate +1
  if(!missing(to)){
    if(!missing(duration)) stop("'duration' and 'to' cannot be both specified.",call.=F)
    to = timeMaster(to, out="sec") * sampRate
  } else if(!missing(duration)) {
    duration = timeMaster(duration, out="sec") * sampRate
    to = min(from + duration, nrow(mea$MEA))
  } else {
    to = nrow(mea$MEA)
  }

  if(to<from) stop("'to' cannot be smaller than 'from'",call.=F)
  #find last nonzero value
  if(to==nrow(mea$MEA)){
    zeroes = apply(mea$MEA,1,sum)
    i=length(zeroes)
    k=zeroes[i]
    while(k == 0){
      i = i-1
      k = zeroes[i]
    }
    end = i + min(10*sampRate-1, 10*sampRate-1) #aggiungi max 10 secondi di zero alla fine
  } else end = to

  # get intervals
  all = from:to
  s1 = mea$MEA[all,1]
  s2 = mea$MEA[all,2]
  if(rescale){
    s1 = rangeRescale(s1,0,1)
    s2 = rangeRescale(s2,0,1)
  }
  myYlim= c(min(s1,s2,na.rm = T),max(s1,s2,na.rm = T))
  myYlim[1] = myYlim[1]-myYlim[2]*0.05
  # colz = c("#2E282A","#EF3E36","red")
  colz = c(mycolz(2), "red")
  origPar = graphics::par(c("xpd","mar","font"))

  #add ccf
  if(!is.logical(ccf)){
    if(is.null(x$ccf)) stop("The 'mea' object does not have a CCF. Please run MEAccf()" )
    if(!ccf %in% names(mea$ccfRes)) stop("'ccf value not recognized. It should be one of ",paste(names(mea$ccfRes),collapse = ", " ),call.=F)
    ccf_lab = ccf
    ccf = list(data.frame(x$ccfRes[[ccf]]))
    if(any(stats::na.omit(ccf[[1]]) < 0)) ABS = F else ABS = T
    ccf = unlist(winInter(ccf, winSec = attr(x,"ccf")$win, incSec =attr(x,"ccf")$inc, sampRate = attr(x,"sampRate")))
    ccf = ccf[all]

    if(!ABS){ ccf = c(ccf,-1,1)
    } else {ccf=c(ccf,0,1)}
    ccf = rangeRescale(ccf,0,myYlim[2])
    ccf = ccf[1:(length(ccf)-2)]
  } else ccf = NULL

  #draw plot
  graphics::par(xpd=TRUE, mar=c(3,3.5,2.5,0.5))

  defPar = list(
    x = s1, type = "l",col=colz ,axes = F, ann=FALSE, ylim = myYlim, lty= c(1,1,3), lwd=c(2,2,2),
    main = paste0("Group: ",attr(mea,"group"),", Id: ",attr(mea,"id"),", Session: ",attr(mea,"session") )
  )
  resPar = dotsList(defPar,...)
  # print(resPar)
  #recycle parameters
  for(k in c("lty","lwd","col")){if(length(resPar[[k]])<3) resPar[[k]] = rep(resPar[[k]], length.out=3) }
  do.call(base::plot, resPar)
  graphics::lines(s1, col=resPar$col[1],lwd = resPar$lwd[1], lty=resPar$lty[1])
  graphics::lines(s2, col=resPar$col[2],lwd = resPar$lwd[2], lty=resPar$lty[2])
  if(!is.null(ccf)){
    graphics::lines(ccf,col=resPar$col[3],lwd = resPar$lwd[3] ,lty=resPar$lty[3])
  }
  if(!is.null(ccf)){
    if(ABS){
      graphics::abline(h=0, lwd=1,lty=3)
      graphics::text( 0,-0.3, "ccf = 0")
    }    else {graphics::abline(h=myYlim[2]/2, lwd=1,lty=3)
      graphics::text(0, myYlim[2]/2, "ccf = 0")}
  }

  graphics::axis(2)
  graphics::title(main = resPar$main)
  graphics::text(x = length(all)/2, myYlim[1]-myYlim[2]*0.01 ,
       labels =  paste(timeMaster(round((from-1)/sampRate),out="min"), "to", timeMaster(round(to/sampRate),out="min"))
       ,cex=1.1 )
  graphics::par(font=2)
  legend_labs = c(attr(mea,"s1Name"),attr(mea,"s2Name"))
  if(!is.null(ccf)) legend_labs = c(legend_labs, paste0("CCF",ccf_lab))
  graphics::legend(length(all)/2,myYlim[1]-myYlim[2]*0.05 ,xjust = 0.5,legend=legend_labs,lwd=resPar$lwd,col = resPar$col, bty = "n",horiz=T,cex =1.3,lty=resPar$lty)
  graphics::par(origPar)
}


#' Adds lines of a \code{MEA} object to a Plot
#'
#'
#'
#' @param x an object of class \code{MEA} (see function \code{\link{readMEA}}).
#' @param from either an integer or a string in the format hh:mm:ss or mm:ss representing the starting second.
#' @param to if \code{duration} is not specified, either an integer or a string in the format hh:mm:ss or mm:ss representing the ending second.
#' @param duration if \code{to} is not specified, either an integer or a string in the format hh:mm:ss or mm:ss representing the amount of seconds to be plotted.
#' @param ccf either FALSE or a string representing the type of ccf to be overlayed. Possible values can be found with the \code{\link{ccfResNames}} function.
#' @param rescale logical. Should the motion energy time-series be rescaled?
#' @param ... further arguments passed to \code{\link[graphics]{lines}}
#'
#' @details Note: if more of than 10s of trailing zeroes are found at the end of both s1 and s2 signals they are truncated.
#' @examples
#' ## read a single file
#' path_normal <- system.file("extdata/normal/200_01.txt", package = "rMEA")
#' mea_normal <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
#'                      s1Name = "Patient", s2Name = "Therapist", skip=1,
#'                      idOrder = c("id","session"), idSep="_")
#' mea_normal <- MEAccf(mea_normal, lagSec = 5, winSec = 30, incSec = 10, ABS = FALSE)
#' mea_smoothed <- MEAsmooth(mea_normal)
#' ## Visual inspection of the data
#' plot(mea_normal[[1]], from = 240, duration=20)
#' lines(mea_smoothed[[1]], from = 240, duration=20, lty=3, col=c(1,2))
#' @export
#'
lines.MEA = function(x, from=0, to = NULL, duration=NULL, ccf=F, rescale =F,... ){
  #### debug
  # mea = mea3$all_38508_8
  # width = 60
  #######
  mea = x
  sampRate = attr(mea,"sampRate")
  from = timeMaster(from, out="sec") * sampRate +1
  if(!missing(to)){
    if(!missing(duration)) stop("'duration' and 'to' cannot be both specified.",call.=F)
    to = timeMaster(to, out="sec") * sampRate
  } else if(!missing(duration)) {
    duration = timeMaster(duration, out="sec") * sampRate
    to = min(from + duration, nrow(mea$MEA))
  } else {
    to = nrow(mea$MEA)
  }

  if(to<from) stop("'to' cannot be smaller than 'from'",call.=F)
  #find last nonzero value
  if(to==nrow(mea$MEA)){
    zeroes = apply(mea$MEA,1,sum)
    i=length(zeroes)
    k=zeroes[i]
    while(k == 0){
      i = i-1
      k = zeroes[i]
    }
    end = i + min(10*sampRate-1, 10*sampRate-1) #aggiungi max 10 secondi di zero alla fine
  } else end = to

  # get intervals
  all = from:to
  s1 = mea$MEA[all,1]
  s2 = mea$MEA[all,2]
  if(rescale){
    s1 = rangeRescale(s1,0,1)
    s2 = rangeRescale(s2,0,1)
  }
  myYlim= c(min(s1,s2,na.rm = T),max(s1,s2,na.rm = T))
  myYlim[1] = myYlim[1]-myYlim[2]*0.05

  colz = c(mycolz(2), "red")

  #add ccf
  if(!is.logical(ccf)){
    if(is.null(x$ccf)) stop("The 'mea' object does not have a CCF. Please run MEAccf()" )
    if(!ccf %in% names(mea$ccfRes)) stop("'ccf value not recognized. It should be one of ",paste(names(mea$ccfRes),collapse = ", " ),call.=F)
    ccf_lab = ccf
    ccf = list(data.frame(x$ccfRes[[ccf]]))
    if(any(stats::na.omit(ccf[[1]]) < 0)) ABS = F else ABS = T
    ccf = unlist(winInter(ccf, winSec = attr(x,"ccf")$win, incSec =attr(x,"ccf")$inc, sampRate = attr(x,"sampRate")))
    ccf = ccf[all]

    if(!ABS){ ccf = c(ccf,-1,1)
    } else {ccf=c(ccf,0,1)}
    ccf = rangeRescale(ccf,0,myYlim[2])
    ccf = ccf[1:(length(ccf)-2)]
  } else ccf = NULL

  defPar = list(x = s1, col = colz, lty = c(1, 1, 3), lwd = c(2, 2, 2))

  resPar = dotsList(defPar,...)
  # print(resPar)
  #recycle parameters
  for(k in c("lty","lwd","col")){if(length(resPar[[k]])<3) resPar[[k]] = rep(resPar[[k]], length.out=3) }
  graphics::lines(s1, col=resPar$col[1],lwd = resPar$lwd[1], lty=resPar$lty[1])
  graphics::lines(s2, col=resPar$col[2],lwd = resPar$lwd[2], lty=resPar$lty[2])
  if(!is.null(ccf)){
    graphics::lines(ccf,col=resPar$col[3],lwd = resPar$lwd[3] ,lty=resPar$lty[3])
  }
  if(!is.null(ccf)){
    if(ABS){
      graphics::abline(h=0, lwd=1,lty=3)
      graphics::text( 0,-0.3, "ccf = 0")
    }    else {graphics::abline(h=myYlim[2]/2, lwd=1,lty=3)
      graphics::text(0, myYlim[2]/2-0.3, "ccf = 0")}
  }
}

# #' Plots all MEA signals contained in an object of class \code{MEAlist}
# #'
# #'
# #' @param x an object of class \code{MEAlist}
# #' @param from integer. The first sample to be plotted. Defaults to the beginning of the signals.
# #' @param to integer. The last sample to be plotted. Defaults to the length of the longest MEA signal
# #' @param ccf either FALSE or a string representing the type of ccf to be overlayed.
# #' One of "all_lags", "s1_lead", "s2_lead", "lag_zero", "s1_lead_0", "s2_lead_0".
# #' @param rescale logical. Should the motion energy time-series be rescaled?
# #' @param ... further graphical parameters passed to  \code{\link[graphics]{plot}}
# #'
# #' @details Note: if more of than 10s of trailing zeroes are found at the end of both s1 and s2 signals they are truncated.
# #'
# #' @export

# plot.MEAlist = function(x,from=1,to=NULL, duration=NULL,ccf = F ,rescale =F,...){
#   mea = x
#   v="y"
#   if(length(mea)>10) {
#     v <- readline(paste("Warning: you are trying to create",length(mea),"plots. Are you sure you want to continue?\r\ny/n: "))
#   }
#   if(tolower(substr(as.character(v),1,1)) == "y"){
#     cat("\r\nDrawing plots:\r\n")
#     i=1
#     for(x in mea){
#       prog(i,length(mea)); i= i+1
#       if(to>nrow(x$MEA)) to = nrow(x$MEA)
#       base::plot(x,from=from,to=to, ccf=ccf, rescale = rescale,...)
#     }
#   } else cat("\r\nOperation aborted.")
# }

mycolz = function(n,demo=F,alpha=1){
  fmod= function(k,m){
    j = floor(k/m)
    a = k-m*j
    return(a)
  }
  colz =numeric(n)
  for (i in 1:n) {
    colz[i] = grDevices::hsv(fmod(i * 0.618033988749895, 1.0),
                  0.8, sqrt(1.0 - fmod(i * 0.618033988749895, 0.5)),alpha = alpha)

  }
  if(demo)graphics::pie(rep(1,n), col=colz)
  return(colz)
}



#' Plot a heatmap of dyadic cross-correlations
#'
#' Graphical representation of the lagged cross-correlations in time. Provides an intuitive description of synchronization dynamics.
#'
#' @param mea an object of class \code{MEA} (see function \code{\link{readMEA}}).
#' @param legendSteps integer. the number of levels used for the color-coding of the legend.
#' @param rescale logical. If TRUE, the color range will represent the minimum and maximum of the data. Otherwise the theoretical correlation range -1 to 1.
#' @param colors a vector of colors defining the plot scale.
#' @param bias a positive number. Allows to skew the color scale. Values larger than 1 give more widely spaced colors at the high end, and vice versa.
#' @param mirror logical. If TRUE, colors are mirrored for negative correlation values. This has effect only if \code{\link{MEAccf}} was run with \code{ABS=FALSE}
#'
#' @details The cross-correlation values are rescaled to be in a range from 0 to 1 before plotting.
#' @export

MEAheatmap  = function(mea, legendSteps = 10, rescale = FALSE, colors = c("#F5FBFF","#86E89E","#FFF83F","#E8A022","#FF3700"), bias = 1, mirror=TRUE){
  if(!is.MEA(mea) || is.null(mea$ccf)) stop("Only MEA objects with ccf analysis can be plotted by this function.",call.=F)
  ABS = !any(stats::na.omit(mea$ccf)<0) #do we have negative numbers?
  mat = mea$ccf
  if(grep("z",attributes(mea)$ccf$filter)) mat = tanh(mat) #revert fisher's z transform to have -1:1 range
  if(rescale){
    if(ABS) mat = rangeRescale(mat, 0,1)
    else    mat = rangeRescale(mat,-1,1)
  }

  sampRate = attr(mea,"sampRate")

  if(ABS){
    bins = seq(1/legendSteps,1, by=1/legendSteps)
    colfunc <- grDevices::colorRampPalette(colors=colors, bias=bias)
  }  else {
    bins = seq(-1+2/legendSteps,1, by=2/legendSteps)
    if(mirror) colfunc <- grDevices::colorRampPalette(colors=c(rev(colors[2:length(colors)]), colors), bias=bias)
    else colfunc <- grDevices::colorRampPalette(colors=colors, bias=bias)
  }

  colz = colfunc(legendSteps)
  colz = c(colz,"#aaaaaa") #colore degli NA

  plotH = ncol(mat)
  plotW = nrow(mat)
  leg_area = 1/5*plotW #legend area
  myYlim = c(1,plotH+1)
  myXlim = c(1,plotW+leg_area)
  oldPar = graphics::par("mar")
  graphics::par(mar=c(5, 6, 4, 2) + 0.1)
  base::plot(0,type="n",xlim=myXlim, ylim=myYlim,bty="n",axes=F,xlab="ccf windows",
       ylab=paste(attr(mea,"s2Name"),"leading    <<< ---   simultaneous  --- >>>    ",attr(mea,"s1Name"),"leading\nseconds"),
       main = paste0("Group: ",attr(mea,"group"),", id: ",attr(mea,"id"),", session: ",attr(mea,"session") ))
  mtext_lab = paste0(attr(mea,"ccf")$filter," with ", attr(mea,"ccf")$win,"s windows and ",attr(mea,"ccf")$inc,"s increments." )
  graphics::mtext(mtext_lab,side = 3)
  #axis
  xleft = rep(1:plotW, each =plotH)
  ybottom = rep(1:plotH,plotW)
  graphics::axis(1,at = 0.5+(1:attr(mea,"ccf")$n_win),labels=(1:attr(mea,"ccf")$n_win) )
  graphics::axis(2,at = 1.5+(-attr(mea,"ccf")$lag :attr(mea,"ccf")$lag )*sampRate + attr(mea,"ccf")$lag*sampRate,
       labels=-attr(mea,"ccf")$lag:attr(mea,"ccf")$lag )


  #to debug
  # mat[,1] = 1
  # mat[20:40,50:70]=NA

  vals = c(t(mat))

  #debug
   # vals[ncol(mat)*(0:nrow(mat))+1] = 1

  iCol = sapply(vals,function(x)sum(x>bins)+1)
  iFill = iCol
  iCol[is.na(iCol)]=length(colz)
  iFill[!is.na(iFill)] = NA
  iFill[is.na(iFill)] = 40

  graphics::rect(xleft,ybottom,xleft+1,ybottom+1,col = colz[iCol],border = NA)

  getNA = which(is.na(mat),arr.ind = T)
  if(length(getNA)>0)
    graphics::rect(getNA[,1],getNA[,2],getNA[,1]+1,getNA[,2]+1,col = 0,border = NA,density = 20)


  #legend

  na_y1 = 1
  na_y2 = plotH/100 * 10
  leg_x1 = plotW+1.5/5*leg_area
  leg_x2 = plotW+3.5/5*leg_area
  leg_y1 = plotH/100*15
  leg_y2 = plotH
  leg_incr = (plotH/100*85 )/legendSteps

  label_0 = plotH/100*20
  label_1 = plotH/100*97
  label_00 =label_0+(label_1-label_0)/2



  graphics::rect(leg_x1, na_y1,leg_x2, na_y2
       ,col = "#aaaaaa",border = NA)
  graphics::rect(leg_x1, na_y1,leg_x2, na_y2
       ,col = 0,border = NA,density = 20)



  legend_border = ifelse(legendSteps>50,NA,0)
  graphics::rect(
    leg_x1,
    seq(leg_y1,leg_y2-leg_incr, by=leg_incr),
    leg_x2,
    seq(leg_y1+leg_incr,leg_y2, by=leg_incr),
    col = colz, border = legend_border
  )

  if(rescale)
    graphics::text(leg_x2, c(na_y2/2,label_0,label_1), labels = c("NA",min(mea$ccf,na.rm=T),max(mea$ccf,na.rm=T)),cex = 1.5,pos=4)
  else
    graphics::text(leg_x2, c(na_y2/2,label_0,label_1), labels = c("NA",ifelse(ABS,0,-1),1),cex = 1.5,pos=4)

}

