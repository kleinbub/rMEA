#' Moving-windows lagged cross-correlation routine for \code{MEA} and \code{MEAlist} objects
#'
#' This function analyzes a bivariate MEA signal represented by two time-series (subject 1 "s1", subject 2 "s2") resulting from a dyadic interaction.
#' MEAccf performs windowed cross-correlations with specified increments. The cross-correlation analysis is repeated for each
#' lag step, with discrete increments of 1 second in both directions.
#'
#' @param mea an object of class \code{MEA} or \code{MEAlist} (see function \code{\link{readMEA}})
#' @param lagSec an integer specifying the maximum number of lags (in seconds) for which the time-series will be shifted forwards and backwards.
#' @param winSec an integer specifying the cross-correlation window size (in seconds).
#' @param incSec an integer specifying the step size (in seconds) between successive windows. Values lower than \code{winSec} result in overlapping windows.
#' @param r2Z logical. The default value TRUE applies Fisher's r to Z transformation to all computed correlations.
#' @param ABS logical. The default value TRUE transforms the (Fisher's Z-transformed) correlations to absolute values.
#'
#' @details The choice of \code{lagSec} depends on the type of synchronization expected from the specific interaction. In the literature, lags of ±5 seconds have been reported by multiple authors.
#' Function \code{\link{MEAlagplot}} can be used for visual inspection of the appropriateness of the chosen lag.
#'
#'  The choice of \code{winSec} represents the temporal resolution of the analysis.
#'  The combination of \code{incSec} and \code{winSec} settings has a big impact on the results. These parameters should be chosen carefully, guided by theoretical and empirical considerations.
#'
#'  If \code{r2Z} is TRUE, values of Fisher's Z are constrained to an upper bound of 10.
#'
#' Using absolute values (\code{ABS}) treats positive and negative cross-correlations as equal. The underlying assumption is that both simultaneous movement (positive correlation) and when
#' one subject accelerates and the other decelerates (negative correlation), are both signs of interrelatedness and should thus contribute equally to overall synchrony.
#' @return The function returns a copy of the \code{mea} object which includes a cross-correlation table
#' @examples ## read a single file
#' path_normal <- system.file("extdata/normal/200_01.txt", package = "rMEA")
#' mea_normal <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
#'                      s1Name = "Patient", s2Name = "Therapist", skip=1,
#'                      idOrder = c("id","session"), idSep="_")
#'
#' ## perform ccf analysis
#' mea_ccf = MEAccf(mea_normal, lagSec = 5, winSec = 60, incSec = 30, r2Z = TRUE, ABS = TRUE)
#' summary(mea_ccf)
#'
#' #visualize the analysis results for the first file
#' MEAheatmap(mea_ccf[[1]])
#'

#' @export
MEAccf = function(mea, lagSec, winSec, incSec, r2Z=T, ABS=T){
  UseMethod("MEAccf",mea)
}

#' @export
MEAccf.MEAlist = function(mea, lagSec, winSec, incSec, r2Z=T, ABS=T) {
  if(any(!sapply(mea, is.MEA)))
    stop("This function accepts only MEA objects (individual or a list of them). Please use readMEA() to import files")
  cat("\r\nComputing CCF:\r\n")
  res = Map(function(k,i) {
    prog(i,length(mea))
    MEAccf(k, lagSec, winSec, incSec, r2Z=r2Z, ABS=ABS)
    },mea, seq_along(mea) )
  res = MEAlist(res)
  return(res)
}

#' @export
MEAccf.MEA = function(mea, lagSec, winSec, incSec, r2Z=T, ABS=T){
  ####debug
  # mea = mearaw[[1]]
  # lagSec = 5
  # winSec = 15
  # incSec = 15
  # ################

  #import C correlation function
  C_cor=get("C_cor", asNamespace("stats"))

  if(!is.MEA(mea)) stop("Only objects of class MEA can be processed by this function")
  sampRate = attr(mea, "sampRate")

  #cat("\r\nMy Setting: lagSec = ",lagSec,"; winSec:",winSec,"; incSec:",incSec,"; sampRate:",sampRate)
  #cat("\r\ndyad CCF function -- v.1.2\r\ninput:",length(fileList),"dyads.\r\nFrequency:",sampRate)
  #   xname = signal$dyadNames[1];yname = signal$dyadNames[2];
  #   cat(paste0("\r\nHigh Sync at positive lags implies that the ",yname, " follows the ",xname,"\r\n"))
  win = winSec*sampRate
  lagSamp = lagSec*sampRate
  #if(lagUnit == "sample")
  ran = ((-lagSamp)): ((lagSamp))
  #else
  #  ran = seq(-lagSamp,lagSamp, sampRate)
  inc = incSec * sampRate
  #if(winSec!=incSec) warning("With overlapping windows, the sync series has to be shifted forward by half window in order to be aligned to the original signal!\r\n")
  iFile = mea$MEA
  n_win = ceiling((nrow(iFile)-win-lagSamp+1)/inc) #calcola il numero di finestre per ciascun file
  lcc=lapply(seq_len(n_win)-1, function(iWin)
  { #-----per ciascuna finestra--------
    ab = (iWin*inc +1):(iWin*inc +win+lagSamp) #calcola il range di sample di ciascuna finestra
    xWin = iFile[ab,1];yWin = iFile[ab,2] #estrai i dati della finestra in un vettore per sogg x e y
    if(max(sum(is.na(xWin)),sum(is.na(yWin))) > length(xWin)/2){ #se ci sono troppi NA, restituisci NA per tutti i lag
      res = rep(NA, length(ran))
    } else {
      res = vector("numeric",length(ran))
      i=0
      for(iLag in ran) { #-------per ciascun lag---------
        i= i+1
        LAG = abs(iLag)
        xRange = 1:(win);
        yRange = (LAG+1) :(win+LAG) #applica il lag spostando indietro il secondo soggetto.
        if(iLag<0){k=xRange;xRange=yRange;yRange=k;} #valori alti a lag positivi implicano che il sogg 2 segue sogg 1
        x = xWin[xRange]
        y = yWin[yRange]
        if(sum(x!=y)<2 ) cor_res = NA #controlla che ci siano almeno 3 punti !=0
        else cor_res = suppressWarnings(.Call(C_cor, x, y, 2L, FALSE))
        # if(!is.na(cor_res) && cor_res==1) cor_res = 0.999
        #if(!is.na(cor_res) && cor_res == 1) warning("In dyad: ", attr(mea,"uid"), ", correlation was 1 in window: ",iWin," and lag: ",iLag,".\r\nX: ",paste(x,collapse = " "),"\r\nY: ",paste(y,collapse = " "),"\r\nPlease check your raw data and report this diagnostic message to developers." ,call.=F)
        res[i] = cor_res
        #################################
      }
      res
    }#fine else
  }) #fine lapply finestre
  ccfmat = data.frame(matrix(unlist(lcc),ncol=length(ran), byrow = T, dimnames=list(paste0("w",seq_len(n_win)),paste0("lag",ran))))
  colnames(ccfmat) = paste0("lag",ran/sampRate)
  if(r2Z) ccfmat = fisher.r2z(ccfmat)
  if(ABS) ccfmat = abs(ccfmat)

  # analytics
  ccfRes = list(
    "all_lags" = apply(ccfmat, 1, mean ,na.rm=T),
    "s1_lead"  = apply(ccfmat[, (lagSec*sampRate+2):(lagSec*sampRate*2+1)], 1, mean ,na.rm=T),
    "s2_lead"  = apply(ccfmat[, 1:lagSec*sampRate], 1, mean ,na.rm=T),
    "lag_zero"  = ccfmat[, lagSec*sampRate+1],
    "s1_lead_0" = apply(ccfmat[, (lagSec*sampRate+1):(lagSec*sampRate*2+1)], 1, mean ,na.rm=T),
    "s2_lead_0" = apply(ccfmat[, 1:(lagSec*sampRate +1)], 1, mean ,na.rm=T),
    "grandAver" = mean(unlist(ccfmat),na.rm=T)
  )
  names(ccfRes$zero) = names(ccfRes$pace)

  filter = "CCF"
  if(r2Z) filter = paste0("z",filter)
  if(ABS) filter = paste0("|",filter,"|")
  attributes(mea)$ccf[["filter"]] = filter
  attributes(mea)$ccf = list(
                          filter = filter,
                          lag   = lagSec,
                          win   = winSec,
                          inc   = incSec,
                          n_win = n_win)
  mea$ccf = ccfmat
  mea$ccfRes = ccfRes

  #debug
  # cat("\r\n",mea$ccf[1:10,10])
  return(mea)
}

fisher.r2z <- function(r) {
  r = 0.5 * (log(1+r) - log(1-r))
  r[r >  10] =  10
  r[r < -10] = -10
  r
}
