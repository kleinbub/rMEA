
cohens_d <- function(x, y) {
  lx <- length(x)- 1
  ly <- length(y)- 1
  md  <- abs(mean(x) - mean(y))        ## mean difference (numerator)
  csd <- lx * stats::var(x) + ly * stats::var(y)
  csd <- csd/(lx + ly)
  csd <- sqrt(csd)                     ## common sd computation

  cd  <- md/csd                        ## cohen's d
}



## binds unequal columns to a same data.frame padding NAs to the end of the shorter
unequalCbind = function(...) {
  dots <- list(...)
  #debug
  #dots = list(my.orig,my.ccf)
  #
  dots = dots[!sapply(dots,is.null)]
  #print(str(dots,max.level=2))
  dots = Map(function(x,i){
    if(is.null(dim(x))) {
      k = data.frame(x)
      colnames(k)=paste0("x",i)
    } else k = x
    k
    },dots, seq_along(dots))

  #print(str(dots))
  if(length(dots)>1){
    dotsNames = unlist(sapply(dots,colnames))
    #print(dotsNames)
    maxlen = max(sapply(dots, nrow))
    fdots = lapply(dots, function(x){
      #deb
      #x= dots[[2]]
      #rm(x,y,fdots,dots,pad,maxlen)
      if(nrow(x)<maxlen){
        pad = maxlen - nrow(x)
        y = data.frame(matrix(rep(NA,ncol(x)*pad),ncol = ncol(x)))
        colnames(y) = colnames(x)
        rownames(y) = paste0("NA",1:pad)
        rbind(x,y)
      } else x
    })
    res = data.frame(do.call("cbind",fdots))
    colnames(res) = dotsNames
    #print(colnames(res))
    return(res)
  } else return(dots[[1]])
}

prog = function(i,n,step=50){
  progStep = c(1,round((n/step)*2:(step-1)),n)
  progStep[progStep==0] = 1
  s = sum(i == progStep)
  if(s) {
    if(i==1) cat0(rep('.',50),"|100%\r\n")
    cat0(rep('.',s))}
  if(i==n) cat0("|Done ;)\r\n")

}

cat0 = function(...) {cat(..., sep="")}

lead0 = function(x, width = 2){
  formatC(x,width = width, format = "d", flag = "0")
}

## timeMaster allows to:
## transform time from different formats to different formats.
## add amounts of time to baseTime expressed in different formats.
## This is verified. And useful. For reasons. :-D

#' Transform time values between different formats
#'
#' This function allows to
#' @param baseTime,add either integer of seconds or a time string in the format h:m:s, m:s, or s, with or without leading zeroes
#' @param out a character string indicating the format of the output. One of "auto" (the default which tries to keep the format of \code{'baseTime'}), "hour", "min", or "sec": can be abbreviated.
#' @param baseSep a character string or a regular expression identifying separators in \code{baseTime}
#' @examples
#' ## Adding seconds to minutes
#' timeMaster(30, add="10:00", out = "min")
#'
#' ## Adding strings to integer seconds and returning a numeric value
#' timeMaster(30, add="10:00")
#'
#' ## Automatic detection of format
#' timeMaster("01:30:55",add="10:00",out = "auto")
#' @export
timeMaster = function(baseTime, add=0, out=c("auto", "hour", "min","sec"), baseSep = "[\\.,:,\\,',-,\"]"){
  #baseTime and add can either be  integers of seconds or a time string in the format h:m:s, m:s, or s, with or without leading zeroes
  #output forces the sum to be reported either as string h:m:s or m:s or as a integer of seconds. auto keeps the 'baseTime' format.
  out = match.arg(out)
  if(is.numeric(baseTime) && sum(baseTime%%1)>0 ) {
    baseTime=as.character(baseTime)
    warning("The . was considered as a ':' in the format mm:ss. timeMaster does not support fractional times yet")
  }
  if(length(baseTime)>1)
  {
    sapply(baseTime,timeMaster,add,out,baseSep,USE.NAMES = F)
  } else {
    #da qui baseTime è contenente un tempo singolo, non un vettore
    if(is.character(baseTime)){
      #è negativo?
      if(substr(baseTime,1,1)=="-"){
        negative = T
        baseTime = substring(baseTime,2)
      } else negative = F
      baseTime = strsplit(baseTime, split=baseSep)
      if      (length(baseTime[[1]])==1) auto = "sec"
      else if (length(baseTime[[1]])==2) auto = "min"
      else if (length(baseTime[[1]])==3) auto = "hour"
      else stop ("baseTime format not recognized. It should either be \"min:sec\" or \"hour:min:sec\" or an integer of seconds")
      sapply(unlist(baseTime), function(k){if(k=='') stop("baseTime contains an empty cell. Please check your separators")})
      #transform to seconds
      baseTime = unlist(lapply(baseTime, function(x){
        if(length(x)==1)
          as.numeric(x[1])
        else if(length(x)==2)
          as.numeric(x[1])*60 + as.numeric(x[2])
        else if(length(x==3))
          as.numeric(x[1])*3600 + as.numeric(x[2])*60 + as.numeric(x[3])
        else stop("Time format must either be \"min:sec\" or \"hour:min:sec\"")
      } ))
      if(negative) baseTime = -baseTime
    } else {
      auto = "sec"
    }
    if(is.character(add)) add = timeMaster(add,0,"sec")
    x = baseTime + add #that's the final amount in seconds

    if(out=="auto") out=auto
    if(out == "sec") {
      return(x)
    } else if(out %in% c("min", "hour")){
      if (x<0) { x = abs(x); negative = T} else negative =F
      mins = floor(x/60)
      secs = (x-mins*60)
      hours = trunc(mins/60)
      mins_h = mins - hours*60
      if(out == "min")
        return(paste0(ifelse(negative,"-",""), lead0(mins),":", lead0(secs)))
      else if (out =="hour"){
        return(paste0(ifelse(negative,"-",""),lead0(hours),":",lead0(mins_h),":", lead0(secs)))
      }
    } else stop("timeMaster failed in a weird way!")
  }
}

#Interpolate windowed data to original samplerate.
#useful to overlay computed indexes on original timeseries
winInter = function(windowsList, winSec, incSec, sampRate){
  #cat("Interpolating",incSec*sampRate,"sample between each HRV datapoint (linear) \r\n")
  inc=incSec*sampRate
  win= winSec*sampRate
  nList=length(windowsList)
  Map(function(daba,i){
    #prog(i,nList)
    res = data.frame(apply(daba,2,function (series){
      stats::approx(x    = seq(1,(length(series)*inc), by=inc)+ceiling(win/2)-1, #windowed datapoints must be at the center of the window
             y    = series, #the exact calculated values
             xout = seq(1,(length(series)-1)*inc + win,by=1) #all samples covered by the windows
      )$y
    }))
    colnames(res) = colnames(daba)
    res
  },windowsList,seq_along(windowsList))
}

#quick and dirty text sanitizing to create safe column names
sanitize = function(string){
  iconv(gsub("\\W","",string), #gsub deletes every non alphanumeric character
        "latin1", "ASCII", sub="")} #iconv deletes every remaining non ascii character


