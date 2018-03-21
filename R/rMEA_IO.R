#' Import MEA raw data
#'
###############################################################################|
#' \code{readMEA} reads the output of MEA software.

#' @param path a character vector of full path names; may point to an individual
#'   file or a directory containing MEA files. Only .txt or .csv file extensions
#'   are considered in directories.
#' @param s1Col,s2Col the index of one or multiple columns in the data,
#'   identifying the two dyad's members (e.g. patient and therapist) motion energy data. If multiple columns
#'   are selected for a subject (e.g. because of multiple regions of interest per subject), their MEA values will be summed.
#' @param sampRate sampling rate at which the data is acquired (usually frames per
#'   second of the original video recording).
#' @param namefilt either NA or a character string specifying a pattern to be matched in the filenames.
#'  Regular expressions can be used.
#' @param s1Name,s2Name the label describing each participant. (e.g. Right/Left, or Patient/Therapist, etc).
#' @param idOrder either NA or a character vector that contains one or more of the three strings: "id",
#'   "session","group" in a given order. These are used to interpret the
#'   filenames and correctly label the cases. The strings can be abbreviated.
#'   If the filenames contains other data the character "x" can be used to skip a position.
#'   If NA, no attempt to identify cases will be done.
#' @param idSep If idOrder is not NA, this character will be used as separator between "id", "session", and "group"
#' information in the filenames.
#' @param removeShortFiles Either NULL or an number ranging from 0 to 1.
#'  Specifies the proportion of the average file length below which a file should be excluded.
#'  (E.g. a value of 0.7 will exclude all files with a duration smaller than 70\% of the mean duration of all other files in the directory.)
#' @param ... Additional arguments passed to \code{\link[utils]{read.table}}. E.g. sep, skip, header, etc.
#'
#' @details For instance if \code{s1Col = c(1,3)} and \code{s2Col = c(2,4)}, the
#'   returned values will be the sum of column 1 and 3 for the first participant and
#'   columns 2 and 4 for the second one.
#' @return an object of class \code{MEAlist}
#' @examples ## read the first sample (intake interviews of patients that carried on therapy)
#' path_normal <- system.file("extdata/normal", package = "rMEA")
#' mea_normal <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
#'                      s1Name = "Patient", s2Name = "Therapist", skip=1,
#'                      idOrder = c("id","session"), idSep="_")
#'mea_normal <- setGroup(mea_normal, "normal")
#'
#' ## read the second sample (intake interviews of patients that dropped out)
#'path_dropout <- system.file("extdata/dropout", package = "rMEA")
#'mea_dropout <- readMEA(path_dropout, sampRate = 25, s1Col = 1, s2Col = 2,
#'                      s1Name = "Patient", s2Name = "Therapist", skip=1,
#'                      idOrder = c("id","session"), idSep="_")
#'mea_dropout <- setGroup(mea_dropout, "dropout")
#'
#'## Combine into a single object
#'mea_all = c(mea_normal, mea_dropout)
#'
#'summary(mea_all)
#'
#'
#'

#' @export
readMEA = function(
  path, s1Col, s2Col, sampRate, namefilt = NA,
  s1Name = "s1", s2Name = "s2",
  idOrder= c("id","session","group"),
  idSep = "_",
  removeShortFiles = NULL,
  ... #additional options passed to read.table (es: skip, header, etc)
){
  # ####debug
  # path = example_data_path
  # path = system.file("data",  package="rMEA")
  # path = "/Users/kleinbub/OneDrive - Università degli Studi di Padova/__Synchro Analysis/__RLab/MEA/wannabe files/data"
  # s1Col = c(1)
  # s2Col = c(2)
  # sampRate = 30
  # namefilt = ""
  # idOrder = c("i","g")
  # idSep = "_"
  # removeShortFiles = 0.7
  # s1Name = "asd";  s2Name = "lol"
  # options=list("skip"=1, "sep"=" ", "header"=F)

  if(length(s1Col)!=length(s2Col)) warning ("s1Col and s2Col have different lengths. Check if this is intended.")

  iStep=0
  if(!file.exists(path)){stop("Selected file or directory does not exists")}
  patt = "\\.(txt|csv)$"
  if (!is.na(namefilt)){
    patt = paste0(namefilt,".*\\.(txt|csv)$")
  }
  nFiles=length(list.files(path, pattern=patt))
  if(nFiles==0){ #no files in the directory. Maybe a single file was provided?
    if(utils::file_test("-f",path)){ #check if file exists and not a directory
      filenames = path
      shortNames = basename(path)
      path=dirname(path)
    } else{
      stop("No [",patt,"] files were found in ",path,"\r\n\r\nFiles in path:\r\n",paste(list.files(path),"\r\n" ), call.=F)
      #stop("The input must be a file or a directory!\r\nTry to use choose.dir() as the first parameter, e.g.: meaMaster(choose.dir(), ....")
    }
  } else {
    filenames = list.files(path, pattern=patt, full.names=T)
    shortNames = list.files(path, pattern=patt, full.names=F)
  }
  nFiles = length(filenames)


  ####now tries to identify cases!

  #first check if separator exists
  if(any(!is.na(idOrder)))
    lapply(shortNames, function(name) {
      if(!grepl(idSep,tools::file_path_sans_ext(name))) stop("idSep: '",idSep, "' could not be found in file: ",name, call.=F)
    })

  #then build lists of groups, sessions, and id for each file
  idOrder = tolower(idOrder) #lowercase
  if(any(substr(idOrder,1,1)=="i",na.rm=T)){
    dyadIds = lapply(seq_along(shortNames), function(i) {
      unlist(strsplit(tools::file_path_sans_ext(shortNames[i]), idSep))[which(substr(idOrder,1,1)=="i")]
    })
  } else dyadIds = as.list(lead0(seq_len(nFiles),width=ifelse(nchar(nFiles)<2 ,2,nchar(nFiles)) ))
  if(any(substr(idOrder,1,1)=="g",na.rm=T)){
    group = lapply(seq_along(shortNames), function(i) {
      unlist(strsplit(tools::file_path_sans_ext(shortNames[i]), idSep))[which(substr(idOrder,1,1)=="g")]
    })
  } else group = as.list(rep("all",nFiles))
  if(any(substr(idOrder,1,1)=="s",na.rm=T)){
    sess = lapply(seq_along(shortNames), function(i) {
      ax = unlist(strsplit(tools::file_path_sans_ext(shortNames[i]), idSep))[which(substr(idOrder,1,1)=="s")]
      x = as.numeric(gsub("[[:alpha:]]","", ax))
      if(is.na(x)){
        warning("No numeric information was found for session identifier '", ax, "' in signal ",shortNames[i],". Please check filenames and idOrder and idSep arguments:\r\n", call.=F)
        ax
      } else x
    })
  } else sess = as.list(rep("01",nFiles))


#if sessions are specified, check their order
  if(any(substr(idOrder,1,1)=="s",na.rm=T) ){
    #now to check if the filenames have repetitions

    nCheck = do.call(paste0, list(group,dyadIds, sess) ) #paste together group,id and session to obtain unique session identifier
    if(length(unique(nCheck)) != length(nCheck))
    warning("Two equal session identifier were found, please check that the session order/names are correct:\r\n",paste0("\r\n",shortNames,"\t",nCheck),call. = F)

    #now check if the session are in progressive order
    deltaSess = sapply(split(unlist(sess),unlist(dyadIds)), function(x) ifelse(length(x)>1,diff(x),1)) #if there are multiple sessions with the same id, check that they are in progressive order
    if(any(deltaSess<=0))
      warning("The sessions may not be in sequential order, please check file names:\r\n",paste0("\r\n",shortNames,"\t",nCheck),call. = F)
  }


  iStep = iStep +1
  cat("\r\nSTEP",iStep,"| Reading",nFiles,"dyads\r\n")
  options = list(...)
  #ugly stuff to set read.table options
  if("skip" %in% names(options)){
    if(length(options$skip) ==1) options$skip=rep(options[["skip"]],nFiles)
    if(length(options$skip) !=nFiles) stop("skip must be provided for each file")
    skipRow = options$skip
    options$skip = NULL
  } else skipRow =rep(0,nFiles)
  lf <- mapply(function(x,iFile) {  prog(iFile,nFiles); do.call(utils::read.table,c(list(x, skip=skipRow[iFile]), options)) },filenames,seq_along(filenames),SIMPLIFY = F )
  if(ncol(lf[[1]])==1) {print(utils::str(lf[[1]]));stop("Import failed. Check 'sep' argument?",call.=F)}
  if(!is.numeric(unlist(lf[[1]][s1Col]))) stop("s1Col column is not numeric. Maybe there is a header? Set 'skip' argument to 1 or more?",call.=F)
  if(!is.numeric(unlist(lf[[1]][s2Col]))) stop("s2Col column is not numeric. Maybe there is a header? Set 'skip' argument to 1 or more?",call.=F)

  len <- lapply(lf, function(x) length(x[[1]]))
  #check if some file are shorter than 70% of the mean length and remove them
  if(!is.null(removeShortFiles)){
    removeFile = unlist(lapply(seq_along(len), function(i){if(len[[i]] / mean(unlist(len)) < removeShortFiles) {
      warning("File ",i,': ',shortNames[i]," is shorter than the ",removeShortFiles*100,"% of the mean length and was removed", call. = F);return(i)}
    } ))
    if(length(removeFile)>0){
      lf[removeFile] = len[removeFile] = dyadIds[removeFile] = group[removeFile] = sess[removeFile] = NULL
      filenames =filenames[-removeFile]; shortNames = shortNames[-removeFile]
      nFiles = length(filenames)
    }
  }


  #sums the ROI and formats the output
  iStep = iStep +1
  cat("\r\nSTEP",iStep,"| Formatting data frames:\r\n")
  lf = Map(function(x,i){
    prog(i,nFiles)
    if(length(s1Col)>1){
      x$s1 = as.numeric(apply(x[s1Col],1,sum))
    } else x$s1 = as.numeric(x[,s1Col])
    if(length(s2Col)>1){
      x$s2 = as.numeric(apply(x[s2Col],1,sum))
    } else x$s2 = as.numeric(x[,s2Col])
    x[c("s1", "s2")]
  }, lf, seq_along(lf))


  #check if any MEA is > 10sd
  Map(function(x,i){
    myDat = c(x$s1,x$s2)
    mySD = stats::sd(myDat)
    if(any(myDat > 10*mySD))
      warning( round(length(myDat[myDat > 10*mySD])/length(myDat)*100,2) ,"% of the data was higher than 10 standard deviations in dyad: ",dyadIds[[i]], ", session: ",sess[[i]],", group:",group[[i]], ". Check the raw data!", call. = F)
  },lf,seq_along(lf))

  # rename columns
  s1Name = sanitize(s1Name)
  s2Name = sanitize(s2Name)
  lf = Map(function(x,i){
    names(x) = c(s1Name,s2Name)
    x
  }, lf, seq_along(lf))

  ##CALCULATE PERCENTAGES BEFORE FILTERING ###########################
  s1Perc = sapply(lf, function(x){length(x[[1]][x[[1]]>0])/nrow(x)})#
  s2Perc = sapply(lf, function(x){length(x[[2]][x[[2]]>0])/nrow(x)})#
  ####################################################################


  ###############################
  ## reading report
  iStep = iStep +1
  cat("\r\nSTEP",iStep,"| ReadMEA report\r\n")
  outtable =data.frame(
    "Filename" = unlist(shortNames),
    "id_dyad" = unlist(dyadIds),
    "session" = unlist(sess),
    "group" = unlist(group),
    "duration_hh.mm.ss" = timeMaster(floor(unlist(len)/sampRate), out="hour"),
    "s1Perc" = round(s1Perc*100,1),
    "s2Perc" = round(s2Perc*100,1),
    row.names = NULL
  )
  colnames(outtable)[colnames(outtable)%in%c("s1Perc","s2Perc")] = c(paste0(s1Name,"_%"),paste0(s2Name,"_%"))
  print(outtable)


  #Populates the objects
  experiment = Map(function(file, i){
    MEA(file, sampRate = sampRate, filter = "raw",
        id = dyadIds[[i]], session = sess[[i]], group = group[[i]],
        s1Name = s1Name, s2Name = s2Name
        )
  }, lf, seq_along(lf))
  names(experiment) = mapply(paste,group,dyadIds,sess, MoreArgs=list(sep="_"))

    # class(experiment) = c("MEAlist",class(experiment))
    # attributes(experiment) = c(attributes(experiment), list(
    #   nId = length(unique(unlist(dyadIds))),
    #   n = nFiles,
    #   groups = unique(unlist(group)),
    #   sampRate = sampRate,
    #   s1Name = s1Name,
    #   s2Name = s2Name
    # ))

  return(MEAlist(experiment))

}


#' Exports analyzed MEA data to .txt files
#'
#' @param mea an object of class \code{MEA} or \code{MEAlist} (see function \code{\link{readMEA}})
#' @param save.directory a character string naming a directory
#' @param mea_data logical. Should the (filtered) MEA data be included in the export
#' @param ccf_data logical. Should the cross correlation results be exported?
#' @param ...  further arguments passed to \code{\link[utils]{write.table}}
#'
#' @details If both \code{mea_data} and \code{ccf_data} are \code{TRUE}, the cross-correlation data will be linearly interpolated to match the sampling rate of MEA data.
#'
#' @export
#' @examples
#' \donttest{
#' ## This example is excluded from test as it takes more than 10s to run
#' ## define a regex filter for the filenames to be read
#' to_read = c("20[123]_*")
#'
#' ## read the first 4 minutes of the files
#' path_normal <- system.file("extdata/normal", package = "rMEA")
#' mea_normal <- readMEA(path_normal, namefilt = to_read, sampRate = 25, s1Col = 1, s2Col = 2,
#'                      s1Name = "Patient", s2Name = "Therapist",
#'                      idOrder = c("id","session"), idSep = "_",  skip = 1, nrow = 6000)
#'
#' ## perform ccf analysis
#' mea_ccf = MEAccf(mea_normal, lagSec = 5, winSec = 60, incSec = 30, r2Z = TRUE, ABS = TRUE)
#'
#' ## export data and analysis
#' save_path = tempdir()
#' writeMEA(mea_ccf, save.directory = save_path, mea_data = TRUE, ccf_data = TRUE)
#'}
writeMEA <- function(mea, save.directory, mea_data, ccf_data, ...){
  UseMethod("writeMEA",mea)
}

#' @export
writeMEA.MEA <-  function(mea, save.directory, mea_data=T, ccf_data=T, ...){
  ###DEBUG
  # mea = mea3$all_38508_1
  # save.directory = "/Users/kleinbub/Desktop"
  ########


  if(!methods::is(mea,"MEA")) stop("Only objects of class MEA can be processed by this function")
  if(!file.exists(save.directory)){stop("save.directory must be an existing directory")}
  if(mea_data){
    Q = mea$MEA
  }
  if(ccf_data){
    if(is.null(mea$ccf)) stop("No CCF calculation was found. Set ccf_data to FALSE, or run MEAccf()",call. = F)
    if(mea_data){
      Qccf = winInter(list(mea$ccf),winSec = attr(mea,"ccf")$win, incSec = attr(mea,"ccf")$inc,sampRate = attr(mea,"sampRate"))[[1]]
    } else {
      Qccf = mea$ccf
    }
    origCols = ncol(Qccf)
    lag0 = which(colnames(Qccf)=="lag0")
    Qccf$grandAverage = apply(Qccf[1:origCols],1,mean,na.rm=T)
    Qccf$leading = apply(Qccf[1:(lag0-1)],1,mean,na.rm=T)
    Qccf$pacing =  apply(Qccf[(lag0+1):origCols],1,mean,na.rm=T)
    Qccf$bestLag = apply(Qccf[1:origCols],1,  function(r){best = which.max(r); ifelse(length(best)!=0,as.numeric(gsub("lag","", names(best))), NA) }) #this may require a smoothing function
    if(mea_data){
      Q = unequalCbind(Q,Qccf)
    } else {
      Q = Qccf
    }
  }
  dots <-  list(...)
  if(!"row.names" %in% names(dots)) dots[["row.names"]]=FALSE
  fileName=paste0(save.directory,ifelse(substr(save.directory, nchar(save.directory), nchar(save.directory))%in%c("/","\\"),"","/"),attr(mea,"uid"),".csv")
  do.call(utils::write.table, c(list(Q,fileName), dots))
}

#' @export
writeMEA.MEAlist <-  function(mea, save.directory, mea_data, ccf_data, ...){
  nFiles= length(mea)
  cat("Exporting",nFiles,"dyads\r\n")
  Map(function(iMEA,i){
    prog(i,nFiles);
    writeMEA(iMEA,save.directory,mea_data,ccf_data,...)
  }, mea, seq_along(mea))
  invisible()
}




