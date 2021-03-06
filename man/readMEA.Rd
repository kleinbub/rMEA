% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rMEA_IO.R
\name{readMEA}
\alias{readMEA}
\title{Import MEA raw data}
\usage{
readMEA(
  path,
  s1Col,
  s2Col,
  sampRate,
  namefilt = NA,
  s1Name = "s1",
  s2Name = "s2",
  idOrder = c("id", "session", "group"),
  idSep = "_",
  removeShortFiles = NULL,
  ...
)
}
\arguments{
\item{path}{a character vector of full path names; may point to an individual
file or a directory containing MEA files. Only .txt or .csv file extensions
are considered in directories.}

\item{s1Col, s2Col}{the index of one or multiple columns in the data,
identifying the two dyad's members (e.g. patient and therapist) motion energy data. If multiple columns
are selected for a subject (e.g. because of multiple regions of interest per subject), their MEA values will be summed.}

\item{sampRate}{sampling rate at which the data is acquired (usually frames per
second of the original video recording).}

\item{namefilt}{either NA or a character string specifying a pattern to be matched in the filenames.
Regular expressions can be used.}

\item{s1Name, s2Name}{the label describing each participant. (e.g. Right/Left, or Patient/Therapist, etc).}

\item{idOrder}{either NA or a character vector that contains one or more of the three strings: "id",
"session","group" in a given order. These are used to interpret the
filenames and correctly label the cases. The strings can be abbreviated.
If the filenames contains other data the character "x" can be used to skip a position.
If NA, no attempt to identify cases will be done.}

\item{idSep}{character vector (or object which can be coerced to such) containing regular expression(s).
If idOrder is not NA, this will be used as separator to split the filenames and identify "id", "session", and "group"
informations.}

\item{removeShortFiles}{Either NULL or an number ranging from 0 to 1.
Specifies the proportion of the average file length below which a file should be excluded.
(E.g. a value of 0.7 will exclude all files with a duration smaller than 70\% of the mean duration of all other files in the directory.)}

\item{...}{Additional arguments passed to \code{\link[utils]{read.table}}. E.g. sep, skip, header, etc.}
}
\value{
an object of class \code{MEAlist}
}
\description{
\code{readMEA} reads the output of MEA software.
}
\details{
For instance if \code{s1Col = c(1,3)} and \code{s2Col = c(2,4)}, the
  returned values will be the sum of column 1 and 3 for the first participant and
  columns 2 and 4 for the second one.
}
\examples{
## read the first sample (intake interviews of patients that carried on therapy)
path_normal <- system.file("extdata/normal", package = "rMEA")
mea_normal <- readMEA(path_normal, sampRate = 25, s1Col = 1, s2Col = 2,
                     s1Name = "Patient", s2Name = "Therapist", skip=1,
                     idOrder = c("id","session"), idSep="_")
mea_normal <- setGroup(mea_normal, "normal")

## read the second sample (intake interviews of patients that dropped out)
path_dropout <- system.file("extdata/dropout", package = "rMEA")
mea_dropout <- readMEA(path_dropout, sampRate = 25, s1Col = 1, s2Col = 2,
                     s1Name = "Patient", s2Name = "Therapist", skip=1,
                     idOrder = c("id","session"), idSep="_")
mea_dropout <- setGroup(mea_dropout, "dropout")

## Combine into a single object
mea_all = c(mea_normal, mea_dropout)

summary(mea_all)



}
