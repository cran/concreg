\name{vcov.concreg}
\alias{vcov.concreg}
\title{
Method to extract the sandwich covariance matrix for the regression coefficients from a \code{concreg} object
}
\description{
 This method returns the sandwich covariance matrix for the estimated regression coefficients from a \code{concreg} object. }
\usage{
\method{vcov}{concreg}(object,  ...)
}
\arguments{
  \item{object}{ 
  a \code{concreg} object
%%     ~~Describe \code{obj} here~~
}
\item{...}{
  additional argument(s) for methods.
}
}
\value{
 A pxp covariance matrix for the p regression coefficients. 
}
\author{
Georg Heinze
}
\seealso{
\code{confint.concreg}
\code{coef.concreg}
}
\examples{
gastric <-
  structure(list(patnr = as.integer(c(46, 1, 2, 3, 4, 5, 47, 6,
                   7, 8, 9, 48, 10, 11, 49, 12, 13, 14, 50, 15, 16, 17, 18, 19,
                   20, 51, 21, 22, 52, 23, 53, 54, 55, 24, 25, 56, 57, 58, 59, 60,
                   61, 62, 63, 64, 26, 65, 27, 66, 28, 29, 67, 68, 69, 70, 30, 71,
                   31, 72, 32, 73, 33, 34, 74, 75, 76, 77, 78, 35, 79, 36, 80, 81,
                   82, 37, 38, 39, 83, 84, 40, 85, 41, 86, 87, 88, 42, 43, 44, 89,
                   90, 45)),
                 treat = as.integer(c(0, 1, 1, 1, 1, 1, 0, 1, 1, 1,
                   1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0,
                   0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0,
                   0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,
                   1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1)),
                 time = as.integer(c(1,
                   17, 42, 44, 48, 60, 63, 72, 74, 95, 103, 105, 108, 122, 125,
                   144, 167, 170, 182, 183, 185, 193, 195, 197, 208, 216, 234, 235,
                   250, 254, 262, 301, 301, 307, 315, 342, 354, 356, 358, 380, 383,
                   383, 388, 394, 401, 408, 445, 460, 464, 484, 489, 499, 523, 524,
                   528, 535, 542, 562, 567, 569, 577, 580, 675, 676, 748, 778, 786,
                   795, 797, 855, 955, 968, 977, 1174, 1214, 1232, 1245, 1271, 1366,
                   1420, 1455, 1460, 1516, 1551, 1585, 1622, 1626, 1690, 1694, 1736
                   )),
                 status = as.integer(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
                   0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0))),
            .Names = c("patnr",
              "treat", "time", "status"), class = "data.frame",
            row.names = c("1",
              "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
              "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24",
              "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35",
              "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46",
              "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57",
              "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68",
              "69", "70", "71", "72", "73", "74", "75", "76", "77", "78", "79",
              "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90"
              ))
  
library(survival)
fit<-concreg(data=gastric, Surv(time,status)~treat)
coef(fit)
vcov(fit)
 }
\keyword{survival}
\keyword{regression}
\keyword{models}
