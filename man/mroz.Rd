% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mroz.R
\docType{data}
\name{mroz}
\alias{mroz}
\title{mroz}
\format{
A data.frame with 428 observations on 22 variables:
\itemize{
 \item \strong{inlf:} =1 if in lab frce, 1975
 \item \strong{hours:} hours worked, 1975
 \item \strong{kidslt6:} # kids < 6 years
 \item \strong{kidsge6:} # kids 6-18
 \item \strong{age:} woman's age in yrs
 \item \strong{educ:} years of schooling
 \item \strong{wage:} est. wage from earn, hrs
 \item \strong{repwage:} rep. wage at interview in 1976
 \item \strong{hushrs:} hours worked by husband, 1975
 \item \strong{husage:} husband's age
 \item \strong{huseduc:} husband's years of schooling
 \item \strong{huswage:} husband's hourly wage, 1975
 \item \strong{faminc:} family income, 1975
 \item \strong{mtr:} fed. marg. tax rte facing woman
 \item \strong{motheduc:} mother's years of schooling
 \item \strong{fatheduc:} father's years of schooling
 \item \strong{unem:} unem. rate in county of resid.
 \item \strong{city:} =1 if live in SMSA
 \item \strong{exper:} actual labor mkt exper
 \item \strong{nwifeinc:} (faminc - wage*hours)/1000
 \item \strong{lwage:} log(wage)
 \item \strong{expersq:} exper^2
}
}
\source{
\url{https://www.cengage.com/cgi-wadsworth/course_products_wp.pl?fid=M20b&product_isbn_issn=9781111531041}
}
\usage{
data(mroz)
}
\description{
Data provided by Ernst R. Berndt, which was used by Mroz (1987) and Wooldridge (2010)
}
\examples{
data(mroz)

}
\references{
{
Mroz, Thomas, (1987), The Sensitivity of an Empirical Model of Married Women's Hours of Work to Economic and Statistical Assumptions, \emph{Econometrica}, 55, issue 4, p. 765-99. \cr
Jeffrey M Wooldridge, 2010. "Econometric Analysis of Cross Section and Panel Data," MIT Press Books, The MIT Press, edition 2, volume 1, number 0262232588, December. \cr
}
}
\keyword{datasets}
