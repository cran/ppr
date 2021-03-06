\name{smooth.spline}
\title{Fit a Smoothing Spline}
\usage{
smooth.spline(x, y, w=rep(1, length(x)), df=5, spar=0, cv=FALSE,
              all.knots=FALSE, df.offset=0, penalty=1)
}
\alias{smooth.spline}
\arguments{
 \item{x}{a vector giving the values of the predictor variable, or  a
   list or a two-column matrix specifying x and y. }
 \item{y}{responses. If \code{y} is missing, the responses are assumed
   to be specified by \code{x}.}
 \item{w}{optional vector of weights}
 \item{df}{the desired equivalent number of degrees of freedom (trace of
 the smoother matrix).}
 \item{spar}{the coeffcient "lambda" of the integral of the squared
   second derivative in the fit criterion.}
 \item{cv}{ordinary (TRUE) or `generalized' (FALSE) cross-validation.}
 \item{all.knots}{if TRUE, all points in \code{x} are uses as knots. If
   FAlSE, a suitably fine grid of knots is used.}
 \item{df.offset}{allows the degrees of freedom to be increased by
   \code{df.offset} in the GCV criterion.}
 \item{penalty}{the coefficient of the penalty for degrees of freedom
   in the GCV criterion.}
}
\description{
 Fits a cubic smoothing spline to the supplied data.
}
\details{
  The \code{x} vextor should contain at lest ten distinct values.

  If \code{spar} is missing or 0, the value of \code{df} is used to
  determine the degree of smoothing. If both are missing, leave-one-out
  cross-validation is used to determine "lambda".
}
\value{
  An object of class \code{"smooth.spline"} with components
  \item{x}{the distinct \code{x} values in increasing order.}
  \item{y}{the fitted values corresponding to \code{x}.}
  \item{w}{the weights used at the unique values of \code{x}.}
  \item{yin}{the y values used at the unique \code{y} values.}
  \item{lev}{leverages, the diagonal values of the smoother matrix.}
  \item{cv.crit}{(generalized) cross-validation score.}
  \item{pen.crit}{penalized criterion}
  \item{df}{equivalent degrees of freedom used.}
  \item{spar}{the value of "lambda" chosen.}
  \item{fit}{list for use by \code{predict.smooth.spline}.}
  \item{call}{}
}
\author{B.D. Ripley}

\seealso{\code{\link{predict.smooth.spline}}}

\examples{
data(cars)
attach(cars)
plot(speed, dist)
cars.spl <- smooth.spline(speed, dist)
print(cars.spl)
lines(cars.spl)
lines(smooth.spline(speed, dist, df=10), lty=2)
detach()
}
\keyword{smooth}
