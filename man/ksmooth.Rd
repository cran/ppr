\name{ksmooth}
\title{Kernel Regression Smoother}
\usage{
ksmooth(x, y, kernel=c("box", "normal"), bandwidth=0.5,
        range.x=range(x), n.points=max(100, length(x)), x.points)
}
\alias{ksmooth}
\arguments{
 \item{x}{input x values}
 \item{y}{input y values}
 \item{kernel}{The kernel to be used.}
 \item{bandwidth}{the bandwidth. The kernels are scaled so that their
   quartiles (viewed as probability densities) are at
   \code{+/-0.25*bandwidth}.}
 \item{range.x}{the range of points to be covered in the output.}
 \item{n.points}{the number of points at which to evaluate the fit.}
 \item{x.points}{points at which to evaluate the smoothed fit. If
   missing, \code{n.points} are chosen uniformly to cover \code{range.x}.}
}
\description{
 The Nadaraya-Watson kernel regression estimate.}
\value{
  A list with components
  \item{x}{values at which the smoothed fit is evaluated. Guaranteed to
    be in increasing order.}
  \item{y}{fitted values corresponding to \code{x}.}
}
\author{B.D. Ripley}
\note{This function is implemented purely for compatibility with S,
  although it is nowhere near as slow as the S function. Better kernel
  smoothers are available in other packages.}
\examples{
data(cars)
attach(cars)
plot(speed, dist)
lines(ksmooth(speed, dist, "normal", bandwidth=2), col=2)
lines(ksmooth(speed, dist, "normal", bandwidth=5), col=3)
lines(ksmooth(speed, dist, "normal", bandwidth=10), col=4)
}
\keyword{smooth}
