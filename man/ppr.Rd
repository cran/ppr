\name{ppr}
\alias{ppr}
\title{
Projection Pursuit Regression
}
\usage{
ppr(formula, data=sys.parent(), weights, subset,
	na.action, contrasts=NULL, ww=rep(1,q),
	nterms, max.terms=nterms, optlevel=2, 
	sm.method=c("supsmu", "spline", "gcvspline"),
	bass=0, span=0, df=5, gcvpen=1)
ppr(x, y, weights=rep(1,n), ww=rep(1,q), nterms,
	max.terms=nterms, optlevel=2,
	sm.method=c("supsmu", "spline", "gcvspline"),
	bass=0, span=0, df=5, gcvpen=1)
}
\arguments{
\item{formula}{
a regression formula specifying one or more response variables and the
explanatory variables.
}
\item{x}{
matrix of explanatory variables.  Rows represent observations, and
columns represent variables.  Missing values are not accepted. 
}
\item{nterms}{
number of terms to include in the final model.
}
\item{data}{
Data frame from which variables specified in \code{formula} are
preferentially to be taken.
}
\item{weights}{
a vector of weights for each case.
}
\item{ww}{
a vector of weights for each response, so the fit criterion is
the sum over case \code{i} and responses \code{j} of
\code{w_i ww_j (y_ij - fit_ij)^2} divided by the sum of \code{w_i}.
}
\item{subset}{
An index vector specifying the cases to be used in the training
sample.  (NOTE: If given, this argument must be named.)
}
\item{na.action}{
A function to specify the action to be taken if \code{NA}s are found. The
default action is for the procedure to fail.  An alternative is
\code{na.omit}, which leads to rejection of cases with missing values on
any required variable.  (NOTE: If given, this argument must be named.)
}
\item{contrasts}{
the contrasts to be used when any factor explanatory variables are coded.
}
\item{max.terms}{
maximum number of terms to choose from when building the model.
}
\item{optlevel}{
integer from 0 to 3 which determines the throughness of an
optimization routine in the SMART program. See the DETAILS section.
}
\item{sm.method}{
the method used for smoothing the ridge functions.  The default is to
use Friedman's super smoother \code{supsmu}.  The alternatives are to use
the smoothing spline code underlying \code{smooth.spline}, either with a
specified (equivalent) degrees of freedom for each ridge functions, or
to allow the smoothness to be chosen by GCV.
}
\item{bass}{
super smoother bass tone control used with automatic span selection
(see \code{supsmu}); the range of values is 0 to 10, with larger values
resulting in increased smoothing.
}
\item{span}{
super smoother span control (see \code{supsmu}).  The default, \code{0},
results in automatic span selection by local cross validation. \code{span}
can also take a value in \code{(0, 1]}.
}
\item{df}{
if \code{sm.method} is \code{"spline"} specifies the smoothness of each ridge
term via the requested equivalent degrees of freedom.
}
\item{gcvpen}{
if \code{sm.method} is \code{"gcvspline"} this is the penalty used in the GCV
selection for each degree of freedom used.
}}
\value{
A list with the following components, many of which are for use by the
method functions.

\item{call}{
the matched call
}
\item{p}{
the number of explanatory variables (after any coding)
}
\item{q}{
the number of response variables
}
\item{ml}{
the argument \code{max.terms}
}
\item{gof}{
the overall residual (weighted) sum of squares for the selected model
}
\item{gofn}{
the overall residual (weighted) sum of squares against the number of
terms, up to \code{max.terms}.  Will be invalid (and zero) for less than
\code{nterms}.
}
\item{df}{
the argument \code{df}
}
\item{edf}{
if \code{sm.method} is \code{"spline"} or \code{"gcvspline"} the equivalent number of
degrees of freedom for each ridge term used.
}
\item{xnames}{
the names of the explanatory variables
}
\item{ynames}{
the names of the response variables
}
\item{alpha}{
a matrix of the projection directions, with a column for each ridge term
}
\item{beta}{
a matrix of the coefficients applied for each response to the ridge
terms: the rows are the responses and the columns the ridge terms
}
\item{yb}{
the weighted means of each response
}
\item{ys}{
the overall scale factor used: internally the responses are divided by
\code{ys} to have unit total weighted sum of squares.
}
\item{fitted.values}{
the fitted values, as a matrix if \code{q > 1}
}
\item{residuals}{
the residuals, as a matrix if \code{q > 1}
}
\item{smod}{
internal work array, which includes the ridge functions evaluated at
the training set points.
}}
\description{
The basic method is given by Friedman (1984), and is essentially the
same code used by \code{ppreg}.  The answers will be very similar on a givn
machine, but this code is extremely sensitive to the compiler used.
The differences are the ability to use spline smoothers and the
interface which should be much easier to use.


The algorithm first adds up to \code{max.terms} ridge terms one at a time;
it will use less if it is unable to find a term to add that makes
sufficient difference.  It then removes the least "important" term at
each step until \code{nterm} terms are left.  The levels of optimization
differ in how thoroughly the models are refitted during this process.
At level 0 the existing ridge terms are not refitted. At level 1 the
projection directions are not refitted, but the ridge functions and
the regression coefficients are.  Levels 2 and 3 refit all the terms
and are equivalent for one response; level 3 is more careful to
re-balance the contributions from each regressor at each step and so
is a little less likely to converge to a saddle point of the sum of
squares criterion.
}
\references{
Friedman, J. H. and Stuetzle, W. (1981)
Projection pursuit regression. 
Journal of the American Statistical Association
\bold{76}, 817-823.


Friedman, J. H. (1984)
SMART User's Guide.
Laboratory for Computational Statistics, Stanford University Technical
Report No. 1.
}
\seealso{
\code{\link{plot.ppr}},  \code{\link{ppreg}},  \code{\link{supsmu}},  \code{\link{smooth.spline}}
}
\examples{
# Note: your numerical values may differ
attach(rock)
area1 <- area/10000; peri1 <- peri/10000
rock.ppr <- ppr(log(perm) ~ area1 + peri1 + shape,
		data=rock, nterms=2, max.terms=5)
rock.ppr
Call:
ppr.formula(formula = log(perm) ~ area1 + peri1 + 
	shape, data = rock, nterms = 2, max.terms = 5)


Goodness of fit:
 2 terms 3 terms 4 terms 5 terms 
 11.2196  7.1895  6.4565  5.8592


summary(rock.ppr)
Call:
ppr.formula(formula = log(perm) ~ area1 + peri1 + 
	shape, data = rock, nterms = 2, max.terms = 5)


Goodness of fit:
 2 terms 3 terms 4 terms 5 terms 
 11.2196  7.1895  6.4565  5.8592


Projection direction vectors:
         term 1    term 2 
 area  0.319492  0.435617
 peri -0.945544 -0.866757
shape  0.062226  0.242839


Coefficients of ridge terms:
  term 1  term 2 
 1.00638 0.72915
par(mfrow=c(1,2), pty="s")
plot(rock.ppr)
plot(update(rock.ppr, bass=5))
plot(update(rock.ppr, sm.method="gcv", gcvpen=2))
}
\keyword{regression}
% Converted by Sd2Rd version 0.2-a5.
