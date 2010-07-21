\name{ps}
\alias{ps}
\title{Specify a smoothing spline fit in a PSgam formula}
\description{A symbolic wrapper to indicate a smooth term in a formula
             argument to PSgam}

\usage{
ps(..., k=50,degree=3,pord=1)
}

\arguments{
     \item{...}{the predictors.}

     \item{k}{an integer giving the number of knots. The number of basis functions is k+degree.}
     
     \item{degree}{an integer giving the degree of the B-splines, e.g. degree=3 gives a cubic spline.}

     \item{pord}{an integer giving the order of difference penalty.}
}

\value{ \code{ps} returns the vectors of predictors, endowed with a number of
attributes. The vector itself is used in the construction of the model
matrix.

Note that \code{ps} does not do the smoothing; it simply sets things up
for \code{PSgam}.
}

\references{
Eilers, P.H.C. and Marx, B.D. (1996) Flexible Smoothing with B-splines and 
  penalties. Statistical Science, 11(2): 89-121.
}

\author{
  Alejandro Jara \email{<atjara@uc.cl>}
 }

\examples{
   # fit Start using a smoothing spline with 4 df.
     y ~ Age + ps(Start, degree=4)
   # fit log(Start) using a smoothing spline with 5 df.
     y ~ Age + ps(log(Start), degree=5)
}

\keyword{models}
\keyword{regression}
\keyword{nonparametric}
\keyword{smooth}
