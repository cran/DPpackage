\name{DPelicit}
\alias{DPelicit}
\alias{DPelicit.default}

\title{Performs a prior elicitation for the precision parameter of a DP prior}
       
\description{
    This function performs a prior elicitation for the precision parameter
    of a DP prior. The function calculates: 
    
    1) the expected value and the standard deviation of the number of clusters, 
       given the values of the parameters of the \code{gamma} 
       prior for the precision parameter, \code{a0} and \code{b0}, or
    
    2) the value of the parameters \code{a0} and \code{b0} of the \code{gamma} prior distribution
       for the precision parameter, \code{alpha}, given the prior expected number and the 
       standard deviation of the number of clusters.
}
  
  
\usage{
DPelicit(n,method='JGL',a0=NULL,b0=NULL,mean=NULL,std=NULL)
}

\arguments{
    \item{n}{        number of observations which distribution follows a
                     DP prior.}
    
    \item{method}{   the method to be used. See \code{details}.}

    \item{a0}{       hyperparameter for the \code{Gamma} prior distribution of 
                     the precision parameter of the Dirichlet process
                     prior, \code{alpha ~ Gamma(a0,b0)}. 
                     }

    \item{b0}{       hyperparameter for the \code{Gamma} prior distribution of 
                     the precision parameter of the Dirichlet process
                     prior, \eqn{alpha ~ Gamma(a0,b0)}. 
                     }

    \item{mean}{     prior expected number of clusters when 
                     \eqn{alpha ~ Gamma(a0,b0)}.}

    \item{std}{      prior standard deviation for the number of clusters 
                     when \eqn{alpha ~ Gamma(a0,b0)}.}
}

\details{

     The methods supported by these functions are based on the fact that a priori
     \code{E(alpha) = a0/b0} and \code{Var(alpha) = a0/b0^2}, and an additional 
     approximation based on Taylor series expansion.

     The default method, \code{"JGL"}, is based on the exact value of the mean and 
     the variance of the number of clusters given the precision parameter 
     alpha (see, Jara, Garcia-Zatera and Lesaffre, 2007).
    
     The Method \code{"KMQ"} is base on the Liu (1996) approximation to 
     the expected value and the variance of the number of clusters given the
     precision parameter alpha (see, Kottas, Muller and Quintana, 2005).
     
     Given the prior judgement for the mean and variance of the number of 
     clusters, the equations are numerically solve for \code{a0} and \code{b0}. 
     With this objective, the Newton-Raphson algorithm and the forward-difference 
     approximation to Jacobian are used.
}

\references{

Jara, A., Garcia-Zattera, M.J., Lesaffre, E. (2007) A Dirichlet Process mixture model for 
  the analysis of correlated binary responses. Computational Statistics and 
  Data Analysis 51: 5402-5415.

Kottas, A., Muller, P., Quintana, F. (2005) Nonparametric Bayesian modeling for multivariate
  ordinal data, Journal of Computational and Graphical Statistics 14: 610-625.

Liu, J.S. (1996) Nonparametric Hierarchical Bayes via Sequential Imputations, The Annals 
  of Statistics, 24: 911-930.

}

\examples{

    # Calculate the expected value and the standard deviation
    # for the number of cluster given alpha ~ Gamma(a0,b0).

      DPelicit(200,a0=2.01,b0=2.01,method="JGL")
      DPelicit(200,a0=2.01,b0=2.01,method="KMQ")

    # Calculate the values of a0 and b0, given the expected value 
    # and the standard deviation of the number of clusters

      DPelicit(200,mean=3.1,std=2.7,method="JGL")
      DPelicit(200,mean=3.1,std=2.7,method="KMQ")
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

}

\keyword{models}
\keyword{nonparametric}