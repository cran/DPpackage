\name{DPbetabinom}
\alias{DPbetabinom}
\alias{DPbetabinom.default}

\title{Bayesian Semiparametric Beta-Binomial Model using a DP prior}
\description{
    This function generates a posterior density sample for a 
    semiparametric version of the Beta-Binomial model using
    a Dirichlet process prior for the mixing distribution.
}
  
\usage{
DPbetabinom(y,ngrid,prior,mcmc,state,status,
           data=sys.frame(sys.parent()),work.dir=NULL)
}

\arguments{
    \item{y}{        a matrix giving the binomial data. The first column must include the 
                     number of sucess and the second column the number of trials.}

    \item{ngrid}{    number of grid points where the predictive density estimate is 
                     evaluated. }
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Dirichlet process
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{a0} is missing, see details
                     below), and \code{a1} and \code{b1} giving the parameters of the
                     beta centering distribution.}

    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \code{nburn} giving the number of burn-in 
                     scans, \code{nskip} giving the thinning interval, \code{nsave} giving
                     the total number of scans to be saved, and \code{ndisplay} giving
                     the number of saved scans to be displayed on screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out).}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}
    
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state}.}

    \item{data}{     data frame.}       
    
     \item{work.dir}{working directory.} 
}

\details{
  This generic function fits a semiparametric version of the Beta-Binomial model (Liu, 1996):
  \deqn{y_i | n_i, p_i \sim Binom(n_i,p_i), i=1,\ldots,n}{yi | ni, pi ~ Binom(ni,pi), i=1,\ldots,n}
  \deqn{p_i | G \sim G}{pi | G ~ G}
  \deqn{G | \alpha, G_0 \sim DP(\alpha G_0)}{G | alpha, G0 ~ DP(alpha G0)}
  
  where, the baseline distribution is the beta distribution,
  \deqn{G_0 = Beta(a_1,b_1)}{G0 = Beta(a1,b1)} 
  
  To complete the model specification, the following hyperprior can be assumed for the total mass parameter:
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}

  Notice that the baseline distribution, \eqn{G_0}{G0}, is a conjugate prior in
  this model specification. Therefore, standard algorihtms for conjugate DP models are used (see, e.g.,
  Escobar and West, 1995; MacEachern, 1998). 
}

\value{
  An object of class \code{DPbetabinom} representing the DP Beta-Binomial model fit. 
  Generic functions such as \code{print}, \code{summary}, and \code{plot} have methods to 
  show the results of the fit. The results include the baseline parameters, \code{alpha}, and the 
  number of clusters.
  
  The MCMC samples of the parameters in the model are stored in the object 
  \code{thetasave}. The object is included in the 
  list \code{save.state} and are matrices which can be analyzed directly by functions 
  provided by the coda package. The subject-specific binomial probabilities are stored
  in the object \code{randsave}.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 
  
  \item{ncluster}{an integer giving the number of clusters.} 
  
  \item{p}{ a vector of dimension (no. observations+1) giving the current value of
  the binomial probabilities.} 
  
  \item{ss}{ an interger vector defining to which of the \code{ncluster} clusters each observation belongs.}

  \item{alpha}{ giving the value of the precision parameter.} 
}

\references{

Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
  Using Mixtures. Journal of the American Statistical Association, 90: 577-588.

Liu, J.S. (1996). Nonparametric Hierarchical Bayes via Sequential Imputations. 
  The Annals of Statistics, 24: 911-930.

MacEachern, S.N. (1998) Computational Methods for Mixture of Dirichlet Process Models,
  in Practical Nonparametric and Semiparametric Bayesian Statistics, 
  eds: D. Dey, P. Muller, D. Sinha, New York: Springer-Verlag, pp. 1-22.
}

\examples{
\dontrun{
    # Data
      data(rolling)
      y <- cbind(rolling$y1,rolling$y2)


    # Prior information

      prior<-list(alpha=1,
                  a1=1,
                  b1=1)

    # Initial state
      state <- NULL

    # MCMC parameters

      mcmc <- list(nburn=5000,
                   nsave=10000,
                   nskip=3,
                   ndisplay=100)

    # Fitting the model

      fit <- DPbetabinom(y=y,ngrid=100, 
                         prior=prior, 
                         mcmc=mcmc, 
                         state=state, 
                         status=TRUE)

      fit
      summary(fit)

    # density estimate
      plot(fit,output="density")

    # parameters
      plot(fit,output="param")
}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

Fernando Quintana \email{<quintana@mat.puc.cl>}

}

\keyword{models}
\keyword{nonparametric}
