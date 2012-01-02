\name{HDPMcdensity}
\alias{HDPMcdensity}
\alias{HDPMcdensity.default}

\title{Bayesian analysis for a hierarchical Dirichlet Process mixture of normals model for
       conditional density estimation}
\description{
    This function generates a posterior density sample for a 
    DP mixture of normals model for related random probability measures.
    Support provided by the NIH/NCI R01CA75981 grant.
}
   
\usage{
HDPMcdensity(formula,study,xpred,ngrid=100,prior,mcmc,
             state,status,data=sys.frame(sys.parent()),
             na.action=na.fail,work.dir=NULL)
}

\arguments{
  \item{formula}{   a two-sided linear formula object describing the
                     model fit, with the response on the
                     left of a \code{~} operator and the terms, separated by \code{+}
                     operators, on the right.} 

   \item{study}{     a (1 by \code{nrec}) vector of study indicators. The i-th index is the study i 
                     that response j belongs to.}

    \item{ngrid}{    integer giving the number of grid points where the density estimate is 
                     evaluated. The default is 100.}

    \item{xpred}{    a matrix giving the covariate values where the predictive 
                     density is evaluated.}
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameters: \code{pe1} and \code{pe0} giving the prior weights for the point mass at 
                     \eqn{\epsilon=1} and at \eqn{\epsilon=1}{epsilon=1}, respectively, \code{ae} and \code{be} 
                     giving the prior parameters for a Beta prior on \eqn{\epsilon}, \code{eps} giving 
                     the value of \eqn{\epsilon} (it must be specified if \code{pe1} is missing),
                     \code{a0} and \code{b0} vectors giving the hyperparameters for
                     prior distribution of the precision parameter of the Dirichlet process
                     prior, \code{alpha} giving the vector of precision parameters (it 
                     must be specified if \code{a0} is missing), \code{m0} and \code{S0} 
                     giving the hyperparameters of the normal prior distribution 
                     for the mean of the normal baseline distribution, \code{mub} giving the mean 
                     of the normal baseline distribution (is must be specified if \code{m0} is missing),
                     \code{nub} and \code{tbinv} giving the hyperparameters of the 
                     inverse Wishart prior distribution for the variance of the normal
                     baseline distribution, \code{sigmab} giving the variance 
                     of the normal baseline distribution (is must be specified if \code{nub} is missing),
                     \code{nu} and \code{tinv} giving the hyperparameters of the 
                     inverse Wishart prior distribution for the variance of the normal
                     kernel, and \code{sigma} giving the covariance matrix of the normal
                     kernel (is must be specified if \code{nu} is missing).}

    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \code{nburn} giving the number of burn-in 
                     scans, \code{nskip} giving the thinning interval, \code{nsave} giving
                     the total number of scans to be saved, \code{ndisplay} giving
                     the number of saved scans to be displayed on screen.}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis (not
                     available yet).}
    
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state} (not available yet).}

    \item{data}{     data frame.}     
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{HDPdensity} to print an error message and terminate if there are any
                     incomplete observations.}       

     \item{work.dir}{working directory.} 
}

\details{
  This generic function fits a hierarchical mixture of DPM of normals model for conditional density
  estimation (Mueller, Quintana and Rosner, 2004). Let \eqn{d_i=(y_i,x_i)}{di=(yi,xi)} be the vector of
  full data, including the responses \eqn{y_i}{yi} and predictors \eqn{x_i}{xi}. The \code{HDPMcdensity} function
  fits the hierarchical mixture of DPM of normals model to the full data and then look at the implied conditional
  distribution of the responses given the predictors. The model is given by:
  \deqn{d_{ij} | F_i \sim F_i}
  where, \eqn{d_{ij}}{yij} denote the j-th 
  observation in the i-th study, i=1,...,I, and \eqn{F_i} is assumed to arise
  as a mixture \eqn{F_i = \epsilon H_0 + (1-\epsilon) H_i} of one common
  distribution \eqn{H_0} and a distribution \eqn{H_i} that is specific or idiosyncratic to
  the i-th study.

  The random probability measures \eqn{H_i} in turn are given a
  Dirichlet process mixture of normal prior. We assume 
  \deqn{
  H_i(d) = \int N(\mu,\Sigma)  d G_i(\mu),~ i=0,1,\ldots,I
  }
  with
  \deqn{G_i | \alpha_i, G_0 \sim DP(\alpha G_0)}
  where, the baseline distribution is \deqn{G_0 = N(\mu| \mu_b,\Sigma_b)}. 
  
  To complete the model specification, independent hyperpriors are assumed (optional),
  \deqn{\Sigma | \nu, T \sim IW(\nu,T)}
  \deqn{\alpha_i | a_{0i}, b_{0i} \sim Gamma(a_{0i},b_{0i})}
  \deqn{\mu_b | m_0, S_0 \sim N(m_0,S_0)}
  \deqn{\Sigma_b | \nu_b, Tb \sim IW(\nu_b,Tb)}
   and
  \deqn{
  p(\epsilon) = \pi_0 \delta_0+ \pi_1 \delta_1+(1- \pi_0 -\ pi_1)
               Be(a_\epsilon,b_\epsilon)} 

  Note that the inverted-Wishart prior is parametrized such that if \eqn{A \sim IW_q(\nu, \psi)}{A ~ IWq(nu, psi)} then
  \eqn{E(A)= \psi^{-1}/(\nu-q-1)}{E(A)= psiinv/(nu-q-1)}.
}

\value{
  An object of class \code{HDPMcdensity} representing the hierarchical DPM of normals model. 
  Generic functions such as \code{print}, \code{plot},
  and \code{summary} have methods to show the results of the fit. The results include 
  \code{sigma}, \code{eps}, the vector of precision parameters 
  \code{alpha}, \code{mub} and \code{sigmab}.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 

  \item{ncluster}{ an integer giving the number of clusters.} 

  \item{ss}{ an interger vector defining to which of the \code{ncluster} clusters each subject belongs.}

  \item{sc}{ an integer vector defining to which DP each cluster belongs. Note that length(sc)=nrec (only 
             the first \code{ncluster} elements are considered to start the chain.}

  \item{alpha}{ giving the vector of dimension nsuties+1 of precision parameters.} 

  \item{muclus}{a matrix of dimension (number of subject + 100) times the 
                total number of variables (responses + predictors), giving the 
                means for each cluster (only the first \code{ncluster} rows are 
                considered to start the chain).}
                  
  \item{mub}{ giving the mean of the normal baseline distributions.}

  \item{sigmab}{ giving the covariance matrix the normal baseline distributions.}
 
  \item{sigma}{giving the normal kernel covariance matrix.}

  \item{eps}{giving the value of \code{eps}.}
}

\seealso{
\code{\link{predict.HDPMcdensity}}
}

\references{
Mueller, P., Quintana, F. and Rosner, G. (2004). A Method for Combining Inference 
  over Related Nonparametric Bayesian Models. Journal of the Royal Statistical 
  Society, Series B, 66: 735-749.
}

\examples{
\dontrun{
    # Data
      data(calgb)
      attach(calgb)
      y <- cbind(Z1,Z2,Z3,T1,T2,B0,B1)
      x <- cbind(CTX,GM,AMOF)
  
      z <- cbind(y,x)

    #  Data for prediction
      data(calgb.pred)
      xpred <- as.matrix(calgb.pred[,8:10])


    # Prior information
      prior <- list(pe1=0.1,
                    pe0=0.1,
                    ae=1,
                    be=1,
                    a0=rep(1,3),
                    b0=rep(1,3),
                    nu=12,
                    tinv=0.25*var(z),
 		  m0=apply(z,2,mean),
                    S0=var(z),
 		  nub=12,
                    tbinv=var(z))		


    # Initial state
      state <- NULL

    # MCMC parameters

      mcmc <- list(nburn=5000,
                   nsave=5000,
                   nskip=3,
                   ndisplay=100)

    # Fitting the model
      fit1 <- HDPMcdensity(formula=y~x,
                          study=~study,
                          xpred=xpred,
                          prior=prior,
                          mcmc=mcmc,
                          state=state,
                          status=TRUE)

    # Posterior inference
      fit1
      summary(fit1)
       
    # Plot the parameters
    # (to see the plots gradually set ask=TRUE)
      plot(fit1,ask=FALSE)

    # Plot the a specific parameters 
    # (to see the plots gradually set ask=TRUE)
      plot(fit1,ask=FALSE,param="eps",nfigr=1,nfigc=2)

    # Plot the measure for each study 
    # under first values for the predictors, xpred[1,]
      predict(fit1,pred=1,i=1,r=1) # pred1, study 1
      predict(fit1,pred=1,i=2,r=1) # pred1, study 2

    # Plot the measure for each study 
    # under second values for the predictors, xpred[2,]
      predict(fit1,pred=2,i=1,r=1) # pred2, study 1
      predict(fit1,pred=2,i=2,r=1) # pred2, study 2

    # Plot the idiosyncratic measure for each study
    # under first values for the predictors, xpred[1,]
      predict(fit1,pred=1,i=1,r=0) # study 1
      predict(fit1,pred=1,i=2,r=0) # study 2

    # Plot the common measure
    # under first values for the predictors, xpred[1,]
      predict(fit1,pred=1,i=0)
}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

Peter Mueller \email{<pmueller@mdanderson.org>}
}

\keyword{models}
\keyword{nonparametric}
