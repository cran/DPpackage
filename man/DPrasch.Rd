\name{DPrasch}
\alias{DPrasch}
\alias{DPrasch.default}

\title{Bayesian analysis for a semiparametric Rasch model}
\description{
    This function generates a posterior density sample for a 
    semiparametric Rasch model, using a DP or a MDP prior
    for the distribution of the random effects. 
}
   
\usage{

DPrasch(y,prior,mcmc,offset,state,status,
        grid=seq(-10,10,length=1000),data=sys.frame(sys.parent()),
        compute.band=FALSE)
     
}

\arguments{
    \item{y}{        a matrix giving the data for which the Rasch Model
                     is to be fitted.}
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Dirichlet process
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{a0} is missing), \code{mub} and \code{Sb} 
                     giving the hyperparameters of the normal prior distribution 
                     for the mean of the normal baseline distribution, \code{mu}
                     giving the mean of the normal baseline distribution
                     (is must be specified if \code{mub} and \code{Sb} are missing),
                     \code{tau1} and \code{tau2} giving the hyperparameters for the 
                     prior distribution of variance of the normal baseline distribution,
                     \code{sigma2} giving the variance of the normal baseline distribution
                     (is must be specified if \code{tau1} and \code{tau2} are missing),
                     and \code{beta0} and \code{Sbeta0} giving the 
                     hyperparameters of the normal prior distribution for the difficulty
                     parameters.}

    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \code{nburn} giving the number of burn-in 
                     scans, \code{nskip} giving the thinning interval, \code{nsave} giving
                     the total number of scans to be saved, and \code{ndisplay} giving
                     the number of saved scans to be displayed on screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out).}   

    \item{offset}{   this can be used to specify an a priori known component to
                     be included in the linear predictor during the fitting.}

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}
    
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state}.}

    \item{grid}{     grid points where the density estimate is 
                     evaluated. The default is seq(-10,10,length=1000).}

    \item{data}{     data frame.}       
    
    \item{compute.band}{  logical variable indicating whether the confidence band for the CDF must be computed.}       
}

\details{
  This generic function fits a semiparametric Rasch model as in
  San Martin et al. (2011), where the linear predictor is modeled as follows:
  \deqn{\eta_{ij} = \theta_i - \beta_j, i=1,\ldots,n, j=1,\ldots,k}{etaij = thetai - betaj, i=1,\ldots,n, j=1,\ldots,k}
  \deqn{\theta_i | G \sim G}{thetai | G ~ G}
  \deqn{G | \alpha, G_0 \sim DP(\alpha G_0)}{G | alpha, G0 ~ DP(alpha G0)}
  \deqn{\beta | \beta_0, S_{\beta_0} \sim N(\beta_0,S_{\beta_0})}{beta | beta0, Sbeta0 ~ N(beta0,Sbeta0)}
  
  where, \eqn{G_0 = N(\theta| \mu, \sigma^2)}{G0 = N(theta| mu, sigma2)}. To
  complete the model specification, independent hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\mu | \mu_b, S_b \sim N(\mu_b,S_b)}{mu | mub, Sb ~ N(mub,Sb)}
  \deqn{\sigma^{-2} | \tau_1, \tau_2 \sim Gamma(\tau_1/2,\tau_2/2)}{sigma2^-1 | tau1, tau2 ~ Gamma(tau1/2,tau2/2)}
  
  The precision or total mass parameter, \eqn{\alpha}{alpha}, of the \code{DP} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. When \eqn{\alpha}{alpha} is random the method described by
  Escobar and West (1995) is used. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value, set \eqn{a_0}{a0} to NULL in the prior specification.

  Each of the parameters of the baseline distribution, \eqn{\mu}{mu} and \eqn{\sigma^2}{sigma2} 
  can be considered as random or fixed at some particular value. In the 
  first case, a Mixture of Dirichlet Process is considered as a prior
  for the distribution of the random effects. To let \eqn{\sigma^2}{sigma2} to be fixed at a particular
  value, set \eqn{\tau_1}{tau1} to NULL in the prior specification. 
  To let \eqn{\mu}{mu} to be fixed at a particular
  value, set \eqn{\mu_b}{mub} to NULL in the prior specification.

  The computational implementation of the model is based on the marginalization of
  the \code{DP} and on the use of MCMC methods for nonconjugate priors 
  (see, MacEachern and Muller, 1998; Neal, 2000).  Specifically, 
  the algorithm 8 with \code{m=1} of Neal (2000), is considered in 
  the \code{DPraschpoisson} function. In this case, the full conditional 
  distributions for the difficulty parameters and in the resampling 
  step of random effects are generated through the Metropolis-Hastings algorithm 
  with a IWLS proposal (see, West, 1985 and Gamerman, 1997).

  The functionals parameters are sampled using
  the \eqn{\epsilon}{epsilon}-DP approximation proposed by Muliere and Tardella (1998), with
  \eqn{\epsilon}{epsilon}=0.01.  
 
}

\value{
  An object of class \code{DPrasch} representing the Rasch
  model fit. Generic functions such as \code{print}, \code{plot},
  and \code{summary} have methods to show the results of the fit. The results include 
  \code{beta}, \code{mu}, \code{sigma2}, the precision parameter 
  \code{alpha}, and the number of clusters.
  
  The function \code{DPrandom} can be used to extract the posterior mean of the 
  random effects.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 
  
  \item{ncluster}{ an integer giving the number of clusters.} 
  
  \item{alpha}{ giving the value of the precision parameter.} 
  
  \item{b}{ a vector of dimension nsubjects giving the value of the random effects
  for each subject.} 

  \item{bclus}{ a vector of dimension nsubjects giving the value of the random
  effects for each clusters (only the first \code{ncluster} are considered to start the chain).} 

  \item{ss}{ an interger vector defining to which of the \code{ncluster} clusters each subject belongs.}
  
  \item{beta}{ giving the value of the difficulty parameters.}
  
  \item{mu}{ giving the mean of the normal baseline distributions.}
  
  \item{sigma2}{ giving the variance of the normal baseline distributions.}
}

\seealso{
\code{\link{DPrandom}}, \code{\link{FPTrasch}}
}

\references{

Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
  Using Mixtures. Journal of the American Statistical Association, 90: 577-588.

Gamerman, D. (1997) Sampling from the posterior distribution in generalized linear
  mixed models. Statistics and Computing, 7: 57-68. 

MacEachern, S. N. and Muller, P. (1998) Estimating mixture of Dirichlet Process
  Models. Journal of Computational and Graphical Statistics, 7 (2): 223-338.

Muliere, P. and Tardella, L. (1998) Approximating distributions of random functionals
  of Ferguson-Dirichlet priors. The Canadian Journal of Statistics, 26(2): 283-297.

Neal, R. M. (2000) Markov Chain sampling methods for Dirichlet process mixture models.
  Journal of Computational and Graphical Statistics, 9:249-265.

San Martin, E., Jara, A., Rolin, J.-M., and Mouchart, M. (2011) On the Bayesian 
  nonparametric generalization of IRT-type models. Psychometrika (To appear)

West, M. (1985) Generalized linear models: outlier accomodation, scale parameter and
  prior distributions. In Bayesian Statistics 2 (eds Bernardo et al.), 531-558, Amsterdam:
  North Holland.

}

\examples{
\dontrun{
    ####################################
    # A simulated Data Set
    ####################################
      nsubject <- 200
      nitem <- 40
      
      y <- matrix(0,nrow=nsubject,ncol=nitem)
      dimnames(y) <- list(paste("id",seq(1:nsubject)), 
                          paste("item",seq(1,nitem)))

      
      ind <- rbinom(nsubject,1,0.5)
      theta <- ind*rnorm(nsubject,1,0.25)+(1-ind)*rnorm(nsubject,3,0.25)
      beta <- c(0,seq(-1,3,length=nitem-1))
      true.cdf <- function(grid)
      {
         0.5*pnorm(grid,1,0.25)+0.5*pnorm(grid,3,0.25) 
      }  
      for(i in 1:nsubject)
      {
         for(j in 1:nitem)
         {
            eta<-theta[i]-beta[j]         
            mean<-exp(eta)/(1+exp(eta))
            y[i,j]<-rbinom(1,1,mean)
         }
      }

    # Prior information

      beta0 <- rep(0,nitem-1)
      Sbeta0 <- diag(1000,nitem-1)

      prior <- list(alpha=1,
                    tau1=6.02,
                    tau2=2.02,
                    mub=0,
                    Sb=100,
                    beta0=beta0,
                    Sbeta0=Sbeta0)

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn <- 5000
      nsave <- 5000
      nskip <- 0
      ndisplay<- 1000
      mcmc <- list(nburn=nburn,
                   nsave=nsave,
                   nskip=nskip,
                   ndisplay=ndisplay)

    # Fit the model
      fit1 <- DPrasch(y=y,prior=prior,mcmc=mcmc, 
                      state=state,status=TRUE,grid=seq(-1,5,0.01),
                      compute.band=TRUE)

    # CDF estimate and truth
      plot(fit1$grid,true.cdf(fit1$grid),type="l",lwd=2,col="red",
           xlab=expression(theta),ylab="CDF")
      lines(fit1$grid,fit1$cdf,lwd=2,col="blue")
      lines(fit1$grid,fit1$cdf.l,lwd=2,col="blue",lty=2)
      lines(fit1$grid,fit1$cdf.u,lwd=2,col="blue",lty=2)

    # Summary with HPD and Credibility intervals
      summary(fit1)
      summary(fit1,hpd=FALSE)

    # Plot model parameters 
    # (to see the plots gradually set ask=TRUE)
      plot(fit1,ask=FALSE)
      plot(fit1,ask=FALSE,nfigr=2,nfigc=2)	

    # Extract random effects
    
      DPrandom(fit1)
      plot(DPrandom(fit1))
      DPcaterpillar(DPrandom(fit1))
}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

}

\keyword{models}
\keyword{nonparametric}
