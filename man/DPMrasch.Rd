\name{DPMrasch}
\alias{DPMrasch}
\alias{DPMrasch.default}

\title{Bayesian analysis for a semiparametric Rasch model}
\description{
    This function generates a posterior density sample for a 
    semiparametric Rasch model, using a DPM of normals prior
    for the distribution of the random effects. 
}
   
\usage{

DPMrasch(y,prior,mcmc,offset=NULL,state,status,
        grid=seq(-10,10,length=1000),data=sys.frame(sys.parent()),
        compute.band=FALSE)
     
}

\arguments{
    \item{y}{        a matrix giving the data for which the Rasch Model
                     is to be fitted.}
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{N} giving the truncation of the Dirichlet 
                     process prior, 
                     \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Dirichlet process
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{a0} is missing), 
                     \code{m0} and \code{s0} 
                     giving the hyperparameters of the normal prior distribution 
                     for the mean, \code{mub}, of the normal baseline distribution, 
                     \code{mub} giving the mean of the baseline distribution (it must be
                     specified if \code{s0} is missing),
                     \code{taub1} and \code{taub2} giving the hyperparameters of the 
                     inverted gamma prior distribution for the variance, \code{sigmab}, 
                     of the baseline distribution, \code{sigmab} giving the variance 
                     of the baseline distribution (is must be specified if \code{taub1} is missing),
                     \code{tauk1} giving the hyperparameter for the 
                     prior distribution of variance of the normal kernel, and
                     \code{taus1} and \code{taus2} giving th hyperparameters of the gamma
                     distribution for \code{tauk2},
                     \code{beta0} and \code{Sbeta0} giving the 
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
    
    \item{compute.band}{  logical variable indicating whether the confidence band for the density and CDF must be computed.}       
}

\details{
  This generic function fits a semiparametric Rasch model as in
  San Martin et al. (2011), where
  \deqn{\eta_{ij} = \theta_i - \beta_j, i=1,\ldots,n, j=1,\ldots,k}{etaij = thetai - betaj, i=1,\ldots,n, j=1,\ldots,k}
  \deqn{\theta_i | G \sim \int N(\mu,\sigma) G(d\mu,\sigma)}{thetai | G ~ int N(mu,sigma)G(d mu, d sigma)}
  \deqn{\beta | \beta_0, S_{\beta_0} \sim N(\beta_0,S_{\beta_0})}{beta | beta0, Sbeta0 ~ N(beta0,Sbeta0)}
  \deqn{G | \alpha, G_0 \sim DP(\alpha G_0)}{G | alpha, G0 ~ DP(alpha G0)}
  where \eqn{G_0 = N(\mu |\mu_b, \sigma_b) IG(\sigma|\tau_{k1},\tau_{k_2})}{G0 = N(mu| mub, sigmab)IG(sigma|tauk1,tauk2)}. 
  To complete the model specification, independent hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\mu_b | m0, s0 \sim N(m0,s0)}{mub | m0, s0 ~ N(m0,s0)}
  \deqn{\sigma_b^{-2} | \tau_{b1}, \tau_{b2} \sim Gamma(\tau_{b1}/2,\tau_{b2}/2)}{sigma_b^-2 | taub1, taub2 ~ Gamma(taub1/2,taub2/2)}
  \deqn{\tau_{k2} | \tau_{s1}, \tau_{s2} \sim Gamma(\tau_{s1}/2,\tau_{s2}/2)}{tauk2 | taus1, taus2 ~ Gamma(taus1/2,taus2/2)}
  
  The precision or total mass parameter, \eqn{\alpha}{alpha}, of the \code{DP} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value, set \eqn{a_0}{a0} to NULL in the prior specification.

  The computational implementation of the model is based on the finite approximation for DP 
  proposed by Ishwaran and James (2002). The full conditional 
  distributions for the difficulty parameters and in the resampling 
  step of random effects are generated through the Metropolis-Hastings algorithm 
  with a IWLS proposal (see, West, 1985 and Gamerman, 1997).
}

\value{
  An object of class \code{DPMrasch} representing the Rasch
  model fit. Generic functions such as \code{print}, \code{plot},
  and \code{summary} have methods to show the results of the fit. The results include 
  \code{beta}, \code{mub}, \code{sigmab}, \code{sigmak2}, the precision parameter 
  \code{alpha}, and the number of clusters.
  
  The function \code{DPrandom} can be used to extract the posterior mean of the 
  random effects.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 

  \item{alpha}{ giving the value of the precision parameter.} 

  \item{b}{ a vector of dimension nsubjects giving the value of the random effects
  for each subject.} 

  \item{beta}{ giving the value of the difficulty parameters.}
  
  \item{ncluster}{ an integer giving the number of clusters.} 
  
  \item{ss}{ an interger vector defining to which of the \code{ncluster} clusters each subject belongs.}
  
  \item{muclus}{ a vector of dimension N giving the value of the normal means.} 

  \item{sigmaclus}{ a vector of dimension N giving the value of the normal variances.} 
  
  \item{mub}{ giving the mean of the normal baseline distributions.}
  
  \item{sigmab}{ giving the variance of the normal baseline distributions.}

  \item{tauk2}{ giving the parameter of the inverse-gamma prior for the normal kernel variances.}

  \item{wdp}{ giving the vector of DP weights.}

  \item{vdp}{ giving the vector of stick-breaking beta random variables used to create the DP weights.}
}

\seealso{
\code{\link{DPrandom}}, \code{\link{DPrasch}}, \code{\link{FPTrasch}}
}

\references{
Gamerman, D. (1997) Sampling from the posterior distribution in generalized linear
  mixed models. Statistics and Computing, 7: 57-68. 

Ishwaran, H. and James, L.F. (2002) Approximate Dirichlet process computing finite normal mixtures:  smoothing and prior information. Journal of Computational and Graphical Statistics, 11: 508-532.

San Martin, E., Jara, A., Rolin, J.-M., and Mouchart, M. (2011) On the Bayesian nonparametric 
  generalization of IRT-type models. Psychometrika (To appear).

West, M. (1985) Generalized linear models: outlier accomodation, scale parameter and
  prior distributions. In Bayesian Statistics 2 (eds Bernardo et al.), 531-558, Amsterdam:
  North Holland.

}

\examples{
\dontrun{
    ####################################
    # A simulated Data Set
    ####################################
      nsubject <- 250
      nitem <- 40
      
      y <- matrix(0,nrow=nsubject,ncol=nitem)
      dimnames(y)<-list(paste("id",seq(1:nsubject)), 
                        paste("item",seq(1,nitem)))

      ind <- rbinom(nsubject,1,0.5)
      theta <- ind*rnorm(nsubject,-1,sqrt(0.25))+
               (1-ind)*rnorm(nsubject,2,sqrt(0.065))
      beta <- c(0,seq(-3,3,length=nitem-1))

      true.density <- function(grid)
      {
            0.5*dnorm(grid,-1,sqrt(0.25))+0.5*dnorm(grid,2,sqrt(0.065))  
      } 

      true.cdf <- function(grid)
      {
            0.5*pnorm(grid,-1,sqrt(0.25))+0.5*pnorm(grid,2,sqrt(0.065))  
      } 

      for(i in 1:nsubject)
      {
         for(j in 1:nitem)
         {
            eta <- theta[i]-beta[j]         
            prob <- exp(eta)/(1+exp(eta))
            y[i,j] <- rbinom(1,1,prob)
         }
      }

    # Prior information

      beta0 <- rep(0,nitem-1)
      Sbeta0 <- diag(100,nitem-1)

      prior <- list(N=50,
                          alpha=1,
                          taub1=6.01,
                          taub2=2.01,
                          taus1=6.01,
                          taus2=2.01,
                          tauk1=6.01,
                          m0=0,
                          s0=100,
                          beta0=beta0,
                          Sbeta0=Sbeta0)

    # Initial state
      state <- NULL      

    # MCMC parameters

      nburn <- 4000
      nsave <- 4000
      nskip <- 0
      ndisplay <- 100
      mcmc <- list(nburn=nburn,
                            nsave=nsave,
                            nskip=nskip,
                            ndisplay=ndisplay)

    # Fit the model
      fit1 <- DPMrasch(y=y,prior=prior,mcmc=mcmc,
                       state=state,status=TRUE,grid=seq(-3,4,0.01))
   
      plot(fit1$grid,fit1$dens.m,type="l",lty=1,col="red",
           xlim=c(-3,4),ylim=c(0,0.8))
      lines(fit1$grid,true.density(fit1$grid),
            lty=2,col="blue")

      plot(fit1$grid,fit1$cdf.m,type="l",lty=1,col="red")
      lines(fit1$grid,true.cdf(fit1$grid),lty=2,col="blue")

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
      DPcaterpillar(DPrandom(fit1))}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

}

\keyword{models}
\keyword{nonparametric}
