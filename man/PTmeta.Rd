\name{PTmeta}
\alias{PTmeta}
\alias{PTmeta.default}

\title{Bayesian analysis for a semiparametric linear mixed effects meta-analysis model using a MPT}
\description{
    This function generates a posterior density sample for a 
    semiparametric linear mixed effects meta-analysis model using a Polya Tree
    or a Mixture of Polya Trees prior for the distribution of the 
    random effects. 
}
  
  
\usage{

PTmeta(formula,prior,mcmc,state,status,data=sys.frame(sys.parent()),
       na.action=na.fail)

}

\arguments{
    \item{formula}{  a two-sided linear formula object describing the
                     fixed-effects part of the model, with the response on the
                     left of a \code{~} operator and the terms, separated by \code{+}
                     operators, on the right. Both effect and variance must be included 
                     in the LHS of the formula object}
    
    \item{prior}{    a list giving the prior information. The list include the following
                     parameter: \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Polya tree
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{a0} and \code{b0} are missing, see details
                     below), \code{tau1} and \code{tau2} giving the hyperparameters for the 
                     prior distribution of the variance of the centering distribution,
                     \code{sigma} giving the value of the variance
                     of the centering distribution (it must be specified if 
                     \code{tau1} and \code{tau2} are missing),
                     \code{mub} and \code{Sb} giving the hyperparameters 
                     of the normal prior distribution for the mean of the normal
                     baseline distribution, \code{mu} giving the value of the mean of the
                     centering distribution (it must be specified if 
                     \code{mub} and \code{Sb} are missing), and
                     \code{beta0} and \code{Sbeta0} giving the 
                     hyperparameters of the normal prior distribution for the fixed effects
                     (must be specified only if fixed effects are considered in the model),
                     \code{M} giving the finite level
                     of the PT prior to be considered, and
                     \code{frstlprob} a logical variable
                     indicating whether the first level probabilities of the PT are fixed
                     or not (the default is FALSE) (see, details).
                     }

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
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{PTmeta} to print an error message and terminate if there are any
                     incomplete observations.}       
}

\details{
  This generic function fits a semiparametric linear mixed effects meta-analysis model using
  a Polya tree prior on the distribution (see, Lavine (1992; 1994) and Hanson (2006) for details about PT)
  on the distribution of the random effects:

  \deqn{y_i \sim N(\theta_i +X_i \beta, \sigma^2_{ei}), i=1,\ldots,n}{yi ~ N(thetai+ Xi beta, sigma2ei), i=1,\ldots,n}
  \deqn{\theta_i | G \sim G}{thetai | G ~ G}
  \deqn{G | \alpha,\mu,\sigma \sim PT(\Pi^{\mu,\sigma},\textit{A})}{G | alpha,mu,sigma ~ PT(Pi^{mu,sigma^2},\textit{A})}

  where the PT prior is centered around a \eqn{N(\mu,\sigma^2)}{N(mu,sigma2)} distribution. 
  If \code{frstlprob} is equal to TRUE, \eqn{\mu=0}{mu=0} and a median zero PT prior 
  is considered (see, Branscum and Hanson, 2008).
 
  To complete the model specification, independent hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\beta | \beta_0, S_{\beta_0} \sim N(\beta_0,S_{\beta_0})}{beta | beta0, Sbeta0 ~ N(beta0,Sbeta0)}
  \deqn{\mu | \mu_b, S_b \sim N(\mu_b,S_b)}{mu | mub, Sb ~ N(mub,Sb)}
  \deqn{\sigma^{-2} | \tau_1, \tau_2 \sim Gamma(\tau_1/2,\tau_2/2)}{sigma^-2 | tau1, tau2 ~ Gamma(tau1/2,tau2/2)}

  The precision parameter, \eqn{\alpha}, of the \code{PT} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. 
  
  The computational implementation of the model is based on the marginalization of
  the \code{PT} and on the MCMC algorihtms described in Hanson (2006) and 
  Jara, Hanson and Lesaffre (2009).
  
  The average effect is sampled using the method of composition described in Jara, Hanson and Lesaffre (2009).  
}

\value{
  An object of class \code{PTmeta} representing the linear
  mixed-effects model fit. Generic functions such as \code{print}, \code{plot},
  \code{summary}, and \code{anova} have methods to show the results of the fit. 
  The results include \code{beta}, \code{mu}, \code{sigma2}, and \code{alpha}.
  
  The function \code{PTrandom} can be used to extract the posterior mean of the 
  random effects.

  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 
  
  \item{alpha}{ giving the value of the precision parameter} 
  
  \item{b}{ a vector of dimension (nsubjects) giving the value of the random effects
   for each subject.} 

  \item{beta}{ giving the value of the fixed effects.}
  
  \item{mu}{ giving the mean of the normal baseline distributions.} 
  
  \item{sigma2}{ giving the variance of the normal baseline distributions.}
}

\seealso{
\code{\link{PTrandom}},
\code{\link{DPMmeta}} , \code{\link{DPMmeta}} ,
\code{\link{DPlmm}}   , \code{\link{DPglmm}}  , \code{\link{DPolmm}} ,
\code{\link{DPMlmm}}  , \code{\link{DPMglmm}} , \code{\link{DPMolmm}}
}

\references{
Branscum, A. and Hanson, T. (2008) Bayesian nonparametric meta-analysis 
   using Polya tree mixture models. Biometrics, 64: 825-833.

Christensen, R., Hanson, T.  Jara, A.. 2008. Parametric Nonparametric Statistics: An 
  Introduction to Mixtures of Finite Polya Trees Models. The American 
  Statistician, 62: 296-306.

Hanson, T. (2006) Inference for Mixtures of Finite Polya Trees.
  Journal of the American Statistical Association, 101: 1548-1565.

Jara, A., Hanson, T., Lesaffre, E. (2009) Robustifying Generalized Linear 
   Mixed Models using a New Class of Mixtures of Multivariate Polya Trees. 
   Journal of Computational and Graphical Statistics, 18(4): 838-860.

Lavine, M. (1992) Some aspects of Polya tree distributions for statistical 
  modelling. The Annals of Statistics, 20: 1222-11235.

Lavine, M. (1994) More aspects of Polya tree distributions for statistical 
  modelling. The Annals of Statistics, 22: 1161-1176.
}

\examples{
\dontrun{

    ##################################################################    
    # Data on the effectiveness of silver sulfadiazine coating
    # on venous catheters for preventing bacterial colonisation of 
    # the catheter and bloodstream infection. 
    # Veenstra D et al (1998) "Efficacy of Antiseptic Impregnated 
    # Central Venous Catheters in Preventing Nosocomial Infections: 
    # A Meta-analysis" JAMA 281:261-267. 
    #
    # Note that -Inf and Inf have been replaced by NA.
    ##################################################################    
    
      studies <- c("Tennenberg","Maki","vanHeerden",
                    "Hannan","Bach(a)","Bach(b)",
                    "Heard","Collins","Ciresi","Ramsay",
                    "Trazzera","George")    

      logOR <- c(-1.5187189,-0.7136877,-1.3217558,-0.1910552,
                  NA,-2.2005195,-0.5057461,-2.3538784,-0.3643810,
                  -0.5371429,-0.7608058,-2.1400662)
       
      varlogOR <- c(0.4157541,0.2632550,0.6739189,0.3727788,NA,
                    0.7623470,0.2306169,0.7477891,0.3645463,0.2291839,
                    0.3561542,0.5190489)^2

      names(logOR) <- studies
      names(varlogOR) <- studies
      y <- cbind(logOR,varlogOR)
      colnames(y) <- c("logOR","varlogOR")


    # Initial state
      state <- NULL

    # MCMC parameters

      nburn<-20000
      nsave<-10000
      nskip<-20
      ndisplay<-100
      mcmc <- list(nburn=nburn,
                   nsave=nsave,
                   nskip=nskip,
                   ndisplay=ndisplay)

    # Prior information 1: non-median zero PT

      prior1<-list(alpha=1,
                   tau1=20,
                   tau2=10,
                   mub=0,
                   Sb=100,
                   M=4)

    # Prior information 2: median zero PT

      prior2<-list(alpha=1,
                   tau1=20,
                   tau2=10,
                   mub=0,
                   Sb=100,
                   M=4,
                   frstlprob=TRUE,
                   Sbeta0=diag(1000,1),
                   beta0=rep(0,1))
                   

    # Fitting the models

      fit1<-PTmeta(formula=y~1,prior=prior1,mcmc=mcmc,
                   state=state,status=TRUE)
      fit1


      fit2<-PTmeta(formula=y~1,prior=prior2,mcmc=mcmc,
                   state=state,status=TRUE)
      fit2

    # Summary with HPD and Credibility intervals
      summary(fit1)
      summary(fit1,hpd=FALSE)

      summary(fit2)
      summary(fit2,hpd=FALSE)

    # Plot model parameters (to see the plots gradually set ask=TRUE)
      plot(fit1,ask=FALSE)
      plot(fit1,ask=FALSE,nfigr=2,nfigc=2)	

      plot(fit2,ask=FALSE)
      plot(fit2,ask=FALSE,nfigr=2,nfigc=2)	
}
}

\author{
Alejandro Jara \email{<atjara@uc.cl>}
}

\keyword{models}
\keyword{nonparametric}
