\name{DPmeta}
\alias{DPmeta}
\alias{DPmeta.default}

\title{Bayesian analysis for a semiparametric linear mixed effects meta-analysis model using a MDP}
\description{
    This function generates a posterior density sample for a 
    semiparametric linear mixed effects meta-analysis model using a Dirichlet process
    or a Mixture of Dirichlet process prior for the distribution of the 
    random effects. 
}
  
  
\usage{

DPmeta(formula,prior,mcmc,state,status,data=sys.frame(sys.parent()),
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
                     prior distribution of the precision parameter of the Dirichlet process
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{a0} and \code{b0} are missing, see details
                     below), \code{tau1} and \code{tau2} giving the hyperparameters for the 
                     prior distribution of the variance of the centering distribution,
                     \code{sigma2} giving the value of the variance
                     of the centering distribution (it must be specified if 
                     \code{tau1} and \code{tau2} are missing),
                     \code{mub} and \code{Sb} giving the hyperparameters 
                     of the normal prior distribution for the mean of the normal
                     baseline distribution, \code{mu} giving the value of the mean of the
                     centering distribution (it must be specified if 
                     \code{mub} and \code{Sb} are missing), and
                     \code{beta0} and \code{Sbeta0} giving the 
                     hyperparameters of the normal prior distribution for the fixed effects
                     (must be specified only if fixed effects are considered in the model). 
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
                     \code{DPmeta} to print an error message and terminate if there are any
                     incomplete observations.}       
}

\details{
  This generic function fits a semiparametric linear mixed effects meta-analysis model: 

  \deqn{y_i \sim N(\theta_i +X_i \beta, \sigma^2_{ei}), i=1,\ldots,n}{yi ~ N(thetai+ Xi beta, sigma2ei), i=1,\ldots,n}
  \deqn{\theta_i | G \sim G}{thetai | G ~ G}
  \deqn{G | \alpha, G_0 \sim DP(\alpha G_0)}{G | alpha, G0 ~ DP(alpha G0)}

  where, \eqn{G_0=N(\theta| \mu, \sigma^2)}{G0=N(theta| mu, sigma2)}.  To complete the model 
  specification, independent hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\beta | \beta_0, S_{\beta_0} \sim N(\beta_0,S_{\beta_0})}{beta | beta0, Sbeta0 ~ N(beta0,Sbeta0)}
  \deqn{\mu | \mu_b, S_b \sim N(\mu_b,S_b)}{mu | mub, Sb ~ N(mub,Sb)}
  \deqn{\sigma^{-2} | \tau_1, \tau_2 \sim Gamma(\tau_1/2,\tau_2/2)}{sigma^-2 | tau1, tau2 ~ Gamma(tau1/2,tau2/2)}

  The precision or total mass parameter, \eqn{\alpha}, of the \code{DP} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. When \eqn{\alpha}{alpha} is random the method described by
  Escobar and West (1995) is used. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value, set \eqn{a_0}{a0} to NULL in the prior specification.
  
  The computational implementation of the model is based on the marginalization of
  the \code{DP} and on the use of MCMC methods for conjugate priors (Escobar, 1994;
  Escobar and West, 1998). 
  
  The average effect is sampled using
  the \eqn{\epsilon}{epsilon}-DP approximation proposed by Muliere and Tardella (1998), with
  \eqn{\epsilon}{epsilon}=0.01.  
}

\value{
  An object of class \code{DPmeta} representing the linear
  mixed-effects model fit. Generic functions such as \code{print}, \code{plot},
  \code{summary}, and \code{anova} have methods to show the results of the fit. 
  The results include \code{beta}, \code{mu}, \code{sigma2}, \code{alpha}, and the 
  number of clusters.
  
  The function \code{DPrandom} can be used to extract the posterior mean of the 
  random effects.

  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 
  
  \item{ncluster}{ an integer giving the number of clusters.} 
  
  \item{alpha}{ giving the value of the precision parameter} 
  
  \item{b}{ a vector of dimension (nsubjects) giving the value of the random effects
   for each subject.} 

  \item{bclus}{ a vector of dimension (nsubjects) giving the value of the random
  effects for each clusters (only the first \code{ncluster}  are considered to start the chain).} 

  \item{ss}{ an interger vector defining to which of the \code{ncluster} clusters each subject belongs.}
  
  \item{beta}{ giving the value of the fixed effects.}
  
  \item{mu}{ giving the mean of the normal baseline distributions.} 
  
  \item{sigma2}{ giving the variance of the normal baseline distributions.}
}

\seealso{
\code{\link{DPrandom}},
\code{\link{DPMmeta}},
\code{\link{DPlmm}}  , \code{\link{DPglmm}} , \code{\link{DPolmm}} ,
\code{\link{DPMlmm}} , \code{\link{DPMglmm}}, \code{\link{DPMolmm}}
}

\references{

Escobar, M.D. (1994) Estimating Normal Means with a Dirichlet Process Prior,
  Journal of the American Statistical Association, 89: 268-277.

Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
  Using Mixtures. Journal of the American Statistical Association, 90: 577-588.

Escobar, M.D. and West, M. (1998) Computing Bayesian Nonparametric Hierarchical Models,
  in Practical Nonparametric and Semiparametric Bayesian Statistics, 
  eds: D. Dey, P. Muller, D. Sinha, New York: Springer-Verlag, pp. 1-22.
  
Muliere, P. and Tardella, L. (1998) Approximating distributions of random functionals
  of Ferguson-Dirichlet priors. The Canadian Journal of Statistics, 26(2): 283-297.
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

    # Prior information

      prior<-list(alpha=1,
                  tau1=20,
                  tau2=10,
                  mub=0,
                  Sb=100)

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

    # Fit the model: First run
    
      fit1<-DPmeta(formula=y~1,prior=prior,mcmc=mcmc,
                   state=state,status=TRUE)
      fit1

    # Summary with HPD and Credibility intervals
      summary(fit1)
      summary(fit1,hpd=FALSE)

    # Plot model parameters (to see the plots gradually set ask=TRUE)
      plot(fit1,ask=FALSE)
      plot(fit1,ask=FALSE,nfigr=2,nfigc=2)	

}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

}

\keyword{models}
\keyword{nonparametric}