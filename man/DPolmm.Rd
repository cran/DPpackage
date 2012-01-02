\name{DPolmm}
\alias{DPolmm}
\alias{DPolmm.default}

\title{Bayesian analysis for a semiparametric ordinal linear mixed model using a MDP}
\description{
    This function generates a posterior density sample for a 
    semiparametric ordinal linear mixed model using a Dirichlet Process
    or a Mixture of Dirichlet process prior for the distribution of the random effects.  
}
  
  
\usage{

DPolmm(fixed,random,prior,mcmc,state,status,
       data=sys.frame(sys.parent()),
       na.action=na.fail)
}

\arguments{
    \item{fixed}{    a two-sided linear formula object describing the
                     fixed-effects part of the model, with the response on the
                     left of a \code{~} operator and the terms, separated by \code{+}
                     operators, on the right.}
    
    \item{random}{   a one-sided formula of the form \code{~z1+...+zn | g}, with 
                     \code{z1+...+zn} specifying the model for the random effects and 
                     \code{g} the grouping variable. The random effects formula will be
                     repeated for all levels of grouping.}

    \item{prior}{    a list giving the prior information. The list include the following
                     parameter: \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Dirichlet process
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{a0} and \code{b0} are missing, see details
                     below), \code{nu0} and \code{Tinv} giving the hyperparameters of the 
                     inverted Wishart prior distribution for the scale matrix of the normal
                     baseline distribution, \code{sigma} giving the value of the covariance
                     matrix of the centering distribution (it must be specified if 
                     \code{nu0} and \code{tinv} are missing),
                     \code{mub} and \code{Sb} giving the hyperparameters 
                     of the normal prior distribution for the mean of the normal
                     baseline distribution, \code{mu} giving the value of the mean of the
                     centering distribution (it must be specified if 
                     \code{mub} and \code{Sb} are missing),
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
                     \code{DPolmm} to print an error message and terminate if there are any
                     incomplete observations.}       
}

\details{
  This generic function fits an ordinal linear mixed-effects model with a probit link (see, e.g.,
  Molenberghs and Verbeke, 2005):

  \deqn{Y_{ij} = k, \mbox{\ if } \gamma_{k-1} \le  W_{ij} < \gamma_{k}, k=1,\ldots,K}{Yij = k, if gammak-1 \le Wij < gammak, k=1,\ldots,K}
  \deqn{W_{ij} \mid \beta_F, \beta_R , b_i  \sim N(X_{ij} \beta_F + Z_{ij} \beta_R + Z_{ij} b_i, 1), i=1,\ldots,N, j=1,\ldots,n_i}{Wij | betaF, betaR , bi ~ N(Xij betaF + Zij betaR + Zij bi, 1), i=1,\ldots,N, j=1,\ldots,ni}
  \deqn{\theta_i | G \sim G}{thetai | G ~ G}
  \deqn{G | \alpha, G_0 \sim DP(\alpha G_0)}{G | alpha, G0 ~ DP(alpha G0)}
  
  where, \eqn{\theta_i = \beta_R + b_i}{thetai = betaR + bi}, \eqn{\beta = \beta_F}{beta = betaF}, and \eqn{G_0=N(\theta| \mu, \Sigma)}{G0=N(theta| mu, Sigma)}. To complete the model 
  specification, independent hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\beta | \beta_0, S_{\beta_0} \sim N(\beta_0,S_{\beta_0})}{beta | beta0, Sbeta0 ~ N(beta0,Sbeta0)}
  \deqn{\mu | \mu_b, S_b \sim N(\mu_b,S_b)}{mu | mub, Sb ~ N(mub,Sb)}
  \deqn{\Sigma | \nu_0, T \sim IW(\nu_0,T)}{Sigma | nu0, T ~ IW(nu0,T)}

  A uniform prior is used for the cutoff points. Note that the inverted-Wishart prior is parametrized such that
  \eqn{E(\Sigma)= T^{-1}/(\nu_0-q-1)}{E(Sigma)= T^{-1}/(nu0-q-1)}.

  The precision or total mass parameter, \eqn{\alpha}, of the \code{DP} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. When \eqn{\alpha}{alpha} is random the method described by
  Escobar and West (1995) is used. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value, set \eqn{a_0}{a0} to NULL in the prior specification.
  
  The computational implementation of the model is based on the marginalization of
  the \code{DP} and on the use of MCMC methods for conjugate priors (Escobar, 1994;
  Escobar and West, 1998). The \eqn{\beta_R}{betaR} parameters are sampled using
  the \eqn{\epsilon}{epsilon}-DP approximation proposed by Muliere and Tardella (1998), with
  \eqn{\epsilon}{epsilon}=0.01.  
}

\value{
  An object of class \code{DPolmm} representing the linear
  mixed-effects model fit. Generic functions such as \code{print}, \code{plot},
  \code{summary}, and \code{anova} have methods to show the results of the fit. 
  The results include \code{betaR}, \code{betaF}, \code{mu}, the elements of 
  \code{Sigma}, the cutoff points, \code{alpha}, and the number of clusters.
  
  The function \code{DPrandom} can be used to extract the posterior mean of the 
  random effects.

  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 
  
  \item{ncluster}{ an integer giving the number of clusters.} 
  
  \item{alpha}{ giving the value of the precision parameter} 
  
  \item{b}{ a matrix of dimension (nsubjects)*(nrandom effects) giving the value of the random effects
  for each subject.} 

  \item{bclus}{ a matrix of dimension (nsubjects)*(nrandom effects) giving the value of the random
  effects for each clusters (only the first \code{ncluster}  are considered to start the chain).} 

  \item{cutoff}{ a real vector defining the cutoff points. Note that the first cutoff must be fixed to
                 0 in this function.}

  \item{ss}{ an interger vector defining to which of the \code{ncluster} clusters each subject belongs.}
  
  \item{beta}{ giving the value of the fixed effects.}
  
  \item{mu}{ giving the mean of the normal baseline distributions.} 
  
  \item{sigma}{ giving the variance matrix of the normal baseline distributions.}
}

\seealso{
\code{\link{DPrandom}},
\code{\link{DPlmm}}  , \code{\link{DPglmm}}, 
\code{\link{DPMglmm}}, \code{\link{DPMlmm}}, \code{\link{DPMolmm}},
\code{\link{PTlmm}}  , \code{\link{PTglmm}}, \code{\link{PTolmm}}
}

\references{

Escobar, M.D. (1994) Estimating Normal Means with a Dirichlet Process Prior,
  Journal of the American Statistical Association, 89: 268-277.

Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
  Using Mixtures. Journal of the American Statistical Association, 90: 577-588.

Escobar, M.D. and West, M. (1998) Computing Bayesian Nonparametric Hierarchical Models,
  in Practical Nonparametric and Semiparametric Bayesian Statistics, 
  eds: D. Dey, P. Muller, D. Sinha, New York: Springer-Verlag, pp. 1-22.

Molenberghs, G. and Verbeke, G. (2005). Models for discrete longitudinal data,
  New York: Springer-Verlag.
  
Muliere, P. and Tardella, L. (1998) Approximating distributions of random functionals
  of Ferguson-Dirichlet priors. The Canadian Journal of Statistics, 26(2): 283-297.
}

\examples{
\dontrun{

    # Schizophrenia Data
      data(psychiatric)
      attach(psychiatric)

    # MCMC parameters

      nburn<-5000
      nsave<-10000
      nskip<-10
      ndisplay<-100
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

    # Initial state
      state <- NULL


    # Prior information

      tinv<-diag(10,1)
      prior<-list(alpha=1,nu0=4.01,tinv=tinv,mub=rep(0,1),Sb=diag(100,1),
                  beta0=rep(0,3),Sbeta0=diag(1000,3))


    # Fitting the model


      fit1<-DPolmm(fixed=imps79o~sweek+tx+sweek*tx,random=~1|id,
                   prior=prior,mcmc=mcmc,state=state,status=TRUE)
      fit1

    # Summary with HPD and Credibility intervals
      summary(fit1)
      summary(fit1,hpd=FALSE)

    # Plot model parameters
      plot(fit1)

    # Plot an specific model parameter
      plot(fit1,ask=FALSE,nfigr=1,nfigc=2,param="sigma-(Intercept)")	
      plot(fit1,ask=FALSE,nfigr=1,nfigc=2,param="ncluster")	

    # Extract random effects
      DPrandom(fit1)
      DPrandom(fit1,centered=TRUE)

    # Extract predictive information of random effects
      DPrandom(fit1,predictive=TRUE)
      DPrandom(fit1,centered=TRUE,predictive=TRUE)
      
      plot(DPrandom(fit1,predictive=TRUE))
      plot(DPrandom(fit1,centered=TRUE,predictive=TRUE))
}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

}

\keyword{models}
\keyword{nonparametric}