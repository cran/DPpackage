\name{DPMolmm}
\alias{DPMolmm}
\alias{DPMolmm.default}

\title{Bayesian analysis for a semiparametric ordinal linear mixed model using a DPM of normals}
\description{
    This function generates a posterior density sample for a 
    semiparametric ordinal linear mixed model using a Dirichlet Process
    Mixture of Normals prior for the distribution of the random effects. 
}
  
  
\usage{

DPMolmm(fixed,random,prior,mcmc,state,status,
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
                     kernel, \code{mb} and \code{Sb} giving the hyperparameters 
                     of the normal prior distribution for the mean of the normal
                     baseline distribution,\code{nub} and \code{Tbinv} giving the hyperparameters of the 
                     inverted Wishart prior distribution for the scale matrix of the normal
                     baseline distribution, and \code{beta0} and \code{Sbeta0} giving the 
                     hyperparameters of the normal prior distribution for the fixed effects
                     (must be specified only if fixed effects are considered in the model).}
                     
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
                     \code{DPMolmm} to print an error message and terminate if there are any
                     incomplete observations.}       
}

\details{
  This generic function fits an ordinal linear mixed-effects model with a probit link (see, e.g.,
  Molenberghs and Verbeke, 2005):

  \deqn{Y_{ij} = k, \mbox{\ if } \gamma_{k-1} \le  W_{ij} < \gamma_{k}, k=1,\ldots,K}{Yij = k, if gammak-1 \le Wij < gammak, k=1,\ldots,K}
  \deqn{W_{ij} \mid \beta_F, \beta_R , b_i  \sim N(X_{ij} \beta_F + Z_{ij} \beta_R + Z_{ij} b_i, 1), i=1,\ldots,N, j=1,\ldots,n_i}{Wij | betaF, betaR , bi ~ N(Xij betaF + Zij betaR + Zij bi, 1), i=1,\ldots,N, j=1,\ldots,ni}
  \deqn{\theta_i | G,\Sigma \sim \int N(m,\Sigma) G(dm)}{thetai | G, Sigma ~ int N(m,Sigma)G(dm)}
  \deqn{G | \alpha, \mu_b,\Sigma_b \sim DP(\alpha N(\mu_b,\Sigma_b))}{G | alpha, mub, Sigmab ~ DP(alpha N(mub, Sigmab))}

  where, \eqn{\theta_i = \beta_R + b_i}{thetai = betaR + bi}, \eqn{\beta = \beta_F}{beta = betaF}, and \eqn{G_0=N(\theta| \mu, \Sigma)}{G0=N(theta| mu, Sigma)}. To complete the model 
  specification, independent hyperpriors are assumed,
  \deqn{\beta | \beta_0, S_{\beta_0} \sim N(\beta_0,S_{\beta_0})}{beta | beta0, Sbeta0 ~ N(beta0,Sbeta0)}
  \deqn{\Sigma | \nu_0, T \sim IW(\nu_0,T)}{Sigma | nu0, T ~ IW(nu0,T)}
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\mu_b | m_b, S_b \sim N(m_b,S_b)}{mub | mb, Sb ~ N(mb,Sb)}
  \deqn{\Sigma_b | \nu_b, Tb \sim IW(\nu_b,Tb)}{Sigma | nub, Tb ~ IW(nub,Tb)}

  A uniform prior is used for the cutoff points. Note that the inverted-Wishart prior is parametrized such that
  \eqn{E(\Sigma)= T^{-1}/(\nu_0-q-1)}{E(Sigma)= T^{-1}/(nu0-q-1)}.

  The precision or total mass parameter, \eqn{\alpha}, of the \code{DP} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. When \eqn{\alpha}{alpha} is random the method described by
  Escobar and West (1995) is used. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value, set \eqn{a_0}{a0} to NULL in the prior specification.

  The computational implementation of the model is based on the marginalization of
  the \code{DP} and on the use of MCMC methods for conjugate priors 
  for a collapsed state of MacEachern (1998). 
  
  The \eqn{\beta_R}{betaR} parameters are sampled using
  the \eqn{\epsilon}{epsilon}-DP approximation proposed by Muliere and Tardella (1998), with
  \eqn{\epsilon}{epsilon}=0.01.  
}

\value{
  An object of class \code{DPMolmm} representing the linear
  mixed-effects model fit. Generic functions such as \code{print}, \code{plot},
  \code{summary}, and \code{anova} have methods to show the results of the fit. 
  The results include \code{betaR}, \code{betaF}, \code{mu}, the elements of 
  \code{Sigma}, \code{mub}, the elements of \code{Sigmab}, 
  the cutoff points, \code{alpha}, and the number of clusters.

  The function \code{DPMrandom} can be used to extract the posterior mean of the 
  random effects.

  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 

  \item{ncluster}{ an integer giving the number of clusters.} 
  
  \item{alpha}{ giving the value of the precision parameter} 
  
  \item{b}{ a matrix of dimension (nsubjects)*(nrandom effects) giving the value of the random effects
  for each subject.} 

  \item{cutoff}{ a real vector defining the cutoff points. Note that the first cutoff must be fixed to
                 0 in this function.}

  \item{mu}{ a matrix of dimension (nsubjects)*(nrandom effects) giving the value of the means
   of the normal kernel for each cluster (only the first \code{ncluster}  are considered to start the chain).} 

  \item{ss}{ an interger vector defining to which of the \code{ncluster} clusters each subject belongs.}
  
  \item{beta}{ giving the value of the fixed effects.}
  
  \item{sigma}{ giving the variance matrix of the normal kernel.}  
  
  \item{mub}{ giving the mean of the normal baseline distributions.} 

  \item{sigmab}{ giving the variance matrix of the normal baseline distributions.}

}

\seealso{
\code{\link{DPMrandom}},
\code{\link{DPMglmm}}, \code{\link{DPMlmm}},
\code{\link{DPlmm}}  , \code{\link{DPglmm}}, \code{\link{DPolmm}},
\code{\link{PTlmm}}  , \code{\link{PTglmm}}, \code{\link{PTolmm}}
}

\references{

Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
  Using Mixtures. Journal of the American Statistical Association, 90: 577-588.

MacEachern, S.N. (1998) Computational Methods for Mixture of Dirichlet Process Models,
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

      prior<-list(alpha=1,
                  tau1=0.01,tau2=0.01,
                  nu0=4.01,
                  tinv=diag(10,1),
                  nub=4.01,
                  tbinv=diag(10,1),
                  mb=rep(0,1),
                  Sb=diag(1000,1),
                  beta0=rep(0,3), 
                  Sbeta0=diag(1000,3))

    # Fitting the model


      fit1 <- DPMolmm(fixed=imps79o~sweek+tx+sweek*tx,random=~1|id,
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
    
      DPMrandom(fit1)
      DPMrandom(fit1,centered=TRUE)


    # Extract predictive information of random effects
    
      pred <- DPMrandom(fit1,predictive=TRUE)
      plot(pred)

}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

}

\keyword{models}
\keyword{nonparametric}