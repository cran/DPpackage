\name{DPMglmm}
\alias{DPMglmm}
\alias{DPMglmm.default}

\title{Bayesian analysis for a semiparametric generalized linear mixed model using a DPM of normals}
\description{
    This function generates a posterior density sample for a 
    semiparametric generalized linear mixed model using a Dirichlet Process
    Mixture of Normals prior for the distribution of the random effects. 
}
  
  
\usage{

DPMglmm(fixed,random,family,offset,n,prior,mcmc,state,status,
      data=sys.frame(sys.parent()),na.action=na.fail)
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

    \item{family}{   a description of the error distribution and link function to
                     be used in the model. This can be a character string naming a
                     family function, a family function or the result of a call to
                     a family function. The families(links) considered by 
                     \code{DPglmm} so far are binomial(logit), binomial(probit),
                      Gamma(log), and poisson(log). The gaussian(identity) case is 
                     implemented separately in the function \code{DPlmm}.}
                     
    \item{offset}{   this can be used to specify an a priori known component to
                     be included in the linear predictor during the fitting (only for
                     poisson and gamma models).}

    \item{n}{        this can be used to indicate the total number of cases in a binomial
                     model (only implemented for the logistic link). If it is not specified
                     the response variable must be binary.}

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
                     baseline distribution, \code{beta0} and \code{Sbeta0} giving the 
                     hyperparameters of the normal prior distribution for the fixed effects
                     (must be specified only if fixed effects are considered in the model)
                     and, \code{tau1} and \code{tau2} giving the hyperparameters for the 
                     prior distribution for the inverse of the dispersion parameter of 
                     the Gamma model 
                     (they must be specified only if the Gamma model is considered).}

    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \code{nburn} giving the number of burn-in 
                     scans, \code{nskip} giving the thinning interval, \code{nsave} giving
                     the total number of scans to be saved, \code{ndisplay} giving
                     the number of saved scans to be displayed on the screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out), \code{tune1} giving the positive Metropolis tuning parameter for the 
                     precision parameter of the Gamma model (the default value is 1.1).}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}
    
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state}.}

    \item{data}{     data frame.}       
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{DPMglmm} to print an error message and terminate if there are any
                     incomplete observations.}       
}

\details{
  This generic function fits a generalized linear mixed-effects model, where the linear
  predictor is modeled as follows:
  \deqn{\eta_i = X_i \beta_F + Z_i \beta_R + Z_i b_i, i=1,\ldots,n}{etai = Xi betaF + Zi betaR + Zi bi, i=1,\ldots,n}
  \deqn{\theta_i | G,\Sigma \sim \int N(\mu,\Sigma) G(d\mu)}{thetai | G, Sigma ~ int N(mu,Sigma)G(d mu)}
  \deqn{G | \alpha, \mu_b,\Sigma_b \sim DP(\alpha N(\mu_b,\Sigma_b))}{G | alpha, mub, Sigmab ~ DP(alpha N(mub, Sigmab))}

  where, \eqn{\theta_i = \beta_R + b_i}{thetai = betaR + bi} , \eqn{\beta = \beta_F}{beta = betaF}, and \eqn{G_0 = N(\theta| \mu, \Sigma)}{G0 = N(theta| mu, Sigma)}. To complete the model 
  specification, independent hyperpriors are assumed,
  \deqn{\beta | \beta_0, S_{\beta_0} \sim N(\beta_0,S_{\beta_0})}{beta | beta0, Sbeta0 ~ N(beta0,Sbeta0)}
  \deqn{\Sigma | \nu_0, T \sim IW(\nu_0,T)}{Sigma | nu0, T ~ IW(nu0,T)}
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\mu_b | m_b, S_b \sim N(m_b,S_b)}{mub | mb, Sb ~ N(mb,Sb)}
  \deqn{\Sigma_b | \nu_b, Tb \sim IW(\nu_b,Tb)}{Sigma | nub, Tb ~ IW(nub,Tb)}

  Note that the inverted-Wishart prior is parametrized such that
  \eqn{E(\Sigma)= T^{-1}/(\nu_0-q-1)}{E(Sigma)= Tinv/(nu0-q-1)}.
  
  The precision or total mass parameter, \eqn{\alpha}{alpha}, of the \code{DP} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. When \eqn{\alpha}{alpha} is random the method described by
  Escobar and West (1995) is used. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value set, \eqn{a_0}{a0} to NULL in the prior specification.

  The inverse of the dispersion parameter of the Gamma model
  is modeled using \code{gamma} distribution, \eqn{\Gamma(\tau_1/2,\tau_2/2)}{Gamma(tau1/2,tau2/2)}. 
  
  The computational implementation of the model is based on the marginalization of
  the \code{DP} and the MCMC is model-specific. 
  
  For the \code{poisson}, \code{Gamma}, and \code{binomial(logit)}, MCMC methods for nonconjugate 
  priors (see, MacEachern and Muller, 1998; Neal, 2000) are used. Specifically, the algorithm 8 
  with \code{m=1} of Neal (2000), is considered in 
  the \code{DPMglmm} function. In this case, the fully conditional distributions for fixed and 
  in the resampling step of random effects are generated through the Metropolis-Hastings algorithm 
  with a IWLS proposal (see, West, 1985 and Gamerman, 1997).
  
  For the \code{binomial(probit)} the following latent variable representation is used:
  \deqn{y_{ij} = I(w_{ij} > 0), j=1,\ldots,n_i}{yij = I(wij > 0), j=1,\ldots,ni}
  \deqn{w_{ij} | \beta, \theta_i, \lambda_i \sim N(X_{ij} \beta + Z_{ij} \theta_i, 1)}{wij | beta, thetai, lambdai ~ N(Xij beta + Zij thetai, 1)}

  In this case, the computational implementation of the model is based on the marginalization of
  the \code{DP} and on the use of MCMC methods for conjugate priors 
  for a collapsed state described by MacEachern (1998).
  
  The \eqn{\beta_R}{betaR} parameters are sampled using
  the \eqn{\epsilon}{epsilon}-DP approximation proposed by Muliere and Tardella (1998), with
  \eqn{\epsilon}{epsilon}=0.01.  
}

\value{
  An object of class \code{DPMglmm} representing the generalized linear
  mixed-effects model fit. Generic functions such as \code{print}, \code{plot},
  \code{summary}, and \code{anova} have methods to show the results of the fit. 
  The results include \code{betaR}, \code{betaF}, \code{sigma2e}, 
  \code{Sigma}, \code{mub}, the elements of \code{Sigmab}, \code{alpha}, and the 
  number of clusters.
  
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

  \item{mu}{ a matrix of dimension (nsubjects)*(nrandom effects) giving the value of the means
   of the normal kernel for each cluster (only the first \code{ncluster}  are considered to start the chain).} 

  \item{ss}{ an interger vector defining to which of the \code{ncluster} clusters each subject belongs.}
  
  \item{beta}{ giving the value of the fixed effects.}
  
  \item{sigma}{ giving the variance matrix of the normal kernel.}  
  
  \item{mub}{ giving the mean of the normal baseline distributions.} 

  \item{sigmab}{ giving the variance matrix of the normal baseline distributions.}
 
  \item{phi}{ giving the dispersion parameter for the Gamma model (if needed).}

}

\seealso{
\code{\link{DPMrandom}},
\code{\link{DPMlmm}}   , \code{\link{DPMolmm}},
\code{\link{DPlmm}}    , \code{\link{DPglmm}} , \code{\link{DPolmm}},
\code{\link{PTlmm}}  , \code{\link{PTglmm}}, \code{\link{PTolmm}}
}

\references{

Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
  Using Mixtures. Journal of the American Statistical Association, 90: 577-588.

Gamerman, D. (1997) Sampling from the posterior distribution in generalized linear
  mixed models. Statistics and Computing, 7: 57-68. 

MacEachern, S.N. (1998) Computational Methods for Mixture of Dirichlet Process Models,
  in Practical Nonparametric and Semiparametric Bayesian Statistics, 
  eds: D. Dey, P. Muller, D. Sinha, New York: Springer-Verlag, pp. 1-22.

MacEachern, S. N. and Muller, P. (1998) Estimating mixture of Dirichlet Process
  Models. Journal of Computational and Graphical Statistics, 7 (2): 223-338.

Muliere, P. and Tardella, L. (1998) Approximating distributions of random functionals
  of Ferguson-Dirichlet priors. The Canadian Journal of Statistics, 26(2): 283-297.

Neal, R. M. (2000) Markov Chain sampling methods for Dirichlet process mixture models.
  Journal of Computational and Graphical Statistics, 9:249-265.

West, M. (1985) Generalized linear models: outlier accomodation, scale parameter and
  prior distributions. In Bayesian Statistics 2 (eds Bernardo et al.), 531-558, Amsterdam:
  North Holland.

}

\examples{
\dontrun{
    # Respiratory Data Example

      data(indon)
      attach(indon)

      baseage2<-baseage**2
      follow<-age-baseage
      follow2<-follow**2 

    # Prior information

      prior<-list(alpha=1,
                  nu0=4.01,
                  tinv=diag(1,1),
                  nub=4.01,
                  tbinv=diag(1,1),
                  mb=rep(0,1),
                  Sb=diag(1000,1),
                  beta0=rep(0,9),
                  Sbeta0=diag(1000,9))

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn<-5000
      nsave<-25000
      nskip<-20
      ndisplay<-1000
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

    # Fit the Probit model
      fit1<-DPMglmm(fixed=infect~gender+height+cosv+sinv+xero+baseage+
                    baseage2+follow+follow2,
                    random=~1|id,family=binomial(probit),
                    prior=prior,mcmc=mcmc,state=state,status=TRUE)

    # Fit the Logit model
      fit2<-DPMglmm(fixed=infect~gender+height+cosv+sinv+xero+baseage+
                    baseage2+follow+follow2,random=~1|id,
                    family=binomial(logit),
                    prior=prior,mcmc=mcmc,state=state,status=TRUE)

    # Summary with HPD and Credibility intervals
      summary(fit1)
      summary(fit1,hpd=FALSE)

      summary(fit2)
      summary(fit2,hpd=FALSE)


    # Plot model parameters (to see the plots gradually set ask=TRUE)
      plot(fit1,ask=FALSE)
      plot(fit1,ask=FALSE,nfigr=2,nfigc=2)	

    # Plot an specific model parameter (to see the plots gradually set ask=TRUE)
      plot(fit1,ask=FALSE,nfigr=1,nfigc=2,param="baseage")	
      plot(fit1,ask=FALSE,nfigr=1,nfigc=2,param="ncluster")	
}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

}

\keyword{models}
\keyword{nonparametric}