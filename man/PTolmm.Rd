\name{PTolmm}
\alias{PTolmm}
\alias{PTolmm.default}

\title{Bayesian analysis for a semiparametric ordinal linear mixed model using a MMPT}
\description{
    This function generates a posterior density sample for a 
    semiparametric ordinal linear mixed model, using a Mixture of
    Multivariate Polya Trees prior for the distribution of the 
    random effects. 
}
  
\usage{

PTolmm(fixed,random,prior,mcmc,state,status,data=sys.frame(sys.parent()),
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
                     prior distribution of the precision parameter of the Polya Tree (PT)
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{a0} and \code{b0} are missing, see details
                     below), \code{nu0} and \code{tinv} giving the hyperparameters of the 
                     inverted Wishart prior distribution for the scale matrix of the normal
                     baseline distribution, \code{sigma} giving the value of the covariance
                     matrix of the centering distribution (it must be specified if 
                     \code{nu0} and \code{tinv} are missing),
                     \code{mub} and \code{Sb} giving the hyperparameters 
                     of the normal prior distribution for the mean of the normal
                     baseline distribution, \code{mu} giving the value of the mean of the
                     centering distribution (it must be specified if 
                     \code{mub} and \code{Sb} are missing), \code{beta0} and \code{Sbeta0} 
                     giving the 
                     hyperparameters of the normal prior distribution for the fixed effects
                     (must be specified only if fixed effects are considered in the model), 
                     \code{M} giving the finite level of the PT prior to be considered,
                     \code{frstlprob} a logical variable
                     indicating whether the first level probabilities of the PT are fixed
                     or not (the default is FALSE), and
                     \code{typepr} indicating whether the type of decomposition of the centering 
                     covariance is random (1) or not (0).} 
                     
    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \code{nburn} giving the number of burn-in 
                     scans, \code{nskip} giving the thinning interval, \code{nsave} giving
                     the total number of scans to be saved, \code{ndisplay} giving
                     the number of saved scans to be displayed on screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out), \code{nbase} giving the number scans to be performed before the 
                     parameters of the centering distribution and the precision parameter are
                     updated (i.e., the update of this parameters is invoked only once in every 
                     \code{nbase} scans) (the default value is 1), \code{tune1}, \code{tune2}, and \code{tune3}, 
                     giving the Metropolis tuning parameter for the baseline mean, 
                     variance, and precision parameter, respectively. If \code{tune1}, \code{tune2},
                     or \code{tune3} are not specified or negative, an adpative Metropolis algorithm is performed.
                     Finally, the integer \code{samplef} indicates whether
                     the functional parameters must be sample (1) or not (0).}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}
    
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state}.}

    \item{data}{     data frame.}       
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{PTolmm} to print an error message and terminate if there are any
                     incomplete observations.}       
}

\details{
  This generic function fits an ordinal linear mixed-effects model with a probit link 
  and a Mixture of Multivariate Polya Trees prior (see, Lavine 1992; 1994, for details about univariate PT) 
  for the distribution of the random effects as described in Jara, Hanson and Lessaffre (2009):

  \deqn{Y_{ij} = k, \mbox{\ if } \gamma_{k-1} \leq  W_{ij} < \gamma_{k}, k=1,\ldots,K}{Yij = k, if gammak-1 \leq Wij < gammak, k=1,\ldots,K}
  \deqn{W_{ij} \mid \beta_F, \beta_R , b_i  \sim N(X_{ij} \beta_F + Z_{ij} \beta_R + Z_{ij} b_i, 1), i=1,\ldots,N, j=1,\ldots,n_i}{Wij | betaF, betaR , bi ~ N(Xij betaF + Zij betaR + Zij bi, 1), i=1,\ldots,N, j=1,\ldots,ni}
  \deqn{\theta_i | G \sim G}{thetai | G ~ G}
  \deqn{G | \alpha,\mu,\Sigma,O \sim PT^M(\Pi^{\mu,\Sigma,O},\mathcal{A})}{G | alpha,mu,Sigma,O ~ PT^M(Pi^{mu,Sigma,O},A)}

  where, \eqn{\theta_i = \beta_R + b_i}{thetai = betaR + bi}, \eqn{\beta = \beta_F}{beta = betaF}, and \eqn{O}{O}
  is an orthogonal matrix defining the decomposition of the centering covariance matrix. As in Hanson (2006), the PT prior is centered around a 
  \eqn{N_d(\mu,\Sigma)}{N_d(mu,Sigma)} distribution. However, we consider the class of partitions \eqn{\Pi^{\mu,\Sigma, O}}{Pi^{mu,Sigma, O}}. The
  partitions starts with base sets that are Cartesian products of intervals obtained as quantiles from the standard normal distribution. A multivariate 
  location-scale transformation, \eqn{\theta=\mu+\Sigma^{1/2} z}{theta=mu+Sigma^{1/2} z}, is applied to each 
  base set yielding the final sets. Here \eqn{\Sigma^{1/2}=T'O'}{Sigma^{1/2}=T'O'} where \eqn{T}{T}
  is the unique upper triangular Cholesky matrix of \eqn{\Sigma}{Sigma}. The family \eqn{\mathcal{A}=\{\alpha_e: e \in E^{*}\}}{A=\{alphae: e \in E*\}}, 
  where \eqn{E^{*}=\bigcup_{m=0}^{M} E_d^m}{E*=U_{m=0}^{M} E_d^m}, 
  with \eqn{E_d}{E_d} and \eqn{E_d^m}{E_m} the \eqn{d}{d}-fold product of \eqn{E=\{0,1\}}
  and the the \eqn{m}-fold product of \eqn{E_d}{E_d}, respectively. The family \eqn{\mathcal{A}}{A}
  was specified 
  as \eqn{\alpha_{e_1 \ldots e_m}=\alpha m^2}{alpha{e1 \ldots em}=\alpha m^2}. 

  To complete the model specification, independent hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\beta | \beta_0, S_{\beta_0} \sim N(\beta_0,S_{\beta_0})}{beta | beta0, Sbeta0 ~ N(beta0,Sbeta0)}
  \deqn{\mu | \mu_b, S_b \sim N(\mu_b,S_b)}{mu | mub, Sb ~ N(mub,Sb)}
  \deqn{\Sigma | \nu_0, T \sim IW(\nu_0,T)}{Sigma | nu0, T ~ IW(nu0,T)}
  \deqn{O \sim Haar(q)}{O ~ Haar(q)}

  A uniform prior is used for the cutoff points. Note that the inverted-Wishart prior is parametrized such that
  \eqn{E(\Sigma)= T^{-1}/(\nu_0-q-1)}{E(Sigma)= T^{-1}/(nu0-q-1)}.

  The precision or total mass parameter, \eqn{\alpha}, of the \code{DP} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. 
  
  The computational implementation of the model is based on the marginalization of
  the \code{PT} as descried in Jara, Hanson and Lessaffre (2009). 
}

\value{
  An object of class \code{PTolmm} representing the linear
  mixed-effects model fit. Generic functions such as \code{print}, \code{plot},
  and \code{summary} have methods to show the results of the fit. 
  The results include \code{betaR}, \code{betaF}, \code{mu}, the elements of 
  \code{Sigma}, \code{alpha}, and \code{ortho}.
  
  The function \code{PTrandom} can be used to extract the posterior mean of the 
  random effects.

  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 
  
  \item{alpha}{ giving the value of the precision parameter} 
  
  \item{b}{ a matrix of dimension (nsubjects)*(nrandom effects) giving the value of the random effects
  for each subject.} 

  \item{cutoff}{ a real vector defining the cutoff points. Note that the first cutoff must be fixed at
                 0 in this function.}

  \item{beta}{ giving the value of the fixed effects.}
  
  \item{mu}{ giving the mean of the normal baseline distributions.} 
  
  \item{sigma}{ giving the variance matrix of the normal baseline distributions.}

  \item{ortho}{ giving the orthogonal matrix \code{H}, used in the decomposition of the covariance matrix.} 
}

\seealso{
\code{\link{PTrandom}},
\code{\link{PTlmm}}  , \code{\link{PTglmm}},
\code{\link{DPMglmm}}, \code{\link{DPMlmm}}, \code{\link{DPMolmm}},
\code{\link{DPlmm}}  , \code{\link{DPglmm}}, \code{\link{DPolmm}}
}

\references{

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

    # Schizophrenia Data
      data(psychiatric)
      attach(psychiatric)

    # Prior information
      prior <- list(M=4,
                    frstlprob=FALSE,
                    alpha=1,
                    nu0=4.01,
                    tinv=diag(1,1),
                    mub=rep(0,1),
                    Sb=diag(100,1),
                    beta0=rep(0,3),
                    Sbeta0=diag(1000,3))

    # MCMC parameters
      mcmc <- list(nburn=10000,
                   nsave=10000,
                   nskip=20,
                   ndisplay=100,
                   samplef=1)

    # Initial state
      state <- NULL

    # Fitting the model
      fit1 <- PTolmm(fixed=imps79o~sweek+tx+sweek*tx,random=~1|id,prior=prior,
                     mcmc=mcmc,state=state,status=TRUE)
      fit1

    # Summary with HPD and Credibility intervals
      summary(fit1)
      summary(fit1,hpd=FALSE)

    # Plot model parameters
      plot(fit1)

    # Plot an specific model parameter
      plot(fit1,ask=FALSE,nfigr=1,nfigc=2,param="sigma-(Intercept)")	

    # Extract random effects
      PTrandom(fit1)

    # Extract predictive information of random effects
      aa<-PTrandom(fit1,predictive=TRUE)
      aa

    # Predictive marginal and joint distributions      
      plot(aa)
}
}

\author{
Alejandro Jara \email{<atjara@uc.cl>}

Tim Hanson \email{<hansont@stat.sc.edu>}
}

\keyword{models}
\keyword{nonparametric}
