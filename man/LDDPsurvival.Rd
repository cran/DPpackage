\name{LDDPsurvival}
\alias{LDDPsurvival}
\alias{LDDPsurvival.default}

\title{Bayesian analysis for a Linear Dependent Dirichlet Process Mixture of
 Survival Models}
\description{
    This function generates a posterior density sample for a 
    Linear Dependent Dirichlet Process Mixture of Survival models.
}
   
\usage{

LDDPsurvival(formula,zpred,prior,mcmc,
             state,status,grid,
             data=sys.frame(sys.parent()),
             na.action=na.fail,work.dir=NULL)

}

\arguments{
   \item{formula}{   a two-sided linear formula object describing the
                     model fit, with the response on the
                     left of a \code{~} operator and the terms, separated by \code{+}
                     operators, on the right. The design matrix is used to model 
                     the distribution of the response in the LDPP mixture of normals model.
                     In the response matrix, the unknown limits should be -999.} 

    \item{zpred}{    a matrix giving the covariate values where the predictive 
                     density is evaluated.}
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Dirichlet process
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{a0} is missing), \code{m0} and \code{Sbeta0} 
                     giving the hyperparameters of the normal prior distribution 
                     for the mean of the normal baseline distribution, \code{mub}
                     giving the mean of the normal baseline distribution of the regression
                     coefficients (is must be specified if \code{m0} is missing),
                     \code{nu} and \code{psiinv} giving the hyperparameters of the 
                     inverted Wishart prior distribution for the scale matrix, \code{sigmab}, 
                     of the baseline distribution, \code{sigmab} giving the variance 
                     of the baseline distribution (is must be specified if \code{nu} is missing),
                     \code{tau1} giving the hyperparameter for the 
                     prior distribution of variance of the normal kernel, and
                     \code{taus1} and \code{taus2} giving th hyperparameters of the gamma
                     distribution for \code{tau2}.}

    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \code{nburn} giving the number of burn-in 
                     scans, \code{nskip} giving the thinning interval, \code{nsave} giving
                     the total number of scans to be saved, \code{ndisplay} giving
                     the number of saved scans to be displayed on screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out).}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}
    
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state}.}

    \item{grid}{     real vector giving the grid points where the density, survival and hazard estimates
                     are evaluated.}

    \item{data}{     data frame.}     
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{LDDPsurvival} to print an error message and terminate if there are any
                     incomplete observations.}       

    \item{work.dir}{working directory.} 

}

\details{
  This generic function fits a Linear Dependent Dirichlet Process Mixture of Survival models
  as in De Iorio et al. (2009):
  \deqn{T_i \in (a_i,b_i]}{Ti in (ai,bi)}
  \deqn{log(T_i)=y_i | f_{X_i} \sim f_{X_i}}{log(Ti)=yi | fXi ~ fXi}
  \deqn{f_{X_i} = \int N(X_i \beta, \sigma^2) G(d \beta d \sigma^2)}{fXi = \int N(Xi beta, sigma2) G(d beta d sigma2)}
  \deqn{G | \alpha, G_0 \sim DP(\alpha G_0)}{G | alpha, G0 ~ DP(alpha G0)}
   where, \eqn{G_0 = N(\beta| \mu_b, s_b)\Gamma(\sigma^2|\tau_1/2,\tau_2/2)}{G0 = N(beta| mub, sb)Gamma(sigma2|tau1/2,tau2/2)}. To
   complete the model specification, independent hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\mu_b | m_0, S_0 \sim N(m_0,S_0)}{mub | m0, S0 ~ N(m0,S0)}
  \deqn{s_b | \nu, \Psi \sim IW(\nu,\Psi)}{sb | nu, psi ~ IW(nu,psi)}
  \deqn{\tau_2 | \tau_{s1}, \tau_{s2} \sim Gamma(\tau_{s1}/2,\tau_{s2}/2)}{tau2 ~ Gamma(tau2 | taus1, taus2 ~ Gamma(taus1/2,taus2/2)}
  
  Note that the inverted-Wishart prior is parametrized such that if \eqn{A \sim IW_q(\nu, \psi)}{A ~ IWq(nu, psi)} then
  \eqn{E(A)= \psi^{-1}/(\nu-q-1)}{E(A)= psiinv/(nu-q-1)}.

  Note also that the LDDP model is a natural and simple extension of the the ANOVA DDP model discussed in 
  in De Iorio et al. (2004). The same model is used in Mueller et al.(2005) as the random effects
  distribution in a repeated measurements model.

  The precision or total mass parameter, \eqn{\alpha}{alpha}, of the \code{DP} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. When \eqn{\alpha}{alpha} is random the method described by
  Escobar and West (1995) is used. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value, set \eqn{a_0}{a0} to NULL in the prior specification.

  The computational implementation of the model is based on the marginalization of
  the \code{DP} and on the use of MCMC methods for non-conjugate DPM models (see, e.g,
  MacEachern and Muller, 1998; Neal, 2000). 
}

\value{
  An object of class \code{LDDPsurvival} representing the LDDP mixture of normals model fit. 
  Generic functions such as \code{print}, \code{plot},
  and \code{summary} have methods to show the results of the fit. The results include 
  \code{mub}, \code{sb}, \code{tau2}, the precision parameter 
  \code{alpha}, and the number of clusters.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 

  \item{betaclus}{a matrix of dimension (number of subject + 100) times the 
                  number of columns in the design matrix, giving the 
                  regression coefficients for each cluster (only the first \code{ncluster} are 
                  considered to start the chain).}
                  
  \item{sigmaclus}{a vector of dimension (number of subjects + 100) giving the variance of the normal kernel for
                   each cluster (only the first \code{ncluster} are 
                  considered to start the chain).}

  \item{alpha}{ giving the value of the precision parameter.} 

  \item{mub}{ giving the mean of the normal baseline distributions.}

  \item{sb}{ giving the covariance matrix the normal baseline distributions.}
 
  \item{ncluster}{ an integer giving the number of clusters.} 
  
  \item{ss}{ an interger vector defining to which of the \code{ncluster} clusters each subject belongs.}

  \item{tau2}{giving the value of the tau2 parameter.}

  \item{y}{giving the value of imputed log-survival times.}

}

\seealso{
\code{\link{LDDPdensity}}
}

\references{
De Iorio, M., Muller, P., Rosner, G., and MacEachern, S. (2004) An ANOVA model for dependent  random measures. Journal of the American Statistical Association, 99(465): 205-215.

De Iorio, M., Johnson, W., Muller, P., and Rosner, G.L. (2009) Bayesian Nonparametric Nonproportional
  Hazards Survival Modeling. Biometrics, To Appear.

Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
  Using Mixtures. Journal of the American Statistical Association, 90: 577-588.

MacEachern, S. N. and Muller, P. (1998) Estimating mixture of Dirichlet Process
  Models. Journal of Computational and Graphical Statistics, 7 (2): 223-338.

Mueller, P., Rosner, G., De Iorio, M., and MacEachern, S. (2005). A Nonparametric 
  Bayesian Model for Inference in Related Studies. Applied Statistics, 54 (3), 611-626.

Neal, R. M. (2000). Markov Chain sampling methods for Dirichlet process mixture models.
  Journal of Computational and Graphical Statistics, 9: 249-265.
}

\examples{
\dontrun{

    #############################################################
    # Time to Cosmetic Deterioration of Breast Cancer Patients
    #############################################################

      data(deterioration)
      attach(deterioration)
      ymat <- cbind(left,right)

    # design matrix for posterior predictive 
      zpred <- rbind(c(1,0),c(1,1))  

    # Prior information

      S0 <- diag(100,2)
      m0 <- rep(0,2)
      psiinv <- diag(1,2)
     
      prior <- list(a0=10,
                    b0=1,
                    nu=4,
                    m0=m0,
                    S0=S0,
                    psiinv=psiinv,
                    tau1=6.01,
                    taus1=6.01,
                    taus2=2.01)

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn <- 5000
      nsave <- 5000
      nskip <- 3
      ndisplay <- 100
      mcmc <- list(nburn=nburn,
                   nsave=nsave,
                   nskip=nskip,
                   ndisplay=ndisplay)

    # Fit the model
      fit1 <- LDDPsurvival(ymat~trt,prior=prior,
                           mcmc=mcmc,state=state,status=TRUE,
                           grid=seq(0.01,70,1),zpred=zpred)


    # Plot posterior density estimates
    # with design vector x0=(1,0) 

      plot(fit1$grid,fit1$densp.h[1,],type="l",
           xlab="time",ylab="density",lty=2,lwd=2)
      lines(fit1$grid,fit1$densp.l[1,],lty=2,lwd=2)
      lines(fit1$grid,fit1$densp.m[1,],lty=1,lwd=3)

    # Plot posterior density estimates
    # with design vector x0=(1,1) 

      plot(fit1$grid,fit1$densp.h[2,],type="l",
           xlab="time",ylab="density",lty=2,lwd=2)
      lines(fit1$grid,fit1$densp.l[2,],lty=2,lwd=2)
      lines(fit1$grid,fit1$densp.m[2,],lty=1,lwd=3)

    # Plot posterior survival estimates
    # with design vector x0=(1,0) 

      plot(fit1$grid,fit1$survp.h[1,],type="l",
           xlab="time",ylab="survival",lty=2,lwd=2,ylim=c(0,1))
      lines(fit1$grid,fit1$survp.l[1,],lty=2,lwd=2)
      lines(fit1$grid,fit1$survp.m[1,],lty=1,lwd=3)

    # Plot posterior survival estimates
    # with design vector x0=(1,1) 

      plot(fit1$grid,fit1$survp.h[2,],type="l",
           xlab="time",ylab="survival",lty=2,lwd=2,ylim=c(0,1))
      lines(fit1$grid,fit1$survp.l[2,],lty=2,lwd=2)
      lines(fit1$grid,fit1$survp.m[2,],lty=1,lwd=3)

    # Plot posterior hazard estimates
    # with design vector x0=(1,0) 

      plot(fit1$grid,fit1$hazp.h[1,],type="l",
           xlab="time",ylab="hazard",lty=2,lwd=2)
      lines(fit1$grid,fit1$hazp.l[1,],lty=2,lwd=2)
      lines(fit1$grid,fit1$hazp.m[1,],lty=1,lwd=3)

    # Plot posterior hazard estimates
    # with design vector x0=(1,1) 

      plot(fit1$grid,fit1$hazp.h[2,],type="l",
           xlab="time",ylab="survival",lty=2,lwd=2)
      lines(fit1$grid,fit1$hazp.l[2,],lty=2,lwd=2)
      lines(fit1$grid,fit1$hazp.m[2,],lty=1,lwd=3)

}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

Peter Mueller \email{<pmueller@math.utexas.edu>}

Gary L. Rosner \email{<grosner@jhmi.edu>}

}

\keyword{models}
\keyword{nonparametric}
