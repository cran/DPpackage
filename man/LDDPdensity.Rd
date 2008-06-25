\name{LDDPdensity}
\alias{LDDPdensity}
\alias{LDDPdensity.default}

\title{Bayesian analysis for a Linear Dependent Dirichlet Process Mixture Model}
\description{
    This function generates a posterior density sample for a 
    Linear Dependent Dirichlet Process Mixture of Normals model. 
    Support provided by the NIH/NCI R01CA75981 grant.
}
   
\usage{

LDDPdensity(formula,prior,mcmc,state,status,grid=seq(-10,10,length=1000),
            xpred,data=sys.frame(sys.parent()),na.action=na.fail)

}

\arguments{
   \item{formula}{   a two-sided linear formula object describing the
                     model fit, with the response on the
                     left of a \code{~} operator and the terms, separated by \code{+}
                     operators, on the right. The design matrix is used to model 
                     the distribution of the response in the LDPP mixture of normals model.} 
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Dirichlet process
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{a0} is missing), code{m} and \code{S} 
                     giving the hyperparameters of the normal prior distribution 
                     for the mean of the normal baseline distribution, \code{mub}
                     giving the mean of the normal baseline distribution of the regression
                     coefficients (is must be specified if \code{m} is missing),
                     \code{nu} and \code{psiinv} giving the hyperparameters of the 
                     inverted Wishart prior distribution for the scale matrix, \code{sigmab}, 
                     of the baseline distribution, code{sigmab} giving the variance 
                     of the baseline distribution (is must be specified if \code{nu} is missing),
                     \code{tau1} and \code{tau2} giving the hyperparameters for the 
                     prior distribution of variance of the normal kernel, and
                     code{sigma2k} giving the variance of the normal kernel
                     (is must be specified if \code{tau1} is missing).}

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

    \item{grid}{     grid points where the density estimate is 
                     evaluated. The default is seq(-10,10,length=1000).}

    \item{xpred}{    a matrix giving the covariate values where the predictive 
                     density is evaluated.}

    \item{data}{     data frame.}     
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{LDDPdensity} to print an error message and terminate if there are any
                     incomplete observations.}       
}

\details{
  This generic function fits a Linear Dependent Dirichlet Process Mixture of Normals model,
  \deqn{y_i | G_{X_i} \sim G_{X_i}}{yi | GXi ~ GXi}
  \deqn{G_{X_i} = \int N(X_i \beta, \sigma^2_k) H(d \beta)}{GXi = \int N(Xi beta, sigma2k) H(d beta)}
  \deqn{H | \alpha, H_0 \sim DP(M H_0)}{H | alpha, H0 ~ DP(M H0)}
   where, \eqn{H_0 = N(\beta| \mu_b, \Sigma_b)}{H0 = N(beta| mub, sigmab)}. To
  complete the model specification, independent hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\mu_b | m, S \sim N(m,S)}{mub | m, S ~ N(m,S)}
  \deqn{\Sigma | \nu, \Psi \sim IW(\nu,\Psi)}{sigma | nu, psi ~ IW(nu,psi)}
  \deqn{\sigma^{-2}_k | \tau_1, \tau_2 \sim Gamma(\tau_1/2,\tau_2/2)}{sigma2k^-1 | tau1, tau2 ~ Gamma(tau1/2,tau2/2)}
  
  Note that the inverted-Wishart prior is parametrized such that if \eqn{A \sim IW_q(\nu, \psi)}{A ~ IWq(nu, psi)} then
  \eqn{E(A)= \psi^{-1}/(\nu-q-1)}{E(A)= psiinv/(nu-q-1)}.

  Note also that the LDDP model is a natural ans simple extension of the the ANOVA DDP model discussed in 
  in De Iorio et al. (2004). The same model is used in Mueller et al. (2005) as the random effects
  distribution in a repeated measurements model.

  The precision or total mass parameter, \eqn{\alpha}{alpha}, of the \code{DP} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. When \eqn{\alpha}{alpha} is random the method described by
  Escobar and West (1995) is used. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value, set \eqn{a_0}{a0} to NULL in the prior specification.

  Each of the parameters of the baseline distribution, \eqn{\mu_b}{mub} and \eqn{\Sigma_b}{sigmab} 
  can be considered as random or fixed at some particular value. In the 
  first case, a Mixture of Dirichlet Process is considered as a prior
  for the distribution of the regression coefficients. To let 
  \eqn{\Sigma_b}{sigmab} to be fixed at a particular
  value, set \eqn{\nu}{nu} to NULL in the prior specification. 
  To let \eqn{\mu_b}{mu} to be fixed at a particular
  value, set \eqn{m}{m} to NULL in the prior specification.

  The variance of the normal kernel, \eqn{\sigma^2_k}{sigma2k} 
  can be considered as random or fixed at some particular value. To let 
  \eqn{\sigma^2_k}{sigma2k} to be fixed at a particular
  value, set \eqn{\tau_1}{tau1} to NULL in the prior specification. 

  The computational implementation of the model is based on the marginalization of
  the \code{DP} and on the use of MCMC methods for conjugate priors 
  for a collapsed state of MacEachern (1998). 

}

\value{
  An object of class \code{LDDPdensity} representing the LDDP mixture of normals model fit. 
  Generic functions such as \code{print}, \code{plot},
  and \code{summary} have methods to show the results of the fit. The results include 
  \code{mub}, \code{sigmab}, \code{sigam2k}, the precision parameter 
  \code{alpha}, and the number of clusters.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 

  \item{b}{ a matrix of dimension number of subject times the 
            number of columns in the design matrix, giving the 
            regression coefficients for each cluster (only the first \code{ncluster} are 
            considered to start the chain).} 

  \item{alpha}{ giving the value of the precision parameter.} 

  \item{mub}{ giving the mean of the normal baseline distributions.}

  \item{sigmab}{ giving the covariance matrix the normal baseline distributions.}

  \item{sigma2k}{ giving the variance of the normal kernel.}
  
  \item{ncluster}{ an integer giving the number of clusters.} 
  
  \item{ss}{ an interger vector defining to which of the \code{ncluster} clusters each subject belongs.}
}

\seealso{
\code{\link{DPdensity}}
}

\references{

De Iorio, M., Muller, P., Rosner, G.L., and MacEachern, S (2004) An ANOVA 
  Model for Dependent Random Measures. Journal of the American Statistical 
  Association, 99: 205-215

Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
  Using Mixtures. Journal of the American Statistical Association, 90: 577-588.

MacEachern, S.N. (1998) Computational Methods for Mixture of Dirichlet Process Models,
  in Practical Nonparametric and Semiparametric Bayesian Statistics, 
  eds: D. Dey, P. Muller, D. Sinha, New York: Springer-Verlag, pp. 23-44.

Mueller, P., Rosner, G., De Iorio, M., and MacEachern, S. (2005). A Nonparametric 
  Bayesian Model for Inference in Related Studies. Applied Statistics, 54 (3), 611-626.

}

\examples{
\dontrun{

    ########################################################## 
    # Simulate data from a mixture of two normal densities
    ##########################################################
      nobs <- 500
      y1   <-rnorm(nobs, 3,.8)

      ## y2 = 0.6
      y21 <- rnorm(nobs,1.5, 0.8)
      y22 <- rnorm(nobs,4.0, 0.6)
      u <- runif(nobs)
      y2 <- ifelse(u<0.6,y21,y22)
      y <- c(y1,y2)

      ## design matrix including a single factor
      trt <- c(rep(0,nobs),rep(1,nobs))

      ## design matrix for posterior predictive 
      xpred <- rbind(c(1,0),c(1,1))  

    # Prior information

      s <- diag(100,2)
      m <- rep(0,2)
      psiinv <- diag(1,2)
     
      prior <- list(a0=1,b0=1/5,
                    nu=4,
                    m=m,
                    s=s,
                    psiinv=psiinv,
                    tau1=0.01,
                    tau2=0.01)

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn <- 500
      nsave <- 1000
      nskip <- 10
      ndisplay <- 100
      mcmc <- list(nburn=nburn,
                   nsave=nsave,
                   nskip=nskip,
                   ndisplay=ndisplay)

    # Fit the model
      fit1 <- LDDPdensity(y~trt,prior=prior,mcmc=mcmc,state=state,status=TRUE,
                          grid=seq(-1,7,length=200),xpred=xpred)

    # Get posterior predictive density for future observations
    # with design vector x0=(1,0) 

      plot(fit1$grid,fit1$dens[1,],type="l",xlab="Y",ylab="p",bty="l")
      p1 <- dnorm(fit1$grid, 3.0, 0.8)
      lines(fit1$grid,p1, lty=2)

    # ... and x0=(1,1) 
      plot(fit1$grid,fit1$dens[2,],type="l",xlab="Y",ylab="p",bty="l")
      p2 <- 0.6*dnorm(fit1$grid, 1.5, 0.8) +
            0.4*dnorm(fit1$grid, 4.0, 0.6) 
      lines(fit1$grid,p2, lty=2)

}
}

\author{

Alejandro Jara \email{<Alejandro.JaraVallejos@med.kuleuven.be>}

Peter Mueller \email{<pmueller@mdanderson.org>}

Gary L. Rosner \email{<glrosner@mdanderson.org>}

}

\keyword{models}
\keyword{nonparametric}
