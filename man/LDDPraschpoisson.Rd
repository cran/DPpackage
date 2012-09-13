\name{LDDPraschpoisson}
\alias{LDDPraschpoisson}
\alias{LDDPraschpoisson.default}

\title{Bayesian analysis for a dependent semiparametric Rasch Poisson model}
\description{
    This function generates a posterior density sample for a 
    semiparametric Rasch Poisson model, using a LDDP mixture of normals prior
    for the distribution of the random effects. 
}
   
\usage{
LDDPraschpoisson(formula,prior,mcmc,
                 offset=NULL,state,status,
                 grid=seq(-10,10,length=1000), 
                 zpred,data=sys.frame(sys.parent()),
                 compute.band=FALSE)
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
                     must be specified if \code{a0} is missing), \code{m0} and \code{S0} 
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
                     distribution for \code{tau2},
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

    \item{zpred}{    a matrix giving the covariate values where the predictive 
                     density is evaluated.}

    \item{data}{     data frame.}       
    
    \item{compute.band}{  logical variable indicating whether the confidence band for the density and CDF must be computed.}       
}

\details{
  This generic function fits a linear dependent semiparametric Rasch Poisson model as in
  Farina et al. (2009), where
  \deqn{\eta_{ij} = \theta_i - \beta_j, i=1,\ldots,n, j=1,\ldots,k}{etaij = thetai - betaj, i=1,...,n, j=1,...,k}
  \deqn{\beta | \beta_0, S_{\beta_0} \sim N(\beta_0,S_{\beta_0})}{\beta | beta0, Sbeta0 ~ N(beta0,Sbeta0)}
  \deqn{\theta_i | f_{X_i} \sim f_{X_i}}{thetai | fXi ~ fXi}
  \deqn{f_{X_i} = \int N(X_i \alpha_c, \sigma^2) G(d \alpha_c d \sigma^2)}{fXi = \int N(Xi alphac, sigma2) G(d alphac d sigma2)}
  \deqn{G | \alpha, G_0 \sim DP(\alpha G_0)}{G | alpha, G0 ~ DP(alpha G0)}
   where, \eqn{G_0 = N(\alpha_c| \mu_b, s_b)\Gamma(\sigma^{-2}|\tau_1/2,\tau_2/2)}{G0 = N(alphac| mub, sb)Gamma(sigma^-2|tau1/2,tau2/2)}. To
   complete the model specification, the following independent hyperpriors are assumed,
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
  An object of class \code{LDDPraschpoisson} representing the LDDP mixture of normals Rasch Poisson model. 
  Generic functions such as \code{print}, \code{plot},
  and \code{summary} have methods to show the results of the fit. The results include 
  \code{beta}, \code{mub}, \code{sb}, \code{tau2}, the precision parameter 
  \code{alpha}, and the number of clusters.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 

  \item{b}{ a vector of dimension nsubjects giving the value of the random effects
  for each subject.} 

  \item{beta}{ giving the value of the difficulty parameters.}

  \item{alphaclus}{a matrix of dimension (number of subject + 100) times the 
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
}

\seealso{
\code{\link{DPrandom}}, \code{\link{DPMraschpoisson}}, \code{\link{DPraschpoisson}}, \code{\link{FPTraschpoisson}}
}

\references{
De Iorio, M., Muller, P., Rosner, G., and MacEachern, S. (2004), An ANOVA model for dependent  random measures," Journal of the American Statistical Association, 99(465): 205-215.

De Iorio, M., Johnson, W., Muller, P., and Rosner, G.L. (2009) Bayesian Nonparametric Nonproportional
  Hazards Survival Modeling. Biometrics, To Appear.

Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
  Using Mixtures. Journal of the American Statistical Association, 90: 577-588.

Farina, P., Quintana, E., San Martin, E., Jara, A. (2009). A Dependent Semiparametric Rasch Model for 
  the Analysis of Chilean Educational Data. In preparation.

MacEachern, S. N. and Muller, P. (1998) Estimating mixture of Dirichlet Process
  Models. Journal of Computational and Graphical Statistics, 7 (2): 223-338.

Mueller, P., Rosner, G., De Iorio, M., and MacEachern, S. (2005). A Nonparametric 
  Bayesian Model for Inference in Related Studies. Applied Statistics, 54 (3), 611-626.

Neal, R. M. (2000). Markov Chain sampling methods for Dirichlet process mixture models.
  Journal of Computational and Graphical Statistics, 9: 249-265.
}

\examples{
\dontrun{
    ####################################
    # A simulated Data Set
    ####################################
     
      grid <- seq(-4,4,0.01)

      dtrue1 <- function(grid)
      {
         0.6*dnorm(grid,-1,0.4)+
         0.3*dnorm(grid,0,0.5)+
         0.1*dnorm(grid,1,0.5)
      }
      
      dtrue2 <- function(grid)
      {
         0.5*dnorm(grid,-1,0.5)+
         0.5*dnorm(grid,1,0.5)
      }

      dtrue3 <- function(grid)
      {
         0.1*dnorm(grid,-1,0.5)+ 
         0.3*dnorm(grid,0,0.5)+
         0.6*dnorm(grid,1,0.4)
      }

      rtrue1 <- function(n)
      {
          ind <- sample(x=c(1,2,3),
                        size=n,replace=TRUE, 
                        prob=c(0.6,0.3,0.1))
          x1 <- rnorm(n,-1,0.4)
          x2 <- rnorm(n, 0,0.5)
          x3 <- rnorm(n, 1,0.5)
          x <- rep(0,n)
          x[ind==1] <- x1[ind==1] 
          x[ind==2] <- x2[ind==2] 
          x[ind==3] <- x3[ind==3]
          return(x)  
      }

      rtrue2 <- function(n)
      {
          ind <- sample(x=c(1,2),
                        size=n,replace=TRUE, 
                        prob=c(0.5,0.5))
          x1 <- rnorm(n,-1,0.5)
          x2 <- rnorm(n, 1,0.5)
          x <- rep(0,n)
          x[ind==1] <- x1[ind==1] 
          x[ind==2] <- x2[ind==2] 
          return(x)  
      }

      rtrue3 <- function(n)
      {
          ind <- sample(x=c(1,2,3),
                        size=n,replace=TRUE, 
                        prob=c(0.1,0.3,0.6))
          x1 <- rnorm(n,-1,0.5)
          x2 <- rnorm(n, 0,0.5)
          x3 <- rnorm(n, 1,0.4)
          x <- rep(0,n)
          x[ind==1] <- x1[ind==1] 
          x[ind==2] <- x2[ind==2] 
          x[ind==3] <- x3[ind==3]
          return(x)  
      }

      b1 <- rtrue1(n=200)
      hist(b1,prob=TRUE,xlim=c(-4,4),ylim=c(0,0.7)) 
      lines(grid,dtrue1(grid))

      b2 <- rtrue2(n=200)
      hist(b2,prob=TRUE,xlim=c(-4,4),ylim=c(0,0.7)) 
      lines(grid,dtrue2(grid))

      b3 <- rtrue3(n=200)
      hist(b3,prob=TRUE,xlim=c(-4,4),ylim=c(0,0.7)) 
      lines(grid,dtrue3(grid))

      nsubject <- 600
      theta <- c(b1,b2,b3)
      trt <- as.factor(c(rep(1,200),rep(2,200),rep(3,200)))
      nitem <- 5
      
      y <- matrix(0,nrow=nsubject,ncol=nitem)
      dimnames(y)<-list(paste("id",seq(1:nsubject)), 
                        paste("item",seq(1,nitem)))

      beta <- c(0,seq(-3,-1,length=nitem-1))

      for(i in 1:nsubject)
      {
         for(j in 1:nitem)
         {
            eta <- theta[i]-beta[j]         
            mm <- exp(eta)
            y[i,j] <- rpois(1,mm)
         }
      }

   ##############################
   # design's prediction matrix
   ##############################

     zpred <- matrix(c(1,0,0,
                       1,1,0,
                       1,0,1),nrow=3,ncol=3,byrow=TRUE)

   ###########################
   # prior
   ###########################

     prior <- list(alpha=1, 
                   beta0=rep(0,nitem-1),
                   Sbeta0=diag(1000,nitem-1),
                   mu0=rep(0,3),
                   S0=diag(100,3),
                   tau1=6.01,
                   taus1=6.01,
                   taus2=2.01,
                   nu=5,
                   psiinv=diag(1,3))

   ###########################
   # mcmc
   ###########################
     mcmc <- list(nburn=5000,
                  nskip=3,
                  ndisplay=100,
                  nsave=5000)

   ###########################
   # fitting the model
   ###########################
 
     fitLDDP <-  LDDPraschpoisson(formula=y ~ trt,
                                  prior=prior,
                                  mcmc=mcmc,
                                  state=NULL,
                                  status=TRUE,
                                  zpred=zpred,
                                  grid=grid,compute.band=TRUE)
  
     fitLDDP

     summary(fitLDDP)

   #########################################
   # plots
   #########################################
     plot(fitLDDP)

     plot(fitLDDP,param="prediction")

   #########################################
   # plot the estimated and true densities
   #########################################

     par(cex=1.5,mar=c(4.1, 4.1, 1, 1))
     plot(fitLDDP$grid,fitLDDP$dens.m[1,],xlim=c(-4,4),ylim=c(0,0.8),
          type="l",lty=1,lwd=3,xlab="Ability",ylab="density",col=1)
     lines(fitLDDP$grid,fitLDDP$dens.u[1,],lty=2,lwd=3,col=1)
     lines(fitLDDP$grid,fitLDDP$dens.l[1,],lty=2,lwd=3,col=1)
     lines(grid,dtrue1(grid),lwd=3,col="red",lty=3)

     par(cex=1.5,mar=c(4.1, 4.1, 1, 1))
     plot(fitLDDP$grid,fitLDDP$dens.m[2,],xlim=c(-4,4),ylim=c(0,0.8),
          type="l",lty=1,lwd=3,xlab="Ability",ylab="density",col=1)
     lines(fitLDDP$grid,fitLDDP$dens.u[2,],lty=2,lwd=3,col=1)
     lines(fitLDDP$grid,fitLDDP$dens.l[2,],lty=2,lwd=3,col=1)
     lines(grid,dtrue2(grid),lwd=3,col="red",lty=3)

     par(cex=1.5,mar=c(4.1, 4.1, 1, 1))
     plot(fitLDDP$grid,fitLDDP$dens.m[3,],xlim=c(-4,4),ylim=c(0,0.8),
          type="l",lty=1,lwd=3,xlab="Ability",ylab="density",col=1)
     lines(fitLDDP$grid,fitLDDP$dens.u[3,],lty=2,lwd=3,col=1)
     lines(fitLDDP$grid,fitLDDP$dens.l[3,],lty=2,lwd=3,col=1)
     lines(grid,dtrue3(grid),lwd=3,col="red",lty=3)

   #########################################
   # Extract random effects
   #########################################
     DPrandom(fitLDDP)
     plot(DPrandom(fitLDDP))
     DPcaterpillar(DPrandom(fitLDPP))

}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

}

\keyword{models}
\keyword{nonparametric}
