\name{LDTFPglmm}
\alias{LDTFPglmm}
\alias{LDTFPglmm.default}

\title{Generalized linear mixed model using a linear dependent tailfree prior}

\description{
    This function generates a posterior density sample for a generalized
    linear mixed model using a linear dependent tail free prior for the
    random intercept distribution.
}
   
\usage{

LDTFPglmm(y,x,roffset=NULL,g,family,
          xtf,prediction,prior,mcmc,
          state,status,ngrid=100,
          grid=NULL,compute.band=FALSE,
          type.band="PD",
          data=sys.frame(sys.parent()),
          na.action=na.fail,
          work.dir=NULL)
}

\arguments{
   \item{y}{         a vector giving the response variables.}

   \item{x}{         a matrix giving the design matrix for the fixed effects. This matrix
                     must include the constant term.} 

   \item{roffset}{   this can be used to specify an a priori known component to
                     be included in the linear predictor during the fitting (only for
                     poisson and gamma models).}
   
   \item{g}{        a vector giving the group indicator for each observation.} 

   \item{family}{   a description of the error distribution and link function to
                    be used in the model. This can be a character string naming a
                    family function, a family function or the result of a call to
                    a family function. The families(links) considered by 
                    \code{LDTFPglmm} so far are binomial(logit), binomial(probit),
                    Gamma(log), and poisson(log). The gaussian(identity) case is 
                    implemented separately in the function \code{LDTFPlmm}.}
                     
   \item{xtf}{       a matrix giving the design matrix for the conditional 
                     probabilities of the random intercepts distributions.} 

   \item{prediction}{a list giving the information used to obtain conditional
                     inferences. The list includes the following
                     elements: \code{xpred} and \code{xtfnpred} giving the 
                     design matrices for the median and conditional probabilities,
                     respectively, used to obtain inferences about the conditional
                     densities of the random effects, and \code{quans} a double precision vector        
                     giving THREE quantiles for which inferences are obtained. If \code{quans} 
                     is not specified, the default is quans=c(0.03,0.50,0.97).}

    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{maxm} an integer giving the truncation of the tailfree
                     process, \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the linear 
                     dependent tailfree prior, \code{alpha} giving the value of the 
                     precision parameter (it must be specified if \code{a0} is missing), 
                     \code{mub} giving the mean of the normal prior of the fixed effects,
                     \code{Sb} giving the (co)variance of the normal prior 
                     distribution for the fixed effects, and
                     \code{taub1} and \code{taub2} giving th hyperparameters of the inv-gamma
                     distribution for the centering variance.}

    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following elements: \code{nburn} an integer giving the number of burn-in 
                     scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an 
                     integer giving the total number of scans to be saved, \code{ndisplay} an 
                     integer giving the number of saved scans to be displayed on screen (the
                     function reports on the screen when every \code{ndisplay} iterations 
                     have been carried out).}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}
    
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state}.}

    \item{ngrid}{    integer giving the number of grid points where the conditional density estimate         
                     is evaluated. The default is 100.}

    \item{grid}{     vector of grid points where the conditional densities are 
                     evaluated. The default value is NULL
                     and the grid is chosen according to the range of the
	             data.}

    \item{compute.band}{  logical variable indicating whether the credible band for the 
                          conditional density and mean function must be computed.}       

    \item{type.band}{  string indication the type of credible band to be computed; if equal to 
                       "HPD" or "PD" then the 95 percent pointwise HPD or PD band is computed,
                        respectively.}       

    \item{data}{     data frame.}     
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{LDTFPdensity} to print an error message and terminate if there are any
                     incomplete observations.}       

    \item{work.dir}{working directory.} 

}

\details{
  This generic function fits a generalized linear mixed effects model using a linear dependent
  tailfree prior for the random intercepts (Jara and Hanson, 2011). The linear
  predictor is modeled as follows:

  \deqn{\eta_{ij} = x_{ij}' \beta + b_i, i=1,\ldots,N, j=1,\ldots,n_i}{yij = xij' beta + bi, i=1,\ldots,N, j=1,\ldots,n_i}
  \deqn{b_i | G_{xtf_i} \sim G_{xtfi}}{bi | Gxtfi ~ Gxtfi}
  \deqn{\{G_{xtf}: xtf \in \mathcal{X} \}| maxm,\alpha, \sigma^2_b \sim LDTFP^{maxm}(h,\Pi^{\sigma^2_b},\textit{A}^{\alpha,\rho})}{{Gxtf: xtf in X} | maxm, alpha, sigma2b ~ LDTFP^maxm(h,Pi^{sigma2b},\textit{A}^{alpha,rhi})}
   where, h is the logistic CDF, and \eqn{G_{xtf}}{Gxtf} is median-zero and centered around an  
   \eqn{N(0,\sigma^2b)}{N(0,sigma2b)} distribution. To complete the model specification, independent hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\sigma^{-2}_b | \tau_{b1}, \tau_{b2} \sim Gamma(\tau_{b1}/2,\tau_{b2}/2)}{sigma^-2b | taub1, taub2 ~ Gamma(taub1/2,taub2/2)}

  The precision parameter, \eqn{\alpha}{alpha}, of the \code{LDTFP} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value, set \eqn{a_0}{a0} to NULL in the prior specification.

  The full conditional distribution for the fixed effects is updated using a MH step based on
  an IWLS proposal (see, e.g.,  Jara, Hanson and Lesaffre, 2009). The remaining parameters are
  sampled using the slice sampling algorithm (Neal, 2003).
}

\value{
  An object of class \code{LDTFPglmm} representing the LDTFP model fit. 
  Generic functions such as \code{print}, \code{plot},
  and \code{summary} have methods to show the results of the fit. The results include 
  \code{beta}, \code{alpha} and \code{sigma^2_b}.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 

                  
  \item{alpha}{a double precision giving the value of the precision parameter.} 

  \item{b}{a double precision giving the value of the random effects.}

  \item{beta}{a vector giving the value of the fixed effects.}

  \item{sigma2b}{a double precision giving the value of the centering variance.}

  \item{betatf}{a matrix giving the 
                regression coefficients for each conditional pribability.}
}

\seealso{
\code{\link{LDTFPdensity}}, \code{\link{LDTFPsurvival}}
}

\references{
Jara, A., Hanson, T. (2011). A class of mixtures of dependent tail-free 
  processes. Biometrika, 98(3): 553 - 566.

Jara, A., Hanson, T., Lesaffre, E. (2009) Robustifying generalized linear 
   mixed models using a new class of mixtures of multivariate Polya trees. 
   Journal of Computational and Graphical Statistics, 18(4): 838-860.

Neal, R. (2003) Slice sampling. Anals of Statistics, 31: 705 - 767. 
}

\examples{
\dontrun{
    ##############################################
    # A simulated data using "perfect"
    # simulation from a mixture of two 
    # normals and normal true models for 
    # the random effects.
    # A Poisson sampling distribution 
    # is considered.
    ##############################################

    # Functions needed to simulate random effects
    # and to evaluate true models

      findq <- function(true.cdf,target,low,
                       upp,epsilon=0.0000001)
      {
         plow <- true.cdf(low)
         pupp <- true.cdf(upp)
         pcenter <- true.cdf((upp+low)/2) 
         err <- abs(pcenter-target)
         i <- 0 
         while(err > epsilon)
         {
               i <- i + 1
               if(target< pcenter)
               {  
                  upp <- (upp+low)/2
                  pupp <- pcenter
                  pcenter <- true.cdf((upp+low)/2) 
                  err <- abs(pcenter-target)
               } 
               if(target>= pcenter)
               {  
                  low <- (upp+low)/2
                  plow <- pcenter
                  pcenter <- true.cdf((upp+low)/2) 
                  err <- abs(pcenter-target)
               } 
           }
           return((upp+low)/2)	
      }	

      true.dens1 <- function(x)
      {
           0.5*dnorm(x,2.,sqrt(0.005))+
           0.5*dnorm(x,2.85,sqrt(0.005))
      }

      true.dens2 <- function(x)
      {
           dnorm(x,2.1,sqrt(0.0324))
      }

      true.cdf1 <- function(x)
      {
           0.5*pnorm(x,2.,sqrt(0.005))+
           0.5*pnorm(x,2.85,sqrt(0.005))
      }

      true.cdf2 <- function(x)
      {
           pnorm(x,2.1,sqrt(0.0324))
      }

    # Simulation of random effects

      nsubject <- 200
      nsim <- nsubject/2 
      qq <- seq(1,nsim)/(nsim+1)
      b1 <- rep(0,nsim)
      for(i in 1:nsim)
      {
          aa <- findq(true.cdf1,qq[i],low=-6,upp=6)
          b1[i] <- aa 
      }	

      b2 <- rep(0,nsim)
      for(i in 1:nsim)
      {
         aa <- findq(true.cdf2,qq[i],low=-6,upp=6)
         b2[i] <- aa 
      }	

      trt <- c(rep(0,nsim),rep(1,nsim))
      b <- c(b1,b2)

      xtf <- cbind(rep(1,nsubject),trt)

    # Simulation of responses

      ni <- 5
      nrec <- nsubject*ni
      y <- NULL
      g <- NULL
      x <- NULL

      z <- rnorm(nrec)

      ll <- 0
      for(i in 1:nsubject)
      {     
          g <- c(g,rep(i,ni))
          for(j in 1:ni)
          {   
              ll <- ll +1
              etaij <- b[i] + 1.2*z[ll]
              ytmp <- rpois(1,exp(etaij))
              y <- c(y,ytmp)
              x <- rbind(x,c(1,trt[i],z[ll]))
          }
     }
     colnames(x) <- c("Intercept","trt","z")

   # Design matrix for prediction 
       
     xpred <- rbind(c(1,0,0),c(1,1,0))
     xtfpred <- rbind(c(1,0),c(1,1)) 

     prediction <- list(xpred=xpred,
                        xtfpred=xtfpred,
                        quans=c(0.03,0.50,0.97))

   # Prior information
     prior <- list(maxm=5,
                   alpha=0.5,
                   mub=rep(0,3),
                   Sb=diag(1000,3),
                   taub1=2.002,
                   taub2=2.002)

   # Initial state
     state <- NULL


   # MCMC parameters
     nburn <- 4000
     nsave <- 4000
     nskip <- 3
     ndisplay <- 500
     mcmc <- list(nburn=nburn,
                  nsave=nsave,
                  nskip=nskip,
                  ndisplay=ndisplay)


   # Fitting the model
     fit1 <- LDTFPglmm(y=y,x=x,g=g,family=poisson(log),
                       xtf=xtf,grid=seq(1.2,3.2,len=200),
                       prediction=prediction,
                       prior=prior, 
                       mcmc=mcmc,     
                       state=state,
                       status=TRUE,
                       compute.band=TRUE)

   # Plotting density estimates and true models
   # for the random intercepts

     par(cex=1.7,mar=c(4.1, 4.1, 1, 1))
     plot(fit1$grid,fit1$densu[1,],type="l",xlab="y",
          ylab="f(y|x)",lty=2,lwd=3,main="trt=0")
     lines(fit1$grid,fit1$densl[1,],lty=2,lwd=3)
     lines(fit1$grid,fit1$densm[1,],lty=1,lwd=3)
     tmp1 <- true.dens1(fit1$grid)
     lines(fit1$grid,tmp1,lty=1,lwd=3,col="red") 

     par(cex=1.7,mar=c(4.1, 4.1, 1, 1))
     plot(fit1$grid,fit1$densu[2,],type="l",xlab="y",
          ylab="f(y|x)",lty=2,lwd=3,main="trt=1")
     lines(fit1$grid,fit1$densl[2,],lty=2,lwd=3)
     lines(fit1$grid,fit1$densm[2,],lty=1,lwd=3)
     tmp1 <- true.dens2(fit1$grid)
     lines(fit1$grid,tmp1,lty=1,lwd=3,col="red") 

}
}

\author{
Alejandro Jara \email{<atjara@uc.cl>}
}

\keyword{models}
\keyword{nonparametric}
