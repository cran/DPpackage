\name{LDTFPdensity}
\alias{LDTFPdensity}
\alias{LDTFPdensity.default}

\title{Density Regression using Linear Dependent Tailfree Processes}
\description{
    This function generates a posterior density sample for a 
    Linear Dependent Tailfree Process model for conditional
    density estimation.
}
   
\usage{

LDTFPdensity(y,x,xtf,prediction,prior,mcmc,
             state,status,ngrid=100,
             grid=NULL,compute.band=FALSE,
             type.band="PD",
             data=sys.frame(sys.parent()),
             na.action=na.fail,
             work.dir=NULL)
}

\arguments{
   \item{y}{         a vector giving the response variables.}

   \item{x}{         a matrix giving the design matrix for the median function.} 

   \item{xtf}{       a matrix giving the design matrix for the conditional 
                     probabilities.} 

   \item{prediction}{a list giving the information used to obtain conditional
                     inferences. The list includes the following
                     elements: \code{xdenpred} and \code{xtfdenpred} giving the 
                     design matrices for the median and conditional probabilities,
                     respectively, used to obtain inferences about the conditional
                     densities and survival functions, \code{xmedpred} and 
                     \code{xtfmedpred} giving the design matrices for the median and conditional 
                     probabilities, respectively, used to obtain inferences about quantiles,
                     and \code{quans} a double precision vector giving THREE quantiles
                     for which inferences are obtained. If \code{quans} is not specified,
                     the default is quans=c(0.03,0.50,0.97).}
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{maxn} an integer giving the truncation of the tailfree
                     process, \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the linear 
                     dependent tailfree prior, \code{alpha} giving the value of the 
                     precision parameter (it must be specified if \code{a0} is missing), 
                     \code{mub} giving the mean of the normal prior of the median regression
                     coefficients, \code{Sb} giving the (co)variance of the normal prior 
                     distribution for the median regression coefficents, and
                     \code{tau1} and \code{tau2} giving th hyperparameters of the inv-gamma
                     distribution for the centering variance.}

    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following elements: \code{nburn} an integer giving the number of burn-in 
                     scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an integer 
                     giving the total number of scans to be saved, \code{ndisplay} an integer giving
                     the number of saved scans to be displayed on screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out).}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}
    
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state}.}

    \item{ngrid}{    integer giving the number of grid points where the conditional density estimate is 
                     evaluated. The default is 100.}

    \item{grid}{     vector of grid points where the conditional densities are 
                     evaluated. The default value is NULL
                     and the grid is chosen according to the range of the
	             data.}

    \item{compute.band}{  logical variable indicating whether the credible band for the conditional density  
                          and mean function must be computed.}       

    \item{type.band}{  string indication the type of credible band to be computed; if equal to "HPD" or
                       "PD" then the 95 percent pointwise HPD or PD band is computed, respectively.}       

    \item{data}{     data frame.}     
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{LDTFPdensity} to print an error message and terminate if there are any
                     incomplete observations.}       

    \item{work.dir}{working directory.} 

}

\details{
  This generic function fits a Linear Dependent Tailfree process (Jara and Hanson, 2011), given by:
  \deqn{y_i = x_i' \beta + v_i, i=1,\ldots,n}{yi = xi' beta + vi, i=1,\ldots,n}
  \deqn{v_i | G_{xtf_i} \sim G_{xtfi}}{vi | Gxtfi ~ Gxtfi}

  \deqn{\{G_{xtf}: xtf \in \mathcal{X} \}| maxm,\alpha, \sigma^2 \sim LDTFP^{maxm}(h,\Pi^{\sigma^2},\textit{A}^{\alpha,\rho})}{{Gxtf: xtf in X} | maxm, alpha, sigma2 ~ LDTFP^maxm(h,Pi^{sigma2},\textit{A}^{alpha,rhi})}
   where, h is the logistic CDF, and \eqn{G_{xtf}}{Gxtf} is median-zero and centered around an  
   \eqn{N(0,\sigma^2)}{N(0,sigma2)} distribution. To complete the model specification, independent hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\sigma^{-2} | \tau_1, \tau_2 \sim Gamma(\tau_1/2,\tau_2/2)}{sigma^-2 | tau1, tau2 ~ Gamma(tau1/2,tau2/2)}

  The precision parameter, \eqn{\alpha}{alpha}, of the \code{LDTFP} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value, set \eqn{a_0}{a0} to NULL in the prior specification.

  The computational implementation of the model is based on Slice sampling (Neal, 2003).
}

\value{
  An object of class \code{LDTFPdensity} representing the LDTFP model fit. 
  Generic functions such as \code{print}, \code{plot},
  and \code{summary} have methods to show the results of the fit. The results include 
  \code{beta}, \code{alpha} and \code{sigma^2}.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 

                  
  \item{alpha}{a double precision giving the value of the precision parameter.} 

  \item{betace}{a vector giving the value of the median regression coefficient.}

  \item{sigma^2}{a double precision giving the value of the centering variance.}

  \item{betatf}{a matrix giving the 
                regression coefficients for each conditional pribability.}
}

\seealso{
\code{\link{DPcdensity}}, \code{\link{LDDPdensity}}
}

\references{
Jara, A., Hanson, T. (2011). A class of mixtures of dependent tail-free 
  processes. Biometrika, 98(3): 553 - 566.

Neal, R. (2003) Slice sampling. Anals of Statistics, 31: 705 - 767. 
}

\examples{
\dontrun{

    ######################## 
    # IgG data
    ########################
      data(igg)
      z <- igg$age
      y <- log(igg$igg)

    # Design matrices
      ages1 <- z^2
      ages2 <- 1/ages1

      x <- cbind(rep(1,length(y)),ages1,ages2)
      xtf <- cbind(rep(1,length(y)),ages1,ages2)
       
      colnames(x) <- c("(Intercept)","age^2","age^{-2}")
      colnames(xtf) <- c("(Intercept)","age^2","age^{-2}")

    # Prediction
      xdpred <- c(11/12,25/12,38/12,52/12,65/12)
      agesp1 <- xdpred^2
      agesp2 <- 1/agesp1
      xdenpred <- cbind(rep(1,length(xdpred)),agesp1,agesp2)
      xtfdenpred <- xdenpred

      xmpred <- seq(0.5,6,len=50)
      agesp1 <- xmpred^2
      agesp2 <- 1/agesp1
      xmedpred <- cbind(rep(1,length(xmpred)),agesp1,agesp2)
      xtfmedpred <- xmedpred
 
      prediction <- list(xdenpred=xdenpred,
			 xtfdenpred=xtfdenpred,
                         xmedpred=xmedpred,
                         xtfmedpred=xtfmedpred,
                         quans=c(0.03,0.50,0.97))

    # Prior information
      Sb <- diag(1000,3)
      mub <- rep(0,3)
 
      prior <- list(maxm=4,
                    a0=1,
                    b0=1,
                    mub=mub,
                    Sb=Sb,
                    tau1=2.02,
                    tau2=2.02)

    # Initial state
      state <- NULL

    # MCMC parameters

      mcmc <- list(nburn=5000,
                   nsave=5000,
                   nskip=4,
                   ndisplay=200)

    # Fitting the model
 
      fit1 <- LDTFPdensity(y=y,
			   x=x,
			   xtf=xtf,
			   prediction=prediction,
			   prior=prior, 
			   mcmc=mcmc,
			   state=state, 
			   status=TRUE,
			   compute.band=TRUE)

      fit1
      summary(fit1)
      plot(fit1)

    # Plots predictions
    # (conditional densities and quantile functions)

      par(mfrow=c(3,2))      
 
      for(i in 1:5)
      {
         plot(fit1$grid,fit1$densm[i,],lwd=2,
              type="l",lty=1,col="black",xlab="log IgG",ylab="density",
              ylim=c(0,2))
         lines(fit1$grid,fit1$densl[i,],lwd=1,lty=2,col="black")    
         lines(fit1$grid,fit1$densu[i,],lwd=1,lty=2,col="black")    
      }
 
      plot(z,y,ylab="log IgG",xlab="Age (years)")
      lines(xmpred,fit1$qmm[,2],lwd=2)
      lines(xmpred,fit1$qml[,2],lwd=1,lty=2)
      lines(xmpred,fit1$qmu[,2],lwd=1,lty=2)

      lines(xmpred,fit1$qmm[,1],lwd=2)
      lines(xmpred,fit1$qml[,1],lwd=1,lty=2)
      lines(xmpred,fit1$qmu[,1],lwd=1,lty=2)

      lines(xmpred,fit1$qmm[,3],lwd=2)
      lines(xmpred,fit1$qml[,3],lwd=1,lty=2)
      lines(xmpred,fit1$qmu[,3],lwd=1,lty=2)


    #######################################
    # A simulated data using "perfect"
    # simulation (mixture of two normals 
    # and normal true models).
    #######################################

    # Functions needed to simulate data
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
            0.5*dnorm(x,2.5,sqrt(0.005))+
            0.5*dnorm(x,2.85,sqrt(0.005))
      }

      true.dens2 <- function(x)
      {
            dnorm(x,2.1,sqrt(0.0324))
      }

      true.cdf1 <- function(x)
      {
            0.5*pnorm(x,2.50,sqrt(0.005))+
            0.5*pnorm(x,2.85,sqrt(0.005))
      }

      true.cdf2 <- function(x)
      {
            pnorm(x,2.1,sqrt(0.0324))
      }

    # Simulation
      nsim <- 500 
      qq <- seq(1,nsim)/(nsim+1)

      y1 <- rep(0,nsim)
      for(i in 1:nsim)
      {
          aa <- findq(true.cdf1,qq[i],low=-6,upp=6)
          y1[i] <- aa 
      }	

      y2 <- rep(0,nsim)
      for(i in 1:nsim)
      {
          aa <- findq(true.cdf2,qq[i],low=-6,upp=6)
          y2[i] <- aa 
      }	

      trt <- c(rep(0,nsim),rep(1,nsim))
      y <- c(y1,y2)

    # Design matrices
      W1 <- cbind(rep(1,2*nsim),trt)
      W2 <- cbind(rep(1,2*nsim),trt)
      colnames(W1) <- c("(Intercept)","trt")
      colnames(W2) <- c("(Intercept)","trt")

    # Design matrix for prediction 
      intp <- rep(1,2)
      trtp <- c(0,1)
      zpred <- cbind(intp,trtp) 

      prediction <- list(xdenpred=zpred,
                         xtfdenpred=zpred,
                         xmedpred=zpred,
                         xtfmedpred=zpred,
                         quans=c(0.03,0.50,0.97))

    # Prior information
      prior <- list(maxm=5,
                    a0=1,
                    b0=1,
                    mub=rep(0,2),
                    Sb=diag(1000,2),
                    tau1=2.002,
                    tau2=2.002)

    # Initial state
      state <- NULL

    # MCMC parameters
      nburn <- 5000
      nsave <- 5000
      nskip <- 4
      ndisplay <- 200
      mcmc <- list(nburn=nburn,
                   nsave=nsave,
                   nskip=nskip,
                   ndisplay=ndisplay)

    # Fitting the model
      fit1 <- LDTFPdensity(y=y,
                           x=W1,
                           xtf=W2,
                           grid=seq(1.2,3.2,len=200),
                           prediction=prediction,
                           prior=prior, 
                           mcmc=mcmc,     
                           state=state,
                           status=TRUE,
                           compute.band=TRUE)

    # Plotting density estimates and true models

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
