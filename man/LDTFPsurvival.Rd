\name{LDTFPsurvival}
\alias{LDTFPsurvival}
\alias{LDTFPsurvival.default}

\title{Survival Regression using Linear Dependent Tailfree Processes}
\description{
    This function generates a posterior density sample for a 
    Linear Dependent Tailfree Process model for conditional
    survival estimation of time-to-event data.
}
   
\usage{
LDTFPsurvival(y,x,xtf,prediction,prior,mcmc,
              state,status,grid=seq(0.01,60,len=100),
              compute.band=FALSE,type.band="PD",
              data=sys.frame(sys.parent()),
              na.action=na.fail,work.dir=NULL)
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

    \item{grid}{     vector of grid points where the conditional survival functions are 
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
                     \code{LDTFPsurvival} to print an error message and terminate if there are any
                     incomplete observations.}       

    \item{work.dir}{working directory.} 

}

\details{
  This generic function fits a Linear Dependent Tailfree process (Jara and Hanson, 2011),
  for (potentially) interval-censored data. Let \eqn{T_i \in R^+}{T_i in R+} be the
  time-to-event for subject i and set \eqn{z_i = \log T_i}{z_i = log T_i}. The model for the
  log time-to.event data is given by:  
  \deqn{z_i = x_i' \beta + v_i, i=1,\ldots,n}{zi = xi' beta + vi, i=1,\ldots,n}
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
  An object of class \code{LDTFPsurvival} representing the LDTFP model fit. 
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
                regression coefficients for each conditional probability.}

  \item{z}{a vector giving the current value of the (imputed) survival times.}

}

\seealso{
\code{\link{LDTFPdensity}},\code{\link{LDDPsurvival}}
}

\references{
Jara, A., Hanson, T. (2011). A class of mixtures of dependent tail-free 
  processes. Biometrika, 98(3): 553 - 566.

Neal, R. (2003) Slice sampling. Anals of Statistics, 31: 705 - 767. 
}

\examples{
\dontrun{

   #############################################################
   # Time to Cosmetic Deterioration of Breast Cancer Patients
   #############################################################

     data(deterioration)
     attach(deterioration)
     y <- cbind(left,right)
      
   # Design matrix

     x <- cbind(rep(1,length(trt)),trt)
     xtf <- cbind(rep(1,length(trt)),trt)
     colnames(x) <- c("(Intercept)","trt")
     colnames(xtf) <- c("(Intercept)","trt")
      
   # Prediction

     xdenpred <- cbind(rep(1,2),c(0,1))
     xtfdenpred <- cbind(rep(1,2),c(0,1))
     xmedpred <- cbind(rep(1,2),c(0,1))
     xtfmedpred <- cbind(rep(1,2),c(0,1))

     prediction <- list(xdenpred=xdenpred,
                        xtfdenpred=xtfdenpred,
                        xmedpred=xmedpred,
                        xtfmedpred=xtfmedpred,
                        quans=c(0.03,0.50,0.97))

   # Prior information

     prior <- list(maxm=5,
                   a0=1,
                   b0=1,
                   mub=rep(0,2),
                   Sb=diag(1000,2),
                   tau1=2,002,
                   tau2=2.002)

   # Initial state
     state <- NULL

   # MCMC parameters

     mcmc <- list(nburn=5000,
                  nsave=5000,
                  nskip=4,
                  ndisplay=200)

   # Fitting the model
 
     fit1 <- LDTFPsurvival(y=y,
                           x=x,
                           xtf=xtf,
                           prediction=prediction,
                           prior=prior, 
                           mcmc=mcmc,     
                           state=state,
                           grid=seq(0.01,70,len=200),
                           status=TRUE,
                           compute.band=TRUE)

     fit1
     summary(fit1)
     plot(fit1)

   # Plotting survival functions estimates
 
     par(cex=1.7,mar=c(4.1, 4.1, 1, 1))
     x1 <- fit1$grid
     y1 <- fit1$survml[1,]
     x2 <- fit1$grid
     y2 <- fit1$survmu[1,]
     aa <- rbind(x2,y2)[, order(-x2, y2)]
     x2 <- aa[1,]
     y2 <- aa[2,]
     plot(fit1$grid,fit1$survmu[1,],type="l",
          xlab="months",ylab="survival",
          lty=1,lwd=2,ylim=c(0,1),col="lightgray")
     polygon(x=c(x1,x2),y=c(y1,y2),border=NA,col="lightgray")
     lines(fit1$grid,fit1$survmm[1,],lty=1,lwd=3)

     par(cex=1.7,mar=c(4.1, 4.1, 1, 1))
     x1 <- fit1$grid
     y1 <- fit1$survml[2,] 
     x2 <- fit1$grid
     y2 <- fit1$survmu[2,]
     aa <- rbind(x2,y2)[, order(-x2, y2)]
     x2 <- aa[1,]
     y2 <- aa[2,]
     plot(fit1$grid,fit1$survmu[2,],type="l",
          xlab="months",ylab="survival",
          lty=1,lwd=2,ylim=c(0,1),col="lightgray")
     polygon(x=c(x1,x2),y=c(y1,y2),border=NA,col="lightgray")
     lines(fit1$grid,fit1$survmm[2,],lty=1,lwd=3)
}
}

\author{
Alejandro Jara \email{<atjara@uc.cl>}

}

\keyword{models}
\keyword{nonparametric}
