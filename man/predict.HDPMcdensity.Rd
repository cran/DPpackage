\name{predict.HDPMcdensity}
\alias{predict.HDPMcdensity}

\title{Predictive Information for the Dependent Random Probability Measures.}
\description{
    Plot the probability measures arising from a HDPM of normals model
    for conditional density estimation.
    Support provided by the NIH/NCI R01CA75981 grant.
}
    
\usage{
\method{predict}{HDPMcdensity}(object,pred,i,r,ask=TRUE,nfigr=2,nfigc=2, ...)
}

\arguments{

    \item{object}{   \code{HDPMcdensity} fitted model object.}
 
    \item{pred}{     indicator for the values of the predictors, given
                     by the row pred in xpred, for which the conditional
                     densities must be drawn.}
 
    \item{i}{        study indicator.}

    \item{r}{        indicator for including (0) or not (1) the
                     common measure.}

    \item{ask}{      logical variable indicating whether the plots must be displayed
                     sequentially or not.}

    \item{nfigr}{    number of rows in the figure.}

    \item{nfigc}{    number of columns in the figure.}

    \item{...}{      further arguments to be passed.}     
}

\details{
  Must run \code{\link{HDPMcdensity}} first to generate posterior
  simulations.
}

\seealso{
\code{\link{HDPMcdensity}}
}

\references{
Mueller, P., Quintana, F. and Rosner, G. (2004). A Method for Combining Inference 
  over Related Nonparametric Bayesian Models. Journal of the Royal Statistical 
  Society, Series B, 66: 735-749.
}


\examples{
\dontrun{
    # Data
      data(calgb)
      attach(calgb)
      y <- cbind(Z1,Z2,Z3,T1,T2,B0,B1)
      x <- cbind(CTX,GM,AMOF)
  
      z <- cbind(y,x)

    #  Data for prediction
      data(calgb.pred)
      xpred <- as.matrix(calgb.pred[,8:10])


    # Prior information
      prior <- list(pe1=0.1,
                    pe0=0.1,
                    ae=1,
                    be=1,
                    a0=rep(1,3),
                    b0=rep(1,3),
                    nu=12,
                    tinv=0.25*var(z),
 		  m0=apply(z,2,mean),
                    S0=var(z),
 		  nub=12,
                    tbinv=var(z))		


    # Initial state
      state <- NULL

    # MCMC parameters

      mcmc <- list(nburn=5000,
                   nsave=5000,
                   nskip=3,
                   ndisplay=100)

    # Fitting the model
      fit1 <- HDPMcdensity(formula=y~x,
                          study=~study,
                          xpred=xpred,
                          prior=prior,
                          mcmc=mcmc,
                          state=state,
                          status=TRUE)

    # Posterior inference
      fit1
      summary(fit1)
       
    # Plot the parameters
    # (to see the plots gradually set ask=TRUE)
      plot(fit1,ask=FALSE)

    # Plot the a specific parameters 
    # (to see the plots gradually set ask=TRUE)
      plot(fit1,ask=FALSE,param="eps",nfigr=1,nfigc=2)

    # Plot the measure for each study 
    # under first values for the predictors, xpred[1,]
      predict(fit1,pred=1,i=1,r=1) # pred1, study 1
      predict(fit1,pred=1,i=2,r=1) # pred1, study 2

    # Plot the measure for each study 
    # under second values for the predictors, xpred[2,]
      predict(fit1,pred=2,i=1,r=1) # pred2, study 1
      predict(fit1,pred=2,i=2,r=1) # pred2, study 2

    # Plot the idiosyncratic measure for each study
    # under first values for the predictors, xpred[1,]
      predict(fit1,pred=1,i=1,r=0) # study 1
      predict(fit1,pred=1,i=2,r=0) # study 2

    # Plot the common measure
    # under first values for the predictors, xpred[1,]
      predict(fit1,pred=1,i=0)
}
}

\author{
Alejandro Jara \email{<atjara@uc.cl>}

Peter Mueller \email{<pmueller@mdanderson.org>}

}

\keyword{models}
\keyword{nonparametric}
