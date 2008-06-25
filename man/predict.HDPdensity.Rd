\name{predict.HDPdensity}
\alias{predict.HDPdensity}
\title{Computes Predictive Information for the Dependent Random Probability Measures.}
\description{
    Generates posterior predictive draws for
    future patients (observations) from the random probability measures.
    Support provided by the NIH/NCI R01CA75981 grant.
}
    
\usage{

\method{predict}{HDPdensity}(object, data.pred=NULL, j=1, r=0, nsim=100, idx.x=NULL, ...)

}

\arguments{

    \item{object}{   \code{HDPdensity} fitted model object.}

    \item{data.pred}{ response (dummy values) and covariates
                      for future patients, (\code{npa} by \code{p}) matrix.}
               
    \item{j}{        study.}

    \item{r}{        indicator for including (0) or not (1) the
                     common measure.}

    \item{nsim}{     number of imputed posterior simulations to use.}

    \item{idx.x}{    vector of size \code{px}, columns (starting to count at 1
                     for the 1st column) that contain the \code{px} covariates.
                     The remaining columns are dummies corresponding to the
                     (\code{p}-\code{px})-dimensional response vector.}

    \item{...}{      further arguments to be passed.}     
}

\details{
  Must run \code{\link{HDPdensity}} first to generate posterior
  simulations.

  The function carries out post-processing of the MCMC posterior
  simulation to generate posterior predictive simulation for
  future observations from the random probability measures
  defined in the model.
  
  For \code{npa} assumed future patients with given covariates
  (specified in \code{xpred}) the function computes posterior 
  predictive inference of future responses.
  The subvector of responses is a dummy to match the dimension. 
}

\value{
  The function returns a matrix \code{zout} with
  \code{p+1} columns of posterior predictive simulations for the
  \code{npa} future patients.
  The first column is a patient index. An index i refers to the
  i-th row in the matrix \code{data.pred} of given patient covariates.
  Columns 2 through \code{p+1} are posterior predictive simulations
  including the (unchanged) covariate vector in the locations
  indicated by \code{idx}.
  
  Scatterplots, density estimates etc. of the posterior
  predictive simulations can be used to evaluate posterior means for
  the RPMs, and to evaluate posterior predictive probabilities
  for events of interest for future subjects (patients).

  See the examples below for examples on how to summarize
  the posterior predictive simulations.
}

\seealso{
\code{\link{HDPdensity}}
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
  
    # Prior information
      Z <- calgb[,1:10]
      mhat <- apply(Z,2,mean)
      v <- diag(var(Z))
     
      prior<-list(a0=1,
                  b0=1,
                  pe1=0.1,
                  pe0=0.1,
                  ae=1,
                  be=1,
                  a=mhat,
                  A=diag(v), 
                  q=15,
                  R=0.25*diag(v),
                  cc=15,
                  C=diag(v))

    # Initial state
      state <- NULL

    # MCMC parameters

      mcmc <- list(nburn=1000,
                   nsave=2000,
                   nskip=0,
                   ndisplay=100,
                   npredupdate=100)


    # Fitting the model
      fit1 <- HDPdensity(formula=cbind(Z1,Z2,Z3,T1,T2,B0,B1)~CTX+GM+AMOF,
                         study=~study,
                         prior=prior,
                         mcmc=mcmc,
                         state=state,
                         data=calgb,  
                         status=TRUE)

    # Load data for future patients (for prediction)
      data(calgb.pred)
      X <- calgb.pred 

    # post-process MCMC output for predictive inference
    # save posterior predictive simulations in z00 ... z30

      z10 <- predict(fit1,data.pred=X,j=1,r=0) # post prediction for study 1
      z20 <- predict(fit1,data.pred=X,j=2,r=0) # .. study 2
      z30 <- predict(fit1,data.pred=X,j=3,r=0) # .. population at large (= study 3)

      z11 <- predict(fit1,data.pred=X,j=1,r=1) # idiosyncratic measures study 1
      z21 <- predict(fit1,data.pred=X,j=2,r=1) # .. study 2
      z00 <- predict(fit1,data.pred=X,j=0,r=0) # common measure

    # covariates (and dummy responses) of future patients
      colnames(z00) <- c("PATIENT",colnames(X))

    # plot estimated density for future patients in study 1, 2 and
    # in population at large
      idx <- which(z10[,1]==1)   ## PATIENT 1
      options(digits=2)
      par(mfrow=c(2,1))          

    # plot prediction fo study 1,2,population
      plot  (density(z10[idx,8]),
             ylim=c(0,1.5),xlim=c(-0.5,2.5),
             xlab="SLOPE OF RECOVERY",bty="l",main="FUTURE PAT 1")
      lines (density(z20[idx,8]),type="l",col=2)
      lines (density(z30[idx,8]),type="l",col=3)
      legend(-0.5,1.5,col=1:3,legend=c("STUDY 1","STUDY 2","POPULATION"),
             lty=c(1,1,1),bty="n")

    # common and idiosyncratic measures
      plot (density(z00[idx,8]),type="l",col=4,lty=1,
            ylim=c(0,1.5),xlim=c(-0.5,2.5),
            xlab="SLOPE OF RECOVERY",bty="l",main="COMMON & IDIOSYNC PARTS")
      lines (density(z11[idx,8]),type="l",col=1,lty=2)
      lines (density(z21[idx,8]),type="l",col=2,lty=2)
      legend(1.5,1.5,col=c(1,2,4),lty=c(2,2,1),
             legend=c("STUDY 1 (idiosyn.)",
                      "STUDY 2 (idiosyn.)",
                      "COMMON"),bty="n")

    # plot estimated density for future patients in study 1, 2 and
    # in population at large
      idx <- which(z10[,1]==2)   ## PATIENT 2
      options(digits=2)
      par(mfrow=c(2,1))

      plot  (density(z10[idx,8]),
             ylim=c(0,1.5),xlim=c(-0.5,2.5),
             xlab="SLOPE OF RECOVERY",bty="l",main="FUTURE PAT 2")
      lines (density(z20[idx,8]),type="l",col=2)
      lines (density(z30[idx,8]),type="l",col=3)
      legend(-0.5,1.5,col=1:3,legend=c("STUDY 1","STUDY 2","POPULATION"),
             lty=c(1,1,1),bty="n")

      plot (density(z00[idx,8]),type="l",col=4,lty=1,
            ylim=c(0,1.5),xlim=c(-0.5,2.5),
            xlab="SLOPE OF RECOVERY",bty="l",main="COMMON & IDIOSYNC PARTS")
      lines (density(z11[idx,8]),type="l",col=1,lty=2)
      lines (density(z21[idx,8]),type="l",col=2,lty=2)
      legend(1.5,1.5,col=c(1,2,4),lty=c(2,2,1),
             legend=c("STUDY 1 (idiosyn.)",
                      "STUDY 2 (idiosyn.)",
                      "COMMON"),bty="n")

    # plot nadir count by covariate, for population 
      z2 <- z30[,3]; ctx <- z30[,9]; gm <- z30[,10]; amf <- z30[,11]
    # fix covariates gm (GM-CSF) and amf (aminofostine)
      idx <- which( (gm==-1.78) & (amf== -0.36) )
      boxplot(split(z2,ctx),
              xlab="CYCLOPHOSPHAMIDE",bty="n",ylab="NADIR COUNT")
}
}

\author{

Peter Mueller \email{<pmueller@mdanderson.org>}
}

\keyword{models}
