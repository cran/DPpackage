\name{HDPdensity}
\alias{HDPdensity}
\alias{HDPdensity.default}

\title{Bayesian analysis for a Hierarchical Mixture of Dirichlet Process Mixture of Normals}
\description{
    This function generates a posterior density sample for a 
    DP mixture of normals model for related random probability measures.
    Support provided by the NIH/NCI R01CA75981 grant.
}
   
\usage{
HDPdensity(formula,study,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail,
           work.dir=NULL)
}

\arguments{
   \item{formula}{   a two-sided linear formula object describing the
                     model fit, with the response on the
                     left of a \code{~} operator and the terms, separated by \code{+}
                     operators, on the right. The design matrix is used to model 
                     the distribution of the responses in the HDP mixture of normals model.} 

   \item{study}{     a (1 by \code{n}) vector of study indicators. The i-th index is the study j 
                     that response i belongs to.}
    
    \item{prior}{    a list giving the prior information. The list includes the following
                     parameters: \code{pe1} and \code{pe0} giving the prior weights for the point mass at 
                     \eqn{\epsilon=1} and at \eqn{\epsilon=1}, respectively, \code{ae} and \code{be} 
                     giving the prior parameters for a Beta prior on \eqn{\epsilon}, \code{eps} giving 
                     the value of \eqn{\epsilon} (it must be specified if \code{pe1} is missing),
                     \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Dirichlet process
                     prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{a0} is missing), \code{a} and \code{A} 
                     giving the hyperparameters of the normal prior distribution 
                     for the mean of the normal baseline distribution, \code{m} giving the mean 
                     of the normal baseline distribution (is must be specified if \code{a} is missing),
                     \code{cc} and \code{C} giving the hyperparameters of the 
                     Wishart prior distribution for the inverse of the scale matrix of the normal
                     baseline distribution, \code{B} giving the covariance matrix of the normal
                     baseline distribution (is must be specified if \code{cc} is missing),
                     \code{q} and \code{R} giving the hyperparameters of the 
                     Wishart prior distribution for the inverse of the scale matrix of the normal
                     kernel, and \code{S} giving the covariance matrix of the normal
                     kernal (is must be specified if \code{q} is missing).}

    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \code{nburn} giving the number of burn-in 
                     scans, \code{nskip} giving the thinning interval, \code{nsave} giving
                     the total number of scans to be saved, \code{ndisplay} giving
                     the number of saved scans to be displayed on screen.}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis (not
                     available yet).}
    
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state} (not available yet).}

    \item{data}{     data frame.}     
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{HDPdensity} to print an error message and terminate if there are any
                     incomplete observations.}       

     \item{work.dir}{working directory.} 
}

\details{
  The function sets up and carries out posterior Markov chain
  Monte Carlo (MCMC) simulation for a hierarchical DP mixture
  model.

  The model is a DP mixture of normals for related random probability
  measures \eqn{H_j}. Each random measure is assumed to arise
  as a mixture \eqn{H_j = \epsilon F_0 + (1-\epsilon) F_j} of one common
  distribution \eqn{F_0} and a distribution \eqn{F_j} that is specific to
  the j-th submodel.

  See \cite{Mueller, Quintana and Rosner (2004)} for details of the
  model. In summary, the implemented model is as follows. Without loss of 
  generality we assume that each submodel corresponds to a different study 
  in a set of related studies. Let \eqn{\theta_{ij}} denote the i-th 
  observation in the j-th study (we use \eqn{\theta}, assuming that the model 
  would typically be used for a random effects distribution). We assume that 
  \eqn{\theta_{ji}, i=1,\ldots,n_j} are samples from a random probability 
  measure for the j-th study, which in turn is a mixture of a measure 
  \eqn{F_0} that is common to all studies, and an idiosyncratic measure 
  \eqn{F_j} that is specific to the j-th study.

  \deqn{
  \theta_{ji} \sim \epsilon F_0 + (1-\epsilon) F_j
  }

  The random probability measures \eqn{F_j} in turn are given a
  Dirichlet process mixture of normal prior. We assume 
  \deqn{
  F_j(\theta) = \int N(\mu,S)  dG_j(\mu),~ j=0,1,\ldots,J
  }
  with \eqn{G_j \sim DP(G^\star(\eta),\alpha)}.
  Here \eqn{\eta} are hyperparameters that index the base measure of the
  DP prior. We use a normal base measure and a conjugate hyperprior
  \deqn{
  G^\star(\mu) = N(m,B), \mbox{ with } m \sim N(a,A),
  \mbox{ and }
  B^{-1}\sim Wishart(cc,(ccC)^{-1})
  }
  The Wishart prior is parametrized such that \eqn{E(B^{-1}=C^{-1})}.
  Let \eqn{\delta_x} denote a point mass at x.
  We complete the model with the hyperpriors
  \deqn{
  S^{-1} \sim W(q,(qR)^{-1}),~
  p(\epsilon) = \pi_0\delta_0+\pi_1\delta_1+(1-\pi_0-\pi_1)
               Be(a_\epsilon,b_\epsilon) 
  }
  Regression on observation-specific covariates \eqn{x_{ji}} can be
  achieved by including \eqn{x_{ji}} with the outcome \eqn{\theta_{ji}},
  and proceeding as if \eqn{(x_{ji},\theta_{ji})} were generated
  as \eqn{\theta_{ji}} in the model described above.
  See \cite{Mueller et al. (2004, section 3.3)} for details.
}

\value{
  The function returns no value. MCMC simulations
  are saved in files in the designated working directory.

  Use \code{\link{predict.HDPdensity}} to plot summaries.

}

\seealso{
\code{\link{predict.HDPdensity}}
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
\keyword{nonparametric}
