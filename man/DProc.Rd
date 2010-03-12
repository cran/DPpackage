\name{DProc}
\alias{DProc}
\alias{DProc.default}

\title{Semiparametric Bayesian ROC curve analysis using DPM of normals}
\description{
    This function performs a ROC curve analysis based on 
    a posterior density sample for Dirichlet process mixture of normals 
    models.
}
  
\usage{
DProc(x,y,fitx=NULL,fity=NULL,ngrid=1000,priorx,priory,
      mcmcx,mcmcy,statex,statey,
      statusx,statusy,data=sys.frame(sys.parent()),
      na.action=na.fail) 
}


\arguments{
    \item{x}{        a vector giving the diagnostic marker measurements 
                     for the healthy subjects.}

    \item{y}{        a vector giving the diagnostic marker measurements 
                     for the diseased subjects.}
                     
    \item{fitx}{     a object containing the results 
                     returned by the \link{DPdensity} model fitting function
                     for the diagnostic marker measurements in the healthy 
                     subjects (Optional).}

    \item{fity}{     a object containing the results 
                     returned by the \link{DPdensity} model fitting function
                     for the diagnostic marker measurements in the diseased
                     subjects (Optional).}

    \item{ngrid}{    number of grid points where the ROC curve is 
                     evaluated. The default value is 1000.}
    
    \item{priorx}{   a list giving the prior information for the 
                     diagnostic marker measurements in the healthy 
                     subjects. See the \link{DPdensity} function for details.
                     (it must be specified if \code{fitx} is missing).}

    \item{priory}{   a list giving the prior information for the 
                     diagnostic marker measurements in the diseased 
                     subjects. See the \link{DPdensity} function for details.
                     (it must be specified if \code{fity} is missing).}

    \item{mcmcx}{    a list giving the MCMC parameters for the 
                     diagnostic marker measurements in the healthy 
                     subjects. See the \link{DPdensity} function for details.
                     (it must be specified if \code{fitx} is missing).} 
    
    \item{mcmcy}{    a list giving the MCMC parameters for the 
                     diagnostic marker measurements in the diseased 
                     subjects. See the \link{DPdensity} function for details.
                     (it must be specified if \code{fity} is missing).} 

    \item{statex}{   a list giving the current value of the parameters. 
                     See the \link{DPdensity} function for details.
                     (it must be specified if \code{fitx} is missing).}

    \item{statey}{   a list giving the current value of the parameters. 
                     See the \link{DPdensity} function for details.
                     (it must be specified if \code{fity} is missing).}

    
    \item{statusx}{  a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{statex}.}

    \item{statusy}{  a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{statex}.}

    \item{data}{     data frame.}       
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{DProc} to print an error message and terminate if there are any
                     incomplete observations.}       
}

\details{
  This generic function performs a ROC curve analysis based on 
  Dirichlet process mixture of normals models for density 
  estimation (Escobar and West, 1995):
  \deqn{x_i | \mu_{x_i}, \Sigma_{x_i} \sim N(\mu_{x_i},\Sigma_{x_i}), i=1,\ldots,n}{xi | muxi, Sigmaxi ~ N(muxi,Sigmaxi), i=1,\ldots,n}
  \deqn{(\mu_{x_i},\Sigma_{x_i}) | G_x \sim G_x}{(muxi,Sigmaxi) | Gx ~ Gx}
  \deqn{G_x | \alpha_x, G_{x_0} \sim DP(\alpha_x G_{x_0})}{G_x | alphax, Gx0 ~ DP(alphax Gx0)}

  \deqn{y_j | \mu_{y_j}, \Sigma_{y_j} \sim N(\mu_{y_j},\Sigma_{y_j}), j=1,\ldots,m}{yj | muyj, Sigmayj ~ N(muyj,Sigmayj), j=1,\ldots,m}
  \deqn{(\mu_{y_j},\Sigma_{y_j}) | G_y \sim G_y}{(muyj,Sigmayj) | Gy ~ Gy}
  \deqn{G_y | \alpha_y, G_{y_0} \sim DP(\alpha_y G_{y_0})}{G_y | alphay, Gy0 ~ DP(alphay Gy0)}

  where, \eqn{x} and \eqn{y} is the vector containing the diagnostic marker 
  measurements in the healthy and diseased subjects, respectively. 
  We refer to the help of \link{DPdensity} functions for details regarding 
  parametrization, prior specification, and implementation.

  The survival and ROC curves are estimated by using a Monte Carlo 
  approximation to the posterior means \eqn{E(G_x | x)}{E(Gx | x)} and
  \eqn{E(G_y | y)}{E(Gy | y)}, which is based on MCMC samples from 
  posterior predictive distribution for a future observation. The optimal
  cut-off point is based on the efficiency test, EFF = TP + TN, and is built 
  on Cohen's kappa as defined in Kraemer (1992).
  
  The ROC curve analysis can be performed from the data directly or from
  the outputs of the \link{DPdensity} function (see, example).
}

\value{
  An object of class \code{DProc} representing the ROC curve analysis 
  based on DP mixture of normals models fit. Generic functions such as 
  \code{print}, and \code{plot} have methods to show the results of 
  the fit. The results include the estimated densities, cdf's, and ROC curve.
}

\seealso{
\code{\link{DPdensity}}
}

\references{
Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
  Using Mixtures. Journal of the American Statistical Association, 90: 577-588.
  
Kraemer, H. C. (1992). Evaluating Medical Tests. Sage Publications.  
}

\examples{
\dontrun{
    ##############################################################
    # Fertility data example:
    # The following are Sperm Deformity Index (SDI) values from 
    # semen samples of men in an infertility study. They are 
    # divided into a "condition" present group defined as those 
    # whose partners achieved pregnancy and "condition" absent 
    # where there was no pregnancy.
    #
    # Aziz et al. (1996) Sperm deformity index: a reliable 
    # predictor of the outcome of fertilization in vitro. 
    # Fertility and Sterility, 66(6):1000-1008.
    #
    ##############################################################
    
     "pregnancy"<- c(165, 140, 154, 139, 134, 154, 120, 133, 
                     150, 146, 140, 114, 128, 131, 116, 128, 
                     122, 129, 145, 117, 140, 149, 116, 147, 
                     125, 149, 129, 157, 144, 123, 107, 129, 
                     152, 164, 134, 120, 148, 151, 149, 138, 
                     159, 169, 137, 151, 141, 145, 135, 135, 
                     153, 125, 159, 148, 142, 130, 111, 140, 
                     136, 142, 139, 137, 187, 154, 151, 149, 
                     148, 157, 159, 143, 124, 141, 114, 136, 
                     110, 129, 145, 132, 125, 149, 146, 138, 
                     151, 147, 154, 147, 158, 156, 156, 128, 
                     151, 138, 193, 131, 127, 129, 120, 159, 
                     147, 159, 156, 143, 149, 160, 126, 136, 
                     150, 136, 151, 140, 145, 140, 134, 140, 
                     138, 144, 140, 140)

     "nopregnancy"<-c(159, 136, 149, 156, 191, 169, 194, 182, 
                      163, 152, 145, 176, 122, 141, 172, 162, 
                      165, 184, 239, 178, 178, 164, 185, 154, 
                      164, 140, 207, 214, 165, 183, 218, 142, 
                      161, 168, 181, 162, 166, 150, 205, 163, 
                      166, 176)


    #########################################################
    # Estimating the ROC curve from the data
    #########################################################

    # Initial state

      statex <- NULL
      statey <- NULL

    # Prior information

      priorx <-list(alpha=10,m2=rep(0,1),
                    s2=diag(100000,1),
                    psiinv2=solve(diag(5,1)),
                    nu1=6,nu2=4,
                    tau1=1,tau2=100)

      priory <-list(alpha=20,m2=rep(0,1),
                    s2=diag(100000,1),
                    psiinv2=solve(diag(2,1)),
                    nu1=6,nu2=4,
                    tau1=1,tau2=100)

    # MCMC parameters

      nburn<-1000
      nsave<-2000
      nskip<-0
      ndisplay<-100

      mcmcx <- list(nburn=nburn,nsave=nsave,nskip=nskip,
                    ndisplay=ndisplay)
      mcmcy <- mcmcx

    # Estimating the ROC

      fit1<-DProc(x=pregnancy,y=nopregnancy,priorx=priorx,priory=priory,
                  mcmcx=mcmcx,mcmcy=mcmcy,statex=statex,statey=statey,
                  statusx=TRUE,statusy=TRUE)
      fit1
      plot(fit1)
      
      
    #########################################################
    # Estimating the ROC curve from DPdensity objects
    #########################################################
     
      fitx<-DPdensity(y=pregnancy,prior=priorx,mcmc=mcmcx,
                      state=statex,status=TRUE)
     
      fity<-DPdensity(y=nopregnancy,prior=priory,mcmc=mcmcy,
                      state=statey,status=TRUE)
 
    # Estimating the ROC

      fit2<-DProc(fitx=fitx,fity=fity)
      
      fit2
      plot(fit2)
      
}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

}

\keyword{models}
\keyword{nonparametric}