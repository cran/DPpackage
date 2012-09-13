\name{LDDProc}
\alias{LDDProc}
\alias{LDDProc.default}

\title{Linear dependent DP model for conditional ROC curve estimation.}
\description{
    This function generates a posterior density sample for a 
    Linear Dependent Dirichlet Process Mixture of Normals model for
    conditional ROC curve estimations.
}
   
\usage{

LDDProc(y.d,z.d,y.nond,z.nond,
        zpred.d,zpred.nond=NULL,prior.d,prior.nond=NULL,
        mcmc,state,status,ngrid=100,
        grid=NULL,compute.band=FALSE,type.band="PD",
        data=sys.frame(sys.parent()),na.action=na.fail,
        work.dir=NULL)

}

\arguments{
   \item{y.d}{       a vector giving the responses for the diseased group.} 

   \item{z.d}{       a matrix giving the design matrix for the diseased group.} 

   \item{y.nond}{    a vector giving the responses for the non-diseased group.} 

   \item{z.nond}{    a matrix giving the design matrix for the non-diseased group.} 

    \item{zpred.d}{    a matrix giving the covariate values where the predictive 
                       density is evaluated for the diseased group.}

    \item{zpred.nond}{ a matrix giving the covariate values where the predictive 
                       density is evaluated for the non-diseased group. By default,
                       zpred.nond=NULL which means that zpred.nond=zpred.d.}
    
    \item{prior.d}{  a list giving the prior information for the diseased group. 
                     The list includes the following
                     parameters: \code{a0} and \code{b0} giving the hyperparameters for
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

    \item{prior.nond}{a list giving the prior information for the non-diseased group. The list
                      includes the same parameters than \code{prior.d}. The default specification 
                      (prior.nond = NULL) uses \code{prior.nond}=\code{prior.d}.} 


    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \code{nburn} giving the number of burn-in 
                     scans, \code{nskip} giving the thinning interval, \code{nsave} giving
                     the total number of scans to be saved, \code{ndisplay} giving
                     the number of saved scans to be displayed on screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out).}   

    \item{state}{    a list giving the current value of the parameters for each 
                     individual model. This list is used
                     if the current analysis is the continuation of a previous analysis.}
    
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state}.}

    \item{ngrid}{    integer giving the number of grid points where the conditional density 
                     estimates are evaluated. The default is 100.}

    \item{grid}{     vector of grid points where the conditional density estimate is 
                     evaluated. The default value is NULL
                     and the grid is chosen according to the range of the
	             data.}

    \item{compute.band}{  logical variable indicating whether the credible band for the conditional  
                          density and mean function must be computed.}       

    \item{type.band}{  string indication the type of credible band to be computed; if equal to "HPD" or
                       "PD" then the 95 percent pointwise HPD or PD band is computed, respectively.}       


    \item{data}{     data frame.}     
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{LDDProc} to print an error message and terminate if there are any
                     incomplete observations.}       

    \item{work.dir}{working directory.} 

}

\details{
  The generic function fits the model described in Inacio et al. (2012) for conditional ROC curve estimation.
  Specifically, the  function fits independent Linear Dependent Dirichlet Process Mixture of Normals models
  for the diseased (i=1) and non-diseased (i=2) groups. The conditional ROC curves are obtained from
  the conditional densities. The model is given by:
  \deqn{y_{ij} | f_{X_{ij}} \sim f_{X_{ij}}}{yij | fXij ~ fXij}
  \deqn{f_{X_{ij}} = \int N(X_{ij} \beta, \sigma^2) G_i(d \beta d \sigma^2)}{fXi = \int N(Xij beta, sigma2) Gi(d beta d sigma2)}
  \deqn{G_i | \alpha_i, G_{0i} \sim DP(\alpha_i G_{0i})}{Gi | alphai, G0i ~ DP(alphai G0i)}
   where, \eqn{G_{0i} = N(\beta | \mu_{bi}, s_{bi})\Gamma(\sigma^2|\tau_{1}/2,\tau_{2i}/2)}{G0 = N(beta| mubi, sbi)Gamma(sigma2|tau1/2,tau2i/2)}.   
   To complete the model specification, independent hyperpriors are assumed,
  \deqn{\alpha_i | a_{0i}, b_{0i} \sim Gamma(a_{0i},b_{0i})}{alphai | a0i, b0i ~ Gamma(a0i,b0i)}
  \deqn{\mu_{bi} | m_{0i}, S_{\beta_{0i}} \sim N(m_{0i},S_{\beta_{0i}})}{mubi | m0i, Sbeta0i ~ N(m0i,Sbeta0i)}
  \deqn{s_{bi} | \nu_i, \Psi_i \sim IW(\nu_i,\Psi_i)}{sbi | nui, psii ~ IW(nui,psii)}
  \deqn{\tau_{2i} | \tau_{s1i}, \tau_{s2i} \sim Gamma(\tau_{s1i}/2,\tau_{s2i}/2)}{tau2i | taus1i, taus2i ~   Gamma(taus1i/2,taus2i/2)}

  The precision or total mass parameters, \eqn{\alpha_i}{alphai}, of the \code{DP} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_{0i},b_{0i})}{Gamma(a0i,b0i)}, 
  or fixed at some particular value. 

  We refer the reader to the help file associated with the \code{\link{LDDPdensity}} function for more
  details about the prior specification, parameterizations and computational strategy.  
}

\value{
  An object of class \code{LDDProc} representing the two LDDP mixture of normals model fits. 
  Generic functions such as \code{print}, \code{plot},
  and \code{summary} have methods to show the results of the fit. The results for each model
  include  \code{mub}, \code{sb}, \code{tau2}, the precision parameter 
  \code{alpha}, and the number of clusters.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. Two different objects are included: \code{state.d} and \code{state.nd}.
  If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case, each of the lists included in \code{state} must include the following objects: 

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

}

\seealso{
\code{\link{LDDPdensity}}
}

\references{
Inacio, V., Jara, A., Hanson, T.E., de Carvalho, M. (2012) Bayesian nonparametric 
   ROC regression modeling with application to diabetes diagnosis. Technical report.
}

\examples{
\dontrun{

    ##############################################################
    # Simulated data example.
    # - Data generated using "perfect" simulation.
    # - one binary predictor.
    # - 250 observations in each 
    #   combination of predictor and 
    #   status.
    ##############################################################

    # Functions required for simulation

      findq <- function(true.cdf,target,low,upp,
                        epsilon=0.0000001)
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


      true.cdf.nond1 <- function(x)
      {
         pnorm(x,2.1,sqrt(0.0324))
      }

      true.cdf.nond2 <- function(x)
      {
	 0.5*pnorm(x,1.85,sqrt(0.005))+
         0.5*pnorm(x,2.25,sqrt(0.005))
      }

      true.cdf.d1 <- function(x)
      {
	 0.5*pnorm(x,1.95,sqrt(0.005))+
         0.5*pnorm(x,2.35,sqrt(0.005))
      }

      true.cdf.d2 <- function(x)
      {
          pnorm(x,2.5,sqrt(0.0324))
      }
	
    # Simulating the data

      nsim <- 250 
      qq <- seq(1,nsim)/(nsim+1)

      y.nond1 <- rep(0,nsim)
      for(i in 1:nsim)
      {
          aa <- findq(true.cdf.nond1,qq[i],
                      low=-6,upp=6,epsilon=0.0000001)
          y.nond1[i] <- aa 
      }	

      y.nond2 <- rep(0,nsim)
      for(i in 1:nsim)
      {
          aa <- findq(true.cdf.nond2,qq[i],
                      low=-6,upp=6,epsilon=0.0000001)
          y.nond2[i] <- aa 
      }	
      y.nond <- c(y.nond1,y.nond2)
      trt.nond <- c(rep(0,nsim),rep(1,nsim))

      y.d1 <- rep(0,nsim)
      for(i in 1:nsim)
      {
          aa <- findq(true.cdf.d1,qq[i],
                      low=-6,upp=6,epsilon=0.0000001)
          y.d1[i] <- aa 
      }	

      y.d2 <- rep(0,nsim)
      for(i in 1:nsim)
      {
          aa <- findq(true.cdf.d2,qq[i],
                      low=-6,upp=6,epsilon=0.0000001)
          y.d2[i] <- aa 
      }	

      y.d <- c(y.d1,y.d2)
      trt.d <- c(rep(0,nsim),rep(1,nsim))

    # Design matrices

      z.d <- cbind(rep(1,2*nsim),trt.d)
      colnames(z.d) <- c("(Intercept)","trt")
      z.nond <- cbind(rep(1,2*nsim),trt.nond)
      colnames(z.nond) <- c("(Intercept)","trt")

    # design matrix for posterior predictive inference 

      zpred <- rbind(c(1,0),c(1,1))  

    # Prior information
      prior <- list(a0=10,
                    b0=1,
                    nu=4,
                    m0=rep(0,2),
                    S0=diag(100,2),
                    psiinv=diag(1,2),
                    tau1=6.01,
                    taus1=6.01,
                    taus2=2.01)

    # Initial state
      state <- NULL

    # MCMC parameters

      nburn <- 5000
      nsave <- 5000
      nskip <- 4
      ndisplay <- 500
      mcmc <- list(nburn=nburn,
                   nsave=nsave,
                   nskip=nskip,
                   ndisplay=ndisplay)

    # Fitting the model

      fit1 <- LDDProc(y.d=y.d,z.d=z.d,
                      y.nond=y.nond,z.nond=z.nond,
                      zpred.d=zpred,
                      prior.d=prior,
                      prior.nond=prior,
                      mcmc=mcmc,
                      state=state, 
                      status=TRUE,
                      compute.band=TRUE)

      fit1
      summary(fit1)
      plot(fit1)


     # Ploting the conditional
     # ROC curve for x=c(1,0),
     # along with the true curve

       par(cex=1.7,mar=c(4.1, 4.1, 1, 1))

       plot(fit1$rocgrid,fit1$rocp.h[1,],type="l",
            lty=2,lwd=2,ylim=c(0,1),xlim=c(0,1),
            xlab="False positive rate",
            ylab="True positive rate")
       lines(fit1$rocgrid,fit1$rocp.l[1,],lty=2,lwd=2)
       lines(fit1$rocgrid,fit1$rocp.m[1,],lty=1,lwd=2)

       nn <- length(fit1$rocgrid)
       tt <- rep(0,nn)
       for(i in 1:nn)
       {
	   tt[i] <- findq(true.cdf.nond1,
                          1-fit1$rocgrid[i],
                          low=-6,upp=6,
                          epsilon=0.0000001)
       }	
       true.roc1 <- 1.0 - true.cdf.d1(tt) 
       lines(fit1$rocgrid,true.roc1,
             lty=1,lwd=3,col="red")

     # Ploting the conditional
     # ROC curve for x=c(1,1),
     # along with the true curve

       par(cex=1.7,mar=c(4.1, 4.1, 1, 1))

       plot(fit1$rocgrid,fit1$rocp.h[2,],type="l",
            lty=2,lwd=2,ylim=c(0,1),xlim=c(0,1),
            xlab="False positive rate", 
            ylab="True positive rate")
       lines(fit1$rocgrid,fit1$rocp.l[2,],lty=2,lwd=2)
       lines(fit1$rocgrid,fit1$rocp.m[2,],lty=1,lwd=2)

       nn <- length(fit1$rocgrid)
       tt <- rep(0,nn)
       for(i in 1:nn)
       {
           tt[i] <- findq(true.cdf.nond2,
                          1-fit1$rocgrid[i],
                          low=-6,upp=6,
                          epsilon=0.0000001)
       }	
       true.roc2 <- 1.0 - true.cdf.d2(tt) 
       lines(fit1$rocgrid,true.roc2,lty=1,lwd=3,col="red")

}
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}

}

\keyword{models}
\keyword{nonparametric}
