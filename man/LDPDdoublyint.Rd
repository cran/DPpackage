\name{LDPDdoublyint}
\alias{LDPDdoublyint}
\alias{LDPDdoublyint.default}

\title{Bayesian analysis for a Linear Dependent Poisson Dirichlet Process Mixture Models for the Analysis of Doubly-Interval-Censored Data}
\description{
    This function generates a posterior density sample for a 
    linear dependent Poisson Dirichlet process mixture model for
    the analysis of doubly-invertval-censored survival data.
}
   
\usage{

LDPDdoublyint(onset,failure,p,xonset,q,xfailure,xpred,
              grid,prior,mcmc,state,status,work.dir=NULL)

}

\arguments{
   \item{onset}{     a matrix giving the interval limits for the onset times.
	             For multiple variables the limits must be included in a sequential manner 
                     for each variable, i.e., (L1,U1,L2,U2,...).} 

    \item{failure}{  a matrix giving the interval limits for the time-to-event times.
	             For multiple variables the limits must be included in a sequential manner 
                     for each variable, i.e. (L1,U1,L2,U2,...).} 

    \item{p}{        an integer giving the number of predictors included for each onset variable.
                     Different design matrices are allowed for the different responses but
                     of the same p-dimension.}
    
    \item{xonset}{   a matrix giving the design matrices for each onset time. For multiple
                     variables the design matrices must be included in order and includes
                     the intercepts, i.e. (XO1,XO2,...).}

    \item{q}{        an integer giving the number of predictors included for each time-to-event variable.
                     Different design matrices are allowed for the different responses but
                     of the same q-dimension.}

    \item{xfailure}{ a matrix giving the design matrices for each time-to-event variable. For multiple
                     variables the design matrices must be included in order and includes
                     the intercepts, i.e. (XT1,XT2,...).}

    \item{xpred}{    a matrix giving the value of the predictors for which survival and
                     hazard functions will be evaluated. The number of columns of \code{xpred} should be
                     (p+q)*nvar/2 where nvar is the number of variables. The design matrices 
                     for the predictions must include onset predictors first and then time-to-event
                     predictors, i.e., (XO1,XO2,...,XT1,XT2,...).}

    \item{grid}{     a matrix including the grids where survival and hazard functions are evaluated.
                     Each row must include the grid points for different variable.}


    \item{prior}{    a list giving the prior information. The list includes the following
                     parameter: \code{q}, \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the \code{a} precision parameter of the Poisson-Dirichlet process
                     prior, \code{mub} and \code{sigmab} giving the hyperparameters for the prior
                     distributions of the \code{b} precision parameter of the Poisson-Dirichlet process
                     prior, \code{nu} and \code{tinv} giving the hyperparameters of the 
                     inverted Wishart prior distribution for the kernel covariance matrix, 
		     \code{mb} and \code{Sb} giving the hyperparameters of the normal prior distribution 
                     for the mean of the normal baseline distribution, 
	             \code{nub} and \code{tbinv} giving the hyperparameters of the 
                     inverted Wishart prior distribution for the for the covariance matrix
                     of the normal baseline distribution, and \code{maxm} giving the finite truncation
                     limit of the sitck-breaking representation of the Poisson-Dirichlet process.}

    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \code{nburn} giving the number of burn-in 
                     scans, \code{nskip} giving the thinning interval, \code{nsave} giving
                     the total number of scans to be saved, \code{ndisplay} giving
                     the number of saved scans to be displayed on screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out), and \code{tune1} and \code{tune2} the parameters needed for the
                     MH candidate generating distribution for the precision parameters of
                     the Poisson-Dirichlet process prior.}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}
    
    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state}.}

    \item{work.dir}{ working directory.} 

}

\details{
  This generic function fits a Linear Dependent Poisson-Dirichlet Process Mixture of Survival models
  as in Jara et al. (2010). The joint distribution of the true chronological onset times and true time-to-events, \eqn{T_i}{Ti},
  is modeled using the mixture model
  \deqn{log(T_i)=y_i | f_{X_i} \sim f_{X_i}}{log(Ti)=yi | fXi ~ fXi}
  \deqn{f_{X_i} = \int N(X_i \beta, \Sigma) dG(\beta)}{fXi = \int N(Xi beta, Sigma) dG(beta)}
  \deqn{G | a, b, G_0 \sim PD(a,b,G_0)}{G | a, b, G0 ~ PD(a, b, G0)}
   where, \eqn{G_0 = N(\beta| m_b, S_b)}{G0 = N(beta | mb, Sb)}.
   To complete the model specification, independent hyperpriors are assumed,
   \deqn{\Sigma | \nu, T \sim IW(\nu,T)}{Sigma | nu, T ~ IW(nu,T)}
   \deqn{a | q, a_0, b_0 \sim q \delta_0 + (1-q)Beta(a_0,b_0)}{a | q, a0, b0 ~ q delta0 + (1-q) Beta(a0,b0)}
   \deqn{b | a, \mu_b, \sigma_b \sim N(\mu_b,\sigma_b) I(-a,\infty)}{b | a, mub, sigmab ~ N(mub,sigmab)I(-a,infty)}
   \deqn{m_b | m_0, S_0 \sim N(m_0,S_0)}{mb | m0, S0 ~ N(m0,S0)}
   \deqn{S_b | \nu_b, Tb \sim IW(\nu_b,T_b)}{Sb | nub, Tb ~ IW(nub,Tb)}
  
   Note that the inverted-Wishart prior is parametrized such that if \eqn{A \sim IW_q(\nu, \psi)}{A ~ IWq(nu, psi)} then
   \eqn{E(A)= \psi^{-1}/(\nu-q-1)}{E(A)= psiinv/(nu-q-1)}.

  The precision parameters are updated using a MH step. The candidate generating distribution is of the form
   \deqn{a' | a, tune_1 \sim 0.5 \delta_0 + 0.5 N(a,tune_1^2)}{a' | a, tune1 ~ 0.5 delta0 + 0.5 N(a,tune1**2)}
   \deqn{b' | b, a', tune_2 \sim N(b,tune_2^2) I(-a',\infty)}{b' | b, a', tune2 ~ N(b,tune2**2)I(-a',infty)}

  The computational implementation of the model is based on the \code{maxm} truncation of stick-breaking representation of 
  the model (see, Jara et al., 2009). 
}

\value{
  An object of class \code{LDPDdoublyint} representing the LDPD mixture of survival models fit. 
  Generic functions such as \code{print}, \code{plot},
  \code{summary}, and \code{predict} have methods to show the results of the fit. The results include 
  \code{mb}, \code{Sb}, \code{sigma}, the precision parameter 
  \code{alpha}, and the number of clusters.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects: 

  \item{alpha}{ a vector giving the value of the precision parameters a and b.} 

  \item{b}{       a matrix of dimension maxm times the 
                  number of columns in the design matrix ((p+q)*nvar/2), giving the 
                  regression coefficients for each cluster.}

 \item{sigma}{ a matrix of dimension nvar times nvar giving the kernel covariance matrix.}

  \item{mb}{ giving the mean of the normal baseline distributions.}

  \item{Sb}{ giving the covariance matrix the normal baseline distributions.}

  \item{ncluster}{ an integer giving the number of clusters.} 

  \item{ss}{ an interger vector defining to which of the \code{ncluster} clusters each subject belongs.}

  \item{y}{ a matrix of dimension nsubject times nvar giving the value of the imputed log-survival times.}
}

\seealso{
\code{\link{LDDPdensity}}, \code{\link{LDDPsurvival}}
}

\references{
Jara, A., Lesaffre, E., De Iorio, M., Quintana, F. (2010). Bayesian semiparametric inference for 
   multivariate doubly-interval-censored data. Annals of Applied Statistics, 4: 2126-2149.
}

\examples{
\dontrun{

###############################################################
# HIV-AIDS Data
###############################################################

  data(hiv)
  attach(hiv)

###############################################################
# Working folder
###############################################################

  work.dir <- NULL

###############################################################
# Response and design matrices
###############################################################

  nsubject <- length(onsetL)
  onset <- cbind(onsetL,onsetU)
  failure <- cbind(failureL,failureU)

  intercept <- rep(1,nsubject)
  p <- 2
  xonset <- cbind(IntO=intercept,trtO=trt)
  q <- 2
  xfailure <- cbind(IntF=intercept,trtF=trt)

###############################################################
# Predictions
###############################################################

   grid <- matrix(c(rep(seq(0,30,len=30),1),
                    rep(seq(0,20,len=30),1)),nrow=2,byrow=T)

   xpred <- matrix(c(rep(c(1,0),1),rep(c(1,0),1),
                     rep(c(1,1),1),rep(c(1,1),1)),
                     nrow=2,byrow=T)

   colnames(xpred) <- colnames(cbind(xonset,xfailure))

###############################################################
# Initial state
###############################################################

  state <- NULL

###############################################################
# Prior
###############################################################
  
  prior<-list(a0=1,
              b0=1,
              q=0.5,
              mub=10,
              sigmab=200,
              nu=4,
              tinv=diag(1,2),
              nub=6,
              tbinv=diag(1,4),
              m0=rep(0,4),
              S0=diag(100,4),
              maxm=40)

###############################################################
# MCMC
###############################################################
  
  mcmc<-list(nburn=5000,
             nskip=9,
             ndisplay=100,
             nsave=5000,
             tune1=0.25,
             tune2=1)

###############################################################
# Fitting the Model
###############################################################

  fit1 <- LDPDdoublyint(onset=onset,failure=failure,
                        p=p,xonset=xonset,
                        q=q,xfailure=xfailure,
                        xpred=xpred,grid=grid,
                        prior=prior,
                        mcmc=mcmc,
                        state=state,
                        status=TRUE,
                        work.dir=work.dir)
 
  fit1
  summary(fit1)

###############################################################
# Getting Information for Predictions
###############################################################
  
# Without CI bands and intervals

  fit1.pred <- predict(fit1)
  fit1.pred
  plot(fit1.pred)

# With CI bands and intervals
 
  fit1.pred <- predict(fit1,compute.band=TRUE)
  fit1.pred
  plot(fit1.pred)
}
}

\author{
Alejandro Jara \email{<atjara@uc.cl>}
}

\keyword{models}
\keyword{nonparametric}
