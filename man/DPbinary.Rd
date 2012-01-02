\name{DPbinary}
\alias{DPbinary}
\alias{DPbinary.default}
\encoding{latin1}

\title{Bayesian analysis for a semiparametric Bernoulli regression model}
\description{
    This function generates a posterior density sample
    for a semiparametric binary regression model. 
}
    
\usage{

DPbinary(formula,inter=TRUE,baseline="logistic",prior,mcmc,state,
         status,misc=NULL,data=sys.frame(sys.parent()),
         na.action=na.fail)
}

\arguments{
   \item{formula}{   a two-sided linear formula object describing the
                     model fit, with the response on the
                     left of a \code{~} operator and the terms, separated by \code{+}
                     operators, on the right.} 
                     
   \item{inter}{     a logical variable indicating whether the intercept should be
                     explicitly included in the linear predictor \code{TRUE} or absorved
                     in the unknown link function \code{FALSE}. The default option
                     is \code{TRUE}. Note that if inter is \code{FALSE} the dimension
                     of the prior distribution for the regression coefficients should
                     not consider the intercept in the model.} 

   \item{baseline}{  a description of the baseline error distribution to
                     be used in the model. The baseline distributions considered by 
                     \code{DPbinary} so far are \emph{logistic} (default), \emph{normal},
                     and \emph{cauchy}.} 

   \item{prior}{     a list giving the prior information. The list includes the following
                     parameters: \emph{a0} and \emph{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Dirichlet process
                     prior, \emph{alpha} giving the value of the precision parameter (it 
                     must be specified if \emph{a0} and \emph{b0} are missing, see the details
                     below), and \emph{beta0} and \emph{Sbeta0} giving the 
                     hyperparameters of the normal prior distribution for the regression
                     coefficients.}
               
    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \emph{nburn} giving the number of burn-in 
                     scans, \emph{nskip} giving the thinning interval, \emph{nsave} giving
                     the total number of scans to be saved, \emph{ndisplay} giving
                     the number of saved scans to be displayed on the screen (the function reports 
                     on the screen when every \emph{ndisplay} iterations have been carried
                     out), and \emph{tune} giving the Metropolis tuning parameter
                     (the default value is 1.1).}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}

    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \emph{state}.}

    \item{misc}{     misclassification information. When used, this list must include
                     two objects, \emph{sens} and \emph{spec}, giving the sensitivity and
                     specificity, respectively. Both can be a vector or a scalar. 
                     This information is used to correct for misclassification in the
                     conditional bernoulli model.}

    \item{data}{     data frame.}       
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{DPbinary} to print an error message and terminate if there are any
                     incomplete observations.}       

}

\details{
  This generic function fits a semiparametric binary regression model using
  a Dirichlet process prior (see, Jara, Garcia-Zattera and Lesaffre, 2006):
  \deqn{ y_i = I(V_i \leq X_i \beta), i=1,\ldots,n}{yi = I(Vi <= Xi \beta),\ i=1,\ldots,n} 
  \deqn{V_1,\ldots,V_n | G \sim G}{V1,\ldots,Vn | G ~ G}
  \deqn{G | \alpha, G_0 \sim DP(\alpha G_0)}{G | alpha, G0 ~ DP(alphaG0)}
  where, \eqn{G_0 = Logistic(V| 0, 1)}{G0 = Logistic(V| 0, 1)} if the baseline is logistic,
  \eqn{G_0 = N(V| 0, 1)}{G0 = N(V| 0, 1)} if the baseline is normal, and 
  \eqn{G_0 = Cauchy(V| 0, 1)}{G0 = Cauchy(V| 0, 1)} if the baseline is cauchy. 
  To complete the model specification, the following
  prior distributions are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\beta | \beta_0, S_{\beta_0} \sim N(\beta_0,S_{\beta_0})}{\beta | beta0, Sbeta0 ~ N(beta0,Sbeta0)}
  The precision or total mass parameter, \eqn{\alpha}{alpha}, of the \emph{DP} prior 
  can be considered as random, having a \emph{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. When \eqn{\alpha}{alpha} is random the method described by
  Escobar and West (1995) is used. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value, set \eqn{a_0} to NULL in the prior specification.

  A Metropolis-Hastings step is used to sample the fully conditional distribution
  of the regression coefficients and errors (see, Jara, Garcia-Zattera and 
  Lesaffre, 2006). In the computational implementation of the model, G is 
  considered as latent data and sampled partially with sufficient accuracy to be 
  able to generate \eqn{V_1,\ldots,V_{n+1}}{V1,\ldots,Vn+1} which are exactly iid G, as proposed by Doss (1994). 
  Both Ferguson's definition of DP and the Sethuraman-Tiwari (1982) representation 
  of the process are used, as described in Jara, Garcia-Zattera and Lesaffre (2006).
  An extra step which moves the clusters in such a way that the posterior distribution
  is still a stationary distribution, is performed in order to improve the rate of
  mixing.
}


\value{
  An object of class \code{DPbinary} representing the semiparametric logistic regression
  model fit. Generic functions such as \code{print}, \code{plot}, \code{predict}, \code{summary},
  and \code{anova} have methods to show the results of the fit. The results include \code{beta}, 
  the precision parameter (\emph{alpha}), the number of clusters (\emph{ncluster}), and the \emph{link} 
  function.
  
  The MCMC samples of the parameters and the errors in the model are stored in the object 
  \emph{thetasave} and \emph{randsave}, respectively. Both objects are included in the 
  list \emph{save.state} and are matrices which can be analyzed directly by functions 
  provided by the coda package.
  
  The list \emph{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \emph{status=TRUE} and create the list state based on 
  this starting values. In this case the list \emph{state} must include the following objects:
  
  \item{beta}{ giving the value of the regression coefficients.} 

  \item{v}{ giving the value of the errors (it must be consistent with \emph{yi = I(Vi < xi beta)}.}, 
  
  \item{y}{ giving the value of the true response binary variable (only if the model
  considers correction for misclassification).} 
  
  \item{alpha}{ giving the value of the precision parameter.}
}

\references{

Doss, H. (1994) Bayesian nonparametric estimation for incomplete data using 
  mixtures of Dirichlet priors. The Annals of Statistics, 22: 1763 - 1786.

Escobar, M.D. and West, M. (1995) Bayesian Density Estimation and Inference 
  Using Mixtures. Journal of the American Statistical Association, 90: 577-588.

Jara, A., Garcia-Zattera, M.J., Lesaffre, E. (2006) Semiparametric Bayesian
  Analysis of Misclassified Binary Data. XXIII International Biometric Conference,
  July 16-21, Montréal, Canada.
 
Sethuraman, J., and Tiwari, R. C. (1982) Convergence of Dirichlet Measures and 
  the Interpretation of their Parameter, in Statistical Decision Theory and Related 
  Topics III (vol. 2), eds. S. S. Gupta and J. O. Berger, New York: Academic Press, 
  pp. 305 - 315.  
}

\examples{
\dontrun{
    # Bioassay Data Example
    # Cox, D.R. and Snell, E.J. (1989). Analysis of Binary Data. 2nd ed. 
    # Chapman and Hall. p. 7  
    # In this example there are 150 subjects at 5 different stimulus levels, 
    # 30 at each level.

      y<-c(rep(0,30-2),rep(1,2),
           rep(0,30-8),rep(1,8),
           rep(0,30-15),rep(1,15),
           rep(0,30-23),rep(1,23),
           rep(0,30-27),rep(1,27))

      x<-c(rep(0,30),
           rep(1,30),
           rep(2,30),
           rep(3,30),
           rep(4,30))

    # Initial state
      state <- NULL

    # MCMC parameters
      nburn<-5000
      nsave<-10000
      nskip<-10
      ndisplay<-100
      mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay,
                   tune=1.1)


    # Prior distribution
      prior <- list(a0=2,b0=1,beta0=rep(0,2), Sbeta0=diag(10000,2))

    # Fit the model
      fit1 <- DPbinary(y~x,prior=prior,mcmc=mcmc,state=state,status=TRUE) 
      fit1

    # Summary with HPD and Credibility intervals
      summary(fit1)
      summary(fit1,hpd=FALSE)

    # Plot model parameters (to see the plots gradually set ask=TRUE)
      plot(fit1)
      plot(fit1,nfigr=2,nfigc=2)	

    # Plot an specific model parameter (to see the plots gradually 
    # set ask=TRUE)
      plot(fit1,ask=FALSE,nfigr=1,nfigc=2,param="x")	
      plot(fit1,ask=FALSE,nfigr=1,nfigc=2,param="ncluster")	
      plot(fit1,ask=FALSE,param="link",nfigc=1,nfigr=1)

    # Table of Pseudo Contour Probabilities
      anova(fit1)

    # Predictive Distribution
      npred<-40  
      xnew<-cbind(rep(1,npred),seq(0,4,length=npred))

      pp<-predict(fit1,xnew)       

      plot(seq(0,4,length=npred),pp$pmean,type='l',ylim=c(0,1),
           xlab="log2(concentration)",ylab="Probability")
      
    # Adding MLE estimates
      points(c(0,1,2,3,4),c(0.067,0.267,0.500,0.767,0.900),col="red")

}      
}

\author{

Alejandro Jara \email{<atjara@uc.cl>}
}

\keyword{models}
\keyword{nonparametric}
