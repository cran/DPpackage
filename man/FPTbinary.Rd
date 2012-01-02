\name{FPTbinary}
\alias{FPTbinary}
\alias{FPTbinary.default}

\title{Bayesian analysis for a Finite Polya Tree Bernoulli regression model}
\description{
    This function generates a posterior density sample
    for a binary regression model using a Finite Polya tree
    prior for the link function. 
}
    
\usage{

FPTbinary(formula,baseline="logistic",prior,mcmc,state,
          status,misc=NULL,data=sys.frame(sys.parent()),
          na.action=na.fail)
}

\arguments{
   \item{formula}{   a two-sided linear formula object describing the
                     model fit, with the response on the
                     left of a \code{~} operator and the terms, separated by \code{+}
                     operators, on the right.} 

   \item{baseline}{  a description of the baseline error distribution to
                     be used in the model. The baseline distributions considered by 
                     \code{FPTbinary} so far is \code{logistic}.} 

   \item{prior}{     a list giving the prior information. The list includes the following
                     parameters: \code{a0} and \code{b0} giving the hyperparameters for
                     prior distribution of the precision parameter of the Finite Polya
                     Tree prior, \code{alpha} giving the value of the precision parameter (it 
                     must be specified if \code{a0} and \code{b0} are missing), 
                     \code{beta0} and \code{Sbeta0} giving the 
                     hyperparameters of the normal prior distribution for the regression
                     coefficients, \code{M} giving the finite level
                     to be considered for the Finite Polya tree.}
               
    \item{mcmc}{     a list giving the MCMC parameters. The list must include
                     the following integers: \code{nburn} giving the number of burn-in 
                     scans, \code{nskip} giving the thinning interval, \code{nsave} giving
                     the total number of scans to be saved, \code{ndisplay} giving
                     the number of saved scans to be displayed on the screen (the function reports 
                     on the screen when every \code{ndisplay} iterations have been carried
                     out), and \code{tune1} and \code{tune2} giving the Metropolis tuning 
                     parameters for the regression coefficients and precision parameter,
                     respectively (the default value is 1.1).}   

    \item{state}{    a list giving the current value of the parameters. This list is used
                     if the current analysis is the continuation of a previous analysis.}

    \item{status}{   a logical variable indicating whether this run is new (\code{TRUE}) or the 
                     continuation of a previous analysis (\code{FALSE}). In the latter case
                     the current value of the parameters must be specified in the 
                     object \code{state}.}

    \item{misc}{     misclassification information. When used, this list must include
                     two objects, \code{sens} and \code{spec}, giving the sensitivity and
                     specificity, respectively. Both can be a vector or a scalar. 
                     This information is used to correct for misclassification in the
                     conditional bernoulli model.}

    \item{data}{     data frame.}       
    
    \item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes 
                     \code{FPTbinary} to print an error message and terminate if there are any
                     incomplete observations.}       

}

\details{
  This generic function fits a semiparametric binary regression model using
  a Finite Polya tree prior (FPT) for the link function 
  (see, Hanson, 2006; Jara, Garcia-Zattera and Lesaffre, 2006):
  \deqn{ y_i = I(V_i \leq X_i \beta), i=1,\ldots,n}{yi = I(Vi <= Xi \beta),\ i=1,\ldots,n} 
  \deqn{V_1,\ldots,V_n | G \sim G}{V1,\ldots,Vn | G ~ G}
  \deqn{G | \alpha \sim FPT^M(\Pi,\textit{A})}{G | alpha ~ FPT^M(Pi,A)}

  where, the FPT is centered around a \eqn{Logistic(0,1)}{Logistic(0,1)} distribution if the baseline is \code{logistic}, by
  taking each \eqn{m} level of the partition \eqn{\Pi}{Pi} to coincide 
  with the \eqn{k/2^m, k=0,\ldots,2^m} quantile of the \eqn{Logistic(0,1}{Logistic(0,1)} distribution.
  The family \eqn{\textit{A}=\{\alpha_e: e \in E^{*}\}}{A=\{alphae: e \in E^{*}\}}, 
  where \eqn{E^{*}=\bigcup_{m=1}^{M} E^m}{E^{*}=\bigcup_{m=1}^{M} E^m} 
  and \eqn{E^m} is the \eqn{m}-fold product of \eqn{E=\{0,1\}},
  is specified as \eqn{\alpha_{e_1 \ldots e_m}=\alpha m^2}{alpha{e1 \ldots em}=\alpha m^2}. 
  To complete the model specification, independent hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim Gamma(a_0,b_0)}{alpha | a0, b0 ~ Gamma(a0,b0)}
  \deqn{\beta | \beta_0, S_{\beta_0} \sim N(\beta_0,S_{\beta_0})}{\beta | beta0, Sbeta0 ~ N(beta0,Sbeta0)}

  The precision parameter, \eqn{\alpha}{alpha}, of the \code{FPT} prior 
  can be considered as random, having a \code{gamma} distribution, \eqn{Gamma(a_0,b_0)}{Gamma(a0,b0)}, 
  or fixed at some particular value. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value, set \eqn{a_0}{a0} to NULL in the prior specification.

  In the computational implementation of the model, Metropolis-Hastings steps are used to 
  sample the posterior distribution of the regression coefficients and the precision
  parameter, as described in Hanson (2006), and Jara, Garcia-Zattera and Lesaffre (2006).
}


\value{
  An object of class \code{FPTbinary} representing the semiparametric logistic regression
  model fit. Generic functions such as \code{print}, \code{plot}, \code{predict}, \code{summary},
  and \code{anova} have methods to show the results of the fit. 
  The results include \code{beta}, the precision parameter (\code{alpha}), and the \code{link} function.
  
  The MCMC samples of the parameters and the errors in the model are stored in the object 
  \code{thetasave} and \code{randsave}, respectively. Both objects are included in the 
  list \code{save.state} and are matrices which can be analyzed directly by functions 
  provided by the coda package.
  
  The list \code{state} in the output object contains the current value of the parameters 
  necessary to restart the analysis. If you want to specify different starting values 
  to run multiple chains set \code{status=TRUE} and create the list state based on 
  this starting values. In this case the list \code{state} must include the following objects:
  
  \item{beta}{ giving the value of the regression coefficients.} 

  \item{v}{ giving the value of the errors (it must be consistent with \code{yi = I(Vi < xi beta)}.}, 

  \item{y}{ giving the value of the true response binary variable (only if the model
            considers correction for misclassification).} 
  
  \item{alpha}{ giving the value of the precision parameter.}
}

\references{

Hanson, T. (2006) Inference for Mixtures of Finite Polya tree models. 
  Journal of the American Statistical Association, 101: 1548-1565. 

Jara, A., Garcia-Zattera, M.J., Lesaffre, E. (2006) Semiparametric Bayesian
  Analysis of Misclassified Binary Data. XXIII International Biometric   
  Conference, July 16-21, Montreal, Canada. 

Lavine, M. (1992) Some aspects of Polya tree distributions for statistical 
  modelling. The Annals of Statistics, 20: 1222-11235.

Lavine, M. (1994) More aspects of Polya tree distributions for statistical 
  modelling. The Annals of Statistics, 22: 1161-1176.

}

\examples{
\dontrun{

    # Prostate cancer data example
      data(nodal)
      attach(nodal)
      lacid<-log(acid)

    # Initial state
      state <- NULL

    # MCMC parameters
      nburn<-20000
      nsave<-10000
      nskip<-10
      ndisplay<-100
      mcmc <- list(nburn=nburn,nsave=nsave,
                   nskip=nskip,ndisplay=ndisplay,
                   tune1=1.1,tune2=1.1)

    # Prior distribution
      prior <- list(alpha=1, beta0=c(0,rep(0.75,5)),
                    Sbeta0=diag(c(100,rep(25,5)),6),M=5)


    # Fitting the Finite Polya tree model
      fit1 <- FPTbinary(ssln~age+lacid+xray+size+grade,
                        prior=prior,mcmc=mcmc,
                        state=state,status=TRUE) 
      fit1

    # Summary with HPD and Credibility intervals
      summary(fit1)
      summary(fit1,hpd=FALSE)

    # Plot model parameters (to see the plots gradually set ask=TRUE)
      plot(fit1)
      plot(fit1,nfigr=2,nfigc=2)	

    # Plot an specific model parameter (to see the plots gradually 
    # set ask=TRUE)
      plot(fit1,ask=FALSE,nfigr=1,nfigc=2,param="xray")	
      plot(fit1,ask=FALSE,param="link",nfigc=1,nfigr=1)

    # Table of Pseudo Contour Probabilities
      anova(fit1)


    # Fitting parametric models

      nburn<-20000
      nsave<-10000
      nskip<-10
      ndisplay<-100
      mcmc <- list(nburn=nburn,nsave=nsave,
                   nskip=nskip,ndisplay=ndisplay,
                   tune=1.1)

      fit2 <- Pbinary(ssln~age+lacid+xray+size+grade,link="probit",
                      prior=prior,mcmc=mcmc,state=state,status=TRUE)  
                    
      fit3 <- Pbinary(ssln~age+lacid+xray+size+grade,link="logit",
                      prior=prior,mcmc=mcmc,state=state,status=TRUE)  

    
    # Model comparison

      DPpsBF(fit1,fit2,fit3)
      
}      
}

\author{
Alejandro Jara \email{<atjara@uc.cl>}

Tim Hanson \email{<hansont@stat.sc.edu>}
}

\keyword{models}
\keyword{nonparametric}