### FPTbinary.R
### Fit a bernoulli regression model using
### a simple Finite Polya tree prior for the link function.
###
### Copyright: Alejandro Jara and Tim Hanson, 2006
### Last modification: 16-08-2006.
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or (at
### your option) any later version.
###
### This program is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
###
### The authors's contact information:
###
###      Alejandro Jara Vallejos
###      Biostatistical Centre
###      Katholieke Universiteit Leuven
###      U.Z. Sint-Rafaël
###      Kapucijnenvoer 35
###      B-3000 Leuven
###      Voice: +32 (0)16 336892  URL  : http://student.kuleuven.be/~s0166452/
###      Fax  : +32 (0)16 337015  Email: Alejandro.JaraVallejos@med.kuleuven.be
###
###      Tim Hanson
###      Division of Biostatistics
###      University of Minnesota
###      School of Public Health
###      A460 Mayo Building, 
###      MMC 303
###      420 Delaware St SE
###      Minneapolis, MN 55455
###      Voice: 612-626-7075  URL  : http://www.biostat.umn.edu/~hanson/
###      Fax  : 612-626-0660  Email: hanson@biostat.umn.edu
###

"FPTbinary"<-
function(formula,baseline="logistic",prior,mcmc,state,status,misc=NULL,
data=sys.frame(sys.parent()),na.action=na.fail) 
UseMethod("FPTbinary")

"FPTbinary.default"<-
function(formula,
         baseline="logistic",
         prior,
         mcmc,
         state,
         status,
         misc=NULL,
         data,
         na.action) 
{
         #########################################################################################
         # call parameters
         #########################################################################################

	 cl <- match.call()
	 mf <- match.call(expand.dots = FALSE)
	 m <- match(c("formula", "data","na.action"), names(mf), 0)
	 mf <- mf[c(1, m)]
	 mf$drop.unused.levels <- TRUE
	 mf[[1]] <- as.name("model.frame")
	 mf <- eval(mf, parent.frame())

         #########################################################################################
         # data structure
         #########################################################################################

  	 yobs<- model.response(mf,"numeric")
	 nrec<-length(yobs)
	 x<-as.matrix(model.matrix(formula))
	 p<-dim(x)[2]
	 
         #########################################################################################
         # misclassification
         #########################################################################################
 
 	 if(is.null(misc))
	 {
	  	sens<-rep(1,nrec)
	  	spec<-rep(1,nrec)
	  	model<-0
	 }
	 else
	 {
 		sens<-misc$sens
	 	spec<-misc$spec
		if(length(sens)==1)sens<-rep(sens,nrec)
		if(length(spec)==1)spec<-rep(spec,nrec)
		model<-1
	 }

         #########################################################################################
         # mcmc specification
         #########################################################################################

         MLElogit<-function(x,y,sens,spec)
         {
   	     fn<-function(theta)
	     {
		eta<-x%*%theta
                p<-plogis(eta)
		like <- sens*p+(1-spec)*(1-p)
		
                if (all(like > 0)) 
                     eval<- -sum(log(like[y==1]))-sum(log(1-like[y==0]))
                else eval<-Inf
		return(eval)
	     }
	     
	     start<-coefficients(glm(y~x-1,family=binomial(logit)))
	
	     foo<-optim(start,fn=fn,method="BFGS",hessian=TRUE)

	     out<-NULL
	     out$beta<-foo$par
	     out$stderr<-sqrt(diag(-solve(-foo$hessian)))
	     out$covb<-(-solve(-foo$hessian))
	     return(out)
         }
         

         MLEprobit<-function(x,y,sens,spec)
         {
   	     fn<-function(theta)
	     {
		eta<-x%*%theta
                p<-pnorm(eta)
		like <- sens*p+(1-spec)*(1-p)
		
                if (all(like > 0)) 
                     eval<- -sum(log(like[y==1]))-sum(log(1-like[y==0]))
                else eval<-Inf
		return(eval)
	     }
	     
	     start<-coefficients(glm(y~x-1,family=binomial(logit)))
	
	     foo<-optim(start,fn=fn,method="BFGS",hessian=TRUE)

	     out<-NULL
	     out$beta<-foo$par
	     out$stderr<-sqrt(diag(-solve(-foo$hessian)))
	     out$covb<-(-solve(-foo$hessian))
	     return(out)
         }


         MLEcauchy<-function(x,y,sens,spec)
         {
   	     fn<-function(theta)
	     {
		eta<-x%*%theta
                p<-pcauchy(eta)
		like <- sens*p+(1-spec)*(1-p)
		
                if (all(like > 0)) 
                     eval<- -sum(log(like[y==1]))-sum(log(1-like[y==0]))
                else eval<-Inf
		return(eval)
	     }
	     
	     start<-coefficients(glm(y~x-1,family=binomial(logit)))
	
	     foo<-optim(start,fn=fn,method="BFGS",hessian=TRUE)

	     out<-NULL
	     out$beta<-foo$par
	     out$stderr<-sqrt(diag(-solve(-foo$hessian)))
	     out$covb<-(-solve(-foo$hessian))
	     return(out)
         }
         
         mcmcvec<-c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
         nsave<-mcmc$nsave
          
         xmatrix<-x

         if(baseline=="logistic")
         {
            fit0<- MLElogit(xmatrix,yobs,sens,spec)
         }
         if(baseline=="normal")
         {
  	    fit0<- MLEprobit(xmatrix,yobs,sens,spec)
	 }
         if(baseline=="cauchy")
         {
	    fit0<- MLEcauchy(xmatrix,yobs,sens,spec)
	 }
	 
         propv<- fit0$covb
         
         if(is.null(mcmc$tune1))
         {
            tune1=1.1
         }
         else
         {
            tune1<-mcmc$tune1
         }
        
         if(is.null(mcmc$tune2))
         {
            tune2=1.1
         }
         else
         {
            tune2<-mcmc$tune2
         }


         resop<-1

         #########################################################################################
         # prior information
         #########################################################################################

	 if(is.null(prior$a0))
	 {
		a0b0<-c(-1,-1)
		alpha<-prior$alpha
	 }
	 else
	 {
	 	a0b0<-c(prior$a0,prior$b0)
	 	alpha<-rgamma(1,prior$a0,prior$b0)
	 }
	 
  	 betapm<-prior$beta0
	 betapv<-prior$Sbeta0
	
	 propv<-diag(tune1,p)%*%solve(solve(betapv)+solve(propv))%*%diag(tune1,p)
	 
	 nlevel<-prior$M

         #########################################################################################
         # parameters depending on status
         #########################################################################################

       	 if(status)
	 {
                beta<-fit0$beta
		eta <- x %*% beta
		v <- rep(0,nrec)
		v[yobs==1]<-eta[yobs==1]-0.5
		v[yobs==0]<-eta[yobs==0]+0.5
		y<-yobs		
	 }	
      	 else
	 {
		beta<-state$beta
		eta <- x %*% beta
		v<-state$v
		y<-state$y		
		alpha<-state$alpha
	 }
	 

         #########################################################################################
         # output
         #########################################################################################
         
         ninter<-2**nlevel
         
         xlink=seq(-6,6,length=34)
         nlink<-length(xlink)
         fsave     <- matrix(0, nrow=mcmc$nsave, ncol=nlink)         
	 thetasave <- matrix(0, nrow=mcmc$nsave, ncol=p+1)
	 randsave  <- matrix(0, nrow=mcmc$nsave, ncol=nrec+1)
	 ppsave    <- matrix(0, nrow=mcmc$nsave, ncol=ninter)


         #########################################################################################
         # working space
         #########################################################################################

	 ninter<-2**nlevel

         acrate<-rep(0,2)
         assign<-matrix(0,nrow=nrec,ncol=nlevel)
         accums<-matrix(0,nrow=nlevel,ncol=ninter)
	 betac<-rep(0,p)

         counter<-matrix(0,nrow=nlevel,ncol=ninter)
	 cpo<-rep(0,nrec)

         endp<-rep(0,ninter-1)
	 etan<-rep(0,nrec)
	 
	 iflag<-rep(0,p)
	 intpn<-rep(0,nrec)
	 intpo<-rep(0,nrec)
	 
	 prob<-rep(0,ninter)

         rvecs<-matrix(0,nrow=nlevel,ncol=ninter)
	 seed<-c(sample(1:29000,1),sample(1:29000,1))
	 
	 workm1<-matrix(0,nrow=p,ncol=p)
	 workm2<-matrix(0,nrow=p,ncol=p)
	 workmh1<-rep(0,(p*(p+1)/2))
	 workv1<-rep(0,p)
	 workv2<-rep(0,p)
	 


         #########################################################################################
         # calling the fortran code
         #########################################################################################

         if(baseline=="logistic")
         {
            foo <- .Fortran("fptbinaryl",
	 	model     =as.integer(model),   	    
	 	nrec      =as.integer(nrec),
	 	p         =as.integer(p),
		sens      =as.double(sens),
		spec      =as.double(spec),
		x         =as.double(x),
		yobs      =as.integer(yobs),
		nlink     =as.integer(nlink),
		xlink     =as.double(xlink),
		a0b0      =as.double(a0b0),
		betapm    =as.double(betapm),		
		betapv    =as.double(betapv),		
		ninter    =as.integer(ninter),
		nlevel    =as.integer(nlevel),
		mcmcvec   =as.integer(mcmcvec),
		nsave     =as.integer(nsave),
		tune2     =as.double(tune2),
		propv     =as.double(propv),
		acrate    =as.double(acrate),
		fsave     =as.double(fsave),
		ppsave    =as.double(ppsave),
		randsave  =as.double(randsave),
		thetasave =as.double(thetasave),		
		cpo       =as.double(cpo),		
		alpha     =as.double(alpha),		
		beta      =as.double(beta),
		v         =as.double(v),
		y         =as.integer(y),
		accums    =as.double(accums),
		assign    =as.integer(assign),
		betac      =as.double(betac),
		counter   =as.integer(counter),
		endp      =as.double(endp),
		eta       =as.double(eta),
		etan      =as.double(etan),
		iflag     =as.integer(iflag),		
		intpn     =as.integer(intpn),		
		intpo     =as.integer(intpo),		
		prob      =as.double(prob),
		rvecs     =as.double(rvecs),
		seed      =as.integer(seed),
		workm1    =as.double(workm1),
		workm2    =as.double(workm2),
		workmh1   =as.double(workmh1),
		workv1    =as.double(workv1),
		workv2    =as.double(workv2),		
		PACKAGE="DPpackage")	
         }

         #########################################################################################
         # save state
         #########################################################################################
	
	 fsave<-matrix(foo$fsave,nrow=mcmc$nsave, ncol=nlink)
 	 thetasave<-matrix(foo$thetasave,nrow=mcmc$nsave, ncol=(p+1))
 	 randsave<-matrix(foo$randsave,nrow=mcmc$nsave, ncol=(nrec+1))
 	 ppsave<-matrix(foo$ppsave,nrow=mcmc$nsave, ncol=ninter)

         model.name<-"Bayesian binary regression model using a FPT prior"		

  	 colnames(thetasave)<-c(dimnames(x)[[2]],"alpha")

         coeff<-apply(thetasave, 2, mean)

	 names(coeff)<-c(dimnames(x)[[2]],"alpha")

         qnames<-NULL
         for(i in 1:nrec){
             idname<-paste("(Subject",i,sep="=")
             idname<-paste(idname,")",sep="")
             qnames<-c(qnames,idname)
         }
         qnames<-c(qnames,"Prediction")
         colnames(randsave)<-qnames


	 if(is.null(prior$a0))
	 {
		acrate<-foo$acrate[1]
	 }
	 else
	 {
	 	acrate<-foo$acrate
	 }


	 state <- list(beta=foo$beta,v=foo$v,alpha=foo$alpha,y=foo$y)
	 save.state <- list(thetasave=thetasave,fsave=fsave,randsave=randsave,
	                    ppsave=ppsave)

	 z<-list(modelname=model.name,coefficients=coeff,acrate=acrate,call=cl,
	         prior=prior,mcmc=mcmc,state=state,save.state=save.state,nrec=foo$nrec,
	         cpo=foo$cpo,p=p,nlink=nlink,xlink=xlink,baseline=baseline,x=x,
	         ninter=ninter,nlevel=nlevel)

	 cat("\n\n")

	 class(z)<-c("FPTbinary")
	 return(z) 
}



###                    
### Estimate the probability curve for a fitted binary 
### regression model using a Finite Polya tree prior.
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 16-08-2006.


"predict.FPTbinary"<-
function(object,xnew=NULL,hpd=TRUE, ...)
{
    stde<-function(x)
    {
    	n<-length(x)
    	return(sd(x)/sqrt(n))
    }

    hpdf<-function(x)
    {
         alpha<-0.05
         vec<-x
         n<-length(x)         
         alow<-rep(0,2)
         aupp<-rep(0,2)
         a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
         return(c(a$alow[1],a$aupp[1]))
    }
    
    pdf<-function(x)
    {
         alpha<-0.05
         vec<-x
         n<-length(x)         
         alow<-rep(0,2)
         aupp<-rep(0,2)
         a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
         return(c(a$alow[2],a$aupp[2]))
    }

   if(is.null(xnew))
   {
      xnew<-object$x
   }

   if(is(object, "FPTbinary"))
   {
 
   	 npred<-dim(xnew)[1]
   	 pnew<-dim(xnew)[2]
   	 nrec<-object$nrec
 
         baseline<-object$baseline
   	 alpha<-object$save.state$thetasave[,(object$p+1)]
   	 nsave<-length(alpha)
   	 ninter<-object$ninter
   	 nlevel<-object$nlevel

 	 ppsave<-matrix(object$save.state$ppsave,nrow=nsave, ncol=ninter)
   	 
   	 if (object$p != pnew) 
   	 {
	   stop("Dimension of xnew is not the same that the design matrix
	    in the original model.\n")
	 }

	 covn<-rep(0,npred) 
	 
	 for(i in 1:npred)
	 {
             covnw<-round(xnew[i,1],3)
             for(j in 2:pnew){
                 covnw<-paste(covnw,round(xnew[i,j],3),sep=";") 
             }
             covn[i]<-covnw  
         }    

	 lp<-xnew%*%t(object$save.state$thetasave[,1:object$p])
	 out<-matrix(0,nrow=npred,ncol=nsave)

         if(baseline=="logistic")
         {
   	    foo <- .Fortran("fptpredl",
		ninter    =as.integer(ninter),
		nlevel    =as.integer(nlevel),
		npred     =as.integer(npred),
		nsave     =as.integer(nsave),
		lp        =as.double(lp),
		ppsave    =as.double(ppsave),
		out       =as.double(out),
		PACKAGE="DPpackage")	
         }
           
         out<-t(matrix(foo$out,nrow=npred,ncol=nsave))

         pm <-apply(out, 2, mean)    
         pmed <-apply(out, 2, median)    
         psd<-apply(out, 2, sd)
         pstd<-apply(out, 2, stde)

         if(hpd){             
            limm<-apply(out, 2, hpdf)
            plinf<-limm[1,]
            plsup<-limm[2,]
         }
         else
         {
            plinf<-apply(out, 2, pdf)
            coef.l<-limm[1,]
            plsup<-limm[2,]
         }

   	 
   	 names(pm)<-covn
   	 names(pmed)<-covn
   	 names(psd)<-covn
   	 names(pstd)<-covn
   	 names(plinf)<-covn
   	 names(plsup)<-covn
   	 
   	 out<-NULL
   	 out$pmean<-pm
   	 out$pmedian<-pmed
   	 out$psd<-psd
   	 out$pstd<-pstd
   	 out$plinf<-plinf
   	 out$plsup<-plsup
   	 out$npred<-npred

   	 out$covn<-covn
   }
  
   out
}   




###
### Tools for FPTbinary: print, summary, plot
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 16-08-2006.
	
	        
"print.FPTbinary"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")

    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(x$cpo)), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
    cat("\nPosterior Inference of Parameters:\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)

    cat("\nAcceptance Rate for Metropolis Step = ",x$acrate,"\n")    

    cat("\nNumber of Observations:",x$nrec)
    cat("\n\n")
    invisible(x)
}


"summary.FPTbinary"<-function(object, hpd=TRUE, ...) 
{
    stde<-function(x)
    {
    	n<-length(x)
    	return(sd(x)/sqrt(n))
    }

    hpdf<-function(x)
    {
         alpha<-0.05
         vec<-x
         n<-length(x)         
         alow<-rep(0,2)
         aupp<-rep(0,2)
         a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
         return(c(a$alow[1],a$aupp[1]))
    }
    
    pdf<-function(x)
    {
         alpha<-0.05
         vec<-x
         n<-length(x)         
         alow<-rep(0,2)
         aupp<-rep(0,2)
         a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
         return(c(a$alow[2],a$aupp[2]))
    }

    thetasave<-object$save.state$thetasave


### Fixed part of the model

    dimen1<-object$p

    if(dimen1==1)
    {
       mat<-matrix(thetasave[,1:dimen1],ncol=1) 
    }
    else
    {
       mat<-thetasave[,1:dimen1]
    }

    coef.p<-object$coefficients[1:dimen1]
    coef.m <-apply(mat, 2, median)    
    coef.sd<-apply(mat, 2, sd)
    coef.se<-apply(mat, 2, stde)

    if(hpd){             
         limm<-apply(mat, 2, hpdf)
         coef.l<-limm[1,]
         coef.u<-limm[2,]
    }
    else
    {
         limm<-apply(mat, 2, pdf)
         coef.l<-limm[1,]
         coef.u<-limm[2,]
    }

    names(coef.m)<-names(object$coefficients[1:dimen1])
    names(coef.sd)<-names(object$coefficients[1:dimen1])
    names(coef.se)<-names(object$coefficients[1:dimen1])
    names(coef.l)<-names(object$coefficients[1:dimen1])
    names(coef.u)<-names(object$coefficients[1:dimen1])

    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.se , coef.l , coef.u)

    if(hpd)
    {
       dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
    }
    else
    {
       dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%CI-Low","95%CI-Upp"))
    }
    
    ans <- c(object[c("call", "modelname")])

    ans$coefficients<-coef.table


### CPO
    ans$cpo<-object$cpo


### Baseline Information


    if(is.null(object$prior$a0))
    {
	dimen2<-0	
    }
    else
    {
	dimen2<-1	    
    }

    if(dimen2==1)
    {
       mat<-matrix(thetasave[,(dimen1+1):(dimen1+dimen2)],ncol=1) 
       coef.p<-object$coefficients[(dimen1+1):(dimen1+dimen2)]
       coef.m <-apply(mat, 2, median)    
       coef.sd<-apply(mat, 2, sd)
       coef.se<-apply(mat, 2, stde)

       if(hpd){             
            limm<-apply(mat, 2, hpdf)
            coef.l<-limm[1,]
            coef.u<-limm[2,]
       }
       else
       {
            limm<-apply(mat, 2, pdf)
            coef.l<-limm[1,]
            coef.u<-limm[2,]
       }

       coef.table <- cbind(coef.p, coef.m, coef.sd, coef.se , coef.l , coef.u)
       if(hpd)
       {
            dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
       }
       else
       {
            dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%CI-Low","95%CI-Upp"))
       }

       ans$prec<-coef.table

    }

    ans$acrate<-object$acrate
    
    ans$nrec<-object$nrec

    class(ans) <- "summaryFPTbinary"
    return(ans)
}


"print.summaryFPTbinary"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(x$cpo)), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
    if (length(x$coefficients)) {
        cat("\nRegression coefficients:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")

    if (length(x$prec)) {
        cat("\nPrecision parameter:\n")
        print.default(format(x$prec, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No precision parameter\n")


    cat("\nAcceptance Rate for Metropolis Step = ",x$acrate,"\n")    
    
    cat("\nNumber of Observations:",x$nrec)
    cat("\n\n")
    invisible(x)
}


"plot.FPTbinary"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...)
{

fancydensplot1<-function(x, hpd=TRUE, npts=200, xlab="", ylab="", main="",col="#bdfcc9", ...)
# Author: AJV, 2006
#
{
	dens <- density(x,n=npts)
	densx <- dens$x
	densy <- dens$y

	meanvar <- mean(x)
	densx1 <- max(densx[densx<=meanvar])
	densx2 <- min(densx[densx>=meanvar])
        densy1 <- densy[densx==densx1]
        densy2 <- densy[densx==densx2]
        ymean <- densy1 + ((densy2-densy1)/(densx2-densx1))*(meanvar-densx1)
        

        if(hpd==TRUE)
	{
		alpha<-0.05
		alow<-rep(0,2)
        	aupp<-rep(0,2)
        	n<-length(x)
		a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(x),
		                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
		xlinf<-a$alow[1]            
		xlsup<-a$aupp[1]            
	}
	else
	{
		xlinf <- quantile(x,0.025)
		xlsup <- quantile(x,0.975)
	}

	densx1 <- max(densx[densx<=xlinf])
	densx2 <- min(densx[densx>=xlinf])
	densy1 <- densy[densx==densx1]
	densy2 <- densy[densx==densx2]
	ylinf <- densy1 + ((densy2-densy1)/(densx2-densx1))*(xlinf-densx1)

	densx1 <- max(densx[densx<=xlsup])
	densx2 <- min(densx[densx>=xlsup])
        densy1 <- densy[densx==densx1]
        densy2 <- densy[densx==densx2]
        ylsup <- densy1 + ((densy2-densy1)/(densx2-densx1))*(xlsup-densx1)

        plot(0.,0.,xlim = c(min(densx), max(densx)), ylim = c(min(densy), max(densy)),
             axes = F,type = "n" , xlab=xlab, ylab=ylab, main=main, cex=1.2)

        
        xpol<-c(xlinf,xlinf,densx[densx>=xlinf & densx <=xlsup],xlsup,xlsup)
        ypol<-c(0,ylinf,densy[densx>=xlinf & densx <=xlsup] ,ylsup,0)
             
        polygon(xpol, ypol, border = FALSE,col=col)
        
        lines(c(min(densx), max(densx)),c(0,0),lwd=1.2)
        
        segments(min(densx),0, min(densx),max(densy),lwd=1.2)
        
        lines(densx,densy,lwd=1.2)
             
        segments(meanvar, 0, meanvar, ymean,lwd=1.2)
        segments(xlinf, 0, xlinf, ylinf,lwd=1.2)
        segments(xlsup, 0, xlsup, ylsup,lwd=1.2)

	axis(1., at = round(c(xlinf, meanvar,xlsup), 2.), labels = T,pos = 0.)
        axis(1., at = round(seq(min(densx),max(densx),length=15), 2.), labels = F,pos = 0.)
        axis(2., at = round(seq(0,max(densy),length=5), 2.), labels = T,pos =min(densx))
}

    stde<-function(x)
    {
    	n<-length(x)
    	return(sd(x)/sqrt(n))
    }

    hpdf<-function(x)
    {
         alpha<-0.05
         vec<-x
         n<-length(x)         
         alow<-rep(0,2)
         aupp<-rep(0,2)
         a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
         return(c(a$alow[1],a$aupp[1]))
    }
    
    pdf<-function(x)
    {
         alpha<-0.05
         vec<-x
         n<-length(x)         
         alow<-rep(0,2)
         aupp<-rep(0,2)
         a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
         return(c(a$alow[2],a$aupp[2]))
    }

   if(is(x, "FPTbinary"))
   {
        if(is.null(param))
        {
           coef.p<-x$coefficients
           n<-length(coef.p)
           pnames<-names(coef.p)
           
           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
           for(i in 1:(n-1))
           {
               title1<-paste("Trace of",pnames[i],sep=" ")
               title2<-paste("Density of",pnames[i],sep=" ")       
               plot(ts(x$save.state$thetasave[,i]),main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot1(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }
           
           if(is.null(x$prior$a0))
           {
               cat("")
           }
           else
           {
               title1<-paste("Trace of",pnames[n],sep=" ")
               title2<-paste("Density of",pnames[n],sep=" ")       
               plot(ts(x$save.state$thetasave[,n]),main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot1(x$save.state$thetasave[,n],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }
           
           title1<-c("Predictive Error Density")
           title2<-c("Link Function")
           fancydensplot1(x$save.state$randsave[,(x$nrec+1)],hpd=hpd,main=title1,xlab="values", ylab="density",col=col)

           pml <-apply(x$save.state$fsave, 2, mean)    
           if(hpd){             
                limm<-apply(x$save.state$fsave, 2, hpdf)
                pll<-limm[1,]
                plu<-limm[2,]
           }
           else
           {
                limm<-apply(x$save.state$fsave, 2, pdf)
                pll<-limm[1,]
                plu<-limm[2,]
           } 
       
            plot(x$xlink,pml,xlab="x",ylab="probability",main=title2,lty=1,type='l',lwd=2,ylim=c(0,1))
            lines(x$xlink,pll,lty=2,lwd=2)
            lines(x$xlink,plu,lty=2,lwd=2)            
        }
   
        else
        {
            coef.p<-x$coefficients
	    n<-length(coef.p)
	    pnames<-names(coef.p)
	    poss<-0 
            for(i in 1:n)
            {
               if(pnames[i]==param)poss=i
            }
            if (poss==0 && param !="link") 
	    {
	      stop("This parameter is not present in the original model.\n")
	    }
	    
	    par(ask = ask)
	    layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
	    
	    if(param !="link")
	    {
               title1<-paste("Trace of",pnames[poss],sep=" ")
               title2<-paste("Density of",pnames[poss],sep=" ")       
               plot(x$save.state$thetasave[,poss],type='l',main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot1(x$save.state$thetasave[,poss],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
            }        
        
            else
            {
               title1<-c("Predictive Error Density")
               title2<-c("Link Function")
               fancydensplot1(x$save.state$randsave[,(x$nrec+1)],hpd=hpd,main=title1,xlab="values", ylab="density",col=col)
           
               pml <-apply(x$save.state$fsave, 2, mean)    
               if(hpd){             
                    limm<-apply(x$save.state$fsave, 2, hpdf)
                    pll<-limm[1,]
                    plu<-limm[2,]
               }
               else
               {
                    limm<-apply(x$save.state$fsave, 2, pdf)
                    pll<-limm[1,]
                    plu<-limm[2,]
               } 
       
               plot(x$xlink,pml,xlab="x",ylab="probability",main=title2,lty=1,type='l',lwd=2,ylim=c(0,1))
               lines(x$xlink,pll,lty=2,lwd=2)
               lines(x$xlink,plu,lty=2,lwd=2)            
            }
        }
   }

}


