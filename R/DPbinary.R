### DPbinary.R
### Fit a semiparametric bernoulli regression model.
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 21-08-2006.
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
### The author's contact information:
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

"DPbinary"<-
function(formula,inter=TRUE,baseline="logistic",prior,mcmc,state,status,misc=NULL,
data=sys.frame(sys.parent()),na.action=na.fail) 
UseMethod("DPbinary")

"DPbinary.default"<-
function(formula,
         inter=TRUE,
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
	 
	 if(inter)
	 {
	    x<-x
	    p<-p
	 }
	 else
	 {
	    x<-x[,2:p]
	    p<-(p-1)
	 }   

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
          
         if(inter)
         {
            xmatrix<-x
         }
         else
         {
            xmatrix<-cbind(rep(1,nrec),x) 
         }   

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

         if(inter)
         {
            propv<- fit0$covb
         }
         else
         {
            propv<- fit0$covb[2:(p+1),2:(p+1)]
         }
         
         if(is.null(mcmc$tune))
         {
            tune=1
         }
         else
         {
            tune<-mcmc$tune
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
	
	 propv<-diag(tune,p)%*%solve(solve(betapv)+solve(propv))%*%diag(tune,p)

         #########################################################################################
         # parameters depending on status
         #########################################################################################

       	 if(status)
	 {
	        if(inter)
	        {
	           beta<-fit0$beta
	        }
	        else
	        {
     	           beta<-fit0$beta[2:(p+1)]
     	        }   
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
         
         xlink=rev(seq(-6,6,length=34))
         nlink<-length(xlink)
         fsave     <- matrix(0, nrow=mcmc$nsave, ncol=nlink)         
	 thetasave <- matrix(0, nrow=mcmc$nsave, ncol=p+2)
	 randsave  <- matrix(0, nrow=mcmc$nsave, ncol=nrec+1)


         #########################################################################################
         # working space
         #########################################################################################

	 ntotal<-nrec
	 maxint<-2*ntotal+1
	 maxend<-2*ntotal
	 maxm<-log(0.0001)/log((ntotal+1)/(ntotal+2))

         acrate<-0	 
	 betac<-rep(0,p)
	 endp<-rep(0,maxend)
	 endp2<-rep(0,maxend)
	 etan<-rep(0,nrec)
	 fsavet<-rep(0,nlink)
	 intcount<-rep(0,maxint)
	 intcount2<-rep(0,maxint)
	 intind<-rep(0,(nrec+1))
	 intind2<-rep(0,(nrec+1))
	 iflag<-rep(0,p)
	 intposso<-matrix(0,nrow=maxint,ncol=(nrec+1))
         intpossn<-matrix(0,nrow=maxint,ncol=(nrec+1))
         lpsav<-rep(0,nrec)
         lpsavc<-rep(0,nrec)
	 mass<-rep(0,maxint)
	 ncluster<-nrec
	 prob<-rep(0,maxint)
	 uvec<-rep(0,maxm)
	 vvec<-rep(0,maxm)
	 vnew<-rep(0,(nrec+1))
	 vnew2<-rep(0,nrec)
	 workm1<-matrix(0,nrow=p,ncol=p)
	 workm2<-matrix(0,nrow=p,ncol=p)
	 workmh1<-rep(0,(p*(p+1)/2))
	 workv1<-rep(0,p)
	 workv2<-rep(0,p)
	 wvec<-rep(0,maxm)
	 seed1<-sample(1:29000,1)
	 seed2<-sample(1:29000,1)
	 cpo<-rep(0,nrec)


         #########################################################################################
         # calling the fortran code
         #########################################################################################

         if(baseline=="logistic")
         {
   	    foo <- .Fortran("dpbinaryl2",
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
		mcmcvec   =as.integer(mcmcvec),
		nsave     =as.integer(nsave),
		propv     =as.double(propv),
		acrate    =as.double(acrate),
		fsave     =as.double(fsave),
		randsave  =as.double(randsave),
		thetasave =as.double(thetasave),		
		cpo       =as.double(cpo),		
		alpha     =as.double(alpha),		
		beta      =as.double(beta),
		ncluster  =as.integer(ncluster),
		y         =as.integer(y),
		v         =as.double(v),
		betac     =as.double(betac),
		endp      =as.double(endp),
		endp2     =as.double(endp2),
		eta       =as.double(eta),
		etan      =as.double(etan),
		fsavet    =as.integer(fsavet),
		intcount  =as.integer(intcount),
		intcount2 =as.integer(intcount2),
		intind    =as.integer(intind),
		intind2   =as.integer(intind2),
		iflag     =as.integer(iflag),
		intposso  =as.integer(intposso),
		intpossn  =as.integer(intpossn),
		lpsav     =as.integer(lpsav),
		lpsavc    =as.integer(lpsavc),
		maxint    =as.integer(maxint),
		maxend    =as.integer(maxend),
		maxm      =as.integer(maxm),
		mass      =as.double(mass),
		prob      =as.double(prob),
		seed1     =as.integer(seed1),
		seed2     =as.integer(seed2),
		uvec      =as.double(uvec),
		vvec      =as.double(vvec),
		vnew      =as.double(vnew),
		vnew2     =as.double(vnew2),
		workm1    =as.double(workm1),
		workm2    =as.double(workm2),
		workmh1   =as.double(workmh1),
		workv1    =as.double(workv1),
		workv2    =as.double(workv2),
		wvec      =as.double(wvec),
		PACKAGE="DPpackage")	
         }


         if(baseline=="normal")
         {
   	    foo <- .Fortran("dpbinaryp2",
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
		mcmcvec   =as.integer(mcmcvec),
		nsave     =as.integer(nsave),
		propv     =as.double(propv),
		acrate    =as.double(acrate),
		fsave     =as.double(fsave),
		randsave  =as.double(randsave),
		thetasave =as.double(thetasave),		
		cpo       =as.double(cpo),		
		alpha     =as.double(alpha),		
		beta      =as.double(beta),
		ncluster  =as.integer(ncluster),
		y         =as.integer(y),
		v         =as.double(v),
		betac     =as.double(betac),
		endp      =as.double(endp),
		endp2     =as.double(endp2),
		eta       =as.double(eta),
		etan      =as.double(etan),
		fsavet    =as.integer(fsavet),
		intcount  =as.integer(intcount),
		intcount2 =as.integer(intcount2),
		intind    =as.integer(intind),
		intind2   =as.integer(intind2),
		iflag     =as.integer(iflag),
		intposso  =as.integer(intposso),
		intpossn  =as.integer(intpossn),
		lpsav     =as.integer(lpsav),
		lpsavc    =as.integer(lpsavc),
		maxint    =as.integer(maxint),
		maxend    =as.integer(maxend),
		maxm      =as.integer(maxm),
		mass      =as.double(mass),
		prob      =as.double(prob),
		seed1     =as.integer(seed1),
		seed2     =as.integer(seed2),
		uvec      =as.double(uvec),
		vvec      =as.double(vvec),
		vnew      =as.double(vnew),
		vnew2     =as.double(vnew2),
		workm1    =as.double(workm1),
		workm2    =as.double(workm2),
		workmh1   =as.double(workmh1),
		workv1    =as.double(workv1),
		workv2    =as.double(workv2),
		wvec      =as.double(wvec),
		PACKAGE="DPpackage")	
         }


         if(baseline=="cauchy")
         {
   	    foo <- .Fortran("dpbinaryc2",
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
		mcmcvec   =as.integer(mcmcvec),
		nsave     =as.integer(nsave),
		propv     =as.double(propv),
		acrate    =as.double(acrate),
		fsave     =as.double(fsave),
		randsave  =as.double(randsave),
		thetasave =as.double(thetasave),		
		cpo       =as.double(cpo),		
		alpha     =as.double(alpha),		
		beta      =as.double(beta),
		ncluster  =as.integer(ncluster),
		y         =as.integer(y),
		v         =as.double(v),
		betac     =as.double(betac),
		endp      =as.double(endp),
		endp2     =as.double(endp2),
		eta       =as.double(eta),
		etan      =as.double(etan),
		fsavet    =as.integer(fsavet),
		intcount  =as.integer(intcount),
		intcount2 =as.integer(intcount2),
		intind    =as.integer(intind),
		intind2   =as.integer(intind2),
		iflag     =as.integer(iflag),
		intposso  =as.integer(intposso),
		intpossn  =as.integer(intpossn),
		lpsav     =as.integer(lpsav),
		lpsavc    =as.integer(lpsavc),
		maxint    =as.integer(maxint),
		maxend    =as.integer(maxend),
		maxm      =as.integer(maxm),
		mass      =as.double(mass),
		prob      =as.double(prob),
		seed1     =as.integer(seed1),
		seed2     =as.integer(seed2),
		uvec      =as.double(uvec),
		vvec      =as.double(vvec),
		vnew      =as.double(vnew),
		vnew2     =as.double(vnew2),
		workm1    =as.double(workm1),
		workm2    =as.double(workm2),
		workmh1   =as.double(workmh1),
		workv1    =as.double(workv1),
		workv2    =as.double(workv2),
		wvec      =as.double(wvec),
		PACKAGE="DPpackage")	
         }



         #########################################################################################
         # save state
         #########################################################################################
	
	 fsave<-matrix(foo$fsave,nrow=mcmc$nsave, ncol=nlink)
 	 thetasave<-matrix(foo$thetasave,nrow=mcmc$nsave, ncol=(p+2))
 	 randsave<-matrix(foo$randsave,nrow=mcmc$nsave, ncol=(nrec+1))

         model.name<-"Bayesian semiparametric binary regression model"		

  	 colnames(thetasave)<-c(dimnames(x)[[2]],"ncluster","alpha")
  	 
         coeff<-apply(thetasave, 2, mean)  	 
  	 
	 names(coeff)<-c(dimnames(x)[[2]],"ncluster","alpha")
	
         qnames<-NULL
         for(i in 1:nrec){
             idname<-paste("(Subject",i,sep="=")
             idname<-paste(idname,")",sep="")
             qnames<-c(qnames,idname)
         }
         qnames<-c(qnames,"Prediction")
         colnames(randsave)<-qnames
	
	 state <- list(beta=foo$beta,v=foo$v,y=foo$y,alpha=foo$alpha)
				  
	 save.state <- list(thetasave=thetasave,fsave=fsave,randsave=randsave)

	 z<-list(modelname=model.name,coefficients=coeff,acrate=foo$acrate,call=cl,
	         prior=prior,mcmc=mcmc,state=state,save.state=save.state,nrec=foo$nrec,
	         cpo=foo$cpo,p=p,nlink=nlink,xlink=xlink,baseline=baseline,x=x)
	         
	 cat("\n\n")

	 class(z)<-c("DPbinary")
	 return(z) 
}



###                    
### Estimate the probability curve for a fitted semiparametric binary 
### regression model.
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 21-08-2006.


"predict.DPbinary"<-
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

   if(is.null(xnew))
   {
      xnew<-object$x
   }

   if(is(object, "DPbinary"))
   {
 
   	 npred<-dim(xnew)[1]
   	 pnew<-dim(xnew)[2]
   	 nrec<-object$nrec
 
         baseline<-object$baseline
   	 alpha<-object$save.state$thetasave[,(object$p+2)]
   	 nsave<-length(alpha)
   	 
   	 v<-matrix(object$save.state$randsave,nsave,nrec+1)
   	 vpred<-v[,(nrec+1)]
   	 v<-v[,1:nrec]

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
   	    foo <- .Fortran("dppredl",
                nsave     =as.integer(nsave),
                nrec      =as.integer(nrec),
                npred     =as.integer(npred),
                alpha     =as.double(alpha),
                v         =as.double(v),
                lp        =as.double(lp),
		out       =as.double(out),
		PACKAGE="DPpackage")	
         }

         if(baseline=="normal")
         {
   	    foo <- .Fortran("dppredp",
                nsave     =as.integer(nsave),
                nrec      =as.integer(nrec),
                npred     =as.integer(npred),
                alpha     =as.double(alpha),
                v         =as.double(v),
                lp        =as.double(lp),
		out       =as.double(out),
		PACKAGE="DPpackage")	
         }
         
         if(baseline=="cauchy")
         {
   	    foo <- .Fortran("dppredc",
                nsave     =as.integer(nsave),
                nrec      =as.integer(nrec),
                npred     =as.integer(npred),
                alpha     =as.double(alpha),
                v         =as.double(v),
                lp        =as.double(lp),
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
   	 out$vpred<-vpred
   	 out$npred<-npred

   	 out$covn<-covn
   }
  
   out
}   


###
### Tools for DPbinary: print, summary, plot
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 21-08-2006.
	
	        
"print.DPbinary"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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


"summary.DPbinary"<-function(object, hpd=TRUE, ...) 
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
	dimen2<-1	
    }
    else
    {
	dimen2<-2	
    }

    if(dimen2==1)
    {
       mat<-matrix(thetasave[,(dimen1+1):(dimen1+dimen2)],ncol=1) 
    }
    else
    {
       mat<-thetasave[,(dimen1+1):(dimen1+dimen2)] 
    }   
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

    ans$acrate<-object$acrate
    
    ans$nrec<-object$nrec

    class(ans) <- "summaryDPbinary"
    return(ans)
}



"print.summaryDPbinary"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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




"plot.DPbinary"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...) 
{

fancydensplot<-function(x, hpd=TRUE, npts=200, xlab="", ylab="", main="",col="#bdfcc9", ...)
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



    if(is(x, "DPbinary")){
        if(is.null(param))
	{
           coef.p<-x$coefficients
           n<-length(coef.p)
           pnames<-names(coef.p)

           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
           for(i in 1:(n-1)){
               title1<-paste("Trace of",pnames[i],sep=" ")
               title2<-paste("Density of",pnames[i],sep=" ")       
               plot(ts(x$save.state$thetasave[,i]),type='l',main=title1,xlab="MCMC scan",ylab=" ")
               if(pnames[i]=="ncluster")
	       {
	          hist(x$save.state$thetasave[,i],main=title2,xlab="values", ylab="probability",probability=TRUE)
	       }
	       else
	       {
                  fancydensplot(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
               }   
           }
           
           if(is.null(x$prior$a0))
           {
              cat("")
           }
           else
           {
               title1<-paste("Trace of",pnames[n],sep=" ")
               title2<-paste("Density of",pnames[n],sep=" ")       
               plot(ts(x$save.state$thetasave[,n]),type='l',main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot(x$save.state$thetasave[,n],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }
           
           title1<-c("Predictive Error Density")
           title2<-c("Link Function")
           fancydensplot(x$save.state$randsave[,(x$nrec+1)],hpd=hpd,main=title1,xlab="values", ylab="density",col=col)

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
               plot(ts(x$save.state$thetasave[,poss]),type='l',main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot(x$save.state$thetasave[,poss],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
            }
            else
            {

               title1<-c("Predictive Error Density")
               title2<-c("Link Function")
               fancydensplot(x$save.state$randsave[,(x$nrec+1)],hpd=hpd,main=title1,xlab="values", ylab="density",col=col)
           
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

