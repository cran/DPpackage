### CSDPbinary.R
### Fit a semiparametric bernoulli regression model.
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 05-07-2006.
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

"CSDPbinary"<-
function(formula,prior,mcmc,state,status,misc=NULL,
data=sys.frame(sys.parent()),na.action=na.fail) 
UseMethod("CSDPbinary")

"CSDPbinary.default"<-
function(formula,
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
		if(length(sens)==1)
		{
		   sens<-rep(misc$sens,nrec)
  	  	   spec<-rep(misc$spec,nrec)
		}
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
         

         mcmcvec<-c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay,mcmc$ntheta,model)
         nsave<-mcmc$nsave
          
         xmatrix<-x

         fit0<- MLElogit(xmatrix,yobs,sens,spec)
         propv<- fit0$covb
 
         if(is.null(mcmc$tune))
         {
            tune=1
         }
         else
         {
            tune<-mcmc$tune
         }

         #########################################################################################
         # prior information
         #########################################################################################

	 if(is.null(prior$a0))
	 {
		a0b0<-c(-1,-1,prior$d,prior$p)
		alpha<-prior$alpha
	 }
	 else
	 {
	 	a0b0<-c(prior$a0,prior$b0,prior$d,prior$p)
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
                beta<-fit0$beta
		eta <- x %*% beta
		v <- rep(0,nrec)
		v[yobs==1]<-eta[yobs==1]-0.5
		v[yobs==0]<-eta[yobs==0]+0.5
		y<-yobs
		theta<-log(3)
	 }	
      	 else
	 {
		beta<-state$beta
		eta <- x %*% beta
		v<-state$v
		y<-state$y
		alpha<-state$alpha
		theta<-state$theta
	 }
	 

         #########################################################################################
         # output
         #########################################################################################
         
         xlink=rev(seq(-6,6,length=34))
         nlink<-length(xlink)
         cpo<-rep(0,nrec)
         fsave     <- matrix(0, nrow=nsave, ncol=nlink)         
	 thetasave <- matrix(0, nrow=nsave, ncol=p+3)
	 randsave  <- matrix(0, nrow=nsave, ncol=nrec+1)

         #########################################################################################
         # working space
         #########################################################################################

	 ntotal<-nrec
	 maxint<-2*(ntotal+3)+1
	 maxend<-2*(ntotal+3)
	 maxm<-log(0.0001)/log((ntotal+1)/(ntotal+2))

         acrate<-rep(0,2)
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
	 massurn1<-rep(0,maxint)
	 massurn2<-rep(0,maxint)
	 massurn3<-rep(0,maxint)
	 massurn4<-rep(0,maxint)
	 ncluster<-nrec
	 prob<-rep(0,maxint)
	 proburn1<-rep(0,maxint)
	 proburn2<-rep(0,maxint)
	 proburn3<-rep(0,maxint)
	 proburn4<-rep(0,maxint)
	 urn<-rep(0,maxint)
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
	 seed<-c(sample(1:29000,1),sample(1:29000,1))

         #########################################################################################
         # calling the fortran code
         #########################################################################################

   	    foo <- .Fortran("csdpbinaryl2",
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
		theta     =as.double(theta),
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
		massurn1  =as.double(massurn1),
		massurn2  =as.double(massurn2),
		massurn3  =as.double(massurn3),
		massurn4  =as.double(massurn4),
		prob      =as.double(prob),
		proburn1  =as.double(proburn1),
		proburn2  =as.double(proburn2),
		proburn3  =as.double(proburn3),
		proburn4  =as.double(proburn4),
		seed      =as.integer(seed),
                urn       =as.integer(urn),
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

         #########################################################################################
         # save state
         #########################################################################################
	
	 fsave<-matrix(foo$fsave,nrow=nsave, ncol=nlink)
 	 thetasave<-matrix(foo$thetasave,nrow=nsave, ncol=(p+3))
 	 randsave<-matrix(foo$randsave,nrow=nsave, ncol=(nrec+1))

         model.name<-"Bayesian semiparametric binary regression model"		

  	 colnames(thetasave)<-c(dimnames(x)[[2]],"theta","ncluster","alpha")
  	 
  	 coeff<-rep(0,(p+3))
	 for(i in 1:(p+3)){
	     coeff[i]<-mean(thetasave[,i])
	 }  
	 names(coeff)<-c(dimnames(x)[[2]],"theta","ncluster","alpha")
	
         qnames<-NULL
         for(i in 1:nrec){
             idname<-paste("(Subject",i,sep="=")
             idname<-paste(idname,")",sep="")
             qnames<-c(qnames,idname)
         }
         qnames<-c(qnames,"Prediction")
         colnames(randsave)<-qnames
	
	 state <- list(beta=foo$beta,v=foo$v,y=foo$y,alpha=foo$alpha,theta=foo$theta)
				  
	 save.state <- list(thetasave=thetasave,fsave=fsave,randsave=randsave)

	 z<-list(modelname=model.name,coefficients=coeff,acrate=foo$acrate,call=cl,
	         prior=prior,mcmc=mcmcvec,state=state,save.state=save.state,nrec=foo$nrec,
	         cpo=foo$cpo,p=p,nlink=nlink,xlink=xlink,x=x,
	         ppar=prior$p,d=prior$d)
	         
	 cat("\n\n")

	 class(z)<-c("CSDPbinary")
	 return(z) 
}

###                    
### Estimate the probability curve for a fitted semiparametric binary 
### regression model.
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 01-07-2006.

predict.CSDPbinary<-function(object,xnew=NULL,hpd=TRUE, ...)
{
   if(is.null(xnew))
   {
      xnew<-object$x
   }
   
   pbaselinel<-function(x,theta,d,p)
   {
       if(x<=(theta-d))
       {
          area<-plogis(theta-d)
          out<-((1-p)/2)*plogis(x)/area
       }

       if(x>(theta-d) && x<=0)
       {
          area<-0.5-plogis(theta-d)
          out<-(1-p)/2+ (p/2)*plogis(x/area)
       }

       if(x>0 && x<=theta)
       {
          area<-plogis(theta)-0.5
          out<-0.5+ (p/2)*plogis(x/area)
       }

       if(x>theta)
       {
          area<-1-plogis(theta)
          out<-(1+p)/2+ ((1-p)/2)*plogis(x/area)
       }
       return(out)
   }
	
   if(is(object, "CSDPbinary"))
   {
   	 npred<-dim(xnew)[1]
   	 pnew<-dim(xnew)[2]
   	 nrec<-object$nrec

   	 alpha<-object$save.state$thetasave[,(object$p+3)]
   	 theta<-object$save.state$thetasave[,(object$p+1)]
   	 ppar<-object$ppar
   	 d<-object$d
   	 
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

 	 pm<-rep(0,npred)
	 pmed<-rep(0,npred)
	 psd<-rep(0,npred)
	 pstd<-rep(0,npred)
	 plinf<-rep(0,npred)
	 plsup<-rep(0,npred)

	 for(i in 1:npred)
         {
	     prob<-rep(0,nsave)
	     for(k in 1:nsave)
	     {
		lsup<-lp[i,k]
		awork<-alpha[k]
		vwork<-v[k,]
		Cdelta<- sum(vwork<=lsup)
		Cparam<- awork*pbaselinel(lsup,theta[k],d,ppar)
		prob[k]<-(Cparam+Cdelta)/(awork+nrec)	                
	      }
	      pm[i]<-mean(prob)
	      pmed[i]<-median(prob)
	      psd[i]<-sqrt(var(prob))
	      pstd[i]<-sqrt(var(prob))/sqrt(nsave)
		 
              if(hpd==TRUE)
              {
                    alow<-rep(0,2)
                    aupp<-rep(0,2)
                    sig<-0.05
                    vec<-prob
                    n<-length(vec)
                    a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(sig),x=as.double(vec),
                                 alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                    plinf[i]<-a$alow[1]            
                    plsup[i]<-a$aupp[1]            
              }
              else
              {
		    plinf[i]<-quantile(surv[j,],0.025)		
		    plsup[i]<-quantile(surv[j,],0.975)			                
              }
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
   return(out)   
}


"print.CSDPbinary"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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



"summary.CSDPbinary"<-function(object, hpd=TRUE, ...) 
{
    dimen<-object$p
    coef.p<-object$coefficients[1:dimen]
    coef.sd<-rep(0,dimen)
    coef.se<-rep(0,dimen)
    coef.l<-rep(0,dimen)
    coef.u<-rep(0,dimen)
    coef.m<-rep(0,dimen)
    names(coef.sd)<-names(object$coefficients[1:dimen])
    names(coef.l)<-names(object$coefficients[1:dimen])
    names(coef.u)<-names(object$coefficients[1:dimen])
    
    alpha<-0.05
    
    for(i in 1:dimen){
        alow<-rep(0,2)
        aupp<-rep(0,2)
        coef.sd[i]<-sqrt(var(object$save.state$thetasave[,i]))
        coef.m[i]<-median(object$save.state$thetasave[,i])
        vec<-object$save.state$thetasave[,i]
        n<-length(vec)
        
        if(hpd==TRUE)
        {
        
                a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                                  alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                coef.l[i]<-a$alow[1]            
                coef.u[i]<-a$aupp[1]            
         }
         else
         {
                coef.l[i]<-quantile(vec,0.025) 
                coef.u[i]<-quantile(vec,0.975) 
         }
    }

    coef.se<-coef.sd/sqrt(n)

    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.se , coef.l , coef.u)
    
    if(hpd==TRUE)
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

    ans$cpo<-object$cpo
    

    if(is.null(object$prior$a0))
    {
	dimen<-2	
    }
    else
    {
	dimen<-3	    
    }

    coef.p<-object$coefficients[(object$p+1):(object$p+dimen)]
    coef.sd<-rep(0,dimen)
    coef.se<-rep(0,dimen)
    coef.l<-rep(0,dimen)
    coef.u<-rep(0,dimen)
    coef.m<-rep(0,dimen)
    names(coef.sd)<-names(object$coefficients[(object$p+1):(object$p+dimen)])
    names(coef.l)<-names(object$coefficients[(object$p+1):(object$p+dimen)])
    names(coef.u)<-names(object$coefficients[(object$p+1):(object$p+dimen)])
    for(i in 1:dimen){
         alow<-rep(0,2)
         aupp<-rep(0,2)
         coef.sd[i]<-sqrt(var(object$save.state$thetasave[,(object$p+i)]))
         coef.m[i]<-median(object$save.state$thetasave[,(object$p+i)])
         vec<-object$save.state$thetasave[,(object$p+i)]
         n<-length(vec)
        
         if(hpd==TRUE)
         {
             a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                              alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
             coef.l[i]<-a$alow[1]            
             coef.u[i]<-a$aupp[1]            
         }
         else
         {
                coef.l[i]<-quantile(vec,0.025) 
                coef.u[i]<-quantile(vec,0.975) 
         }
    }

    coef.se<-coef.sd/sqrt(n)

    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.se , coef.l , coef.u)

    if(hpd==TRUE)
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



"print.summaryCSDPbinary"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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




"plot.CSDPbinary"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...) 
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


    if(is(x, "CSDPbinary")){
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
               plot(x$save.state$thetasave[,i],type='l',main=title1,xlab="MCMC scan",ylab=" ")
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
               plot(x$save.state$thetasave[,n],type='l',main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot(x$save.state$thetasave[,n],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }
           
           title1<-c("Predictive Error Density")
           title2<-c("Link Function")
           fancydensplot(x$save.state$randsave[,(x$nrec+1)],hpd=hpd,main=title1,xlab="values", ylab="density",col=col)
           
           pml<-rep(0,x$nlink)
           pll<-rep(0,x$nlink)
           plu<-rep(0,x$nlink)
           alpha<-0.05
           
           for(i in 1:x$nlink)
           {
               pml[i]<-mean(x$save.state$fsave[,i])
               alow<-rep(0,2)
	       aupp<-rep(0,2)
	       vec<-x$save.state$fsave[,i]
	       n<-length(vec)
	       a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
		                 alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
	       if(hpd)
	       {
	           pll[i]<-a$alow[1]  
	           plu[i]<-a$aupp[1]
	       }
	       else
	       {
	           pll[i]<-a$alow[2]  
	           plu[i]<-a$aupp[2]
	       }
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
               fancydensplot(x$save.state$thetasave[,poss],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
            }
            else
            {
               title1<-c("Predictive Error Density")
               title2<-c("Link Function")
               fancydensplot(x$save.state$randsave[,(x$nrec+1)],hpd=hpd,main=title1,xlab="values", ylab="density",col=col)
           
               pml<-rep(0,x$nlink)
               pll<-rep(0,x$nlink)
               plu<-rep(0,x$nlink)
               alpha<-0.05
           
               for(i in 1:x$nlink)
               {
                   pml[i]<-mean(x$save.state$fsave[,i])
                   alow<-rep(0,2)
	           aupp<-rep(0,2)
	           vec<-x$save.state$fsave[,i]
	           n<-length(vec)
	           a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
		                 alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
	           if(hpd)
	           {
	               pll[i]<-a$alow[1]  
	               plu[i]<-a$aupp[1]
	           }
	           else
	           {
	              pll[i]<-a$alow[2]  
	              plu[i]<-a$aupp[2]
	           }
                }
                plot(x$xlink,pml,xlab="x",ylab="probability",main=title2,lty=1,type='l',lwd=2,ylim=c(0,1))
                lines(x$xlink,pll,lty=2,lwd=2)
                lines(x$xlink,plu,lty=2,lwd=2)
            }
        }
   }
}

