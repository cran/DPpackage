### PTlm.R                   
### Fit a semiparametric linear model.
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 08-10-2006.
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

"PTlm"<-
function(formula,ngrid=200,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("PTlm")

"PTlm.default"<-
function(formula,
         ngrid=200,
         prior,
         mcmc,
         state,
         status, 
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

  	 y<- model.response(mf,"numeric")
  	 nrec<-length(y)
  	 x<-model.matrix(formula)
  	 p<-dim(x)[2]


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
	 tau<-c(prior$tau1,prior$tau2)
	
         #########################################################################################
         # mcmc specification
         #########################################################################################

  	 if(is.null(mcmc$tune1))
  	 {
  	    tune1<-1.1
  	 }
  	 else
  	 {
  	    tune1<-mcmc$tune1
  	 }

  	 if(is.null(mcmc$tune2))
  	 {
  	    tune2<-1.1
  	 }
  	 else
  	 {
  	    tune2<-mcmc$tune2
  	 }

  	 if(is.null(mcmc$tune3))
  	 {
  	    tune3<-1.1
  	 }
  	 else
  	 {
  	    tune3<-mcmc$tune3
  	 }
  	 
         mcmcvec<-c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
         nsave<-mcmc$nsave

	 fit0<- lm(y~x[,2:p])
	 propv<- vcov(fit0)
	 
	 propv<-diag(tune1,p)%*%propv%*%diag(tune1,p)

         beta<-coefficients(fit0)
         
         s2<-sum((resid(fit0))**2)/(nrec-p)
	 

         #########################################################################################
         # output
         #########################################################################################
         left<-min(resid(fit0))-0.5*sqrt(var(resid(fit0)))
         right<-max(resid(fit0))+0.5*sqrt(var(resid(fit0)))

	 grid<-seq(left,right,length=ngrid)
	 f<-rep(0,ngrid)

	 thetasave <- matrix(0, nrow=nsave, ncol=p+2)
	 randsave  <- matrix(0, nrow=nsave, ncol=nrec+1)
	 
	
         #########################################################################################
         # parameters depending on status
         #########################################################################################

	 if(status==TRUE)
	 {
	 	beta <- beta
		v <- resid(fit0)
		sigma2<-s2
	 }	
	 if(status==FALSE)
	 {
	        alpha<-state$alpha 
		beta<-state$beta
		v<-state$v
		sigma2<-state$sigma2
	 }

         #########################################################################################
         # working space
         #########################################################################################

         mu<-0

         seed<-c(sample(1:29000,1),sample(1:29000,1),sample(1:29000,1))

 	 acrate<-rep(0,3)
	 cpo<-rep(0,nrec)
	 betac<-rep(0,p)
 	 iflag<-rep(0,p)
 	 vc<-rep(0,nrec)
	 workm1<-matrix(0,nrow=p,ncol=p)
	 workm2<-matrix(0,nrow=p,ncol=p)
	 workmh1<-rep(0,(p*(p+1)/2))
	 workv1<-rep(0,p)
	 workv2<-rep(0,p)
	
         #########################################################################################
         # calling the fortran code
         #########################################################################################
 
         if(is.null(prior$M))
         {
 	        foo <- .Fortran("ptlm",
  	 	ngrid      =as.integer(ngrid),
	 	nrec      =as.integer(nrec),
	 	p         =as.integer(p),
		x         =as.double(x),	 	
		y         =as.double(y),
		a0b0      =as.double(a0b0),
		betapm    =as.double(betapm),		
		betapv    =as.double(betapv),		
		tau       =as.double(tau),
		mcmc      =as.integer(mcmcvec),
		nsave     =as.integer(nsave),
		propv     =as.double(propv),
		tune2     =as.double(tune2),
		tune3     =as.double(tune3),
                seed      =as.integer(seed),
		acrate    =as.double(acrate),
		randsave  =as.double(randsave),
		thetasave =as.double(thetasave),
                cpo       =as.double(cpo),
                f          =as.double(f),
		alpha     =as.double(alpha),		
		beta      =as.double(beta),
		mu        =as.double(mu),
		sigma2    =as.double(sigma2),
                v         =as.double(v),
		betac     =as.double(betac),
		iflag     =as.integer(iflag),
		vc        =as.double(vc),
		workm1    =as.double(workm1),
		workm2    =as.double(workm2),
		workmh1   =as.double(workmh1),
		workv1    =as.double(workv1),
		workv2    =as.double(workv2),
		grid       =as.double(grid),
		PACKAGE="DPpackage")	         
         }
         else
         {
             	maxm<-prior$M
 	        foo <- .Fortran("ptlmp",
	 	maxm      =as.integer(maxm),
  	 	ngrid      =as.integer(ngrid),
	 	nrec      =as.integer(nrec),
	 	p         =as.integer(p),
		x         =as.double(x),	 	
		y         =as.double(y),
		a0b0      =as.double(a0b0),
		betapm    =as.double(betapm),		
		betapv    =as.double(betapv),		
		tau       =as.double(tau),
		mcmc      =as.integer(mcmcvec),
		nsave     =as.integer(nsave),
		propv     =as.double(propv),
		tune2     =as.double(tune2),
		tune3     =as.double(tune3),
		seed      =as.integer(seed),
		acrate    =as.double(acrate),
		randsave  =as.double(randsave),
		thetasave =as.double(thetasave),
                cpo       =as.double(cpo),
 		f          =as.double(f),                
		alpha     =as.double(alpha),		
		beta      =as.double(beta),
		mu        =as.double(mu),
		sigma2    =as.double(sigma2),
                v         =as.double(v),
		betac     =as.double(betac),
		iflag     =as.integer(iflag),
		vc        =as.double(vc),
		workm1    =as.double(workm1),
		workm2    =as.double(workm2),
		workmh1   =as.double(workmh1),
		workv1    =as.double(workv1),
		workv2    =as.double(workv2),
 		grid       =as.double(grid),
		PACKAGE="DPpackage")	
         }
 

         #########################################################################################
         # save state
         #########################################################################################
	
 	 thetasave<-matrix(foo$thetasave,nrow=nsave, ncol=(p+2))
 	 randsave<-matrix(foo$randsave,nrow=nsave, ncol=(nrec+1))
 	 
	 colnames(thetasave)<-c(dimnames(x)[[2]],"sigma2","alpha")

         qnames<-NULL
         for(i in 1:nrec){
             idname<-paste("(Subject",i,sep="=")
             idname<-paste(idname,")",sep="")
             qnames<-c(qnames,idname)
         }
         qnames<-c(qnames,"Prediction")
         
         colnames(randsave)<-qnames

	 model.name<-"Bayesian Semiparametric Median Regression Model"
	 coeff<-rep(0,(p+2))
	 for(i in 1:(p+2)){
	 	coeff[i]<-mean(thetasave[,i])
	 }
	 names(coeff)<-c(dimnames(x)[[2]],"sigma2","alpha")
	
	 state <- list(alpha=foo$alpha,beta=foo$beta,v=foo$v,sigma2=foo$sigma2)	
		      
	 save.state <- list(thetasave=thetasave,randsave=randsave)

         if(is.null(prior$a0))
	 {
            acrate<-foo$acrate[1:2]
         }
         else
         {
            acrate<-foo$acrate
         }   

	 z<-list(modelname=model.name,coefficients=coeff,acrate=acrate,call=cl,
	         prior=prior,mcmc=mcmc,state=state,save.state=save.state,cpo=foo$cpo,
	         nrec=nrec,p=p,dens=foo$f,grid=grid)
	
	 cat("\n\n")
	 class(z)<-c("PTlm")
	 z 
}


###
### Tools for PTlm: print, summary, plot
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 08-10-2006.
	

"print.PTlm"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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


"plot.PTlm"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...) 
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


    if(is(x, "PTlm")){
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
           plot(x$grid,x$dens,ylab="density",main=title1,lty=1,type='l',lwd=2,xlab="values")
           
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

            if (poss==0 && param !="predictive") 
	    {
	      stop("This parameter is not present in the original model.\n")
	    }

            if (param !="predictive") 
            {
	       par(ask = ask)
	       layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
               title1<-paste("Trace of",pnames[poss],sep=" ")
               title2<-paste("Density of",pnames[poss],sep=" ")       
               plot(x$save.state$thetasave[,poss],type='l',main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot(x$save.state$thetasave[,poss],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
            }
            else
            {
               title1<-c("Predictive Error Density")
               plot(x$grid,x$dens,ylab="density",main=title1,lty=1,type='l',lwd=2,xlab="values")
            }
        }
   }
}




"summary.PTlm"<-function(object, hpd=TRUE, ...) 
{
    dimen<-object$p+1
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
	
    }
    else
    {
         dimen<-1
         coef.p<-object$coefficients[(object$p+2)]
         coef.sd<-rep(0,dimen)
         coef.se<-rep(0,dimen)
         coef.l<-rep(0,dimen)
         coef.u<-rep(0,dimen)
         coef.m<-rep(0,dimen)
         names(coef.sd)<-names(object$coefficients[(object$p+2)])
         names(coef.l)<-names(object$coefficients[(object$p+2)])
         names(coef.u)<-names(object$coefficients[(object$p+2)])
         for(i in 1:dimen){
             alow<-rep(0,2)
             aupp<-rep(0,2)
             coef.sd[i]<-sqrt(var(object$save.state$thetasave[,(object$p+1+i)]))
             coef.m[i]<-median(object$save.state$thetasave[,(object$p+1+i)])
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
    }


    ans$acrate<-object$acrate
    
    ans$nrec<-object$nrec

    class(ans) <- "summaryPTlm"
    return(ans)
}



"print.summaryPTlm"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(x$cpo)), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
    if (length(x$coefficients)) 
    {
        cat("\nRegression coefficients:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")

    if (length(x$prec)) 
    {
        cat("\nPrecision parameter:\n")
        print.default(format(x$prec, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat(" \n")

    cat("\nAcceptance Rate for Metropolis Step = ",x$acrate,"\n")    
    
    cat("\nNumber of Observations:",x$nrec)
    cat("\n\n")
    invisible(x)
}


