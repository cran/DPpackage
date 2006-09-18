### PTdensity.R                   
### Fit a Mixture of Polya trees for density estimation
###
### Copyright:Alejandro Jara and Tim Hanson, 2006
###
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


PTdensity<-function(y,ngrid=1000,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("PTdensity")

PTdensity.default<-function(y,ngrid=1000,prior,mcmc,state,status,data,na.action=na.fail)
{
         #########################################################################################
         # call parameters
         #########################################################################################
         cl <- match.call()
	 y<-na.action(as.matrix(y))	
	  
         #########################################################################################
         # data structure
         #########################################################################################
         nrec<-dim(y)[1]
         nvar<-dim(y)[2]
          

         if(nvar>1) stop("So far, this function is only for univariate density estimation")

         left<-min(y)-0.5*sqrt(var(y))
         right<-max(y)+0.5*sqrt(var(y))

         #########################################################################################
         # prior information
         #########################################################################################

  	 if(is.null(prior$a0))
  	 {
  	    ca<--1
  	    cb<--1 
  	    cpar<-prior$alpha
  	    crand<-0
  	 }
         else
         {
            ca<-prior$a0
  	    cb<-prior$b0
  	    cpar<-rgamma(1,shape=ca,scale=cb)
  	    crand<-1
  	 }
  	 ab<-c(ca,cb)
  	 
  	 
         #########################################################################################
         # mcmc specification
         #########################################################################################
         mcmcvec<-c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
         nsave<-mcmc$nsave

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


         #########################################################################################
         # output
         #########################################################################################
         acrate<-rep(0,3)
         thetasave<-matrix(0,nrow=nsave,ncol=3)
         f<-rep(0,ngrid)
         
         #########################################################################################
         # parameters depending on status
         #########################################################################################
         
    	 if(status==TRUE)
	 {
                mu<-mean(y)
                sigma<-sqrt(var(y))
   	 }
	 
      	 if(status==FALSE)
	 {
	        cpar<-state$cpar
                mu<-state$mu 
	        sigma<-state$sigma
	 }    


         #########################################################################################
         # working space
         #########################################################################################

         acrate<-rep(0,3)
	 cpo<-rep(0,nrec)
         grid<-seq(left,right,length=ngrid)  
	 seed<-c(sample(1:29000,1),sample(1:29000,1))

         #########################################################################################
         # calling the fortran code
         #########################################################################################

         if(is.null(prior$M))
         {
             foo <- .Fortran("ptdensityu",
           	ngrid      =as.integer(ngrid),
         	nrec       =as.integer(nrec),
 		y          =as.double(y),
 		ab         =as.double(ab),
		mcmcvec    =as.integer(mcmcvec),
		nsave      =as.integer(nsave),
                tune1      =as.double(tune1),
                tune2      =as.double(tune2),
                tune3      =as.double(tune3),
 		acrate     =as.double(acrate),
 		f          =as.double(f),
		thetasave  =as.double(thetasave),		
		cpo        =as.double(cpo),		
		cpar       =as.double(cpar),		
		mu         =as.double(mu),		
		sigma      =as.double(sigma),		
		grid       =as.double(grid),		
		seed      =as.integer(seed),
		PACKAGE    ="DPpackage")
	 }	
         else
         {
             nlevel<-prior$M
             ninter<-2**nlevel
             assign<-matrix(0,nrow=nrec,ncol=nlevel)
	     accums<-matrix(0,nrow=nlevel,ncol=ninter)
             counter<-matrix(0,nrow=nlevel,ncol=ninter)
             endp<-rep(0,ninter-1)
 	     intpn<-rep(0,nrec)
	     intpo<-rep(0,nrec)
	     prob<-rep(0,ninter)
             rvecs<-matrix(0,nrow=nlevel,ncol=ninter)

             foo <- .Fortran("ptdensityup",
         	ngrid      =as.integer(ngrid),
         	nrec       =as.integer(nrec),
 		y          =as.double(y),
 		ab         =as.double(ab),
 		nlevel     =as.integer(nlevel),
 		ninter     =as.integer(ninter),
		mcmcvec    =as.integer(mcmcvec),
		nsave      =as.integer(nsave),
                tune1      =as.double(tune1),
                tune2      =as.double(tune2),
                tune3      =as.double(tune3),
 		acrate     =as.double(acrate),
 		f          =as.double(f),
		thetasave  =as.double(thetasave),		
		cpo        =as.double(cpo),		
		cpar       =as.double(cpar),		
		mu         =as.double(mu),		
		sigma      =as.double(sigma),		
		grid       =as.double(grid),		
		intpn     =as.integer(intpn),		
		intpo     =as.integer(intpo),		
		accums    =as.double(accums),
		assign    =as.integer(assign),
		counter   =as.integer(counter),
		endp      =as.double(endp),
		prob      =as.double(prob),
		rvecs     =as.double(rvecs),
		seed      =as.integer(seed),
		PACKAGE    ="DPpackage")
	 }

         #########################################################################################
         # save state
         #########################################################################################
         model.name<-"Bayesian semiparametric density estimation"		
                
         varnames<-colnames(y)
         if(is.null(varnames))
         {
               varnames<-all.vars(cl)[1:nvar]
         }
		
         state <- list(
                       cpar=foo$cpar,
	               mu=foo$mu,
	               sigma=foo$sigma
                      )
         thetasave<-matrix(foo$thetasave,nrow=nsave,ncol=3)
 
         coeff<-apply(thetasave,2,mean) 
         names(coeff)<-c("mu","sigma","alpha")
         colnames(thetasave)<-c("mu","sigma","alpha")
             
         save.state <- list(thetasave=thetasave)
         
         if(crand==0)
         {
            acrate<-foo$acrate[1:2]
         }
         else
         {
            acrate<-foo$acrate
         }
         
	 z<-list(call=cl,y=y,varnames=varnames,modelname=model.name,cpo=foo$cpo,
                 prior=prior,mcmc=mcmc,state=state,save.state=save.state,nrec=foo$nrec,
                 nvar=foo$nvar,crand=crand,coefficients=coeff,f=foo$f,
                 grid=foo$grid,acrate=acrate,nvar=nvar,x1=foo$grid,dens=foo$f)
                 
         cat("\n\n")
 	 class(z)<-"PTdensity"
  	 return(z)
}


###                    
### Tools
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 16-09-2006.
###


"print.PTdensity"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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


"summary.PTdensity"<-function(object, hpd=TRUE, ...) 
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


    if(object$crand==0)
    {
       mat<-matrix(thetasave[,1:2],ncol=2) 
       dimen<-2
    }
    else
    {
       mat<-thetasave[,1:3]
       dimen<-3       
    }

    coef.p<-object$coefficients[1:dimen]
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

    names(coef.m)<-names(object$coefficients[1:dimen])
    names(coef.sd)<-names(object$coefficients[1:dimen])
    names(coef.se)<-names(object$coefficients[1:dimen])
    names(coef.l)<-names(object$coefficients[1:dimen])
    names(coef.u)<-names(object$coefficients[1:dimen])

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

    ans$acrate<-object$acrate
    
    ans$nrec<-object$nrec

    class(ans) <- "summaryPTdensity"
    return(ans)
}


"print.summaryPTdensity"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(x$cpo)), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
    if (length(x$coefficients)) {
        cat("\nBaseline parameters:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")

    cat("\nAcceptance Rate for Metropolis Step = ",x$acrate,"\n")    
    
    cat("\nNumber of Observations:",x$nrec)
    cat("\n\n")
    invisible(x)
}


"plot.PTdensity"<-function(x, ask=TRUE, output="density", param=NULL, hpd=TRUE, nfigr=1, nfigc=1, col="#bdfcc9", ...) 
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


   if(is(x, "PTdensity"))
   {

      if(output=="density")
      {

      # Density estimation
	
	par(ask = ask)
	layout(matrix(seq(1,nfigr*nfigc,1),nrow=nfigr,ncol=nfigc,byrow=TRUE))
	
        title1<-paste("Density of",x$varnames,sep=' ')
    
	aa<-hist(x$y,plot=F)
	maxx<-max(aa$intensities+aa$density)+0.1*max(aa$intensities+aa$density)
	miny<-min(x$y)
	maxy<-max(x$y)
	deltay<-(maxy-miny)*0.2
	miny<-miny-deltay
	maxy<-maxy+deltay
	hist(x$y,probability=T,xlim=c(min(x$grid),max(x$grid)),ylim=c(0,maxx),nclas=25,main=title1,xlab="values", ylab="density")	
        lines(x$grid,x$f,lwd=2)
	
	
      }
      else
      {

        if(is.null(param))
        {
           pnames<-colnames(x$save.state$thetasave)
           n<-dim(x$save.state$thetasave)[2]
           cnames<-names(x$coefficients)

           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
           for(i in 1:(n-1))
           {
               title1<-paste("Trace of",pnames[i],sep=" ")
               title2<-paste("Density of",pnames[i],sep=" ")       
               plot(x$save.state$thetasave[,i],type='l',main=title1,xlab="MCMC scan",ylab=" ")
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
               plot(x$save.state$thetasave[,n],type='l',main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot1(x$save.state$thetasave[,n],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }
        }
   
        else
        {

            pnames<-colnames(x$save.state$thetasave)
            n<-dim(x$save.state$thetasave)[2]
	    poss<-0 
            for(i in 1:n)
            {
               if(pnames[i]==param)poss=i
            }
            if (poss==0) 
	    {
	      stop("This parameter is not present in the original model.\n")
	    }
	    
	    par(ask = ask)
	    layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
            title1<-paste("Trace of",pnames[poss],sep=" ")
            title2<-paste("Density of",pnames[poss],sep=" ")       
            plot(x$save.state$thetasave[,poss],type='l',main=title1,xlab="MCMC scan",ylab=" ")
            fancydensplot1(x$save.state$thetasave[,poss],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
        }


      }	
   }
}





