### DPrasch.R                   
### Fit a Rasch model with a mixture of Dirichlet process prior
### for the random effect distribution
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 24-09-2006.
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


DPrasch<-function(y,prior,mcmc,state,status,delete=FALSE,data=sys.frame(sys.parent()))
UseMethod("DPrasch")

DPrasch.default<-function(y,prior,mcmc,state,status,delete=FALSE,data)
{
         #########################################################################################
         # call parameters
         #########################################################################################
         cl <- match.call()
	 y<-as.matrix(y)
	  
         #########################################################################################
         # data structure
         #########################################################################################
     	 nsubject<-dim(y)[1]
	 p<-dim(y)[2]

         if(delete==TRUE)
         {
              # delete missing values
              
                y<-as.matrix(na.omit(y))
            	nsubject<-dim(y)[1]
	        p<-dim(y)[2]

              # delete according to criteria

                ok<-0
                while(ok==0)  
                {
                   yfin<-NULL
                   nused<-0
                   for(i in 1:nsubject)
                   {
                       if(sum(y[i,])<p && sum(y[i,])>0)
                       {
                           nused<-nused+1
                           yfin<-rbind(yfin,y[i,])
                       }
                   } 

                   ans<-apply(yfin, 1, sum) 
                   yfin2<-NULL
                   iused<-0  
                   for(i in 1:p)
                   {
                       if(ans[i]<nused && ans[i]>0)
                       {
                          iused<-iused+1
                          yfin2<-cbind(yfin2,yfin[,i])
                       }
                   }          
            
                   if(nused==nsubject && iused==p)
                   {
                      ok<-1
                      nsubject<-nused
                      p<-iused
                      y<-yfin2
                   }  
                   else
                   {
                      nsubject<-nused
                      p<-iused
                      y<-yfin2
                   } 
                }  
                
                ywork<-y
                imiss<-0
                nmissing<-1
                datastrm<-matrix(0,nrow=1,ncol=2)
                nrec<-nsubject*p-nmissing
         }
         else
         {
               ywork<-y
               datastrm<-NULL
               nmissing<-0
               total<-0
               for(i in 1:nsubject)
               {
                  for(j in 1:p)
                  {
                     if(is.na(y[i,j]))
                     {
                         nmissing<-nmissing+1
                         datastrm<-rbind(datastrm,c(i,j))   
                     }
                     else
                     {
                         total<-total+y[i,j]            
                     }
                  
                  }
               }
               if(nmissing>0)
               {
                  imiss<-1 
                  for(i in 1:nmissing)
                  {
                      ywork[datastrm[i,1],datastrm[i,2]]<-rbinom(1,1,total/nrec)               
                  }
                  nrec<-nsubject*p-nmissing
               }
               else
               {
                  imiss<-0
                  nmissing<-1
                  datastrm<-matrix(0,nrow=1,ncol=2)
                  nrec<-nsubject*p
               }
         }
         

         #########################################################################################
         # prior information
         #########################################################################################

  	 if(is.null(prior$a0))
  	 {
  	    a0<--1
  	    b0<--1 
  	    alpha<-prior$alpha
  	    alpharand<-0
  	 }
         else
         {
            a0<-prior$a0
  	    b0<-prior$b0
  	    alpha<-rgamma(1,shape=a0,scale=b0)
  	    alpharand<-1
  	 }
  	 a0b0<-c(a0,b0)
  	 

  	 if(is.null(prior$tau1))
  	 {
              tau1<--1
              tau2<--1
              sigma<-prior$sigma
              sigmainv<-1/sigma              
              sigmarand<-0
  	 }
  	 else
  	 {
              tau1<-prior$tau1
              tau2<-prior$tau2
              sigma<-1
              sigmainv<-1/sigma
  	      sigmarand<-1
  	 }

  	 if(is.null(prior$mub))
  	 {
  	      psiinv<--1
  	      smu<-0
  	      mu<-prior$mu
  	      murand<-0
  	 }
  	 else
  	 {
  	      psiinv<-(1/prior$Sb)
  	      smu<-psiinv*prior$mub
  	      mu<-rnorm(1,prior$mub,sqrt(prior$Sb))
              murand<-1
         }     

         b0<-prior$beta0
         prec<-solve(prior$Sbeta0)
         
         if(dim(prec)[1] != (p-1)) stop("the dimension of beta0 and Sbeta0 must be p-1")


         #########################################################################################
         # mcmc specification
         #########################################################################################
         mcmcvec<-c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
         nsave<-mcmc$nsave

         #########################################################################################
         # output
         #########################################################################################
         acrate<-rep(0,2)
         cpo<-matrix(0,nrow=nsubject,ncol=p)
         ngrid<-500
	 grid<-rev(seq(-6,6,length=ngrid))
         faccum<-rep(0,ngrid)
         thetasave<-matrix(0,nrow=nsave,ncol=p+3)
         randsave<-matrix(0,nrow=nsave,ncol=nsubject+1)


         #########################################################################################
         # MLE estimation
         #########################################################################################
         
         RaschMLE<-function(y,p,nsubject)
         {
             ywork2<-rep(0,nsubject*p)
             x<-matrix(0,nrow=nsubject*p,ncol=p-1)
             count<-0
             for(i in 1:nsubject)
             {
                 for(j in 1:p)
                 {
                    count<-count+1
                    ywork2[count]<-y[i,j]
                    if(j>1)
                    {
                      x[count,j-1]<--1
                    }
                  }
             }
             out<-NULL

             fit0<-glm.fit(x=x, y=ywork2, family = binomial(logit))
             out$beta<-fit0$coefficients
             p1 <- 1:(p-1)
             Qr <- fit0$qr
             out$cov <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
             return(out)
	 }

         fit0<-RaschMLE(ywork,p,nsubject)

  	 if(is.null(mcmc$tune1))
  	 {
  	    tune1<-1.1
  	 }
  	 else
  	 {
  	    tune1<-mcmc$tune1
  	 }

         propvf = diag(tune1,p-1)%*% solve(prec + solve(fit0$cov))%*% diag(tune1,p-1)


         cohen<-function(y,nsubject,p) 
         # Cohen's (1976) approximation for asymptotic variance
         # COHEN, L. A modified logistic response model for item analysis. Unpublished manuscript, 1976.
         {
                ywork2<-y
                ok<-0
                while(ok==0)  
                {
                   yfin<-NULL
                   nused<-0
                   for(i in 1:nsubject)
                   {
                       if(sum(ywork2[i,])<p && sum(ywork2[i,])>0)
                       {
                           nused<-nused+1
                           yfin<-rbind(yfin,ywork2[i,])
                       }
                   } 

                   ans<-apply(yfin, 1, sum) 
                   yfin2<-NULL
                   iused<-0  
                   for(i in 1:p)
                   {
                       if(ans[i]<nused && ans[i]>0)
                       {
                          iused<-iused+1
                          yfin2<-cbind(yfin2,yfin[,i])
                       }
                   }          
            
                   if(nused==nsubject && iused==p)
                   {
                      ok<-1
                      nsubject<-nused
                      p<-iused
                      ywork2<-yfin2
                   }  
                   else
                   {
                      nsubject<-nused
                      p<-iused
                      ywork2<-yfin2
                   } 
                }  

                ss1<-apply(ywork2, 2, sum)
                nrss<-apply(ywork2, 1, sum) 
                nr<-table(nrss)
                xi<-log((nsubject-ss1)/ss1)
                U<-var(xi)
                yr<-log(nrss/(p-nrss))
                V<-var(yr)
                X <- sqrt((1+U/2.89)/(1-U*V/8.35))
                Y <- sqrt((1+V/2.89)/(1-U*V/8.35))
         
                VV<-( (X**2)*( p/ (nrss*(p-nrss)) ) )
                VV2<-(1+(1/VV))
                return(VV2)
         }


  	 if(is.null(mcmc$tune2))
  	 {
  	    tune2<-1.1
  	 }
  	 else
  	 {
  	    tune2<-mcmc$tune2
  	 }
  	 
  	 VV2<-cohen(ywork,nsubject,p) 
         
         propvr<-mean(tune2*(1/VV2)*tune2)
         
         #########################################################################################
         # parameters depending on status
         #########################################################################################

	 if(status==TRUE)
	 {
                beta<-fit0$beta
                b<-rnorm(nsubject,mean=mu,sd=sigma)
                ss<-seq(1,nsubject)
                sigmainv<-1/sigma
                bclus<-b
                ncluster<-nsubject
   	 }
	 
      	 if(status==FALSE)
	 {
	        alpha<-state$alpha
	        b<-state$b
                bclus<-state$bclus 
                beta<-state$beta
                mu<-state$mu
	        ncluster<-state$ncluster
	        sigma<-state$sigma
	        sigmainv<-1/sigma
	        ss<-state$ss
	 }    

         #########################################################################################
         # working space
         #########################################################################################
        
         betac<-rep(0,p-1)
         ccluster<-rep(0,nsubject)
         fsavet<-rep(0,ngrid)
         prob<-rep(0,nsubject+2)
         seed1<-sample(1:29000,1)
         seed2<-sample(1:29000,1)

         workmh1<-rep(0,(p-1)*p/2)
         workv1<-rep(0,p-1)


         #########################################################################################
         # calling the fortran code
         #########################################################################################

         foo <- .Fortran("sprasch",
  	 	datastr    =as.integer(datastrm),
  	 	imiss      =as.integer(imiss),
  	 	ngrid      =as.integer(ngrid),
  	 	nmissing   =as.integer(nmissing),
  	 	nsubject   =as.integer(nsubject),
	 	p          =as.integer(p),
  	 	y          =as.integer(ywork),
  	 	a0b0       =as.double(a0b0),
  	 	b0         =as.double(b0),
  	 	prec       =as.double(prec),
  	 	psiinv     =as.double(psiinv),
 		smu        =as.double(smu),
  	 	tau1       =as.double(tau1),
  	 	tau2       =as.double(tau2),
 		mcmc       =as.integer(mcmcvec),
 		nsave      =as.integer(nsave),
                acrate     =as.double(acrate),
 		cpo        =as.double(cpo),
 		faccum     =as.double(faccum),
 		randsave   =as.double(randsave),
 		thetasave  =as.double(thetasave),
 		alpha      =as.double(alpha),
 		b          =as.double(b),
		bclus      =as.double(bclus),
 		beta       =as.double(beta),
 		mu         =as.double(mu),
                ncluster   =as.integer(ncluster),
 		sigma      =as.double(sigma),
 		sigmainv   =as.double(sigmainv),
                ss         =as.integer(ss),
 		betac      =as.double(betac),
                ccluster   =as.integer(ccluster),
                fsavet     =as.integer(fsavet),
 		propvf     =as.double(propvf),
 		propvr     =as.double(propvr),
 		prob       =as.double(prob),
 		seed1      =as.integer(seed1),
 		seed2      =as.integer(seed2),
 		workmh1    =as.double(workmh1),
 		workv1     =as.double(workv1),
 		grid       =as.double(grid),
		PACKAGE    ="DPpackage")


         #########################################################################################
         # save state
         #########################################################################################

         model.name<-"Bayesian Semiparametric Rasch Model"		
                
         state <- list(alpha=foo$alpha,
	               b=foo$b,
	               bclus=foo$bclus,
	               beta=foo$beta,
	               mu=foo$mu,
	               ncluster=foo$ncluster,
	               sigma=foo$sigma,
                       ss=foo$ss
                       )

         cpo<-matrix(foo$cpo,nrow=nsubject,ncol=p)
         randsave<-matrix(foo$randsave,nrow=nsave,ncol=nsubject+1)
         thetasave<-matrix(foo$thetasave,nrow=nsave,ncol=p+3)
 
         pnames<-NULL
         for(i in 2:p)
         {
             pnames<-c(pnames,paste("beta",i,sep=""))
         }
         pnames<-c(pnames,"mu","sigma","ncluster","alpha")
         colnames(thetasave)<-pnames
         
         qnames<-NULL
         for(i in 1:nsubject)
         {
             temp<-paste("theta(ID=",i,sep="")
             temp<-paste(temp,")",sep="")
             qnames<-c(qnames,temp)
         }
         qnames<-c(qnames,"theta(Prediction)")
         dimnames(randsave)<-list(NULL,qnames)
         
         coeff<-apply(thetasave, 2, mean)
         
         save.state <- list(thetasave=thetasave,randsave=randsave)

         
         acrate<-foo$acrate
         
	 z<-list(call=cl,y=ywork,modelname=model.name,cpo=cpo,
                 prior=prior,mcmc=mcmc,state=state,save.state=save.state,nrec=nrec,
                 nsubject=nsubject,p=p,alpharand=alpharand,sigmarand=sigmarand,murand=murand,acrate=acrate,
                 coefficients=coeff,cdf=foo$faccum,grid=grid)
                 
         cat("\n\n")
 	 class(z)<-c("DPrasch")
  	 return(z)
}



###                    
### Tools
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 24-09-2006.
###


"print.DPrasch"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")

    if (length(x$coefficients)) {
        cat("Posterior Inference of Parameters:\n")
        if(x$alpharand==1){
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)}
        if(x$alpharand==0){
        print.default(format(x$coefficients[1:(length(x$coefficients)-1)], digits = digits), print.gap = 2, 
            quote = FALSE)}

    }

    cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate,"\n")    
   
    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Groups:",x$nsubject,"\n")    
    cat("\n\n")
    invisible(x)
}



"summary.DPrasch"<-function(object, hpd=TRUE, ...) 
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

    #nsave<-object$nsave
    #dimen<-length(object$coefficients)
    #thetasave<-matrix(object$save.state$thetasave,nrow=nsave, ncol=dimen)
    thetasave<-object$save.state$thetasave


### Difficulty parameters

    dimen1<-object$p-1

    mat<-thetasave[,1:dimen1]

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

    dimen2<-object$murand+object$sigmarand
 
    if(dimen2>0){
    
       if(dimen2==1)
       {
          if(object$murand==1)
          {
             mat<-matrix(thetasave[,dimen1+1],ncol=1) 
             coef.p<-object$coefficients[dimen1+1]
          }
          else
          {
             mat<-matrix(thetasave[,dimen1+2],ncol=1) 
             coef.p<-object$coefficients[dimen1+2]
          }
       }
       else
       {
          mat<-thetasave[,(dimen1+1):(dimen1+2)]
          coef.p<-object$coefficients[(dimen1+1):(dimen1+2)]
       }
    
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

       ans$base<-coef.table
    }

### Precision parameter

    if(is.null(object$prior$a0))
    {
      dimen3<-1
      coef.p<-object$coefficients[(dimen1+2+1)]
      mat<-matrix(thetasave[,(dimen1+2+1)],ncol=1)
    }
    else
    {
      dimen3<-2
      coef.p<-object$coefficients[(dimen1+2+1):(dimen1+2+2)]
      mat<-thetasave[,(dimen1+2+1):(dimen1+2+2)]

    }  

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
    ans$nrec<-object$nrec
    ans$nsubject<-object$nsubject
    ans$acrate<-object$acrate

    class(ans) <- "summaryDPraschpoisson"
    return(ans)
}


"print.summaryDPrasch"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(as.vector(x$cpo))), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
    if (length(x$coefficients)) {
        cat("\nDifficulty parameters:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

    if (length(x$base)) {
        cat("\nBaseline distribution:\n")
        print.default(format(x$base, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No baseline parameters\n")

    if (length(x$prec)) {
        cat("\nPrecision parameter:\n")
        print.default(format(x$prec, digits = digits), print.gap = 2, 
            quote = FALSE)
    }

    cat("\nAcceptance Rate for Metropolis Steps = ",x$acrate,"\n")    

    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Groups:",x$nsubject,"\n")
    cat("\n\n")
    invisible(x)
}


"plot.DPrasch"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...)
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


   if(is(x, "DPrasch"))
   {
        if(is.null(param))
        {
           coef.p<-x$coefficients[1:(x$p-1)]
           n<-length(coef.p)
           pnames<-names(coef.p)
           
           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
           for(i in 1:n)
           {
               title1<-paste("Trace of",pnames[i],sep=" ")
               title2<-paste("Density of",pnames[i],sep=" ")       
               plot(ts(x$save.state$thetasave[,i]),main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot1(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }

           if(x$murand==1)
           {
               title1<-paste("Trace of","mu",sep=" ")
               title2<-paste("Density of","mu",sep=" ")       
               plot(ts(x$save.state$thetasave[,x$p]),main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot1(x$save.state$thetasave[,x$p],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }

           if(x$sigmarand==1)
           {
               title1<-paste("Trace of","sigma",sep=" ")
               title2<-paste("Density of","sigma",sep=" ")       
               plot(ts(x$save.state$thetasave[,x$p+1]),main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot1(x$save.state$thetasave[,x$p+1],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }

           title1<-paste("Trace of","ncluster",sep=" ")
           title2<-paste("Density of","ncluster",sep=" ")       
           plot(ts(x$save.state$thetasave[,x$p+2]),main=title1,xlab="MCMC scan",ylab=" ")
           hist(x$save.state$thetasave[,x$p+2],main=title2,xlab="values", ylab="probability",probability=TRUE)

           if(x$alpharand==1)
           {
               title1<-paste("Trace of","alpha",sep=" ")
               title2<-paste("Density of","alpha",sep=" ")       
               plot(ts(x$save.state$thetasave[,x$p+3]),main=title1,xlab="MCMC scan",ylab=" ")
               fancydensplot1(x$save.state$thetasave[,x$p+3],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
           }
           
           title1<-c("Predictive Density")
           title2<-c("Predictive CDF")
           #fancydensplot1(x$save.state$randsave[,(x$nsubject+1)],hpd=hpd,main=title1,xlab="theta", col=col)
           plot(x$grid,x$cdf,ylab="probability",main=title2,lty=1,type='l',lwd=2,ylim=c(0,1),xlab="theta")
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
            if(poss==0 && param !="predictive")             
	    {
	      stop("This parameter is not present in the original model.\n")
	    }
	    
	    par(ask = ask)
	    layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))

	    if(param !="predictive")
	    {
                title1<-paste("Trace of",pnames[poss],sep=" ")
                title2<-paste("Density of",pnames[poss],sep=" ")       
                plot(ts(x$save.state$thetasave[,poss]),main=title1,xlab="MCMC scan",ylab=" ")
                if(param=="ncluster")
                {
                   hist(x$save.state$thetasave[,poss],main=title2,xlab="values", ylab="probability",probability=TRUE)
                }
                else
                {
                  fancydensplot1(x$save.state$thetasave[,poss],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
                }   
                
            }
            else
            {
               title1<-c("Predictive Error Density")
               title2<-c("Link Function")
               #fancydensplot1(x$save.state$randsave[,(x$nsubject+1)],hpd=hpd,main=title1,xlab="theta", col=col)
               plot(x$grid,x$cdf,ylab="probability",main=title2,lty=1,type='l',lwd=2,ylim=c(0,1),xlab="theta")

            }                
        }
   }

}


