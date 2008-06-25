### LDDPdensity.R                   
### Fit a Linear Dependent Dirichlet Process Mixture of Normals Model
###
### Copyright: Alejandro Jara, Peter Mueller and Gary L. Rosner, 2008
### Last modification: 02-06-2008.
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
###      Alejandro Jara
###      Biostatistical Centre
###      Katholieke Universiteit Leuven
###      U.Z. Sint-Rafaël
###      Kapucijnenvoer 35
###      B-3000 Leuven
###      Voice: +32 (0)16 336892  URL  : http://student.kuleuven.be/~s0166452/
###      Fax  : +32 (0)16 337015  Email: Alejandro.JaraVallejos@med.kuleuven.be
###
###      Peter Mueller
###      Department of Biostatistics
###      The University of Texas MD Anderson Cancer Center
###      1515 Holcombe Blvd, Unit 447 
###      Houston TX 77030-4009, USA
###      Voice: (713) 563-4296  URL  : http://www.mdanderson.org/departments/biostats
###      Fax  : (713) 563-4243  Email: pmueller@mdanderson.org
###
###      Gary L. Rosner
###      Department of Biostatistics
###      The University of Texas MD Anderson Cancer Center
###      1515 Holcombe Blvd, Unit 447 
###      Houston TX 77030-4009, USA
###      Voice: (713) 563-4285  URL  : http://www.mdanderson.org/departments/biostats
###      Fax  : (713) 563-4243  Email: glrosner@mdanderson.org
###

"LDDPdensity"<-
function(formula,prior,mcmc,state,status,grid=seq(-10,10,length=1000),xpred,data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("LDDPdensity")

"LDDPdensity.default"<-
function(formula, 
         prior,
         mcmc,
         state,
         status, 
         grid=seq(-10,10,length=1000),
         xpred,
         data=sys.frame(sys.parent()),
         na.action=na.fail)
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

  	 y <- model.response(mf,"numeric")
	 nrec <- length(y)
	 x <- as.matrix(model.matrix(formula))
	 p <- dim(x)[2]

       #########################################################################################
       # prediction
       #########################################################################################
         if(dim(xpred)[2] != p) stop("the dimension of xpred must be npred*p")
         npred <- dim(xpred)[1]
	 ngrid <- length(grid)


       #########################################################################################
       # prior information
       #########################################################################################

  	 if(is.null(prior$a0))
  	 {
  	    a0 <--1
  	    b0 <--1 
  	    alpha <- prior$alpha
  	    alpharand <- 0
  	 }
         else
         {
            a0 <- prior$a0
  	    b0 <- prior$b0
  	    alpha <-rgamma(1,shape=a0,rate=b0)
  	    alpharand<-1
  	 }
  	 a0b0 <- c(a0,b0)
  	 
  	 if(is.null(prior$tau1))
  	 {
              tau1 <--1
              tau2 <--1
              sigma2k <- prior$sigma2k
              s2krand <- 0
  	 }
  	 else
  	 {
              tau1 <- prior$tau1
              tau2 <- prior$tau2
              sigma2k <- var(y)
              s2krand <-1
  	 }
  	 tau <- c(tau1,tau2)
  	 

  	 if(is.null(prior$nu))
  	 {
              nu <--1
              psiinv <- diag(1,p)
              sigmab <- prior$sigmab
              sigmabrand <- 0
  	 }
  	 else
  	 {
              nu <- prior$nu
              psiinv <- prior$psiinv
              sigmab <- diag(1,p)
              sigmabrand <- 1
  	 }
  	   	 
  	 if(is.null(prior$m))
  	 {
              sm <- rep(0,p)
              sinv <- matrix(0,nrow=p,ncol=p) 
              mub <- prior$mub
  	      murand <- 0
  	 }
  	 else
  	 {
              m <- prior$m
              sinv <-solve(prior$s)
              sm <- sinv%*%m
              mub <- rep(0,p)
              murand <- 1
         }     

       #########################################################################################
       # mcmc specification
       #########################################################################################
         mcmcvec<-c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
         nsave<-mcmc$nsave

       #########################################################################################
       # output
       #########################################################################################
         cpo <- matrix(0,nrow=nrec,ncol=2)
         densr <- matrix(0,nrow=npred,ncol=ngrid)
         thetasave <- matrix(0,nrow=nsave,ncol=(p+(p*(p+1)/2)+3))
         randsave <- matrix(0,nrow=nsave,ncol=nrec*p)

       #########################################################################################
       # parameters depending on status
       #########################################################################################

	 if(status==TRUE)
	 {
                b <- matrix(0,nrow=nrec,ncol=p)
                ncluster <- 1
                ss <- rep(1,nrec)
   	 }

      	 if(status==FALSE)
	 {
                ncluster <- state$ncluster	        
                ss <- state$ss
	        alpha <- state$alpha
	        b <- state$b
	        mub <- state$mub
	        sigmab <- state$sigmab
	        sigma2k <- state$sigma2k
	 }

       #########################################################################################
       # working space
       #########################################################################################

         iflagc <- rep(0,p)
         seed <- c(sample(1:29000,1),sample(1:29000,1))
         sigmabinv <- matrix(0,nrow=p,ncol=p)
         theta <- rep(0,p)
         workmc <- matrix(0,nrow=p,ncol=p)
         workmc2 <- matrix(0,nrow=p,ncol=p)
         workvc <- rep(0,p)
         workmhc <- rep(0,p*(p+1)/2)
         workmhc2 <- rep(0,p*(p+1)/2)
         ztz <- matrix(0,nrow=p,ncol=p)
         zty <- rep(0,p)

         prob <- rep(0,nrec)
         ccluster <- rep(0,nrec)
         cstrt <- matrix(0,nrow=nrec,ncol=nrec) 
          
       #########################################################################################
       # calling the fortran code
       #########################################################################################

         foo <- .Fortran("lddpmnormals",
 	 	nrec       =as.integer(nrec),
 	 	p          =as.integer(p),
  	 	y          =as.double(y),
  	 	x          =as.double(x),
                npred      =as.integer(npred),
                xpred      =as.double(xpred),
                ngrid      =as.integer(ngrid),
                grid       =as.double(grid),
                a0b0       =as.double(a0b0),
                tau        =as.double(tau),
                murand     =as.integer(murand),
                sm         =as.double(sm),
                sinv       =as.double(sinv),
                nu         =as.integer(nu),
                psiinv     =as.double(psiinv),
                mcmc       =as.integer(mcmcvec),
                nsave      =as.integer(nsave),
                cpo        =as.double(cpo),
                densr      =as.double(densr),
                randsave   =as.double(randsave),
                thetasave  =as.double(thetasave),
                alpha      =as.double(alpha),
                b          =as.double(b),
                mub        =as.double(mub),
                sigmab     =as.double(sigmab),
                sigma2k    =as.double(sigma2k),
                ncluster   =as.integer(ncluster),
                ss         =as.integer(ss),
                ccluster   =as.integer(ccluster),
                cstrt      =as.integer(cstrt),
                iflagc     =as.integer(iflagc),
 		prob       =as.double(prob),
 		seed       =as.integer(seed),
                theta      =as.double(theta),
                ztz        =as.double(ztz),
                zty        =as.double(zty),
                sigmabinv  =as.double(sigmabinv),
 		workmc     =as.double(workmc),
 		workmc2    =as.double(workmc2),
 		workvc     =as.double(workvc),
 		workmhc    =as.double(workmhc),
 		workmhc2   =as.double(workmhc2),
		PACKAGE    ="DPpackage")

       #########################################################################################
       # save state
       #########################################################################################

         model.name <- "Linear Dependent Dirichlet Process Mixture Model"		
                
         state <- list(ncluster=foo$ncluster,
                       ss=foo$ss,
                       alpha=foo$alpha,
                       b=matrix(foo$b,nrow=nrec,ncol=p),
	               mub=foo$mub,
	               sigmab=matrix(foo$sigmab,nrow=p,ncol=p),
	               sigma2k=foo$sigma2k
                       )

         densr <- matrix(foo$densr,nrow=npred,ncol=ngrid)
         cpom <- matrix(foo$cpo,nrow=nrec,ncol=2)
         cpo <- cpom[,1]         
         fso <- cpom[,2]
         randsave <- matrix(foo$randsave,nrow=nsave,ncol=nrec*p)
         thetasave <- matrix(foo$thetasave,nrow=nsave,ncol=(p+(p*(p+1)/2)+3))

         coeffname <- dimnames(x)[[2]]

         pnames1 <- NULL
         for(i in 1:p)
         {
             pnames1 <- c(pnames1,paste("mub",coeffname[i],sep=""))
         }
         
         pnames2 <- NULL
         for(i in 1:p)
         {
            for(j in i:p)
            {
                tmp <- paste("sigmab",coeffname[i],sep="")
                tmp <- paste(tmp,coeffname[j],sep=":")
                pnames2 <- c(pnames2,tmp)
            }    
         }

         pnames <- c(pnames1,pnames2,"sigma2k","ncluster","alpha")
         colnames(thetasave) <- pnames


         coeff <- apply(thetasave, 2, mean)

         renames <- NULL
         for(i in 1:nrec)
         {
             tmp <- paste(coeffname,i,sep=":")
             renames <- c(renames,tmp)
         }
         colnames(randsave) <- renames

         
         save.state <- list(densr=densr,
                            thetasave=thetasave,
                            randsave=randsave)

	 z<-list(call=cl,
	         y=y,
	         modelname=model.name,
	         cpo=cpo,
                 fso=fso, 
                 prior=prior,
                 mcmc=mcmc,
                 state=state,
                 save.state=save.state,
                 nrec=nrec,
                 npred=npred,
                 p=p,
                 grid=grid,
                 dens=densr,
                 coefficients=coeff)
                 
         cat("\n\n")
 	 class(z)<-c("LDDPdensity")
  	 return(z)
}


###                    
### Tools
###
### Copyright: Alejandro Jara, 2008
### Last modification: 02-06-2008.
###


"print.LDDPdensity" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
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

    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Predictors:",x$p,"\n")    
    cat("\n\n")
    invisible(x)
}


"summary.LDDPdensity"<-function(object, hpd=TRUE, ...) 
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

    ans <- c(object[c("call", "modelname")])

### CPO
    ans$cpo<-object$cpo

### Baseline Information

    mat<-NULL
    coef.p<-NULL
    
    dimen1 <- object$p+object$p*(object$p+1)/2+1
    
    coef.p <- object$coefficients[1:dimen1]
    mat <- thetasave[,1:dimen1]

    if(dimen1>0)
    {
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

    dimen1<-object$p+object$p*(object$p+1)/2
    
    if(is.null(object$prior$a0))
    {
      dimen2<-1
      coef.p<-object$coefficients[(dimen1+1)]
      mat<-matrix(thetasave[,(dimen1+1)],ncol=1)
    }
    else
    {
      dimen2<-2
      coef.p<-object$coefficients[(dimen1+1):(dimen1+2)]
      mat<-thetasave[,(dimen1+1):(dimen1+2)]

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
    ans$nvar<-object$p

    class(ans) <- "summaryLDDPdensity"
    return(ans)
}


"print.summaryLDDPdensity"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")
    	     
    cat("Posterior Predictive Distributions (log):\n")	     
    print.default(format(summary(log(as.vector(x$cpo))), digits = digits), print.gap = 2, 
            quote = FALSE) 
            
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

    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Predictors:",x$p,"\n")        
    cat("\n\n")
    invisible(x)
}



"plot.LDDPdensity"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...)
{

fancydensplot1<-function(x, hpd=TRUE, npts=200, xlab="", ylab="", main="",col="#bdfcc9", ...)
# Author: AJV, 2007
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


   if(is(x, "LDDPdensity"))
   {
        if(is.null(param))
        {
           coef.p <- x$coefficients
           n <- length(coef.p)
           pnames <- names(coef.p)
           
           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
           for(i in 1:length(coef.p))
           {
               title1<-paste("Trace of",pnames[i],sep=" ")
               title2<-paste("Density of",pnames[i],sep=" ")       
               plot(ts(x$save.state$thetasave[,i]),main=title1,xlab="MCMC scan",ylab=" ")
               if(pnames[i]=="ncluster")
               {
                  hist(x$save.state$thetasave[,i],main=title2,xlab="values", ylab="probability",probability=TRUE)
               }
               else
               {
                  fancydensplot1(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
               }
           }

           for(i in 1:x$npred)
           {
               title1<-c("Predictive Density")           
               plot(x$grid,x$save.state$densr[i,],main=title1,lty=1,type='l',lwd=2,xlab="y",ylab="density")
           }
           
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
               for(i in 1:x$npred)
               {
                   title1<-c("Predictive Density")           
                   plot(x$grid,x$save.state$densr[i,],main=title1,lty=1,type='l',lwd=2,xlab="y",ylab="density")
               }
            }                
        }
   }

}





