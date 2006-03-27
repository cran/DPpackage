###                    
### Estimate the survival cur base on a fitted semiparametric aft model 
### for interval censored data.
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 05-05-2006.

"DPsurvpred"<-
function(object,grid,xnew=NULL,hpd=TRUE)
UseMethod("DPsurvpred")

"DPsurvpred.default"<-
function(object,grid,xnew=NULL,hpd=TRUE)
{
   if(is(object, "DPsurvint"))
   {
   	if(is.null(xnew))
   	{
   	     base<-1	
   	     npred<-1
   	     xnew<-rep(0,object$p)
   	     ngrid<-length(grid)
   	     nrec<-object$nrec
   	     
   	     mulog<-object$save.state$thetasave[,(object$p+1)]
	     sdlog<-sqrt(object$save.state$thetasave[,(object$p+2)])
	     alpha<-sqrt(object$save.state$thetasave[,(object$p+4)])
	      
	     v<-matrix(object$save.state$randsave,nsave,nrec+1)
	     vpred<-v[,(nrec+1)]
	     v<-v[,1:nrec]
	     covn<-rep(0,1)
	     covn[1]<-"Baseline"
   	}
   	else
   	{
   	     base<-0
   	     npred<-dim(xnew)[1]
   	     pnew<-dim(xnew)[2]
   	     ngrid<-length(grid)
   	     nrec<-object$nrec
   	     
   	     mulog<-object$save.state$thetasave[,(object$p+1)]
   	     sdlog<-sqrt(object$save.state$thetasave[,(object$p+2)])
   	     alpha<-sqrt(object$save.state$thetasave[,(object$p+4)])
   	     
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
	 }    

 	 pm<-matrix(0,npred,ngrid)
	 pmed<-matrix(0,npred,ngrid)
	 psd<-matrix(0,npred,ngrid)
	 pstd<-matrix(0,npred,ngrid)
	 plinf<-matrix(0,npred,ngrid)
	 plsup<-matrix(0,npred,ngrid)
	     
	 expxb<-exp(xnew%*%t(object$save.state$thetasave[,1:object$p]))
	        
	 for(i in 1:npred)
         {
	     surv<-matrix(0,nrow=ngrid,ncol=nsave)
	     for(j in 1:ngrid)
	     {
	         for(k in 1:nsave)
	         {
			linf<-grid[j]*expxb[i,k]
			mwork<-mulog[k]
			swork<-sdlog[k]
			awork<-alpha[k]
			vwork<-v[k,]
			Cdelta<- sum(vwork>linf)
			Cparam<- awork*plnorm(linf,meanlog=mwork,sdlog=swork,
			                      lower.tail=FALSE,log.p=FALSE)
			surv[j,k]<-(Cparam+Cdelta)/(awork+nrec)	                
	          }
	          pm[i,j]<-mean(surv[j,])
		  pmed[i,j]<-median(surv[j,])
		  psd[i,j]<-sqrt(var(surv[j,]))
		  pstd[i,j]<-sqrt(var(surv[j,]))/sqrt(nsave)
		  
                  if(hpd==TRUE)
                  {
                        alow<-rep(0,2)
                        aupp<-rep(0,2)
                        sig<-0.05
                        vec<-surv[j,]
                        n<-length(vec)
                        a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(sig),x=as.double(vec),
                                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                        plinf[i,j]<-a$alow[1]            
                        plsup[i,j]<-a$aupp[1]            
                  }
                  else
                  {
		        plinf[i,j]<-quantile(surv[j,],0.025)		
		        plsup[i,j]<-quantile(surv[j,],0.975)			                
                  }
	      }
   	 }
   	 
   	 dimnames(pm)<-list(covn,grid)
   	 dimnames(pmed)<-list(covn,grid)
   	 dimnames(psd)<-list(covn,grid)
   	 dimnames(pstd)<-list(covn,grid)
   	 dimnames(plinf)<-list(covn,grid)
   	 dimnames(plsup)<-list(covn,grid)
   	 
   	 out<-NULL
   	 out$pmean<-pm
   	 out$pmedian<-pmed
   	 out$psd<-psd
   	 out$pstd<-pstd
   	 out$plinf<-plinf
   	 out$plsup<-plsup
   	 out$xnew<-base
   	 out$vpred<-vpred
   	 out$grid<-grid
   	 out$npred<-npred
   	 out$base<-base
   	 out$covn<-covn

	 class(out)<-c("DPsurvpred")
	 out
   }
}   
