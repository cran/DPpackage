###                    
### Fit a semiparametric logistic regression model.
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 09-05-2006.

"DPbinary"<-
function(formula,prior,mcmc,state,status,misc=NULL,data=sys.frame(sys.parent()),na.action=na.fail) 
UseMethod("DPbinary")

"DPbinary.default"<-
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
         # mcmc specification
         #########################################################################################

         mcmcvec<-c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
         nsave<-mcmc$nsave

	 fit0<- glm(yobs~x[,2:p],family=binomial(logit))
	 propv<- vcov(fit0)
         if(is.null(mcmc$tune))
         {
            tune=1
         }
         else
         {
            tune<-mcmc$tune
         }
         propv<-propv*tune

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
	

         #########################################################################################
         # parameters depending on status
         #########################################################################################

       	 if(status)
	 {
		beta<-coefficients(fit0)
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
	 maxint<-2*ntotal+1+nlink
	 maxend<-2*ntotal+nlink
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
	 mass<-rep(0,maxint)
	 ncluster<-nrec
	 prob<-rep(0,maxint)
	 prob2<-rep(0,maxint)
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
	 seed3<-sample(1:29000,1)
	 cpo<-rep(0,nrec)


         #########################################################################################
         # calling the fortran code
         #########################################################################################

  	 foo <- .Fortran("dpbinaryl",
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
		mcmc      =as.integer(mcmcvec),
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
		maxint    =as.integer(maxint),
		maxend    =as.integer(maxend),
		maxm      =as.integer(maxm),
		mass      =as.double(mass),
		prob      =as.double(prob),
		prob2     =as.double(prob2),
		seed1     =as.integer(seed1),
		seed2     =as.integer(seed2),
		seed3     =as.integer(seed3),
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
	
	 fsave<-matrix(foo$fsave,nrow=mcmc$nsave, ncol=nlink)
 	 thetasave<-matrix(foo$thetasave,nrow=mcmc$nsave, ncol=(p+2))
 	 randsave<-matrix(foo$randsave,nrow=mcmc$nsave, ncol=(nrec+1))

         model.name<-"Bayesian semiparametric binary regression model"		

  	 colnames(thetasave)<-c(dimnames(x)[[2]],"ncluster","alpha")
  	 
  	 coeff<-rep(0,(p+2))
	 for(i in 1:(p+2)){
	     coeff[i]<-mean(thetasave[,i])
	 }  
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
	         cpo=foo$cpo,p=p,nlink=nlink,xlink=xlink)
	         
	 cat("\n\n")

	 class(z)<-c("DPbinary")
	 return(z) 
}




	
	        
