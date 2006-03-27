###                    
### Fit a semiparametric aft model for interval censored data.
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 03-05-2006.

"DPsurvint"<-
function(formula,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("DPsurvint")

"DPsurvint.default"<-
function(formula,
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
  	 nrec<-length(y[,1])
  	 x<-model.matrix(formula)
  	 p<-dim(x)[2]
  	 x<-x[,2:p]
  	 p<-(p-1)

         type<-rep(2,nrec)
         for(i in 1:nrec){
	    type[y[,1]==-999]<-1
	    type[y[,2]==-999]<-3
         }

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
	 m0<-prior$m0
	 s0<-prior$s0
	 tau<-c(prior$tau1,prior$tau2)
	
         #########################################################################################
         # mcmc specification
         #########################################################################################

         mcmcvec<-c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
         nsave<-mcmc$nsave
         
  	 tmpy<- y[type==2,]
         mid<-matrix((tmpy[,1]+tmpy[,2])/2,ncol=1)
         nmid<-dim(mid)[1]

         tmpx<- x[type==2,]
         tmpx<- cbind(rep(1,nmid),tmpx)
         beta<-solve(t(tmpx)%*%tmpx)%*%t(tmpx)%*%log(mid)

         s2<-sum((log(mid)-tmpx%*%beta)**2)/(nmid-p-1)
	 propv<-solve(t(tmpx)%*%tmpx)*s2
	 propv<-propv*mcmc$tune
	 
	 beta<--beta[2:(p+1)]
	 propv<-propv[2:(p+1),2:(p+1)]

         #########################################################################################
         # output
         #########################################################################################

	 thetasave <- matrix(0, nrow=mcmc$nsave, ncol=p+4)
	 randsave  <- matrix(0, nrow=mcmc$nsave, ncol=nrec+1)
	 
	
         #########################################################################################
         # parameters depending on status
         #########################################################################################

	 if(status==TRUE)
	 {
	 	beta<-beta
	 	eta <- x %*% beta
		v <- rep(0,nrec)
		v[type==1]<- y[type==1,2]*exp(eta[type==1])/2
		v[type==2]<-(y[type==2,1]*exp(eta[type==2])+y[type==2,2]*exp(eta[type==2]))/2
		v[type==3]<-(y[type==3,1]*exp(eta[type==3])+max(y)*exp(eta[type==3]))/2
		mu<-0
		sigma<-s2
	 }	
	 if(status==FALSE)
	 {
	        alpha<-state$alpha 
		beta<-state$beta
		v<-state$v
		eta <- x %*% beta
		mu<-state$mu
		sigma<-state$sigma
	 }

         #########################################################################################
         # working space
         #########################################################################################

 	 acrate<-0
	 cpo<-rep(0,nrec)
	 ncluster<-0
	 ntotal<-nrec+1
	 maxint<-4*nrec+1
	 maxend<-4*nrec
	 maxm<-log(0.0001)/log((ntotal+1)/(ntotal+2))
 	 iflag<-rep(0,p)
	 intcount<-rep(0,maxint)
	 intcount2<-rep(0,maxint)
	 intind<-rep(0,(nrec+1))
	 intind2<-rep(0,(nrec+1))
	 seed<-c(sample(1:29000,1),sample(1:29000,1),sample(1:29000,1))
	 betac<-rep(0,p)
         etan<-rep(0,nrec)
	 endp<-rep(0,maxend)
	 endp2<-rep(0,maxend)
	 mass<-rep(0,maxint)
	 prob<-rep(0,maxint)
	 prob2<-rep(0,maxint)
	 uvec<-rep(0,maxm)
	 vvec<-rep(0,maxm)
	 vnew<-rep(0,(nrec+1))
	 vnew2<-rep(0,nrec)
         wvec<-rep(0,maxm)
	 workm1<-matrix(0,nrow=p,ncol=p)
	 workm2<-matrix(0,nrow=p,ncol=p)
	 workmh1<-rep(0,(p*(p+1)/2))
	 workv1<-rep(0,p)
	 workv2<-rep(0,p)
	

         #########################################################################################
         # calling the fortran code
         #########################################################################################
 
 	 foo <- .Fortran("dpsurvint",
	 	nrec      =as.integer(nrec),
	 	p         =as.integer(p),
		x         =as.double(x),	 	
		y         =as.double(y),
		interind  =as.integer(type),
		mcmc      =as.integer(mcmcvec),
		nsave     =as.integer(mcmc$nsave),
		propv     =as.double(propv),
		a0b0      =as.integer(a0b0),
		betapm    =as.double(betapm),		
		betapv    =as.double(betapv),		
		m0        =as.double(m0),		
		s0        =as.double(s0),		
		tau       =as.double(tau),
		acrate    =as.double(acrate),
		thetasave =as.double(thetasave),
		randsave  =as.double(randsave),
                cpo       =as.double(cpo),
                ncluster  =as.integer(ncluster),
		alpha     =as.double(alpha),		
		beta      =as.double(beta),
		mu        =as.double(mu),
		sigma2    =as.double(sigma),
                v         =as.double(v),
 		maxint    =as.integer(maxint),
		maxend    =as.integer(maxend),
		maxm      =as.integer(maxm),
		iflag     =as.integer(iflag),
		intcount  =as.integer(intcount),
		intcount2 =as.integer(intcount2),
		intind    =as.integer(intind),
		intind2   =as.integer(intind2),
                seed      =as.integer(seed),
		betac     =as.double(betac),
		eta       =as.double(eta),
		etan      =as.double(etan),
		endp      =as.double(endp),
		endp2     =as.double(endp2),
		mass      =as.double(mass),
		prob      =as.double(prob),
		prob2     =as.double(prob2),		
                uvec      =as.double(uvec),
                vvec      =as.double(vvec),
		vnew      =as.double(vnew),
		vnew2     =as.double(vnew2),
		wvec      =as.double(wvec),
		workm1    =as.double(workm1),
		workm2    =as.double(workm2),
		workmh1   =as.double(workmh1),
		workv1    =as.double(workv1),
		workv2    =as.double(workv2),
		PACKAGE="DPpackage")	

         #########################################################################################
         # save state
         #########################################################################################
	
 	 thetasave<-matrix(foo$thetasave,nrow=mcmc$nsave, ncol=(p+4))
 	 randsave<-matrix(foo$randsave,nrow=mcmc$nsave, ncol=(nrec+1))
 	 
	 colnames(thetasave)<-c(dimnames(x)[[2]],"mu","sigma2","ncluster","alpha")

         qnames<-NULL
         for(i in 1:nrec){
             idname<-paste("(Subject",i,sep="=")
             idname<-paste(idname,")",sep="")
             qnames<-c(qnames,idname)
         }
         qnames<-c(qnames,"Prediction")
         
         colnames(randsave)<-qnames

	 model.name<-"Bayesian Semiparametric AFT Regression Model"
	 coeff<-rep(0,(p+4))
	 for(i in 1:(p+4)){
	 	coeff[i]<-mean(thetasave[,i])
	 }
	 names(coeff)<-c(dimnames(x)[[2]],"mu","sigma2","ncluster","alpha")
	
	 state <- list(alpha=foo$alpha,beta=foo$beta,v=foo$v,mu=foo$mu,sigma=foo$sigma)	
		      
	 save.state <- list(thetasave=thetasave,randsave=randsave)

	 z<-list(modelname=model.name,coefficients=coeff,acrate=foo$acrate,call=cl,
	         prior=prior,mcmc=mcmc,state=state,save.state=save.state,cpo=foo$cpo,
	         nrec=nrec,p=p)
	
	 cat("\n\n")
	 class(z)<-c("DPsurvint")
	 z 
}


