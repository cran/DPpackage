###                    
### Fit a linear Dirichlet Process mixture of normal model.
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 24-03-2006.

"DPlmm"<-
function(fixed,random,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("DPlmm")


"DPlmm.default"<-
function(fixed,
         random,
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
         m <- mcall <- cl <- match.call()
         nm <- names(m)[-1]
         keep <- is.element(nm, c("data", "na.action"))
         for (i in nm[!keep]) m[[i]] <- NULL
         cfixed<-all.vars(fixed)
         crandom<-all.vars(random)
         alist <- lapply(as.list(c(crandom,cfixed)), as.name)
         osubject<-do.call('data.frame', alist[length(crandom)])
         namessub<-names(table(osubject))
         nosubject<-length(namessub)
         idsubject<-rep(0,dim(osubject)[1])
 
         for(i in 1:nosubject){
             idsubject[osubject==namessub[i]]<-i
         }
 
         allvars<-NULL
         if(length(cfixed) >1) allvars <-c(cfixed[2:length(cfixed)])
         if(length(crandom)>1) allvars <-c(allvars,crandom[1:(length(crandom)-1)],"idsubject")
         if(length(crandom)==1) allvars <-c(allvars,"idsubject")
         
         m$formula <- as.formula(paste(cfixed[1],"~", paste(allvars, collapse = "+")))
         environment(m$formula) <- environment(fixed)
         m$drop.unused.levels <- TRUE
         Mdata<-do.call('data.frame', alist)
         m$data<-data.frame(idsubject,Mdata)
         m[[1]] <- as.name("model.frame")
         mf <- eval.parent(m)
         
         #########################################################################################
         # data structure
         #########################################################################################
     	 resp<- model.response(mf,"numeric")
	 nrec<-length(resp)
	 mat<-model.matrix(m$formula,mf)
	 pp<-dim(mat)[2]
        
         #random effects
         idrec<-seq(1,nrec)
         namessub2<-names(table(mf$idsubject))
         nsubject<-length(namessub2)
         namesre<-rep(0,nsubject)
         newid<-rep(0,nrec)
         for(i in 1:nsubject){
             newid[mf$idsubject==namessub2[i]]<-i
             namesre[i]<-namessub[as.integer(namessub2[i])]    
         }
         mf$idsubject<-newid

         freqsub<-table(mf$idsubject)
         maxni<-max(freqsub)
         datastr<-matrix(0,nrow=nsubject,ncol=maxni+1)
         datastr[,1]<-freqsub
         for(i in 1:nsubject){
             for(j in 1:freqsub[i]){
                datastr[i,(j+1)]<-idrec[idsubject==i][j] 
             }
         }

         #########################################################################################
         # design matrix for random and fixed effects
         #########################################################################################
         q<-1
         z<-mat[,1]
         possr<-1
         if(length(crandom)>1){
            for(i in 1:(length(crandom)-1)){
               for(j in 1:(pp-1)){
                  if(crandom[i]==colnames(mat)[j]){
                     z<-cbind(z,mat[,j])
                     q<-q+1
                     possr<-c(possr,j)
                  }   
               }
            }
         }
         if(q==1)z<-matrix(z,ncol=1)
         colnames(z)<-colnames(mat)[possr]

         x<-NULL
         possf<-NULL
         nfixed<-0
         for(i in 1:(pp-1)){
            type<-0
            for(j in 1:q){
               if(i==possr[j])type<-1
            }
            if(type==0){
              x<-cbind(x,mat[,i])   
              nfixed<-nfixed+1
              possf<-c(possf,i)
            }  
         }
         
         #########################################################
         #### NOTE: Is there is no fixed effects in the model ####
         ####       nfixed=0 and p=1.                         ####
         #########################################################

         if(nfixed==0){
            p<-1
            x<-matrix(0,nrow=nrec,ncol=1)
         }

         if(nfixed>0){
            p<-nfixed
            colnames(x)<-colnames(mat)[possf]
         }
        
         xtx<-t(x)%*%x

         #########################################################################################
         # prior information
         #########################################################################################
  	 
  	 if(is.null(prior$a0))
  	 {
  	    a0b0<-c(-1,-1)
  	    alpha<-prior$alpha
  	    alphapr<-0
  	 }
         else
         {
            a0b0<-c(prior$a0,prior$b0)
  	    alpha<-rgamma(1,shape=prior$a0,scale=prior$b0)
  	    alphapr<-1
  	 }
  	 
  	 nu0<-prior$nu0
  	 tau1<-prior$tau1
  	 tau2<-prior$tau2
  	 tau<-c(tau1,tau2)
  	 
	 tinv<-prior$tinv
	 
	 if(nfixed==0){
	    prec<-matrix(0,nrow=1,ncol=1)
	    sb<-matrix(0,nrow=1,ncol=1)
	 }

	 if(nfixed>0){
            prec<-solve(prior$Sbeta0)
            sb<-prec%*%prior$beta0
         }

	 psiinv<-solve(prior$Sb)
	 smu<-psiinv%*%prior$mub
	 
	 
         #########################################################################################
         # mcmc specification
         #########################################################################################
         mcmcvec<-c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
         nsave<-mcmc$nsave


         #########################################################################################
         # output
         #########################################################################################
         nuniq<-(q*(q+1)/2)
         
         randsave<-matrix(0,nrow=nsave,ncol=q*(nsubject+1))
         thetasave<-matrix(0,nrow=nsave,ncol=q+nfixed+1+q+nuniq+2)
         
         cpo<-rep(0,nrec)

         #########################################################################################
         # parameters depending on status
         #########################################################################################
    	 if(status==TRUE)
	 {
	        if(nfixed==0){
	           beta<-matrix(0,nrow=1,ncol=1)
	           fit0<- glm.fit(z, resp, family= gaussian(link = "identity"))   
	           sigma2e<-mean(fit0$residuals**2)
	           b<-matrix(0,nrow=nsubject,ncol=q)
		   bclus<-matrix(0,nrow=nsubject,ncol=q)
		   
	           bclus[1,]<-coefficients(fit0)
	           for(i in 1:nsubject){
	               b[i,]<-coefficients(fit0)
	           }
	           mu<-coefficients(fit0)
                   q1 <- 1:q
	           Qr <- fit0$qr
	           sigma<- chol2inv(Qr$qr[q1, q1, drop = FALSE])
	           sigmainv<-solve(sigma)
	        }

	        if(nfixed>0){
	           fit0<- glm.fit(cbind(x,z), resp, family= gaussian(link = "identity"))   
	           sigma2e<-mean(fit0$residuals**2)
	           b<-matrix(0,nrow=nsubject,ncol=q)
		   bclus<-matrix(0,nrow=nsubject,ncol=q)
                   beta<-coefficients(fit0)[1:p]
		   bclus[1,]<-coefficients(fit0)[(p+1):(p+q)]

	           for(i in 1:nsubject){
	               b[i,]<-coefficients(fit0)[(p+1):(p+q)]
	           }
	           mu<-coefficients(fit0)[(p+1):(p+q)]
	           
	           q1 <- 1:(p+q)
		   Qr <- fit0$qr
		   sigma<- chol2inv(Qr$qr[q1, q1, drop = FALSE])[(p+1):(p+q),(p+1):(p+q)]
		   sigmainv<-solve(sigma)
	        }

                betar<-rep(0,q)
                ncluster<-1
                ss<-rep(1,nsubject)
	 }	
      	 if(status==FALSE)
	 {
	        alpha<-state$alpha
                b<-state$b 
                bclus<-state$bclus 
	        beta<-state$beta
	        betar<-state$betar
	        mu<-state$mu
	        ncluster<-state$ncluster
	        sigma<-state$sigma
	        sigmainv<-solve(sigma)
	        sigma2e<-state$sigma2e         
	        ss<-state$ss
	 }
         
         #########################################################################################
         # working space
         #########################################################################################
         ccluster<-rep(0,nsubject) 
         iflag<-rep(0,p) 
         iflag2<-rep(0,maxni)
         iflagb<-rep(0,q) 
         prob<-rep(0,nsubject+1)
         quadf<-matrix(0,nrow=q,ncol=q)
         res<-rep(0,nrec)
         seed1<-sample(1:29000,1)
         seed2<-sample(1:29000,1)
         seed3<-sample(1:29000,1)
         seed<-c(seed1,seed2,seed2)
         theta<-rep(0,q)
         work1<-matrix(0,nrow=p,ncol=p)
         work2<-matrix(0,nrow=p,ncol=p)
         workb1<-matrix(0,nrow=q,ncol=q)
         workb2<-matrix(0,nrow=q,ncol=q)
         workmh1<-rep(0,p*(p+1)/2) 
         workmh2<-rep(0,q*(q+1)/2) 
         workmh3<-rep(0,q*(q+1)/2) 
         workk1<-matrix(0,nrow=maxni,ncol=q) 
         workkv1<-rep(0,maxni) 
         workkm1<-matrix(0,nrow=maxni,ncol=maxni) 
         workkm2<-matrix(0,nrow=maxni,ncol=maxni) 
         workv1<-rep(0,p) 
         workv2<-rep(0,p) 
         workvb1<-rep(0,q) 
         workvb2<-rep(0,q) 
         xty<-rep(0,p) 
         ywork<-rep(0,maxni) 
         zty<-rep(0,q) 
         ztz<-matrix(0,nrow=q,ncol=q) 
         ztzinv<-matrix(0,nrow=q,ncol=q) 
 
         #########################################################################################
         # calling the fortran code
         #########################################################################################
         foo <- .Fortran("splme",
         	datastr    =as.integer(datastr),
 	 	maxni      =as.integer(maxni),         
 	 	nrec       =as.integer(nrec),
 	 	nsubject   =as.integer(nsubject),
 	 	nfixed     =as.integer(nfixed),
 	 	p          =as.integer(p),
 	 	q          =as.integer(q),
 	 	subject    =as.integer(newid),
 		x          =as.double(x),	 	
 		xtx        =as.double(xtx),	 	
 		y          =as.double(resp),
 		z          =as.double(z),	 
 		a0b0       =as.double(a0b0),
 		nu0        =as.integer(nu0),
 		prec       =as.double(prec),	 
 		psiinv     =as.double(psiinv),	  		
 		sb         =as.double(sb),	  		
 		smu        =as.double(smu),	  		
 		tau        =as.double(tau),	  		
 		tinv       =as.double(tinv),	  		 		
 		mcmc       =as.integer(mcmcvec),
 		nsave      =as.integer(nsave),
 		randsave   =as.double(randsave),
 		thetasave  =as.double(thetasave),
 		cpo        =as.double(cpo),
 		alpha      =as.double(alpha),		
 		b          =as.double(b),		
                bclus      =as.double(bclus),		
 		beta       =as.double(beta),
 		betar      =as.double(betar),
 		mu         =as.double(mu),
 		ncluster   =as.integer(ncluster),
 		sigma      =as.double(sigma),
 		sigma2e    =as.double(sigma2e),
 		ss         =as.integer(ss),
 		ccluster   =as.integer(ccluster),
 		iflag      =as.integer(iflag),
 		iflag2     =as.integer(iflag2),
 		iflagb     =as.integer(iflagb),
 		prob       =as.double(prob),
 		quadf      =as.double(quadf),
 		res        =as.double(res),
 		seed       =as.integer(seed),
 		sigmainv   =as.double(sigmainv),
 		theta      =as.double(theta),
 		work1      =as.double(work1),
 		work2      =as.double(work2),
 		workb1     =as.double(workb1),
 		workb2     =as.double(workb2),
                workmh1    =as.double(workmh1),
                workmh2    =as.double(workmh2),
                workmh3    =as.double(workmh3),
                workk1     =as.double(workk1),
 		workkv1    =as.double(workkv1),
 		workkm1    =as.double(workkm1),
 		workkm2    =as.double(workkm2),
 		workv1     =as.double(workv1),
 		workv2     =as.double(workv2),
 		workvb1    =as.double(workvb1),
 		workvb2    =as.double(workvb2),
 		xty        =as.double(xty),
 		ywork      =as.double(ywork),
 		zty        =as.double(zty), 		
 		ztz        =as.double(ztz), 		
 		ztzinv     =as.double(ztzinv), 		
		PACKAGE    ="DPpackage")	


         #########################################################################################
         # save state
         #########################################################################################

         dimen<-q+nfixed+1+q+nuniq+2
         thetasave<-matrix(foo$thetasave,nrow=nsave, ncol=dimen)
         randsave<-matrix(foo$randsave,nrow=nsave, ncol=q*(nsubject+1))
         cpo<-foo$cpo
 
 	 if(nfixed==0)pnames1<- c(colnames(mat)[possr])
 	 if(nfixed >0)pnames1<- c(colnames(mat)[possr],colnames(mat)[possf])
 	 
 	 pnames2<-"residual"
 	 pnames3<- paste("mu",colnames(mat)[possr],sep="-")
         pnames4<-NULL
         
         for(i in 1:q){
            for(j in i:q){
               if(i==j)aa<-paste("sigma",colnames(mat)[possr][i],sep="-")
               if(i!=j)aa<-paste("sigma",colnames(mat)[possr][i],colnames(mat)[possr][j],sep="-")
               pnames4<-c(pnames4,aa)            
            }
         }
         
         pnames5<- c("ncluster","alpha")

         colnames(thetasave)<-c(pnames1,pnames2,pnames3,pnames4,pnames5)
         
         qnames<-NULL
         for(i in 1:nsubject){
             for(j in 1:q){
                 idname<-paste("(Subject",namesre[i],sep="=")
                 idname<-paste(idname,")")
                 qnamestemp<-paste(colnames(mat)[possr][j],idname,sep=" ")
                 qnames<-c(qnames,qnamestemp)
             }
         }
         for(j in 1:q){
             qnamestemp<-paste("pred",colnames(mat)[possr][j],sep="-")
             qnames<-c(qnames,qnamestemp)
         }
         
         colnames(randsave)<-qnames
         
         
	 model.name<-"Bayesian semiparametric linear mixed effect model"		

         coeff<-rep(0,dimen)
         for(i in 1:dimen){
             coeff[i]<-mean(thetasave[,i])
         }
		
	 names(coeff)<-c(pnames1,pnames2,pnames3,pnames4,pnames5)


	 state <- list(alpha=foo$alpha,b=matrix(foo$b,nrow=nsubject,ncol=q),
	               bclus=matrix(foo$bclus,nrow=nsubject,ncol=q),
	               beta=foo$beta,betar=foo$betar,
	               mu=foo$mu,ncluster=foo$ncluster,
	               sigma=matrix(foo$sigma,nrow=q,ncol=q),
	               sigma2e=foo$sigma2e,ss=foo$ss)

	 save.state <- list(thetasave=thetasave,randsave=randsave)

	 z<-list(modelname=model.name,coefficients=coeff,call=cl,
                 prior=prior,mcmc=mcmc,state=state,save.state=save.state,nrec=foo$nrec,
                 nsubject=foo$nsubject,nfixed=foo$nfixed,nrandom=foo$q,
                 cpo=foo$cpo,alphapr=alphapr,namesre1=namesre,namesre2=colnames(mat)[possr],
                 z=z,x=x,mf=mf,dimen=dimen)
                 
         cat("\n\n")        

         class(z)<-c("DPlmm")
         return(z) 
}
