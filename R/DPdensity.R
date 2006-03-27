###                    
### Fit a linear Dirichlet Process mixture of normal model for
### density estimation
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 13-04-2006.



DPdensity<-function(y,prior,mcmc,state,status,method="neal",data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("DPdensity")

DPdensity.default<-function(y,prior,mcmc,state,status,method="neal",data,na.action=na.fail)
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
  	 
  	 
  	 if(is.null(prior$nu2))
  	 {
              psiinv1<-matrix(prior$psiinv1,nvar,nvar)
              psiinv2<-psiinv1
              psi1<-matrix(solve(psiinv1),nvar,nvar)
              nuvec<-c(prior$nu1,-1)
              psi1rand<-0
  	 }
  	 else
  	 {
  	      psiinv1<-matrix(var(y),nvar,nvar)
              psi1<-matrix(solve(psiinv1),nvar,nvar)
              psiinv2<-matrix(prior$psiinv2,nvar,nvar)
  	      nuvec<-c(prior$nu1,prior$nu2)
  	      psi1rand<-1
  	 }
  	 
  	 
  	 if(is.null(prior$m2) && is.null(prior$s2))
  	 {
  	      s2inv<-matrix(0,nrow=nvar,ncol=nvar) 
  	      s2invm2<-matrix(0,nrow=nvar,ncol=1)
  	      m1<-prior$m1
  	      m1rand<-0
  	 }
  	 else
  	 {
              s2inv<-solve(prior$s2)
              s2invm2<-s2inv%*%prior$m2
  	      m1<-rep(0,nvar)
              for(i in 1:nvar)
      	      {
       		  m1[i]<-mean(y[,i])+rnorm(1,0,100)
       	      }
              m1rand<-1
         }     

         
         
         if(is.null(prior$tau1) && is.null(prior$tau2)) 
         {
             tau<-c(-2,-2)
             k0<-prior$k0
             k0rand<-0
         }
         else
         {
             tau<-c(prior$tau1,prior$tau2)
             k0<-rgamma(1,shape=prior$tau1,scale=prior$tau2)
             k0rand<-1
         }



         #########################################################################################
         # mcmc specification
         #########################################################################################
         mcmcvec<-c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
         nsave<-mcmc$nsave


         #########################################################################################
         # output
         #########################################################################################
         
         thetasave<-matrix(0,nrow=nsave,ncol=nvar+nvar*(nvar+1)/2+2)
         randsave<-matrix(0,nrow=nsave,ncol=(nrec+2)*nvar+(nrec+1)*nvar*(nvar+1)/2)
         

         #########################################################################################
         # parameters depending on status
         #########################################################################################
         nuniq<-nvar*(nvar+1)/2
         
    	 if(status==TRUE)
	 {
          	muclus<-matrix(0,nrow=nrec+2,ncol=nvar)
          	sigmaclus<-matrix(0,nrow=nrec+2,ncol=nuniq)
          	for(i in 1:1)
          	{
          		counter<-0
          		for(j in 1:nvar)
          		{
				muclus[i,j]<-m1[j]
          			for(k in j:nvar)
          			{
          				counter<-counter+1
          				sigmaclus[i,counter]<-psiinv1[j,k]
          			}
          		}
          	}
                ncluster<-1

                ss<-rep(1,nrec)
   	 }
	 
      	 if(status==FALSE)
	 {
	        alpha<-state$alpha
                m1<-state$m1
                muclus<-state$muclus 
	        ncluster<-state$ncluster
	        psi1<-state$psi1
	        psiinv1<-solve(psi1)
	        k0<-state$k0
	        sigmaclus<-state$sigmaclus 
	        ss<-state$ss
	 }    


         #########################################################################################
         # working space
         #########################################################################################
         ccluster<-rep(0,nrec) 
         cpo<-rep(0,nrec) 
         iflag<-rep(0,nvar) 
         muwork<-rep(0,nvar) 
         prob<-rep(0,(nrec+2))
         s1<-matrix(0,nvar,nvar)
         seed1<-sample(1:29000,1)
         seed2<-sample(1:29000,1)
         seed3<-sample(1:29000,1)
         seed<-c(seed1,seed2,seed2)
         sigmawork<-matrix(0,nrow=nvar,ncol=nvar)
         sigworkinv<-matrix(0,nrow=nvar,ncol=nvar)
         theta<-rep(0,nvar)
         workm1<-matrix(0,nrow=nvar,ncol=nvar)
         workm2<-matrix(0,nrow=nvar,ncol=nvar)
         workm3<-matrix(0,nrow=nvar,ncol=nvar)
         workmh1<-rep(0,nvar*(nvar+1)/2) 
         workmh2<-rep(0,nvar*(nvar+1)/2) 
         workv1<-rep(0,nvar) 
         workv2<-rep(0,nvar) 
         workv3<-rep(0,nvar) 
	 ywork<-rep(0,nvar)


         #########################################################################################
         # calling the fortran code
         #########################################################################################
         
         if(method=="no-gaps")
         {
                foo <- .Fortran("spdeng",
  	 	nrec       =as.integer(nrec),
  	 	nvar       =as.integer(nvar),
  	 	y          =as.double(y),
  	 	a0b0       =as.double(a0b0),
  	 	k0         =as.double(k0),
  	 	nuvec      =as.integer(nuvec),
  	 	m1rand     =as.integer(m1rand),
  	 	s2inv      =as.double(s2inv),
  	 	s2invm2    =as.double(s2invm2),
  	 	psiinv2    =as.double(psiinv2),
  	 	tau        =as.double(tau),
 		mcmc       =as.integer(mcmcvec),
 		nsave      =as.integer(nsave),
 		cpo        =as.double(cpo),
 		randsave   =as.double(randsave),
 		thetasave  =as.double(thetasave),
 		alpha      =as.double(alpha),		
 		m1         =as.double(m1),		
                muclus     =as.double(muclus),		 		
 		ncluster   =as.integer(ncluster),
 		psi1       =as.double(psi1),
 		psiinv1    =as.double(psiinv1),
 		s1         =as.double(s1),
 		sigmaclus  =as.double(sigmaclus),
 		ss         =as.integer(ss),
 		ccluster   =as.integer(ccluster),
 		iflag      =as.integer(iflag),
 		muwork     =as.double(muwork),
 		prob       =as.double(prob),
 		seed       =as.integer(seed),
 		sigmawork  =as.double(sigmawork),
 		sigworkinv =as.double(sigworkinv),
 		theta      =as.double(theta),
 		workm1     =as.double(workm1),
 		workm2     =as.double(workm2),
 		workm3     =as.double(workm3),
 		workmh1    =as.double(workmh1),
 		workmh2    =as.double(workmh2),
 		workv1     =as.double(workv1),
 		workv2     =as.double(workv2),
 		workv3     =as.double(workv3),
		ywork      =as.double(ywork),
		PACKAGE    ="DPpackage")
         }		

         if(method=="neal")
         {
                foo <- .Fortran("spdenn",
  	 	nrec       =as.integer(nrec),
  	 	nvar       =as.integer(nvar),
  	 	y          =as.double(y),
  	 	a0b0       =as.double(a0b0),
  	 	k0         =as.double(k0),
  	 	nuvec      =as.integer(nuvec),
  	 	m1rand     =as.integer(m1rand),
  	 	s2inv      =as.double(s2inv),
  	 	s2invm2    =as.double(s2invm2),
  	 	psiinv2    =as.double(psiinv2),
  	 	tau        =as.double(tau),
 		mcmc       =as.integer(mcmcvec),
 		nsave      =as.integer(nsave),
 		cpo        =as.double(cpo),
 		randsave   =as.double(randsave),
 		thetasave  =as.double(thetasave),
 		alpha      =as.double(alpha),		
 		m1         =as.double(m1),		
                muclus     =as.double(muclus),		 		
 		ncluster   =as.integer(ncluster),
 		psi1       =as.double(psi1),
 		psiinv1    =as.double(psiinv1),
 		s1         =as.double(s1),
 		sigmaclus  =as.double(sigmaclus),
 		ss         =as.integer(ss),
 		ccluster   =as.integer(ccluster),
 		iflag      =as.integer(iflag),
 		muwork     =as.double(muwork),
 		prob       =as.double(prob),
 		seed       =as.integer(seed),
 		sigmawork  =as.double(sigmawork),
 		sigworkinv =as.double(sigworkinv),
 		theta      =as.double(theta),
 		workm1     =as.double(workm1),
 		workm2     =as.double(workm2),
 		workm3     =as.double(workm3),
 		workmh1    =as.double(workmh1),
 		workmh2    =as.double(workmh2),
 		workv1     =as.double(workv1),
 		workv2     =as.double(workv2),
 		workv3     =as.double(workv3),
		ywork      =as.double(ywork),
		PACKAGE    ="DPpackage")
         }		

         #########################################################################################
         # save state
         #########################################################################################

                model.name<-"Bayesian semiparametric density estimation"		
                
		varnames<-colnames(y)
		
                state <- list(
                       alpha=foo$alpha,
	               m1=matrix(foo$m1,nrow=nvar,ncol=1),
	               muclus=matrix(foo$muclus,nrow=nrec+2,ncol=nvar),
	               ncluster=foo$ncluster,
	               psi1=matrix(foo$psi1,nrow=nvar,ncol=nvar),
	               k0=foo$k0,
	               sigmaclus=matrix(foo$sigmaclus,nrow=nrec+2,ncol=nuniq),
                       ss=foo$ss
                             )
                randsave<-matrix(foo$randsave,nrow=nsave,ncol=(nrec+2)*nvar+(nrec+1)*nvar*(nvar+1)/2)
                thetasave<-matrix(foo$thetasave,nrow=nsave,ncol=nvar+nvar*(nvar+1)/2+2)
 
                indip<-rep(0,(nvar+nvar*(nvar+1)/2+2))
 
 
                if(is.null(varnames))
                {
                      varnames<-all.vars(cl)[1:nvar]
                }
 
 
                coeff<-NULL
                pnames1<-NULL

                if(is.null(prior$m2) && is.null(prior$s2))
		{
		
		}
		else
		{
		    for(i in 1:nvar)
		    {
		       coeff<-c(coeff,mean(thetasave[,i]))
		       pnames1<-c(pnames1,paste("m1",varnames[i],sep="-"))
		       indip[i]<-1
		    }
		}

           	if(is.null(prior$nu2))
        	{
        	
                }
                else
                {
		    for(i in 1:(nvar*(nvar+1)/2))
		    {
		       coeff<-c(coeff,mean(thetasave[,(nvar+i)]))
		       indip[(nvar+i)]<-1
		    }
		    
                    for(i in 1:nvar)
                    {
                       for(j in i:nvar)
                       {
                         if(i==j)pnames1<-c(pnames1,paste("psi1",varnames[i],sep="-"))
                         if(i!=j)
                         {
                            tempname<-paste(varnames[i],varnames[j],sep="-")
                            pnames1<-c(pnames1,paste("psi1",tempname,sep="-"))
                         }    
                       }  
                    }
                }

                coeff<-c(coeff,mean(thetasave[,(nvar+nvar*(nvar+1)/2+1)]))
                pnames2<-c("ncluster")
		indip[(nvar+nvar*(nvar+1)/2+1)]<-1

                if(alpharand==1)
                {
                   coeff<-c(coeff,mean(thetasave[,(nvar+nvar*(nvar+1)/2+2)]))
                   pnames2<-c(pnames2,"alpha")
		   indip[(nvar+nvar*(nvar+1)/2+2)]<-1
                }


                names(coeff)<-c(pnames1,pnames2)


                pnames1<-NULL
	        for(i in 1:nvar)
	        {
            	       pnames1<-c(pnames1,paste("m1",varnames[i],sep="-"))
                }

                for(i in 1:nvar)
                {
                    for(j in i:nvar)
                    {
                        if(i==j)pnames1<-c(pnames1,paste("psi1",varnames[i],sep="-"))
                        if(i!=j)
                        {
                            tempname<-paste(varnames[i],varnames[j],sep="-")
                            pnames1<-c(pnames1,paste("psi1",tempname,sep="-"))
                        }    
                    }  
                }

                pnames2<-c("ncluster","alpha")
                dimnames(thetasave)<-list(NULL,c(pnames1,pnames2))

                pnamesre<-NULL              
                for(i in 1:nrec)
                {
                    for(j in 1:nvar)
                    {
                          tmpn<-paste("mu",varnames[j],sep="-")
                          tmpn<-paste(tmpn," (Subject=",sep="")
                          tmpn<-paste(tmpn,i,sep="")
                          tmpn<-paste(tmpn,")",sep="")
                          pnamesre<-c(pnamesre,tmpn)
                    }
                    
                    for(j in 1:nvar)
                    {
                    	  for(k in j:nvar)
                    	  {
                               tmpn<-paste("sigma",varnames[j],sep="-")
                               tmpn<-paste(tmpn,varnames[k],sep="-")
                               tmpn<-paste(tmpn," (Subject=",sep="")
                               tmpn<-paste(tmpn,i,sep="")
                               tmpn<-paste(tmpn,")",sep="")
                               pnamesre<-c(pnamesre,tmpn)
                           }    
                    }
                }


                for(j in 1:nvar)
                {
                      tmpn<-paste("mu",varnames[j],sep="-")
                      tmpn<-paste(tmpn," (Prediction)",sep="")
                      tmpn<-paste(tmpn,")",sep="")
                      pnamesre<-c(pnamesre,tmpn)
                }
                
                for(j in 1:nvar)
                {
                	  for(k in j:nvar)
                	  {
                           tmpn<-paste("sigma",varnames[j],sep="-")
                           tmpn<-paste(tmpn,varnames[k],sep="-")
                           tmpn<-paste(tmpn," (Prediction)",sep="")
                           tmpn<-paste(tmpn,")",sep="")
                           pnamesre<-c(pnamesre,tmpn)
                       }    
                }


                for(j in 1:nvar)
                {
                      tmpn<-paste(varnames[j]," (Prediction)",sep="")
                      pnamesre<-c(pnamesre,tmpn)
                }

                dimnames(randsave)<-list(NULL,pnamesre)
              
                save.state <- list(thetasave=thetasave,randsave=randsave)


	 z<-list(call=cl,y=y,varnames=varnames,modelname=model.name,cpo=foo$cpo,
                 prior=prior,mcmc=mcmc,state=state,save.state=save.state,nrec=foo$nrec,
                 nvar=foo$nvar,alpharand=alpharand,psi1rand=psi1rand,m1rand=m1rand,
                 k0rand=k0rand,coefficients=coeff,indip=indip)
                 
         cat("\n\n")
 	 class(z)<-"DPdensity"
  	 return(z)
}


