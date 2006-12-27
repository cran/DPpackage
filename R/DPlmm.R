### DPlmm.R                   
### Fit a linear Dirichlet Process mixture of normal model.
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 24-03-2006.
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
         seed<-c(seed1,seed2)
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




"print.DPlmm"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", sep = "")
    print(x$call)
    cat("\n")

    if (length(x$coefficients)) {
        cat("Posterior Inference of Parameters:\n")
        if(x$alphapr==1){
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)}
        if(x$alphapr==0){
        print.default(format(x$coefficients[1:(length(x$coefficients)-1)], digits = digits), print.gap = 2, 
            quote = FALSE)}

    }
    else cat("No coefficients\n")
    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Groups:",x$nsubject,"\n")    
    cat("\n\n")
    invisible(x)
}


"summary.DPlmm"<-function(object, hpd=TRUE, ...) 
{
    dimen1<-object$nrandom+object$nfixed+1
    coef.p<-object$coefficients[1:(dimen1-1)]
    
    coef.sd<-rep(0,dimen1-1)
    coef.se<-rep(0,dimen1-1)
    coef.l<-rep(0,dimen1-1)
    coef.u<-rep(0,dimen1-1)
    coef.m<-rep(0,dimen1-1)
    names(coef.sd)<-names(object$coefficients[1:(dimen1-1)])
    names(coef.l)<-names(object$coefficients[1:(dimen1-1)])
    names(coef.u)<-names(object$coefficients[1:(dimen1-1)])
    
    alpha<-0.05
    
    for(i in 1:(dimen1-1)){
        alow<-rep(0,2)
        aupp<-rep(0,2)
        coef.sd[i]<-sqrt(var(object$save.state$thetasave[,i]))
        coef.m[i]<-median(object$save.state$thetasave[,i])
        vec<-object$save.state$thetasave[,i]
        n<-length(vec)
        
        a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                     
        if(hpd){             
           coef.l[i]<-a$alow[1]            
           coef.u[i]<-a$aupp[1]
        }   
        else
        {
           coef.l[i]<-a$alow[2]            
           coef.u[i]<-a$aupp[2]
        }
    }

    coef.se<-coef.sd/sqrt(n)

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


    coef.p<-object$coefficients[dimen1]
    coef.sd<-rep(0,1)
    coef.se<-rep(0,1)
    coef.l<-rep(0,1)
    coef.u<-rep(0,1)
    coef.m<-rep(0,1)
    names(coef.sd)<-names(object$coefficients[dimen1])
    names(coef.l)<-names(object$coefficients[dimen1])
    names(coef.u)<-names(object$coefficients[dimen1])
    for(i in 1:1){
        alow<-rep(0,2)
        aupp<-rep(0,2)
        coef.sd[i]<-sqrt(var(object$save.state$thetasave[,dimen1]))
        coef.m[i]<-median(object$save.state$thetasave[,dimen1])
        vec<-object$save.state$thetasave[,dimen1]
        n<-length(vec)
        
        a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                     
        if(hpd)
        {
            coef.l[i]<-a$alow[1]            
            coef.u[i]<-a$aupp[1]            
        }
        else
        {
            coef.l[i]<-a$alow[2]            
            coef.u[i]<-a$aupp[2]            
        }
    }
    coef.se<-coef.sd/sqrt(n)
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
    ans$resvar<-coef.table


    ans$cpo<-object$cpo


    dimen2<-object$nrandom+object$nrandom*(object$nrandom+1)/2
    coef.p<-object$coefficients[(dimen1+1):(dimen1+dimen2)]
    coef.sd<-rep(0,dimen2)
    coef.se<-rep(0,dimen2)
    coef.l<-rep(0,dimen2)
    coef.u<-rep(0,dimen2)
    coef.m<-rep(0,dimen2)
    names(coef.sd)<-names(object$coefficients[(dimen1+1):(dimen1+dimen2)])
    names(coef.l)<-names(object$coefficients[(dimen1+1):(dimen1+dimen2)])
    names(coef.u)<-names(object$coefficients[(dimen1+1):(dimen1+dimen2)])
    for(i in 1:dimen2){
        alow<-rep(0,2)
        aupp<-rep(0,2)
        coef.sd[i]<-sqrt(var(object$save.state$thetasave[,dimen1+i]))
        coef.m[i]<-median(object$save.state$thetasave[,dimen1+i])
        vec<-object$save.state$thetasave[,dimen1+i]
        n<-length(vec)
        
        a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
        
        if(hpd)
        {
            coef.l[i]<-a$alow[1]            
            coef.u[i]<-a$aupp[1]            
        }
        else
        {
            coef.l[i]<-a$alow[2]            
            coef.u[i]<-a$aupp[2]            
        }
    }

    coef.se<-coef.sd/sqrt(n)

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

    if(is.null(object$prior$a0))
    {
      dimen3<-1
      coef.p<-object$coefficients[(dimen1+dimen2+1)]
    }
    else
    {
      dimen3<-length(object$coefficients)-(dimen1+dimen2)
      coef.p<-object$coefficients[(dimen1+dimen2+1):length(object$coefficients)]
    }  

    coef.sd<-rep(0,dimen3)
    coef.se<-rep(0,dimen3)
    coef.l<-rep(0,dimen3)
    coef.u<-rep(0,dimen3)
    coef.m<-rep(0,dimen3)
    names(coef.sd)<-names(coef.p)
    names(coef.l)<-names(coef.p)
    names(coef.u)<-names(coef.p)
    for(i in 1:dimen3){
        alow<-rep(0,2)
        aupp<-rep(0,2)
        coef.sd[i]<-sqrt(var(object$save.state$thetasave[,dimen1+dimen2+i]))
        coef.m[i]<-median(object$save.state$thetasave[,dimen1+dimen2+i])
        vec<-object$save.state$thetasave[,dimen1+dimen2+i]
        n<-length(vec)
        
        a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
        if(hpd)
        {
            coef.l[i]<-a$alow[1]            
            coef.u[i]<-a$aupp[1]            
        }
        else
        {
            coef.l[i]<-a$alow[2]            
            coef.u[i]<-a$aupp[2]            
        }
    }

    coef.se<-coef.sd/sqrt(n)

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

    class(ans) <- "summaryDPlmm"
    return(ans)
}


"print.summaryDPlmm"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
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

    cat("\nResidual variance:\n")
    print.default(format(x$resvar, digits = digits), print.gap = 2, 
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
    else cat("No precision parameter\n")

    cat("\nNumber of Observations:",x$nrec)
    cat("\nNumber of Groups:",x$nsubject,"\n")
    cat("\n\n")
    invisible(x)
}


"plot.DPlmm"<-function(x, hpd=TRUE, ask=TRUE, nfigr=2, nfigc=2, param=NULL, col="#bdfcc9", ...)
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


   if(is(x, "DPlmm"))
   {
        if(is.null(param))
        {
           coef.p<-x$coefficients
           n<-length(coef.p)
           pnames<-names(coef.p)
           
           par(ask = ask)
           layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr , ncol=nfigc ,byrow=TRUE))
           for(i in 1:(n-1))
           {
               title1<-paste("Trace of",pnames[i],sep=" ")
               title2<-paste("Density of",pnames[i],sep=" ")       
               plot(x$save.state$thetasave[,i],type='l',main=title1,xlab="MCMC scan",ylab=" ")
               if(pnames[i]=="ncluster")
               {
                  hist(x$save.state$thetasave[,i],main=title2,xlab="values", ylab="probability",probability=TRUE)
               }
               else
               {
                  fancydensplot1(x$save.state$thetasave[,i],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
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
               fancydensplot1(x$save.state$thetasave[,n],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
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
            if (poss==0) 
	    {
	      stop("This parameter is not present in the original model.\n")
	    }
	    
	    par(ask = ask)
	    layout(matrix(seq(1,nfigr*nfigc,1), nrow=nfigr, ncol=nfigc, byrow = TRUE))
            title1<-paste("Trace of",pnames[poss],sep=" ")
            title2<-paste("Density of",pnames[poss],sep=" ")       
            plot(x$save.state$thetasave[,poss],type='l',main=title1,xlab="MCMC scan",ylab=" ")
            if(param=="ncluster")
            {
               hist(x$save.state$thetasave[,poss],main=title2,xlab="values", ylab="probability",probability=TRUE)
            }
            else
            {
               fancydensplot1(x$save.state$thetasave[,poss],hpd=hpd,main=title2,xlab="values", ylab="density",col=col)
            }   
        }
   }

}

