    
c=======================================================================                      
      subroutine dpsurvint(nrec,p,x,y,interind,
     &                     mcmc,nsave,propv,
     &                     a0b0,betapm,betapv,m0,s0,tau,
     &                     acrate,thetasave,randsave,cpo,
     &                     ncluster,alpha,beta,mu,sigma2,v,
     &                     maxint,maxend,maxm,
     &                     iflag,intcount,intcount2,    
     &                     intind,intind2,seed,
     &                     betac,eta,etan,endp,endp2,mass,
     &                     prob,prob2,
     &                     uvec,vvec,vnew,vnew2,wvec,
     &                     workm1,workm2,
     &                     workmh1,workv1,workv2)      
c=======================================================================                      
c
c     Version 1.0: 
c     Last modification: 02-05-2006.
c
c     Subroutine `dpsurvint' to run a Markov chain in the  
c     semiparametric atf model for interval censored data. 
c
c     Author: A.J.V.
c
c---- Data -------------------------------------------------------------
c
c        interind    :  integer vector giving the type of interval for
c                       each observation, 1=left, 2=inter, and 3=right,
c                       interind(nrec)
c        nrec        :  integer giving the number of observations.
c        p           :  integer giving the number of fixed coefficients.
c        x           :  real matrix giving the design matrix for the 
c                       fixed effects, x(nrec,p). 
c        y           :  real matrix giving the limits of the intervals,
c                       for the unobserved times, y(nrec,2).
c-----------------------------------------------------------------------
c
c---- Prior information ------------------------------------------------
c 
c        aa0, ab0    :  real giving the hyperparameters of the prior
c                       distribution for the precision parameter,
c                       alpha ~ Gamma(aa0,ab0). If aa0<0 the precision 
c                       parameter is considered as a constant.
c        betapm      :  real vector giving the prior mean of regression
c                       coefficients, betapm(p).
c        betapv      :  real matrix giving the prior covariance of 
c                       regression coefficients, betapv(p,p).
c        m0          :  real giving the prior mean of the mean of the
c                       lognormal baseline distribution.
c        s0          :  real giving the prior variance of the mean of 
c                       the lognormal baseline distribution.
c        tau1, tau2  :  reals giving the hyperparameters of the prior 
c                       distribution for the inverse of the 
c                       variance of the lognormal baseline distribution,
c                       1/sigma ~ Gamma(tau1/2,tau2/2).
c
c-----------------------------------------------------------------------
c
c---- MCMC parameters --------------------------------------------------
c
c        nburn       :  integer giving the number of burn-in scans.
c        ndisplay    :  integer giving the number of saved scans to be
c                       displayed on screen.
c        nskip       :  integer giving the thinning interval.
c        nsave       :  integer giving the number of scans to be saved.
c        propv       :  real matrix giving the variance of the normal
c                       proposal for the mh algorithm, propv(p,p).
c        
c-----------------------------------------------------------------------
c
c---- Output ----------------------------------------------------------- 
c
c        acrate      :  real giving the MH acceptance rate. 
c        randsave    :  real matrix containing the mcmc samples for
c                       the errors and prediction,
c                       randsave(nsave,nrec+1).
c        thetasave   :  real vector containing the mcmc sample for the
c                       regression parameters, betsave(nsave,p+4). 
c        cpo         :  real giving the cpo, cpo(nrec).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        alpha       :  real giving the current value of the precision
c                       parameter of the Dirichlet process.
c        beta        :  real vector giving the current value of the 
c                       regression coefficients, beta(p).
c        mu          :  real giving the mean of the baseline 
c                       distribution.
c        ncluster    :  integer giving the number of clusters in the
c                       random effects.
c        sigma2      :  real giving the current value of the variance
c                       of the baseline distribution.
c        v           :  real vector giving the current value of the 
c                       errors, v(nrec).
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        area        :  real giving the area of a region.
c        ainf        :  logical working variable.
c        binf        :  logical working variable.
c        betac       :  real vector giving the current value of the 
c                       candidate for regression parameters, betac(p).
c        cdflnorm    :  cdf of a lognormal distribution.
c        counter     :  index.
c        dispcount   :  index. 
c        endp        :  real vector giving the end of the intervals,
c                       endp(maxend).
c        endp2       :  real vector giving the end of the intervals,
c                       endp2(maxend).
c        eta         :  real vector giving the linear predictor, 
c                       eta(nrec).
c        etan        :  real vector giving the linear predictor, 
c                       etan(nrec).
c        evali       :  integer indicator used in updating the state.
c        i           :  index.
c        indi        :  index.
c        intcount    :  integer vector used to count the number of 
c                       observations for each interval, 
c                       intcount(maxint).
c        intcount2   :  integer vector used to count the number of 
c                       observations for each interval, 
c                       intcount2(maxint).
c        intind      :  integer vector giving the interval where 
c                       each obbservation belong, intind(nrec+1).
c        intind2     :  integer vector giving the interval where 
c                       each obbservation belong, intind2(nrec+1).
c        iflag       :  integer vector used to evaluate the prior
c                       distribution for the regression coefficients, 
c                       iflag(p).
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        j1          :  index.         
c        j2          :  index.         
c        k           :  index.   
c        keepbeta    :  real working variable.
c        liminf      :  real working variable.
c        limsup      :  real working variable.
c        logliko     :  real working variable.
c        loglikn     :  real working variable.
c        logprioro   :  real working variable.
c        logpriorn   :  real working variable.
c        maxend      :  integer giving the maximum number of intervals.
c        maxm        :  integer giving the maximum number of points in
c                       the finite approximation of the DP.
c        mass        :  real vector giving the mass of the DP,
c                       mass(maxint).
c        maxint      :  integer giving the maximum number of intervals.
c        maxu        :  real working variable.
c        mrand       :  index.
c        npoints     :  index.
c        ns          :  index.
c        nso         :  index.
c        nscan       :  index.
c        ok          :  integer indicator.
c        prob        :  real vector used to update the cluster 
c                       structure, prob(maxint).
c        prob2       :  real vector used to update the cluster 
c                       structure, prob2(maxint).
c        ratio       :  real working variable.
c        rbeta       :  real beta random number generator.
c        rgamma      :  real gamma random number generator.
c        rnorm       :  real normal random number generator.
c        runif       :  real uniform random number generator.
c        rtlnorm     :  real truncated lnormal random number generator.
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        seed3       :  seed for random number generation.
c        sel         :  index.
c        skipcount   :  index. 
c        swork       :  real working variable.
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
c        uvec        :  real working vector, uvec(maxm).
c        vbar        :  real working variable.
c        vvec        :  real working vector, vvec(maxm).
c        vnew        :  real working vector, vnew(nrec+1).
c        vnew2       :  real working vector, vnew2(nrec).
c        workm1      :  real matrix used to update the fixed effects,
c                       workm1(p,p).
c        workm2      :  real matrix used to update the fixed effects,
c                       workm2(p,p).
c        workmh1     :  real vector used to update the fixed effects,
c                       workmh1(p*(p+1)/2).
c        workv1      :  real vector used to update the fixed effects,
c                       workv1(p).
c        workv2      :  real vector used to update the fixed effects,
c                       workv2(p).
c        wvec        :  real working vector, wvec(maxm).
c
c=======================================================================                  

      implicit none 

c+++++Data
      integer nrec,p,interind(nrec)  
      real*8 x(nrec,p),y(nrec,2)

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay
      real*8 propv(p,p)

c+++++Prior information
      real*8 a0b0(2),aa0,ab0,betapm(p),betapv(p,p),m0,s0
      real*8 tau(2),tau1,tau2

c+++++Output
      real*8 acrate,randsave(nsave,nrec+1),thetasave(nsave,p+4)
      real*8 cpo(nrec)
      
c+++++Current values of the parameters
      integer ncluster
      real*8 alpha,beta(p),mu,sigma2,v(nrec)

c+++++Working space
      integer maxint 
      integer counter,dispcount,evali,i,indi,intcount(maxint),iflag(p)
      integer intcount2(maxint),intind(nrec+1),intind2(nrec+1)
      integer isave,iscan,j,j1,j2,k
      integer maxend,maxm
      integer mrand,npoints,ns,nso,nscan,ok
      integer seed(3),seed1,seed2,seed3
      integer sel,skipcount
      real*8 area,betac(p),cdflnorm
      real*8 eta(nrec),etan(nrec),endp(maxend),endp2(maxend)
      real*8 keepbeta,liminf,limsup,logliko,loglikn,logprioro,logpriorn
      real*8 mass(maxint),maxu,prob(maxint),prob2(maxint)
      real*8 ratio,rbeta,rgamma,rnorm,rtlnorm,sigma,swork,tmp1,tmp2
      real*8 uvec(maxm),workm1(p,p),workm2(p,p),workmh1(p*(p+1)/2)
      real*8 workv1(p),workv2(p),wvec(maxm),vbar,vvec(maxm),vnew(nrec+1)
      real*8 vnew2(nrec),vpred
      integer sprint
      logical ainf,binf
      real runif

c++++ initialize variables

      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      
      tau1=tau(1)
      tau2=tau(2)
      aa0=a0b0(1)
      ab0=a0b0(2)
      
      sigma=sqrt(sigma2)
      
      call rdisc(1,nrec,evali)
      vpred=v(evali)
      sel=1

c++++ set random number generator

      seed1=seed(1)
      seed2=seed(2)
      seed3=seed(3)

      call setrand(seed1,seed2,seed3)
      

c++++ evaluate log-prior for current value of parameters

      call dmvn(p,beta,betapm,betapv,logprioro,workv1,workm1,
     &          workm2,workv2,iflag)  


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*(nsave)

      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()
	
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ MH to updating regression coefficients and G +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ sample candidates

         call rmvnorm(p,beta,propv,workmh1,workv1,betac)
      
c+++++++ evaluate log-prior for candidate value of parameters

         call dmvn(p,betac,betapm,betapv,logpriorn,workv1,workm1,
     &             workm2,workv2,iflag)  

c+++++++ define the end points base on data
         
         counter=0
         do i=1,nrec
            tmp1=0.d0
            do j=1,p
               tmp1=tmp1+x(i,j)*betac(j)                   
            end do
            
            etan(i)=tmp1
            
            if(interind(i).eq.1)then
               counter=counter+1
               endp(counter)=y(i,2)*exp(eta(i))   
               counter=counter+1
               endp(counter)=y(i,2)*exp(etan(i))   
            end if
            
            if(interind(i).eq.3)then
               counter=counter+1
               endp(counter)=y(i,1)*exp(eta(i))   
               counter=counter+1
               endp(counter)=y(i,1)*exp(etan(i))   
            end if
            
            if(interind(i).eq.2)then
               counter=counter+1
               endp(counter)=y(i,1)*exp(eta(i))   
               counter=counter+1
               endp(counter)=y(i,1)*exp(etan(i))   
               counter=counter+1
               endp(counter)=y(i,2)*exp(eta(i))   
               counter=counter+1
               endp(counter)=y(i,2)*exp(etan(i))   
            end if   
         end do
         
c+++++++ sort end points

         do i=1,counter-1
            do j=i+1,counter
               if(endp(j).lt.endp(i))then
                  tmp1=endp(i)
                  endp(i)=endp(j)
                  endp(j)=tmp1
               end if
            end do
         end do   
         
         npoints=1
         endp2(1)=endp(1)
         do i=2,counter
            if(endp(i).ne.endp2(npoints))then
               npoints=npoints+1
               endp2(npoints)=endp(i)
            end if
         end do

c         call intpr("npoints",-1,npoints,1)

c+++++++ updating G base on the partition defined by beta and betac


         do i=1,maxint
            mass(i)=0.d0
            intcount2(i)=0
         end do

         do i=1,nrec
            indi=npoints+1
            do j=npoints,1,-1
               if(v(i).le.endp2(j))indi=indi-1
            end do
            mass(indi)=mass(indi)+1.d0
            intind2(i)=indi
            intcount2(indi)=intcount2(indi)+1
         end do

c+++++++ pred
         indi=npoints+1
         do j=npoints,1,-1
            if(vpred.le.endp2(j))indi=indi-1
         end do
         intcount2(indi)=intcount2(indi)+1
c+++++++ end pred

         do i=1,npoints+1
            area=0.d0
            if(i.eq.1)then
               area=alpha*cdflnorm(endp2(1),mu,sigma,1,0)
             else if(i.gt.npoints)then
               area=alpha*(1.d0-cdflnorm(endp2(npoints),mu,sigma,1,0))
             else
               area=alpha*(cdflnorm(endp2(i  ),mu,sigma,1,0)-  
     &                     cdflnorm(endp2(i-1),mu,sigma,1,0))
            end if
            mass(i)=mass(i)+area
         end do


         call dirichlet(mass,maxint,npoints+1,prob)
         

c         call dblepr("indi",-1,endp2,npoints)
         
c         return

c+++++++ evaluate log-likelihood for current and candidate value 
c+++++++ of parameters

         logliko=0.d0
         loglikn=0.d0
         do i=1,nrec
            tmp1=0.d0
            tmp2=0.d0

            if(interind(i).eq.1)then
               j=1 
               do while(y(i,2)*exp(eta(i)).ge.endp2(j).and.j.le.npoints)   
                  tmp1=tmp1+prob(j)
                  j=j+1
               end do

               j=1 
               do while(y(i,2)*exp(etan(i)).ge.endp2(j)
     &                  .and.j.le.npoints)   
                  tmp2=tmp2+prob(j)
                  j=j+1
               end do
            end if            

            if(interind(i).eq.3)then
               j=npoints 
               do while(y(i,1)*exp(eta(i)).le.endp2(j).and.j.ge.1)   
                  tmp1=tmp1+prob(j+1)
                  j=j-1
               end do

               j=npoints 
               do while(y(i,1)*exp(etan(i)).le.endp2(j).and.j.ge.1)   
                  tmp2=tmp2+prob(j+1)
                  j=j-1
               end do
            end if            

            if(interind(i).eq.2)then
               j1=npoints 
               do while(y(i,1)*exp(eta(i)).le.endp2(j1).and.j1.ge.1)   
                  j1=j1-1
               end do
               j1=j1+1
               
               j2=1 
               do while(y(i,2)*exp(eta(i)).ge.endp2(j2)
     &                 .and.j2.le.npoints)   
                  j2=j2+1
               end do
               j2=j2-1
               
               do j=j1,j2
                  tmp1=tmp1+prob(j)
               end do
               
               j1=npoints 
               do while(y(i,1)*exp(etan(i)).le.endp2(j1).and.j1.ge.1)   
                  j1=j1-1
               end do
               j1=j1+1
               
               j2=1 
               do while(y(i,2)*exp(etan(i)).ge.endp2(j2)
     &                  .and.j2.le.npoints)   
                  j2=j2+1
               end do
               j2=j2-1
               
               do j=j1,j2
                  tmp2=tmp2+prob(j)
               end do
            end if            

            logliko=logliko+log(tmp1)
            loglikn=loglikn+log(tmp2)
            
         end do

c+++++++ aceptation step

         ok=0
         ratio=dexp(loglikn+logpriorn-logliko-logprioro)

         if(dble(runif()).lt.ratio)then
            do j=1,p
               beta(j)=betac(j)
            end do
            logprioro=logpriorn
            do i=1,nrec
               eta(i)=etan(i)
            end do
            acrate=acrate+1.d0
           else 
            ok=1            
         end if

         if(ok.eq.1)go to 100 

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ updating errors using the Sethuraman-Tiwari (1982) +++
c+++++++ representation                                     +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do i=1,maxint
            intcount(i)=0
            prob2(i)=0.d0
         end do

         do i=1,nrec
         
            if(interind(i).eq.1)then
               do j=1,npoints+1
                  prob2(j)=0.d0
               end do
               
               j=1 
               do while(y(i,2)*exp(eta(i)).ge.endp2(j).and.j.le.npoints)   
                  prob2(j)=prob(j)
                  j=j+1
               end do
               
               call simdisc(prob2,maxint,npoints+1,evali)
            end if                              


            if(interind(i).eq.3)then
               do j=1,npoints+1
                  prob2(j)=0.d0
               end do

               j=npoints 
               do while(y(i,1)*exp(eta(i)).le.endp2(j).and.j.ge.1)   
                  prob2(j+1)=prob(j+1)
                  j=j-1
               end do

               call simdisc(prob2,maxint,npoints+1,evali)
            end if            


            if(interind(i).eq.2)then
               do j=1,npoints+1
                  prob2(j)=0.d0
               end do
               
               j1=npoints 
               do while(y(i,1)*exp(eta(i)).le.endp2(j1).and.j1.ge.1)   
                  j1=j1-1
               end do
               j1=j1+1
               
               j2=1 
               do while(y(i,2)*exp(eta(i)).ge.endp2(j2)
     &                 .and.j2.le.npoints)   
                  j2=j2+1
               end do
               j2=j2-1
               
               do j=j1,j2
                  prob2(j)=prob(j)
               end do
               
               call simdisc(prob2,maxint,npoints+1,evali)
            end if            

            if(evali.lt.1.and.evali.gt.npoints+1)return
            intcount(evali)=intcount(evali)+1
            intind(i)=evali
         end do

c+++++++ pred         
         do j=1,npoints+1
	    prob2(j)=prob(j)
         end do
         call simdisc(prob2,maxint,npoints+1,evali)
         intcount(evali)=intcount(evali)+1
         intind(nrec+1)=evali
c+++++++ end pred         

         do i=1,npoints+1
            ns=intcount(i)
            nso=intcount2(i)
            area=0.d0
            
            if(i.eq.1)then
               area=alpha*cdflnorm(endp2(1),mu,sigma,1,0)
               ainf=.true.
               binf=.false.
               liminf=0.d0
               limsup=endp2(1)
             else if(i.gt.npoints)then
               area=alpha*(1.d0-cdflnorm(endp2(npoints),mu,sigma,1,0))
               ainf=.false.
               binf=.true.
               liminf=endp2(npoints)
               limsup=0.d0
             else
               area=alpha*(cdflnorm(endp2(i  ),mu,sigma,1,0)-  
     &                     cdflnorm(endp2(i-1),mu,sigma,1,0))
               ainf=.false.
               binf=.false.
               liminf=endp2(i-1)
               limsup=endp2(i)
            end if
            
            if(ns.gt.0)then
	       maxu=0.d0
	       do j=1,ns
	          uvec(j)=dble(runif())
	          if(uvec(j).gt.maxu)maxu=uvec(j)
	       end do

               mrand=1
               keepbeta=rbeta(1.d0,mass(i))
               wvec(mrand)=keepbeta
               tmp1=1.d0-keepbeta
               tmp2=keepbeta
               
               do while(maxu.ge.tmp2.and.mrand.lt.maxm)
                  mrand=mrand+1
                  if(mrand.gt.maxm)then
                    call intpr('Increase maxn in R function',-1,0,1)   
                    return
                  end if  
                  keepbeta=rbeta(1.d0,mass(i))
                  wvec(mrand)=keepbeta*tmp1
                  tmp1=tmp1*(1.d0-keepbeta)
                  tmp2=tmp2+wvec(mrand)
               end do
               
               do j=1,mrand
                  if(dble(runif()).le.(1.d0-(area/mass(i))))then
                     if(nso.eq.0)return
                     call rdisc(1,nso,evali)
                     counter=0
                     k=1
                     do while(k.le.(nrec+1).and.counter.lt.evali)
                        if(intind2(k).eq.i)then
                           counter=counter+1 
                           sel=k
                        end if
                        k=k+1
                     end do
                     vvec(j)=v(sel)  
                   else
                     vvec(j)=rtlnorm(mu,sigma,ainf,binf,liminf,limsup)
                  end if
               end do
               
c               call dblepr("vvec",-1,vvec,mrand)

               j=1
               counter=0
               do while(j.le.(nrec+1).and.counter.lt.ns)
                  if(intind(j).eq.i)then
	            counter=counter+1 
	            tmp1=0.d0
	            k=1
	            do while(uvec(counter).gt.tmp1.and.k.lt.mrand)
	               tmp1=tmp1+wvec(k)
	               k=k+1  
	            end do
	            vnew(j)=vvec(k)
                  end if
                  j=j+1
               end do
               
            end if 
         end do

         do i=1,nrec
            v(i)=vnew(i)
         end do
         vpred=vnew(nrec+1)

100      continue


c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ updating parameters of the baseline                +++
c+++++++ distribution                                       +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


c+++++++ sort errors

c+++++++ check if the user has requested an interrupt
         call rchkusr()
 
         do i=1,nrec
            vnew(i)=v(i) 
         end do

         do i=1,nrec-1
            do j=i+1,nrec
               if(vnew(j).lt.vnew(i))then
                  tmp1=vnew(i)
                  vnew(i)=vnew(j)
                  vnew(j)=tmp1
               end if
            end do
         end do   
         
c+++++++ clusters and updating step

c+++++++ check if the user has requested an interrupt
         call rchkusr()
         
         swork=0.d0
         ncluster=1
         vnew2(1)=vnew(1)
         
         do i=2,nrec
            if(vnew(i).ne.vnew2(ncluster))then
               ncluster=ncluster+1
               vnew2(ncluster)=vnew(i)
               swork=swork+log(vnew(i))
            end if
         end do

         
         tmp1=(1.d0/s0)+(dble(ncluster)/sigma)
         
         vbar=(m0/s0)+(swork/sigma)
         
         vbar=vbar/tmp1
         tmp1=1.d0/tmp1

         mu=rnorm(vbar,tmp1)

         swork=0.d0
         do i=1,ncluster
            swork=swork+(log(vnew2(i))-mu)**2            
         end do
         
         sigma2=1.d0/
     &           rgamma(0.5d0*(dble(ncluster)+tau1),0.5d0*(swork+tau2))
         
         sigma=sqrt(sigma2)

c++++++++++++++++++++++++++++++++++         
c+++++++ Precision parameter
c++++++++++++++++++++++++++++++++++
         if(aa0.gt.0.d0)then
            call samalph(alpha,aa0,ab0,ncluster,nrec)
         end if 

c+++++++ save samples
         
         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1
               
c+++++++++++++ regression coefficient information
               do j=1,p
                  thetasave(isave,j)=beta(j)
               end do

c+++++++++++++ base line information
               thetasave(isave,p+1)=mu
               thetasave(isave,p+2)=sigma2

c+++++++++++++ cluster information
               thetasave(isave,p+3)=ncluster  
               thetasave(isave,p+4)=alpha

c+++++++++++++ cpo, errors and predictive information
               do j=1,nrec

                  randsave(isave,j)=v(j)
                  
                  tmp1=0.d0
                  if(interind(j).eq.1)then
                     k=1 
                     do while(y(j,2)*exp(eta(j)).ge.endp2(k)
     &                        .and.k.le.npoints)   
                        tmp1=tmp1+prob(k)
                        k=k+1
                     end do
                  end if      

                  if(interind(j).eq.3)then
                     k=npoints 
                     do while(y(j,1)*exp(eta(j)).le.endp2(k).and.k.ge.1)   
                        tmp1=tmp1+prob(k+1)
                        k=k-1
                     end do
                  end if            

                  if(interind(j).eq.2)then
                     j1=npoints 
                     do while(y(j,1)*exp(eta(j)).le.endp2(j1)
     &                        .and.j1.ge.1)   
                        j1=j1-1
                     end do
                     j1=j1+1
               
                     j2=1 
                     do while(y(j,2)*exp(eta(j)).ge.endp2(j2)
     &                        .and.j2.le.npoints)   
                        j2=j2+1
                     end do
                     j2=j2-1
               
                     do k=j1,j2
                        tmp1=tmp1+prob(k)
                     end do
                  end if     

                  cpo(j)=cpo(j)+1.0d0/tmp1 
               end do
               
               randsave(isave,nrec+1)=vpred

c+++++++++++++ print
               skipcount = 0
               if(dispcount.ge.ndisplay)then
c                  call intpr("isave",5,isave,1)
                  tmp1=sprint(isave,nsave)
                  dispcount=0
               end if   

            end if         
         end if
      end do


c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ post chain analysis                                +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      acrate=acrate/dble(nscan)      
      
      do i=1,nrec
         cpo(i)=dble(nsave)/cpo(i)
      end do

      return
      end


