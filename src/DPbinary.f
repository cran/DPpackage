    
c=======================================================================                      
      subroutine dpbinaryl(model,nrec,p,sens,spec,x,yobs,
     &                     nlink,xlink,
     &                     a0b0,betapm,betapv,
     &                     mcmc,nsave,propv,
     &                     acrate,fsave,randsave,thetasave,cpo,
     &                     alpha,beta,ncluster,y,v,
     &                     betac,endp,endp2,eta,etan,fsavet,
     &                     intcount,intcount2,intind,intind2,iflag,
     &                     maxint,maxend,maxm,mass,
     &                     prob,prob2,
     &                     seed1,seed2,seed3,
     &                     uvec,vvec,vnew,vnew2,
     &                     workm1,workm2,
     &                     workmh1,workv1,workv2,wvec)
c=======================================================================                      
c
c     Version 1.0: 
c     Last modification: 09-05-2006.
c
c     Subroutine `dpbinaryl' to run a Markov chain in the  
c     semiparametric logistic regression model. 
c
c     Author: A.J.V.
c
c---- Data -------------------------------------------------------------
c
c        model       :  integer indicating if the model correct for
c                       missclasification (1) or not (0).
c        nrec        :  integer giving the number of observations.
c        p           :  integer giving the number of fixed coefficients.
c        sens        :  real vector of sensitivity, sens(nrec).
c        spec        :  real vector of specificity, spec(nrec).
c        x           :  real matrix giving the design matrix for the 
c                       fixed effects, x(nrec,p). 
c        yobs        :  integer vector giving the oberved binary 
c                       response, yobs(nrec).
c
c-----------------------------------------------------------------------
c
c---- Link information -------------------------------------------------
c       
c        nlink       :  integer giving the number of grid points to 
c                       evaluate the link.
c        xlink       :  real vector giving the grid points, xlink(nlink)
c
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
c        fsave       :  real matrix containing the mcmc samples for
c                       the link.
c        randsave    :  real matrix containing the mcmc samples for
c                       the errors and prediction,
c                       randsave(nsave,nrec+1).
c        thetasave   :  real vector containing the mcmc sample for the
c                       regression parameters, betsave(nsave,p+2). 
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
c        ncluster    :  integer giving the number of clusters in the
c                       random effects.
c        y           :  integer vector giving the current value of the
c                       true binary responses, y(nrec).
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
c        cdflogis    :  cdf of a logistic distribution.
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
c        fsavet      :  integer vector used to store the number of
c                       observations to calculate the predictive
c                       distribution of the link, fsavet(nlink).
c        i           :  index.
c        imax        :  index.
c        imin        :  index.
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
c        k           :  index.   
c        keepbeta    :  real working variable.
c        liminf      :  real working variable.
c        limsup      :  real working variable.
c        logliko     :  real working variable.
c        loglikn     :  real working variable.
c        logprioro   :  real working variable.
c        logpriorn   :  real working variable.
c        maxint      :  integer giving the maximum number of intervals.
c        maxend      :  integer giving the maximum number of intervals.
c        maxm        :  integer giving the maximum number of points in
c                       the finite approximation of the DP.
c        mass        :  real vector giving the mass of the DP,
c                       mass(maxint).
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
c        runif       :  real uniform random number generator.
c        rtslogistic2:  real truncated logistic random number generator.
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        seed3       :  seed for random number generation.
c        sel         :  index.
c        skipcount   :  index. 
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
c        uvec        :  real working vector, uvec(maxm).
c        vvec        :  real working vector, vvec(maxm).
c        vnew        :  real working vector, vnew(nrec+1).
c        vnew2       :  real working vector, vnew2(nrec).
c        vpred       :  real working variable. 
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
      
c+++++Observed variables
      integer model,nrec,p,yobs(nrec)
      real*8 sens(nrec),spec(nrec)
      real*8 x(nrec,p)

c+++++Link information
      integer nlink 
      real*8 xlink(nlink)
      
c+++++Prior information
      real*8 a0b0(2),aa0,ab0,betapm(p),betapv(p,p)

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay
      real*8 propv(p,p)

c+++++Stored output
      real*8 acrate,fsave(nsave,nlink)
      real*8 randsave(nsave,nrec+1),thetasave(nsave,p+2)
      real*8 cpo(nrec)
      
c+++++Current values of the parameters
      integer ncluster,y(nrec)
      real*8 alpha,beta(p),v(nrec)
      
c+++++Working space
      integer sprint
      integer maxint
      integer counter,dispcount
      integer evali
      integer fsavet(nlink)
      integer i,imax,imin,indi,intcount(maxint),intcount2(maxint)
      integer intind(nrec+1),intind2(nrec+1)
      integer iflag(p),isave,iscan
      integer j
      integer k
      integer maxend,maxm
      integer mrand,npoints,ns,nso,nscan,ok
      integer seed1,seed2,seed3
      integer sel,skipcount

      real*8 area
      real*8 betac(p),cdflogis
      real*8 endp(maxend),endp2(maxend)
      real*8 eta(nrec),etan(nrec)
      real*8 keepbeta
      real*8 liminf,limsup,logliko,loglikn,logprioro,logpriorn
      real*8 mass(maxint),maxu
      real*8 prob(maxint),prob2(maxint)
      real*8 ratio,rbeta,rtslogistic2
      real*8 tmp1,tmp2
      real*8 uvec(maxm)
      real*8 vpred
      real*8 vvec(maxm),vnew(nrec+1),vnew2(nrec)
      real*8 workm1(p,p),workm2(p,p),workmh1(p*(p+1)/2)
      real*8 workv1(p),workv2(p),wvec(maxm)

      real runif
      
      logical ainf,binf

      
c++++ initialize variables

      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      
      aa0=a0b0(1)
      ab0=a0b0(2)
      
      call rdisc(1,nrec,evali)
      vpred=v(evali)
      sel=1

c++++ set random number generator

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

c+++++++ define the end points
         
         counter=0
         do i=1,nrec
            tmp1=0.d0
            do j=1,p
               tmp1=tmp1+x(i,j)*betac(j)                   
            end do
            etan(i)=tmp1
            counter=counter+1
            endp(counter)=eta(i)
            counter=counter+1
            endp(counter)=etan(i)
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
               area=alpha*cdflogis(endp2(1),0.d0,1.d0,1,0)
             else if(i.gt.npoints)then
               area=alpha*cdflogis(endp2(npoints),0.d0,1.d0,0,0)
             else
               area=alpha*(cdflogis(endp2(i  ),0.d0,1.d0,1,0)-  
     &                     cdflogis(endp2(i-1),0.d0,1.d0,1,0))
            end if
            mass(i)=mass(i)+area
         end do

         call dirichlet(mass,maxint,npoints+1,prob)
         

c+++++++ evaluate log-likelihood for current and candidate value 
c+++++++ of parameters

        
         logliko=0.d0
         loglikn=0.d0
         do i=1,nrec
            tmp1=0.d0
            j=1 
            do while(eta(i).ge.endp2(j).and.j.le.npoints)
               tmp1=tmp1+prob(j)
               j=j+1
            end do
            tmp2=0.d0
            j=1 
            do while(etan(i).ge.endp2(j).and.j.le.npoints)
               tmp2=tmp2+prob(j)
               j=j+1
            end do
            
            if(yobs(i).eq.1)then
               tmp1=tmp1*sens(i)+(1-spec(i))*(1-tmp1) 
               if(tmp1.lt.0.d0)tmp1=0.d0
               if(tmp1.gt.1.d0)tmp1=1.d0
               tmp2=tmp2*sens(i)+(1-spec(i))*(1-tmp2) 
               if(tmp2.lt.0.d0)tmp2=0.d0
               if(tmp2.gt.1.d0)tmp2=1.d0
               logliko=logliko+log(tmp1)
               loglikn=loglikn+log(tmp2)
            end if   
            if(yobs(i).eq.0)then
               tmp1=tmp1*sens(i)+(1-spec(i))*(1-tmp1) 
               if(tmp1.lt.0.d0)tmp1=0.d0
               if(tmp1.gt.1.d0)tmp1=1.d0
               tmp2=tmp2*sens(i)+(1-spec(i))*(1-tmp2) 
               if(tmp2.lt.0.d0)tmp2=0.d0
               if(tmp2.gt.1.d0)tmp2=1.d0
               logliko=logliko+log(1.d0-tmp1)
               loglikn=loglikn+log(1.d0-tmp2)
            end if   
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


           if(model.eq.1)then
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++ updating the true binary variables               +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             do i=1,nrec
                tmp1=0.d0
                j=1 
                do while(eta(i).ge.endp2(j).and.j.le.npoints)
                   tmp1=tmp1+prob(j)
                   j=j+1
                end do
            
                if(yobs(i).eq.1)then
                   tmp1=sens(i)*tmp1/
     &                   (sens(i)*tmp1+(1.d0-spec(i))*(1.d0-tmp1))
                   call rbinom(1,tmp1,evali)
                end if   
                if(yobs(i).eq.0)then
                   tmp1=(1.d0-sens(i))*tmp1/
     &                   ((1.d0-sens(i))*tmp1+spec(i)*(1.d0-tmp1))
                   call rbinom(1,tmp1,evali)
                end if            
                y(i)=evali
             end do
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++ end updating the true binary variables           +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
           end if


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
            if(y(i).eq.1)then
              prob2(1)=prob(1)
              imin=npoints+1
              imax=1
              tmp1=0.d0
              do j=1,npoints
                 prob2(j+1)=0.d0
                 if(eta(i).gt.endp2(j))then
                    if(j.lt.imin)imin=j
                    if(j.gt.imax)imax=j
                    tmp1=tmp1+prob(j+1)
                    prob2(j+1)=prob(j+1)
                 end if  
              end do
              if(tmp1.eq.0.d0)then
                 call rdisc(imin,imax,evali)               
               else
                 call simdisc(prob2,maxint,npoints+1,evali)
              end if   
            end if

            if(y(i).eq.0)then
              prob2(npoints+1)=prob(npoints+1) 
              imin=npoints+1
              imax=1
              tmp1=0.d0
              do j=npoints,1,-1
                 prob2(j)=0.d0
                 if(eta(i).lt.endp2(j))then
                    if(j.lt.imin)imin=j
                    if(j.gt.imax)imax=j
                    tmp1=tmp1+prob(j)
                    prob2(j)=prob(j)
                 end if   
              end do
              if(tmp1.eq.0.d0)then
                 call rdisc(imin,imax,evali)               
               else
                 call simdisc(prob2,maxint,npoints+1,evali)
              end if   
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
               area=alpha*cdflogis(endp2(1),0.d0,1.d0,1,0)
               ainf=.true.
               binf=.false.
               liminf=0.d0
               limsup=endp2(1)
             else if(i.gt.npoints)then
               area=alpha*cdflogis(endp2(npoints),0.d0,1.d0,0,0)
               ainf=.false.
               binf=.true.
               liminf=endp2(npoints)
               limsup=0.d0
             else
               area=alpha*(cdflogis(endp2(i  ),0.d0,1.d0,1,0)-  
     &                     cdflogis(endp2(i-1),0.d0,1.d0,1,0))
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
                     vvec(j)=rtslogistic2(ainf,binf,liminf,limsup)
                  end if
               end do

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
c+++++++ measure                                            +++
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
            
c+++++++ clusters and updating step (not now)
c+++++++ check if the user has requested an interrupt
         call rchkusr()
         
         ncluster=1
         vnew2(1)=vnew(1)
         
         do i=2,nrec
            if(vnew(i).ne.vnew2(ncluster))then
               ncluster=ncluster+1
               vnew2(ncluster)=vnew(i)
            end if
         end do

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

c+++++++++++++ cluster information
               thetasave(isave,p+1)=ncluster  
               thetasave(isave,p+2)=alpha

c+++++++++++++ link information
c               do j=1,nlink
c                  tmp1=0.d0
c                  k=1 
c                  do while(xlink(j).ge.endp2(k).and.k.le.npoints)
c                     tmp1=tmp1+prob(k)
c                     k=k+1
c                  end do
c                  fsave(isave,j)=tmp1
c               end do            

               do i=1,nlink
                  fsavet(i)=0
               end do
                
               do i=1,nrec
                  j=1
                  do while(v(i).le.xlink(j).and.j.le.nlink)
                     fsavet(j)=fsavet(j)+1.d0
                     j = j + 1
                  end do
               end do

               do j=1,nlink
                  fsave(isave,j)=
     &             ( alpha*cdflogis(xlink(j),0.d0,1.d0,1,0) +
     &               dble(fsavet(j))  )/
     &             (alpha+dble(nrec))  
               end do

c+++++++++++++ cpo, errors and predictive information
               do j=1,nrec

                  randsave(isave,j)=v(j)

                  tmp1=0.d0
                  k=1 
                  do while(eta(j).ge.endp2(k).and.k.le.npoints)
                     tmp1=tmp1+prob(k)
                     k=k+1
                  end do

                  tmp1=sens(j)*tmp1+(1-spec(j))*(1-tmp1)
            
                  if(yobs(j).eq.1)then
                    cpo(j)=cpo(j)+1.0d0/tmp1 
                  end if   
                  if(yobs(j).eq.0)then
                    tmp1=1.0d0-tmp1 
                    cpo(j)=cpo(j)+1.0d0/tmp1 
                  end if            
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
      
      acrate=acrate/dble(nscan)      
      
      do i=1,nrec
         cpo(i)=dble(nsave)/cpo(i)
      end do
      
      return
      end


