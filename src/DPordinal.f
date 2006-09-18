    
c=======================================================================                      
      subroutine dpordinalp(nrec,p,x,y,interind,ncateg,
     &                      mcmc,nsave,propv,
     &                      a0b0,betapm,betapv,m0,s0,tau,
     &                      acrate,thetasave,randsave,cpo,
     &                      ncluster,alpha,beta,mu,sigma2,v,
     &                      maxint,maxend,maxm,
     &                      iflag,
     &                      imaxs,imaxsc,imins,iminsc,
     &                      intcount,intcount2,    
     &                      intind,intind2,
     &                      intposso,intpossn,
     &                      seed,
     &                      betac,eta,etan,endp,endp2,
     &                      linfs,linfsc,lsups,lsupsc,
     &                      mass,
     &                      prob,
     &                      uvec,vvec,vnew,vnew2,wvec,
     &                      workm1,workm2,
     &                      workmh1,workv1,workv2)      
c=======================================================================                      
c
c     Subroutine `dpordinalp' to run a Markov chain in the  
c     semiparametric regression model for ordinal data. 
c
c     Copyright: Alejandro Jara Vallejos, 2006
c
c     Version 1.0: 
c
c     Last modification: 02-07-2006.
c
c     This program is free software; you can redistribute it and/or modify
c     it under the terms of the GNU General Public License as published by
c     the Free Software Foundation; either version 2 of the License, or (at
c     your option) any later version.
c
c     This program is distributed in the hope that it will be useful, but
c     WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c     General Public License for more details.
c
c     You should have received a copy of the GNU General Public License
c     along with this program; if not, write to the Free Software
c     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
c
c     The author's contact information:
c
c     Alejandro Jara Vallejos
c     Biostatistical Centre
c     Katholieke Universiteit Leuven
c     U.Z. Sint-Rafaël
c     Kapucijnenvoer 35
c     B-3000 Leuven
c     Voice: +32 (0)16 336892 
c     Fax  : +32 (0)16 337015 
c     URL  : http://student.kuleuven.be/~s0166452/
c     Email: Alejandro.JaraVallejos@med.kuleuven.be
c
c---- Data -------------------------------------------------------------
c
c        interind    :  integer vector giving the type of interval for
c                       each observation, 1=left, 2=inter, and 3=right,
c                       interind(nrec)
c        ncateg      :  integer giving the number of levels in the
c                       ordinal response. 
c        nrec        :  integer giving the number of observations.
c        p           :  integer giving the number of fixed coefficients.
c        x           :  real matrix giving the design matrix for the 
c                       fixed effects, x(nrec,p). 
c        y           :  integer vector giving the ordinal response,
c                       y(nrec).
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
c        cutoff      :  real vector giving the current value of the
c                       cutoff points, cutoff(ncateg-1)
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
c        cdfnorm     :  cdf of a normal distribution.
c        counter     :  index.
c        dispcount   :  index. 
c        efind       :  integer function to find the intervals where the
c                       errors should be sampled from.
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
c        iflag       :  integer vector used to evaluate the prior
c                       distribution for the regression coefficients, 
c                       iflag(p).
c        imaxs       :  integer vector used to store the highest
c                       partition where each error can be sampled from,
c                       imaxs(nrec). 
c        imaxsc      :  integer vector used to store the highest
c                       partition defined for the candidate coefficients
c                       where each error can be sampled from,
c                       imaxs(nrec). 
c        imins       :  integer vector used to store the lowest
c                       partition where each error can be sampled from,
c                       imins(nrec). 
c        iminsc      :  integer vector used to store the lowest
c                       partition defined for the candidate coefficients
c                       where each error can be sampled from,
c                       imins(nrec). 
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
c        intposso    :  integer matrix giving the observations that
c                       belong to each cluster, intposso(maxint,nrec+1).
c        intpossn    :  integer matrix giving the observations that
c                       belong to each cluster, intpossn(maxint,nrec+1).
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        k           :  index.   
c        keepbeta    :  real working variable.
c        liminf      :  real working variable.
c        limsup      :  real working variable.
c        linfs       :  real vector used to store the lim inf of the
c                       interval where the error must be sampled,
c                       linfs(nrec).
c        linfsc      :  real vector used to store the lim inf of the
c                       interval defined for the candidate value of the 
c                       coefficients where the error must be sampled, 
c                       linfsc(nrec).
c        logliko     :  real working variable.
c        loglikn     :  real working variable.
c        logprioro   :  real working variable.
c        logpriorn   :  real working variable.
c        lsups       :  real vector used to store the lim sup of the
c                       interval where the error must be sampled,
c                       lsups(nrec).
c        lsupsc      :  real vector used to store the lim sup of the
c                       interval defined for the candidate value of the 
c                       coefficients where the error must be sampled, 
c                       lsupsc(nrec).
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
c        ratio       :  real working variable.
c        rbeta       :  real beta random number generator.
c        rgamma      :  real gamma random number generator.
c        rnorm       :  real normal random number generator.
c        runif       :  real uniform random number generator.
c        rtnorm      :  real truncated normal random number generator.
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        seed3       :  seed for random number generation.
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
      integer ncateg,nrec,p,interind(nrec),y(nrec)  
      real*8 x(nrec,p)

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
      real*8 alpha,beta(p)
      real*8 cutoff(ncateg-1)
      real*8 mu,sigma2,v(nrec)

c+++++Working space
      integer maxint 
      integer counter,dispcount
      integer efind
      integer evali
      integer i,iflag(p)
      integer imax,imaxs(nrec),imaxsc(nrec)
      integer imin,imins(nrec),iminsc(nrec)
      integer indi,intcount(maxint)
      integer intcount2(maxint),intind(nrec+1),intind2(nrec+1)
      integer intposso(maxint,nrec+1),intpossn(maxint,nrec+1)
      integer isave,iscan,j,k
      integer maxend,maxm
      integer mrand,npoints,ns,nso,nscan,ok
      integer seed(3),seed1,seed2,seed3
      integer skipcount
      
      real*8 area,betac(p),cdfnorm
      real*8 eta(nrec),etan(nrec)
      real*8 endp(maxend),endp2(maxend)
      real*8 keepbeta
      real*8 liminf,limsup,linfs(nrec),linfsc(nrec)
      real*8 logliko,loglikn,logprioro,logpriorn
      real*8 lsups(nrec),lsupsc(nrec)
      real*8 mass(maxint),maxu,prob(maxint)
      real*8 ratio,rbeta,rgamma,rnorm,rtnorm,sigma,swork,tmp1,tmp2
      real*8 uvec(maxm),workm1(p,p),workm2(p,p),workmh1(p*(p+1)/2)
      real*8 workv1(p),workv2(p),wvec(maxm),vbar,vvec(maxm),vnew(nrec+1)
      real*8 vnew2(nrec),vpred
      integer sprint
      logical ainf,binf
      real runif

c+++++CPU time
      real*8 sec00,sec0,sec1,sec

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


c++++ set random number generator

      seed1=seed(1)
      seed2=seed(2)
      seed3=seed(3)

      call setall(seed1,seed2)

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ Check consistency with the data before starting 
c++++ the algorithm
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c++++ check if the user has requested an interrupt
      call rchkusr()
 
      do i=1,nrec
         vnew(i)=v(i)
         if(interind(i).eq.1)then
           linfs(i)=0.0
           lsups(i)=cutoff(1)-eta(i)

           if(v(i).gt.lsups(i))then
              call intpr("i",-1,i,1)
              call intpr("type",-1,interind(i),1)
              call dblepr("u",-1,lsups(i),1)
              call dblepr("eta",-1,eta(i),1)
              call dblepr("v",-1,v(i),1)
              call rexit("Errors not consistent with data in S0")
           end if   
         end if
         
         if(interind(i).eq.3)then
           linfs(i)=cutoff(ncateg-1)-eta(i)
           lsups(i)=0.0
           if(v(i).le.linfs(i))then
              call intpr("i",-1,i,1)
              call intpr("type",-1,interind(i),1)
              call dblepr("l",-1,linfs(i),1)
              call dblepr("eta",-1,eta(i),1)
              call dblepr("v",-1,v(i),1)
              call rexit("Errors not consistent with data in S0")
           end if   
         end if
            
         if(interind(i).eq.2)then
           do j=2,ncateg-1
              if(y(i).eq.j)then
                 linfs(i)=cutoff(j-1)-eta(i)
                 lsups(i)=cutoff(j)-eta(i)
              end if
           end do
           if(v(i).le.linfs(i))then
              call intpr("i",-1,i,1)
              call intpr("type",-1,interind(i),1)
              call dblepr("l",-1,linfs(i),1)
              call dblepr("u",-1,lsups(i),1)
              call dblepr("eta",-1,eta(i),1)
              call dblepr("v",-1,v(i),1)
              call rexit("Errors not consistent with data in S0")
           end if   
           if(v(i).gt.lsups(i))then
              call intpr("i",-1,i,1)
              call intpr("type",-1,interind(i),1)
              call dblepr("l",-1,linfs(i),1)
              call dblepr("u",-1,lsups(i),1)
              call dblepr("eta",-1,eta(i),1)
              call dblepr("v",-1,v(i),1)
              call rexit("Errors not consistent with data in S0")
           end if   
         end if   
      end do

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ Sort errors and compute the number of clusters 
c++++ before starting the algorithm
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      call sortvec(nrec+1,vnew,1,nrec)

      ncluster=1
      vnew2(1)=vnew(1)
     
      do i=2,nrec
         call rchkusr()
         if(vnew(i).ne.vnew2(ncluster))then
            ncluster=ncluster+1
            vnew2(ncluster)=vnew(i)
         end if
      end do


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*(nsave)

      call cpu_time(sec0)
      sec00=0.d0

c++++ evaluate log-prior for current value of parameters

      call dmvn(p,beta,betapm,betapv,logprioro,workv1,workm1,
     &          workm2,workv2,iflag)  

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
               linfs(i)=0.0
               lsups(i)=cutoff(1)-eta(i)
               linfsc(i)=0.0
               lsupsc(i)=cutoff(1)-etan(i)
               counter=counter+1
               endp(counter)=cutoff(1)-eta(i)
               counter=counter+1
               endp(counter)=cutoff(1)-etan(i)
            end if
            
            if(interind(i).eq.3)then
               linfs(i)=cutoff(ncateg-1)-eta(i)
               lsups(i)=0.0
               linfsc(i)=cutoff(ncateg-1)-etan(i)
               lsupsc(i)=0.0
               counter=counter+1
               endp(counter)=cutoff(ncateg-1)-eta(i)
               counter=counter+1
               endp(counter)=cutoff(ncateg-1)-etan(i)
            end if
            
            if(interind(i).eq.2)then
               do j=2,ncateg-1
                  if(y(i).eq.j)then
                    linfs(i)=cutoff(j-1)-eta(i)
                    lsups(i)=cutoff(j)-eta(i)
                    linfsc(i)=cutoff(j-1)-etan(i)
                    lsupsc(i)=cutoff(j)-etan(i)
                  end if
               end do
               counter=counter+1
               endp(counter)=linfs(i)
               counter=counter+1
               endp(counter)=linfsc(i)
               counter=counter+1
               endp(counter)=lsups(i)
               counter=counter+1
               endp(counter)=lsupsc(i)
            end if   
         end do
         
c+++++++ sort end points

         call sortvec(maxend,endp,1,counter)
         
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
            
            intposso(indi,intcount2(indi))=i
         end do

         do i=1,npoints+1
            area=0.d0
            if(i.eq.1)then
               area=alpha*cdfnorm(endp2(1),mu,sigma,1,0)
             else if(i.gt.npoints)then
               area=alpha*(1.d0-cdfnorm(endp2(npoints),mu,sigma,1,0))
             else
               area=alpha*(cdfnorm(endp2(i  ),mu,sigma,1,0)-  
     &                     cdfnorm(endp2(i-1),mu,sigma,1,0))
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
            tmp2=0.d0

            if(interind(i).eq.1)then
               imin=1
               imax=efind(lsups(i),maxend,endp2,npoints)
               tmp1=0.d0
               do j=imin,imax
                  tmp1=tmp1+prob(j)
               end do
               imins(i)=imin
               imaxs(i)=imax

               imin=1
               imax=efind(lsupsc(i),maxend,endp2,npoints)
               tmp2=0.d0
               do j=imin,imax
                  tmp2=tmp2+prob(j)
               end do
               iminsc(i)=imin
               imaxsc(i)=imax
            end if   
            
            if(interind(i).eq.3)then
               imin=efind(linfs(i),maxend,endp2,npoints)
               imin=imin+1
               imax=npoints+1 
               tmp1=0.d0
               do j=imin,imax
                  tmp1=tmp1+prob(j)
               end do
               imins(i)=imin
               imaxs(i)=imax

               imin=efind(linfsc(i),maxend,endp2,npoints)
               imin=imin+1
               imax=npoints+1 
               tmp2=0.d0
               do j=imin,imax
                  tmp2=tmp2+prob(j)
               end do
               iminsc(i)=imin
               imaxsc(i)=imax
            end if   

            if(interind(i).eq.2)then
               imin=efind(linfs(i),maxend,endp2,npoints)
               imin=imin+1
               imax=efind(lsups(i),maxend,endp2,npoints)
               tmp1=0.d0
               do j=imin,imax
                  tmp1=tmp1+prob(j)
               end do
               imins(i)=imin
               imaxs(i)=imax
               
               imin=efind(linfsc(i),maxend,endp2,npoints)
               imin=imin+1
               imax=efind(lsupsc(i),maxend,endp2,npoints)
               tmp2=0.d0
               do j=imin,imax
                  tmp2=tmp2+prob(j)
               end do
               iminsc(i)=imin
               imaxsc(i)=imax
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
               linfs(i)=linfsc(i)
               lsups(i)=lsupsc(i)
               imins(i)=iminsc(i)
               imaxs(i)=imaxsc(i)
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
         end do

         do i=1,nrec
            imin=imins(i) 
            imax=imaxs(i)
            tmp1=0.d0
            do j=imin,imax
               tmp1=tmp1+prob(j)
            end do

            if(tmp1.eq.0.d0)then
               call intpr("Imin",-1,imin,1)
               call intpr("Imax",-1,imax,1)
               call rexit("Zero Probability")
              else
               call simdiscint(prob,maxint,imin,imax,evali)               
            end if   

            if(evali.lt.imin)then
               call intpr("Evali",-1,evali,1)
               call intpr("Imin",-1,imin,1)
               call intpr("Imax",-1,imax,1)
               call rexit("Error in the intervals")
            end if  
            if(evali.gt.imax)then
               call intpr("Evali",-1,evali,1)
               call intpr("Imin",-1,imin,1)
               call intpr("Imax",-1,imax,1)
               call rexit("Error in the intervals")
            end if  
            
            if(evali.lt.1.and.evali.gt.(npoints+1))then
               call intpr("Evali",-1,evali,1)
               call intpr("Imin",-1,imin,1)
               call intpr("Imax",-1,imax,1)
               call rexit("Error in the intervals")
            end if  

            intcount(evali)=intcount(evali)+1
            intind(i)=evali

            intpossn(evali,intcount(evali))=i
            
            imins(i)=imin
            imaxs(i)=imax
            
         end do


c+++++++ predictions         

         call simdiscint(prob,maxint,1,npoints+1,evali)                        
         intcount(evali)=intcount(evali)+1
         intind(nrec+1)=evali
         intpossn(evali,intcount(evali))=nrec+1

c+++++++ end predictions         

         do i=1,npoints+1
            ns=intcount(i)
            nso=intcount2(i)
            area=0.d0
            
            if(i.eq.1)then
               area=alpha*cdfnorm(endp2(1),mu,sigma,1,0)
               ainf=.true.
               binf=.false.
               liminf=0.d0
               limsup=endp2(1)
             else if(i.gt.npoints)then
               area=alpha*cdfnorm(endp2(npoints),mu,sigma,0,0)
               ainf=.false.
               binf=.true.
               liminf=endp2(npoints)
               limsup=0.d0
             else
               area=alpha*(cdfnorm(endp2(i  ),mu,sigma,1,0)-  
     &                     cdfnorm(endp2(i-1),mu,sigma,1,0))
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
                    call rexit("Increase maxn in R function")
                    return
                  end if  
                  keepbeta=rbeta(1.d0,mass(i))
                  wvec(mrand)=keepbeta*tmp1
                  tmp1=tmp1*(1.d0-keepbeta)
                  tmp2=tmp2+wvec(mrand)
               end do

               do j=1,mrand
                  if(dble(runif()).le.(1.d0-(area/mass(i))))then
                     call rdisc(1,nso,evali)
                     if(evali.lt.1.or.evali.gt.nso)return
                     vvec(j)=v(intposso(i,evali))
                   else
                     tmp1=rtnorm(mu,sigma,liminf,limsup,ainf,binf)
                     vvec(j)=tmp1
                  end if
               end do

               do j=1,ns
                  tmp1=0.d0
	          k=1
	          do while(uvec(j).gt.tmp1.and.k.lt.mrand)
	             tmp1=tmp1+wvec(k)
	             k=k+1  
	          end do
	          vnew(intpossn(i,j))=vvec(k)
               end do

            end if 
         end do

         do i=1,nrec
            v(i)=vnew(i)
         end do
         vpred=vnew(nrec+1)


c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Check consistency with the data    
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()
 
          do i=1,nrec
             if(interind(i).eq.1)then
               linfs(i)=0.0
               lsups(i)=cutoff(1)-eta(i)
 
              if(v(i).gt.lsups(i))then
                  call intpr("i",-1,i,1)
                  call intpr("type",-1,interind(i),1)
                  call dblepr("u",-1,lsups(i),1)
                  call dblepr("eta",-1,eta(i),1)
                  call dblepr("v",-1,v(i),1)
                  call rexit("Errors not consistent with data in S1")
               end if   
             end if
          
             if(interind(i).eq.3)then
               linfs(i)=cutoff(ncateg-1)-eta(i)
               lsups(i)=0.0
               if(v(i).le.linfs(i))then
                  call intpr("i",-1,i,1)
                  call intpr("type",-1,interind(i),1)
                  call dblepr("l",-1,linfs(i),1)
                  call dblepr("eta",-1,eta(i),1)
                  call dblepr("v",-1,v(i),1)
                  call rexit("Errors not consistent with data in S1")
               end if   
             end if
             
             if(interind(i).eq.2)then
               do j=2,ncateg-1
                  if(y(i).eq.j)then
                     linfs(i)=cutoff(j-1)-eta(i)
                     lsups(i)=cutoff(j)-eta(i)
                  end if
               end do
               if(v(i).le.linfs(i))then
                  call intpr("i",-1,i,1)
                  call intpr("type",-1,interind(i),1)
                  call dblepr("l",-1,linfs(i),1)
                  call dblepr("u",-1,lsups(i),1)
                  call dblepr("eta",-1,eta(i),1)
                  call dblepr("v",-1,v(i),1)
                  call rexit("Errors not consistent with data in S1")
               end if   
               if(v(i).gt.lsups(i))then
                  call intpr("i",-1,i,1)
                  call intpr("type",-1,interind(i),1)
                  call dblepr("l",-1,linfs(i),1)
                  call dblepr("u",-1,lsups(i),1)
                  call dblepr("eta",-1,eta(i),1)
                  call dblepr("v",-1,v(i),1)
                  call rexit("Errors not consistent with data in S1")
               end if   
             end if   
         end do
 
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Sort errors and compute the number of clusters 
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         call sortvec(nrec+1,vnew,1,nrec)

         ncluster=1
         vnew2(1)=vnew(1)
     
         do i=2,nrec
            call rchkusr()
            if(vnew(i).ne.vnew2(ncluster))then
               ncluster=ncluster+1
               vnew2(ncluster)=vnew(i)
            end if
         end do

100      continue


c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ updating parameters of the baseline                +++
c+++++++ distribution                                       +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ clusters and updating step

c+++++++ check if the user has requested an interrupt
         call rchkusr()
         
         swork=0.d0
         
         do i=1,ncluster
            swork=swork+log(vnew2(i))
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
               do i=1,nrec

                  randsave(isave,i)=v(i)

                  imin=imins(i)
                  imax=imaxs(i)
                  tmp1=0.d0
                  do j=imin,imax
                     tmp1=tmp1+prob(j)
                  end do

                  cpo(i)=cpo(i)+1.0d0/tmp1 
               end do
               
               randsave(isave,nrec+1)=vpred

c+++++++++++++ print
               skipcount = 0
               if(dispcount.ge.ndisplay)then
                  call cpu_time(sec1)
                  sec00=sec00+(sec1-sec0)
                  sec=sec00
                  sec0=sec1
                  tmp1=sprint(isave,nsave,sec)
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

