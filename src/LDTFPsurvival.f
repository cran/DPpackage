
c=======================================================================
       subroutine ldtfpsurvival(nrec,ptf,pce,interind,xtf,xce,y,
     &                          ngrid,npredden,npredmed,grid,
     &                          xtfpred,xcepred,
     &                          xtfpredm,xcepredm,quans,
     &                          maxm,ntprob,ntlr,a0b0,betacepm,
     &                          gprior,precce,
     &                          tau,alpha,betace,betatf,sigma2,z,
     &                          mcmc,nsave,seed,
     &                          cpo,densm,densl,densu,
     &                          qmm,qml,qmu,survmm,survml,survmu,
     &                          thetasave,randsave,
     &                          iflag,iflagtf,
     &                          nobsbc,obsbc,c0,
     &                          workm1,workvh1,workv1,
     &                          worksam,worksam2,worksam3,fs,
     &                          workm2,workvh2,
     &                          workv2,workv3,workv4,k,
     &                          prob,probc)
c=======================================================================
c     # 64 arguments
c
c     Subroutine `ldtfpsurvival' to run a Markov chain for a 
c     Linear Dependent Tailfree Process prior for 
c     conditional density estimation.
c
c     Copyright: Alejandro Jara, 2011.
c
c     Version 1.0: 
c
c     Last modification: 11-11-2011.
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
c      Alejandro Jara
c      Department of Statistics
c      Facultad de Matematicas
c      Pontificia Universidad Catolica de Chile
c      Casilla 306, Correo 22 
c      Santiago
c      Chile
c      Voice: +56-2-3544506  URL  : http://www.mat.puc.cl/~ajara
c      Fax  : +56-2-3547729  Email: atjara@uc.cl
c
c---- Data -------------------------------------------------------------
c 
c        interind    :  integer matrix giving the type of data,
c                       interind(nrec).
c        nrec        :  integer giving the number of data points. 
c        ptf         :  integer giving the number of predictors for
c                       the tailfree conditional probabilities.
c        ptf         :  integer giving the number of predictors for
c                       the median function.
c        xtf         :  real matrix giving the design matrix for
c                       the conditional talfree probabilities, 
c                       xtf(nrec,ptf).
c        xce         :  real matrix giving the design matrix for
c                       the median function, 
c                       xce(nrec,ptf).
c        y           :  real matrix giving the limits of the 
c                       intervals, y(nrec,2). 
c
c-----------------------------------------------------------------------
c
c---- Prediction -------------------------------------------------------
c 
c        cband       :  integer value indicating whether the 
c                       credible bands need to be computed or not.
c        ngrid       :  integer giving the number of grid points where 
c                       the density estimates are evaluated.
c        npredden    :  integer giving the number of predictions for
c                       the conditional densities.
c        npredmed    :  integer giving the number of predictions for
c                       the quantile functions.
c        tband       :  integer indicating the type of credible 
c                       band that need to be computed.
c        grid        :  real vector giving the grid, grid(ngrid). 
c        xtfpred     :  real matrix giving the design matrix, for the 
c                       conditional probabilities, used in the
c                       estimation of the conditional densities, 
c                       xtfpred(npreden,ptf).
c        xcepred     :  real matrix giving the design matrix, for the 
c                       median function, used in the
c                       estimation of the conditional densities, 
c                       xcepred(npreden,ptf).
c        xtfpredm    :  real matrix giving the design matrix, for the 
c                       conditional probabilities, used in the
c                       estimation of the quantile functions, 
c                       xtfpredm(npredmed,ptf).
c        xcepredm    :  real matrix giving the design matrix, for the 
c                       median function, used in the
c                       estimation of the quntile functions, 
c                       xcepred(npredmed,ptf).
c        quans       :  real vector giving the quantile functions
c                       to be evaluated, quans(3).
c
c-----------------------------------------------------------------------
c
c---- Prior information ------------------------------------------------
c
c        a0b0        :  real vector giving the hyperparameters of the
c                       prior distribution for the precision parameter,
c                       alpha ~ Gamma(a0b0(1),a0b0(2)). If a0b0(1)<0 
c                       the precision parameter is considered constant.
c        betapm      :  double precision vector giving the mean of the
c                       normal prior distribution for the median
c                       regression coefficients, betapm(pce).
c        maxm        :  integer giving the finite level of the
c                       dependent tailfree process.
c        precce      :  double precision matrix giving the precision
c                       matrix of the normal prior distribution for
c                       the median regression coefficients,
c                       precce(pce,pce).
c        ntprob      :  integer giving the number of final sets
c                       involved in the partition.
c        ntlr        :  integer giving the number of logistic
c                       regressions involved. 
c        tau         :  real vector giving the values of hyperparameters 
c                       of the inv-Gamma prior for sigma^{-2}, tau(2).
c        gprior      :  real matrix giving the common component of the
c                       covariance matrix for the tailfree logistic 
c                       regression coefficients, gprior(ptf,ptf).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        alpha       :  real giving the current value of the precision
c                       parameter of the DTFP.
c        betace      :  real vector giving the value of the median 
c                       regression coefficients, betace(pce).
c        betatf      :  real matrix giving the value of the tailfree
c                       regression coefficients, betatf(ntlr,ptf).
c        sigma2      :  real giving the value of the centering 
c                       variance.
c        z           :  real vector giving the imputed values
c                       of the log-survival times, z(nrec).
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
c        nbase       :  integer giving the the number of scans where 
c                       the baseline distribution and the precision
c                       parameter are sampled.
c        seed        :  integer vector giving the seeds for the random
c                       number generator, seed(2).
c        
c-----------------------------------------------------------------------
c
c---- Output -----------------------------------------------------------
c
c        cpo         :  real giving the cpos and fsos, cpo(nrec,2). 
c        densm       :  real matrix giving the posterior mean of the 
c                       density, densm(npredden,ngrid).
c        densl       :  real matrix giving the lower limit of the 
c                       HPD of the density, densu(npredden,ngrid).
c        densu       :  real matrix giving the upper limit of the  
c                       HPD of the density, densu(npredden,ngrid).
c        qmm         :  real matrix giving the posterior mean of the 
c                       quantile functions, qmm(npredmed,3).
c        qml         :  real matrix giving the lower limit of the 
c                       HPD of the quantile functions, 
c                       qml(npredmed,3).
c        qmu         :  real matrix giving the upper limit of the  
c                       HPD of the quantile functions, 
c                       qmu(npredmed,3).
c        randsave    :  real matrix containing the mcmc samples for
c                       the tailfree regression coefficients, 
c                       randsave(nsave,(ntlr-1)*ptf).
c        survmm      :  real matrix giving the posterior mean of the 
c                       survival functions, survmm(npredden,ngrid).
c        survml      :  real matrix giving the lower limit of the 
c                       HPD of the survival functions, 
c                       survml(npredden,ngrid).
c        survmu      :  real matrix giving the upper limit of the  
c                       HPD of the survival functions, 
c                       survmu(npredden,ngrid).
c        thetasave   :  real matrix containing the mcmc samples for
c                       the parameters, thetasave(nsave,pce+2)
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        c0          :  real matrix used to update the reg. coeff.,
c                       c0(ptf,ptf).
c        fs          :  real vector used to evaluate the conditional
c                       densities, fs(ngrid).
c        iflag       :  integer vector used to update reg. coeff.,
c                       iflag(pce).
c        iflagtf     :  integer vector used to update reg. coeff.,
c                       iflag(ptf).
c        k           :  integer vector, k(maxm).
c        nobsbc      :  integer vector giving the number of
c                       observations in each parition set,
c                       nobsbc(ntprob).
c        obsbc       :  integer matrix giving the 
c                       observations in each partition set,
c                       obsbc(ntprob,nrec).
c        prob        :  real working vector, prob(2**maxm).
c        probc       :  real working vector, prob(2**maxm).
c        workm1      :  real matrix used to update the reg. coeff.,
c                       workm1(pce,pce).
c        workm2      :  real matrix used to update the reg. coeff.,
c                       workm2(ptf,ptf).
c        workmh1     :  real vector used to update the reg. coeff.,
c                       workmh1(pce*(pce+1)/2).
c        workmh2     :  real vector used to update the reg. coeff.,
c                       workmh2(ptf*(ptf+1)/2).
c        worksam     :  real vector used to comput HPD bands,
c                       worksam(nsave).
c        worksam2    :  real matrix used to comput HPD bands,
c                       worksam2(nsave,ngrid).
c        worksam3    :  real matrix used to comput HPD bands,
c                       worksam3(nsave,npredden).
c        workv1      :  real vector used to update the reg. coeff.,
c                       workv1(pce).
c        workv2      :  real vector used to update the reg. coeff.,
c                       workv2(ptf).
c        workv3      :  real vector used to update the reg. coeff.,
c                       workv3(ptf).
c        workv4      :  real vector used to update the reg. coeff.,
c                       workv4(ptf).
c
c=======================================================================                  

      implicit none
c++++ input
      integer nrec,ptf,pce
      integer interind(nrec)
      double precision xtf(nrec,ptf)
      double precision xce(nrec,pce)
      double precision y(nrec,2)

c++++ prediction
      integer ngrid,npredden,npredmed
      double precision grid(ngrid)
      double precision xtfpred(npredden,ptf)
      double precision xcepred(npredden,pce)
      double precision xtfpredm(npredmed,ptf)
      double precision xcepredm(npredmed,pce)
      double precision quans(3)

c++++ prior
      integer maxm,ntprob,ntlr
      double precision a0b0(2)
      double precision betacepm(pce)
      double precision gprior(ptf,ptf)
      double precision precce(pce,pce)
      double precision tau(2),tau1,tau2

c++++ current value
      double precision alpha
      double precision betace(pce)
      double precision betatf(ntlr,ptf)
      double precision sigma2
      double precision z(nrec)

c++++ mcmc
      integer mcmc(5),nburn,nskip,nsave,ndisplay

c++++ bands
      integer cband,tband

c++++ seeds
      integer seed1,seed2,seed(2)

c++++ output
      double precision cpo(nrec,2)
      double precision densm(npredden,ngrid)
      double precision densl(npredden,ngrid)
      double precision densu(npredden,ngrid)
      double precision survmm(npredden,ngrid)
      double precision survml(npredden,ngrid)
      double precision survmu(npredden,ngrid)

      double precision qmm(npredmed,3)
      double precision qml(npredmed,3)
      double precision qmu(npredmed,3)

      double precision thetasave(nsave,pce+2)
      double precision randsave(nsave,(ntlr-1)*ptf)

c++++ external working space - vector and matrices
      integer iflag(pce)
      integer iflagtf(ptf)
      integer nobsbc(ntprob)
      integer obsbc(ntprob,nrec)
      double precision c0(ptf,ptf)
      double precision workm1(pce,pce)
      double precision workvh1(pce*(pce+1)/2)
      double precision workv1(pce)

      double precision worksam(nsave)
      double precision worksam2(nsave,ngrid)
      double precision worksam3(nsave,npredden)
      double precision fs(ngrid)

      double precision workm2(ptf,ptf)
      double precision workvh2(ptf*(ptf+1)/2)
      double precision workv2(ptf)
      double precision workv3(ptf)
      double precision workv4(ptf)

      integer k(maxm)

      double precision prob(2**maxm)
      double precision probc(2**maxm)
 
c++++ internal working space
c++++ internal working space
      integer accept,acceptt
      integer i,ii,i1,j,jj,j1,k1,ll,mm
      integer d1,d2
      integer dispcount,iscan,isave
      integer ntot,n1,n2
      integer nscan,skipcount
      integer sprint
      double precision liminf
      double precision limsup
      double precision loglikn,logliko
      double precision logpriorn
      double precision rgamma
      real runif
      double precision tmp1,tmp2,tmp3
      double precision uni
  
c++++ CPU time
      double precision sec00,sec0,sec1,sec

c++++ Working space slice sampling
      integer evali
      double precision re,uwork
      double precision logy,xx0,xx1,llim,rlim
      double precision grlim,gllim,gxx1

c++++ initializing variables

      open(unit=1,file='dppackage_dftp1.out',status='unknown',
     &     form='unformatted')
      open(unit=2,file='dppackage_dftp2.out',status='unknown',
     &     form='unformatted')
      open(unit=3,file='dppackage_dftp3.out',status='unknown',
     &     form='unformatted')
      open(unit=4,file='dppackage_dftp4.out',status='unknown',
     &     form='unformatted')
      open(unit=5,file='dppackage_dftp5.out',status='unknown',
     &     form='unformatted')
 
      tau1=tau(1)
      tau2=tau(2)  
 
      call loglikldtfp(maxm,ntlr,ntprob,nrec,pce,ptf,betace,betatf,z,
     &                 xce,xtf,sigma2,
     &                 nobsbc,obsbc,logliko,k)
      i1=1
      do i=2,maxm
         j1=2**(i-1)

         do j=1,j1
            k1=i1+j

            c0=gprior/(alpha*dble(i**2))

            ii=(k1-1)*2+1
            jj=(k1-1)*2+2
            n1=nobsbc(ii)
            n2=nobsbc(jj)
            ntot=n1+n2

            if(n1.gt.0.and.n2.gt.0)then

               call startlrcoefldtfp(5,k1,ii,jj,n1,n2,maxm,ntlr,
     &                               ntprob,nrec,ptf,obsbc,betatf,
     &                               xtf,c0,iflagtf,workv2,workm2,
     &                               workv3)
            end if
         end do
         i1=i1+j1
      end do
  
      call loglikldtfp(maxm,ntlr,ntprob,nrec,pce,ptf,betace,
     &                 betatf,z,xce,xtf,sigma2,
     &                 nobsbc,obsbc,logliko,k)

      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)

      cband=mcmc(4)
      tband=mcmc(5)
   
      seed1=seed(1)
      seed2=seed(2)
      call setall(seed1,seed2)

c++++ MCMC
      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*nsave

      call cpu_time(sec0)
      sec00=0.d0

      do iscan=1,nscan
 
c++++++++++++++++++++++++++++++++++++++
c+++++++ imputing censored data     +++
c++++++++++++++++++++++++++++++++++++++

         do i=1,nrec

c++++++++++ check if the user has requested an interrupt
            call rchkusr()
            
            if(interind(i).eq.1)then
               liminf=-999.d0
               limsup=dlog(y(i,2))
            end if

            if(interind(i).eq.2)then
               liminf=dlog(y(i,1))
               limsup=dlog(y(i,2))
            end if

            if(interind(i).eq.3)then
               liminf=dlog(y(i,1))
               limsup=999.d0
            end if

            if(interind(i).ne.4)then
 
               evali=1
               xx0=z(i)

               call ldensldtfp(xx0,i,maxm,nrec,ntlr,
     &                         ntprob,pce,betace,xce,
     &                         ptf,betatf,xtf,
     &                         sigma2,tmp1,k)

               uni=runif()
               re=-log(uni)
               logy=tmp1-re

               uwork=dble(runif())*0.05d0  
               llim=xx0-uwork
               rlim=xx0+(0.05d0-uwork)
           
               if(llim.lt.liminf)llim=liminf 
               if(rlim.gt.limsup)rlim=limsup
 
               evali=evali+1
               call ldensldtfp(llim,i,maxm,nrec,ntlr,
     &                         ntprob,pce,betace,xce,
     &                         ptf,betatf,xtf,
     &                         sigma2,gllim,k)

               evali=evali+1
               call ldensldtfp(rlim,i,maxm,nrec,ntlr,
     &                         ntprob,pce,betace,xce,
     &                         ptf,betatf,xtf,
     &                         sigma2,grlim,k)

               do while(gllim.gt.logy)
                  llim=llim-0.05d0

                  if(llim.lt.liminf)then
                     llim=liminf
                     gllim=logy-1.d0
                    else   
                     evali=evali+1
                     call ldensldtfp(llim,i,maxm,
     &                               nrec,ntlr,ntprob,
     &                               pce,betace,xce,
     &                               ptf,betatf,xtf,
     &                               sigma2,gllim,k)
                  end if
               end do 

               do while(grlim.gt.logy)
                  rlim=rlim+0.05d0

                  if(rlim.gt.limsup)then
                     rlim=limsup
                     grlim=logy-1.d0
                    else   
                     evali=evali+1
                     call ldensldtfp(rlim,i,maxm,nrec,
     &                               ntlr,ntprob,pce,
     &                               betace,xce,
     &                               ptf,betatf,xtf,
     &                               sigma2,grlim,k)
                  end if
               end do 

               xx1=llim+(rlim-llim)*dble(runif())

               evali=evali+1
               call ldensldtfp(xx1,i,maxm,nrec,ntlr,
     &                         ntprob,pce,betace,xce,
     &                         ptf,betatf,xtf,
     &                         sigma2,gxx1,k)

               do while(gxx1.lt.logy)
                  if(xx1.gt.xx0)rlim=xx1
                  if(xx1.lt.xx0)llim=xx1

                  if(llim.lt.liminf)llim=liminf 
                  if(rlim.gt.limsup)rlim=limsup 

                  xx1=llim+(rlim-llim)*dble(runif())

                  evali=evali+1
                  call ldensldtfp(xx1,i,maxm,nrec,
     &                            ntlr,ntprob,pce,
     &                            betace,xce,
     &                            ptf,betatf,xtf,
     &                            sigma2,gxx1,k)
               end do
               z(i)=xx1
 
c             call dblepr("z",-1,z(i),1)
            end if 
         end do

c        call dblepr("z",-1,z,nrec)

c++++++++++++++++++++++++++++++++++++++
c+++++++ baseline reg coefficients  +++
c++++++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do i=1,pce

c++++++++++ check if the user has requested an interrupt
            call rchkusr()
       
            evali=1
            xx0=betace(i)

            call logposldtfp(maxm,ntlr,ntprob,nrec,pce,
     &                       ptf,betace,betatf,z,xce,xtf,
     &                       sigma2,betacepm,precce,
     &                       nobsbc,obsbc,tmp1,k)

            uni=runif()
            re=-log(uni)
            logy=tmp1-re

            uwork=dble(runif())*0.25d0  
            llim=xx0-uwork
            rlim=xx0+(0.25d0-uwork)


            evali=evali+1
            betace(i)=llim
            call logposldtfp(maxm,ntlr,ntprob,nrec,pce,
     &                       ptf,betace,betatf,z,xce,xtf,
     &                       sigma2,betacepm,precce,
     &                       nobsbc,obsbc,gllim,k)

            evali=evali+1
            betace(i)=rlim
            call logposldtfp(maxm,ntlr,ntprob,nrec,pce,
     &                       ptf,betace,betatf,z,xce,xtf,
     &                       sigma2,betacepm,precce,
     &                       nobsbc,obsbc,grlim,k)

            do while(gllim.gt.logy)
               llim=llim-0.25d0
               evali=evali+1
               betace(i)=llim
               call logposldtfp(maxm,ntlr,ntprob,nrec,pce,
     &                       ptf,betace,betatf,z,xce,xtf,
     &                       sigma2,betacepm,precce,
     &                       nobsbc,obsbc,gllim,k)
            end do 

            do while(grlim.gt.logy)
               rlim=rlim+0.25d0
               betace(i)=rlim
               call logposldtfp(maxm,ntlr,ntprob,nrec,pce,
     &                       ptf,betace,betatf,z,xce,xtf,
     &                       sigma2,betacepm,precce,
     &                       nobsbc,obsbc,grlim,k)
            end do 

            xx1=llim+(rlim-llim)*dble(runif())
            betace(i)=xx1
            call logposldtfp(maxm,ntlr,ntprob,nrec,pce,
     &                       ptf,betace,betatf,z,xce,xtf,
     &                       sigma2,betacepm,precce,
     &                       nobsbc,obsbc,gxx1,k)

            do while(gxx1.lt.logy)
               if(xx1.gt.xx0)rlim=xx1
               if(xx1.lt.xx0)llim=xx1

               xx1=llim+(rlim-llim)*dble(runif())
               betace(i)=xx1
               call logposldtfp(maxm,ntlr,ntprob,nrec,pce,
     &                       ptf,betace,betatf,z,xce,xtf,
     &                       sigma2,betacepm,precce,
     &                       nobsbc,obsbc,gxx1,k)
            end do
            betace(i)=xx1

         end do

c++++++++++++++++++++++++++++++++++++++
c+++++++ baseline variance          +++
c++++++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         evali=1
         xx0=sigma2

         call loglikldtfp(maxm,ntlr,ntprob,nrec,pce,ptf,betace,
     &                    betatf,z,xce,xtf,xx0,
     &                    nobsbc,obsbc,loglikn,k)
         if(tau1>0)then
            logpriorn=-(0.5d0*tau1+1.d0)*log(xx0)-0.5d0*tau2/xx0
          else
            logpriorn=-log(xx0)
         end if      
         tmp1=loglikn+logpriorn

         uni=runif()
         re=-log(uni)
         logy=tmp1-re

         uwork=dble(runif())*0.05d0  
         llim=xx0-uwork
         rlim=xx0+(0.05d0-uwork)
           
         if(llim.lt.0.00001d0)llim=0.00001d0 

         evali=evali+1
         call loglikldtfp(maxm,ntlr,ntprob,nrec,pce,ptf,betace,
     &                    betatf,z,xce,xtf,llim,
     &                    nobsbc,obsbc,loglikn,k)
         if(tau1>0)then
            logpriorn=-(0.5d0*tau1+1.d0)*log(llim)-
     &                  0.5d0*tau2/llim
          else
            logpriorn=-log(llim)
         end if      
         gllim=loglikn+logpriorn

         evali=evali+1
         call loglikldtfp(maxm,ntlr,ntprob,nrec,pce,ptf,betace,
     &                    betatf,z,xce,xtf,rlim,
     &                    nobsbc,obsbc,loglikn,k)
         if(tau1>0)then
            logpriorn=-(0.5d0*tau1+1.d0)*log(rlim)-
     &                  0.5d0*tau2/rlim
          else
            logpriorn=-log(rlim)
         end if      
         grlim=loglikn+logpriorn

         do while(gllim.gt.logy)
            llim=llim-0.05d0

            if(llim.lt.0.00001d0)then
               llim=0.00001d0 
               gllim=logy-1.d0
              else   
               evali=evali+1

               call loglikldtfp(maxm,ntlr,ntprob,nrec,pce,ptf,
     &                          betace,betatf,z,
     &                          xce,xtf,llim,
     &                          nobsbc,obsbc,loglikn,k)
               if(tau1>0)then
                  logpriorn=-(0.5d0*tau1+1.d0)*log(llim)-
     &                        0.5d0*tau2/llim
                 else
                  logpriorn=-log(llim)
               end if      
               gllim=loglikn+logpriorn
             end if
         end do 

         do while(grlim.gt.logy)
            rlim=rlim+0.05d0

            evali=evali+1
            call loglikldtfp(maxm,ntlr,ntprob,nrec,pce,ptf,
     &                       betace,betatf,z,
     &                       xce,xtf,rlim,
     &                       nobsbc,obsbc,loglikn,k)
            if(tau1>0)then
               logpriorn=-(0.5d0*tau1+1.d0)*log(rlim)-
     &                     0.5d0*tau2/rlim
             else
               logpriorn=-log(rlim)
            end if      
            grlim=loglikn+logpriorn
         end do 

         xx1=llim+(rlim-llim)*dble(runif())
 
         evali=evali+1
         call loglikldtfp(maxm,ntlr,ntprob,nrec,pce,ptf,
     &                    betace,betatf,z,
     &                    xce,xtf,xx1,
     &                    nobsbc,obsbc,loglikn,k)
         if(tau1>0)then
           logpriorn=-(0.5d0*tau1+1.d0)*log(xx1)-
     &                 0.5d0*tau2/xx1
         else
           logpriorn=-log(xx1)
        end if      
        gxx1=loglikn+logpriorn

        do while(gxx1.lt.logy)
           if(xx1.gt.xx0)rlim=xx1
           if(xx1.lt.xx0)llim=xx1

           if(llim.lt.0.00001d0)llim=0.00001d0 

           xx1=llim+(rlim-llim)*dble(runif())

           evali=evali+1

           call loglikldtfp(maxm,ntlr,ntprob,nrec,pce,ptf,
     &                      betace,betatf,z,
     &                      xce,xtf,xx1,
     &                      nobsbc,obsbc,loglikn,k)
           if(tau1>0)then
              logpriorn=-(0.5d0*tau1+1.d0)*log(xx1)-
     &                    0.5d0*tau2/xx1
            else
              logpriorn=-log(xx1)
           end if
           gxx1=loglikn+logpriorn
        end do

        sigma2=xx1


c+++++++++++++++++++++++++++++++++++++
c++++++ tf logistic regressions    +++
c+++++++++++++++++++++++++++++++++++++

c+++++++++++++++++++++++++++++++++++++
c++++++ tf logistic regressions    +++
c+++++++++++++++++++++++++++++++++++++

        acceptt=0 

        i1=1
        do i=2,maxm
           j1=2**(i-1)

           do j=1,j1

c++++++++++++ check if the user has requested an interrupt
              call rchkusr()

              k1=i1+j

              do d1=1,ptf
                 do d2=1,ptf 
                    c0(d1,d2)=gprior(d1,d2)/
     &                        (alpha*real(i**2))
                 end do
              end do
              ii=(k1-1)*2+1
              jj=(k1-1)*2+2
              n1=nobsbc(ii)
              n2=nobsbc(jj)
              ntot=n1+n2

              if(ntot.gt.0)then

                 call updatelrcoefldtfpss(k1,ii,jj,n1,n2,maxm,
     &                               ntlr,ntprob,
     &                               nrec,ptf,obsbc,betatf,xtf,c0,
     &                               accept,iflagtf,workv2,workm2)

               else
                 accept=1
                 call updatelrcoefldtfp0(k1,maxm,ntlr,ntprob,
     &                                   ptf,betatf,c0,
     &                                   iflagtf,workv2,workv3,
     &                                   workm2,workvh2,workv4)
              end if
 
              acceptt=acceptt+accept

           end do
           i1=i1+j1
        end do

c+++++++++++++++++++++++++++++++++++++
c++++++ alpha                      +++
c+++++++++++++++++++++++++++++++++++++
     
        if(a0b0(1)>0)then

c+++++++++ check if the user has requested an interrupt
           call rchkusr()

           i1=1
           tmp1=0.d0

           do i=2,maxm
              j1=2**(i-1)
        
              do j=1,j1
                 k1=i1+j

                 do d1=1,ptf
                    do d2=1,ptf 
                       c0(d1,d2)=gprior(d1,d2)/
     &                        (dble(i**2))
                    end do
                 end do

                 call inverse(c0,ptf,iflagtf)
                 tmp2=0.d0 
                 do ll=1,ptf
                    do mm=1,ptf
                       tmp2=tmp2+betatf(k1,ll)*c0(ll,mm)*betatf(k1,mm)
                    end do
                 end do
                 tmp1=tmp1+tmp2
              end do
              i1=i1+j1
           end do
           tmp1=0.5d0*tmp1
           alpha=rgamma(a0b0(1)+dble(ptf)*dble(ntlr-1),a0b0(2)+tmp1)

c           call dblepr("alpha",-1,alpha,1)     
        end if

c++++++++++++++++++++++++++++++++         
c++++++ save samples
c++++++++++++++++++++++++++++++++         

        if(iscan.gt.nburn)then
           skipcount=skipcount+1
           if(skipcount.gt.nskip)then
              isave=isave+1
              dispcount=dispcount+1

c++++++++++++ regression coefficient information

              do i=1,pce
                 thetasave(isave,i)=betace(i)
              end do

c++++++++++++ centering variance

              thetasave(isave,pce+1)=sigma2

c++++++++++++ precision parameter

              thetasave(isave,pce+2)=alpha

c++++++++++++ tailfree regression coefficients

              mm=0
              do i=2,ntlr
                 do j=1,ptf
                    mm=mm+1 
                    randsave(isave,mm)=betatf(i,j)
                 end do
              end do

c++++++++++++ density

              call densldtfpaft(maxm,ngrid,ntlr,ntprob,npredden,
     &                          pce,ptf,betace,
     &                          betatf,grid,xcepredm,xtfpredm,
     &                          sigma2,densl,densm,k)

              do i=1,npredden
                 write(1) (densl(i,j),j=1,ngrid)
              end do   

c++++++++++++ survival function

              do i=1,npredden
                 do j=1,ngrid
                    tmp1=log(grid(j))
                    call cdfldtfp(i,tmp1,maxm,ntlr,npredden,
     &                            pce,ptf,betace,betatf,
     &                            sigma2,xtfpred,xcepred,
     &                            tmp3,k,prob,probc)
                    survmm(i,j)=survmm(i,j)+1.d0-tmp3
                    survml(i,j)=1.d0-tmp3
                 end do
                 write(2) (survml(i,j),j=1,ngrid) 
              end do   

c++++++++++++ cpo

              do i=1,nrec

c+++++++++++++++ check if the user has requested an interrupt
                 call rchkusr()

                 if(interind(i).eq.1)then
                    tmp1=log(y(i,2))

                    call cdfldtfp(i,tmp1,maxm,ntlr,nrec,pce,
     &                            ptf,betace,betatf,
     &                            sigma2,xtf,xce,
     &                            tmp3,k,prob,probc)
                    cpo(i,1)=cpo(i,1)+1.d0/tmp3
                    cpo(i,2)=cpo(i,2)+tmp3  
                 end if

                 if(interind(i).eq.2)then
                    tmp1=log(y(i,1))
                    call cdfldtfp(i,tmp1,maxm,ntlr,nrec,pce,
     &                            ptf,betace,betatf,
     &                            sigma2,xtf,xce,
     &                            tmp2,k,prob,probc)

                    tmp1=log(y(i,2))
                    call cdfldtfp(i,tmp1,maxm,ntlr,nrec,pce,
     &                            ptf,betace,betatf,
     &                            sigma2,xtf,xce,
     &                            tmp3,k,prob,probc)
                    cpo(i,1)=cpo(i,1)+1.d0/(tmp3-tmp2)
                    cpo(i,2)=cpo(i,2)+(tmp3-tmp2)  
                 end if

                 if(interind(i).eq.3)then
                    tmp1=log(y(i,1))

                    call cdfldtfp(i,tmp1,maxm,ntlr,nrec,pce,
     &                            ptf,betace,betatf,
     &                            sigma2,xtf,xce,
     &                            tmp2,k,prob,probc)
                    cpo(i,1)=cpo(i,1)+1.d0/(1.d0-tmp2)
                    cpo(i,2)=cpo(i,2)+(1.d0-tmp2)  
                 end if

 
                 if(interind(i).eq.4)then

                    call ldensldtfp(z(i),i,maxm,nrec,ntlr,ntprob,
     &                              pce,betace,xce,
     &                              ptf,betatf,xtf,
     &                              sigma2,tmp1,k)
                    tmp2=exp(tmp1-z(i)) 

                    cpo(i,1)=cpo(i,1)+1.d0/tmp2
                    cpo(i,2)=cpo(i,2)+tmp2  
                 end if
              end do

c++++++++++++ quantiles

              call quantileldtfp(maxm,ntlr,npredmed,pce,ptf,
     &                           betace,betatf,
     &                           sigma2,xtfpredm,xcepredm,
     &                           qmm,qml,k,prob,probc,quans)

              write(3) (exp(qml(j,1)),j=1,npredmed)
              write(4) (exp(qml(j,2)),j=1,npredmed)
              write(5) (exp(qml(j,3)),j=1,npredmed)

c++++++++++++ print
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

      do i=1,npredden
         do j=1,ngrid
            densm(i,j)=densm(i,j)/dble(nsave)
            survmm(i,j)=survmm(i,j)/dble(nsave)
         end do
      end do

      do i=1,npredmed
         do j=1,3
            qmm(i,j)=qmm(i,j)/dble(nsave)
         end do
      end do

      do i=1,nrec
         cpo(i,1)=dble(nsave)/cpo(i,1)
         cpo(i,2)=cpo(i,2)/dble(nsave)                                    
      end do

      close(unit=1)
      close(unit=2)
      close(unit=3)
      close(unit=4)
      close(unit=5)

      if(cband.eq.1)then

      call hpddensregdtfp(nsave,npredden,ngrid,worksam,
     &                    worksam2,fs,densl,densu,tband)

      call hpddensregdtfps(nsave,npredden,ngrid,worksam,
     &                     worksam2,fs,survml,survmu,tband)

      call hpdldtfpq1(nsave,npredmed,worksam,worksam3,
     &                qml,qmu,tband)

      call hpdldtfpq2(nsave,npredmed,worksam,worksam3,
     &                qml,qmu,tband)


      call hpdldtfpq3(nsave,npredmed,worksam,worksam3,
     &                qml,qmu,tband)

      end if
 
      return
      end


