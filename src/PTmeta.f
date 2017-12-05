c=======================================================================                      
      subroutine ptmeta(nrec,nfixed,p,                                  
     &                  x,y,sigma2e,                                    
     &                  a0b0,prec,sb,tau,m0,s0,                         
     &                  maxm,mdzero,                                    
     &                  mcmc,nsave,mcmcad,                              
     &                  acrate,randsave,thetasave,cpo,                  
     &                  alpha,b,                                         
     &                  beta,mu,sigma,                                  
     &                  mc,                                             
     &                  iflagp,workmp,workmhp,workvp,xty,               
     &                  seed,                                           
     &                  whicho,whichn,                                                         
     &                  betasave,bsave)                                 
c=======================================================================                      
c     # of arguments = 37.
c
c     Subroutine `ptmeta' to run a Markov chain in the semiparametric 
c     meta-analytic linear mixed model using a Polya tree prior 
c     for the distributions of the random effecs. 
c
c     Copyright: Alejandro Jara, 2007-2010.
c
c     Version 1.0: 
c
c     Last modification: 24-08-2007.
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
c        nrec        :  integer giving the number of observations.
c        nfixed      :  integer giving the number of fixed effects,
c                       if nfixed is 0 then p=1.
c        p           :  integer giving the number of fixed coefficients.
c        x           :  real matrix giving the design matrix for the 
c                       fixed effects, x(nrec,p). 
c        y           :  real vector giving the response variable,
c                       y(nrec).
c        sigma2e     :  real vector giving the value of the residual
c                       variances, sigma2e(nrec).
c-----------------------------------------------------------------------
c
c---- Prior information ------------------------------------------------
c 
c        a0, b0      :  real giving the hyperparameters of the prior
c                       distribution for the precision parameter,
c                       alpha ~ Gamma(a0,b0). If a0<0 the precision 
c                       parameter is considered as a constant.
c        prec        :  real matrix giving the prior precision matrix
c                       for the fixed effects, prec(p,p).
c        sb          :  real vector giving the product of the prior 
c                       precision and prior mean for the fixed effects,
c                       sb(p).
c        m0          :  real giving the prior mean for the baseline mean.
c        s0          :  real giving the prior variance for the baseline
c                       mean.
c        tau1, tau2  :  reals giving the hyperparameters of the prior 
c                       distribution for the inverse of the baseline
c                       variance, 1/sigma ~ Gamma(tau1/2,tau2/2).
c        maxm        :  integer giving the number of binary partitions
c                       in the Polya tree prior.
c        mdzero      :  integer indicating whether a median zero (1) 
c                       PT prior is used.
c-----------------------------------------------------------------------
c
c---- MCMC parameters --------------------------------------------------
c
c        nburn       :  integer giving the number of burn-in scans.
c        ndisplay    :  integer giving the number of saved scans to be
c                       displayed on screen.
c        nskip       :  integer giving the thinning interval.
c        nsave       :  integer giving the number of scans to be saved.
c        mcmcad      :  real vector used to save the parameters of
c                       the adaptive MH, mcmcad(16).
c-----------------------------------------------------------------------
c
c---- Output -----------------------------------------------------------
c
c        acrate      :  real giving the acceptance rates of MH steps, 
c                       acrate(4).
c        cpo         :  real giving the cpo and fso, cpo(nrec,2). 
c        randsave    :  real matrix containing the mcmc samples for
c                       the random effects and prediction,
c                       randsave(nsave,nrec+1).
c        thetasave   :  real matrix containing the mcmc samples for
c                       the averaged random effects, fixed effects, 
c                       mean and covariance of the baseline distribution, 
c                       thetsave(nsave,nfixed+4).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        alpha       :  real giving the current value of the precision
c                       parameter of the PT prior.
c        b           :  real vector giving the current value of the 
c                       random effects, b(nrec).
c        beta        :  real vector giving the current value of the 
c                       fixed effects, beta(p).
c        mu          :  real giving the mean of the normal 
c                       base line distribution for the random effects.
c        sigma       :  real giving the current value of the
c                       variance for normal base line 
c                       distribution for the random effects.
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        dispcount   :  index. 
c        i           :  index. 
c        ii          :  index. 
c        iflagp      :  integer vector used to invert the of the lhs
c                       least square solution for the fixed effects,
c                       iflagp(p).
c        isave       :  index. 
c        iscan       :  index. 
c        mc          :  real vector used to save model comparison 
c                       information, mc(5).
c        nscan       :  integer indicating the total number of MCMC
c                       scans.
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        skipcount   :  index. 
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
c        tmp3        :  real used to accumulate quantities.
c        whicho      :  integer vector giving the rand. eff. in each
c                       partition, whicho(nrec).
c        whichn      :  integer vector giving the rand. eff. in each
c                       partition, whichn(nrec).
c        workmp      :  real matrix used to update the fixed effects,
c                       workmp(p,p).
c        workmhp     :  real vector used to update the fixed effects,
c                       workmhp(p*(p+1)/2)
c        workvp      :  real vector used to update the fixed effects,
c                       workvp(p)
c        xty         :  real vector used to save the product 
c                       Xt(Y-Zb), xty(p).
c
c=======================================================================                  
      implicit none 

c+++++Data
      integer nrec,nfixed,p
      double precision y(nrec),x(nrec,p),sigma2e(nrec)

c+++++Prior 
      integer maxm,mdzero
      double precision aa0,ab0,a0b0(2),prec(p,p)
      double precision sb(p)
      double precision tau1,tau2,tau(2)
      double precision m0,s0
      
c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay
      double precision mcmcad(16)

c+++++Output
      double precision acrate(4)
      double precision cpo(nrec,2)
      double precision randsave(nsave,nrec+1)
      double precision thetasave(nsave,nfixed+4)

c+++++Current values of the parameters
      double precision alpha,beta(p),b(nrec)
      double precision mu,sigma

c+++++Seeds
      integer seed(2),seed1,seed2

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++External working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++fixed effects
      integer iflagp(p)
      double precision workmp(p,p)
      double precision workmhp(p*(p+1)/2)
      double precision workvp(p)
      double precision xty(p)

c+++++PT
      integer whicho(nrec),whichn(nrec)

c++++ model's performance
      double precision mc(5)
      double precision betasave(p),bsave(nrec)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Internal working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++General
      integer ii,i,j,k
      integer sprint 

      double precision acrate2 
      double precision alphac
      double precision betar      
      double precision logcgkn,logcgko
      double precision loglikn,logliko
      double precision logpriorn,logprioro
      double precision muc
      double precision ratio
      double precision sdc
      double precision sigmainv
      double precision theta,thetac
      double precision tmp1,tmp2
      double precision ztz,zty

c+++++MCMC
      integer dispcount,isave,iscan,nscan,skipcount 

c+++++RNG and distributions
      double precision dnrm,dlnrm,rnorm
      real runif

c++++ model's performance
      double precision dbarc,dbar,dhat,pd,lpml

c+++++CPU time
      double precision sec00,sec0,sec1,sec

c+++++Adaptive MH
      integer acceptb
      double precision logval
      double precision maxa,maxb
      double precision nscanp

c+++++Adaptive MH for mu
      integer skipm
      double precision aratemu(5),counterm
      double precision pilogestmu(2)
      double precision tune1(2)

c+++++Adaptive MH for sigma2
      integer skips
      double precision aratesig(5),counters
      double precision pilogestsig(2)
      double precision tune2(2)

c+++++Adaptive MH for alpha
      integer skipa
      double precision aratea(5),countera
      double precision pilogesta(2)
      double precision tune3(2)

c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      
      tau1=tau(1)
      tau2=tau(2)
      aa0=a0b0(1)
      ab0=a0b0(2)
      
c++++ set random number generator
      seed1=seed(1)
      seed2=seed(2)

      call setall(seed1,seed2)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      sigmainv=1.d0/sigma

      dbar=0.d0
      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*(nsave)

      maxa=100.d0
      maxb=100.d0

      aratemu(1)=0.d0
      aratemu(2)=0.d0
      aratemu(3)=0.d0
      aratemu(4)=0.d0
      aratemu(5)=0.d0
      skipm=0
      tune1(1)=mcmcad(1)
      tune1(2)=mcmcad(2)
      counterm=mcmcad(3)
      pilogestmu(1)=mcmcad(4)
      pilogestmu(2)=mcmcad(5)
     
      aratesig(1)=0.d0
      aratesig(2)=0.d0
      aratesig(3)=0.d0
      aratesig(4)=0.d0
      aratesig(5)=0.d0      
      skips=0
      tune2(1)=mcmcad(6)
      tune2(2)=mcmcad(7)
      counters=mcmcad(8)
      pilogestsig(1)=mcmcad(9)
      pilogestsig(2)=mcmcad(10)

      aratea(1)=0.d0
      aratea(2)=0.d0
      aratea(3)=0.d0
      aratea(4)=0.d0
      aratea(5)=0.d0      
      skipa=0
      tune3(1)=mcmcad(11)
      tune3(2)=mcmcad(12)
      countera=mcmcad(13)
      pilogesta(1)=mcmcad(14)
      pilogesta(2)=mcmcad(15)

      nscanp=mcmcad(16)

      call cpu_time(sec0)
      sec00=0.d0

      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++++
c+++++++ fixed effects
c++++++++++++++++++++++++++++++++++

         if(nfixed.gt.0)then
            do i=1,p
               xty(i)=sb(i)
               workvp(i)=0.d0
               do j=1,p
                  workmp(i,j)=prec(i,j)
               end do
            end do

            do i=1,nrec
               tmp1=y(i)-b(i) 

               do j=1,p
                  xty(j)=xty(j)+x(i,j)*(tmp1/sigma2e(i))
               end do
               
               do j=1,p
                  do k=1,p 
                     workmp(j,k)=workmp(j,k)+x(i,j)*x(i,k)/sigma2e(i)
                  end do   
               end do
            end do

            call inverse(workmp,p,iflagp) 

            do i=1,p
               tmp1=0.d0
               do j=1,p
                  tmp1=tmp1+workmp(i,j)*xty(j) 
               end do
               workvp(i)=tmp1
            end do

            call rmvnorm(p,workvp,workmp,workmhp,xty,beta)
         end if

c+++++++++++++++++++++++++++++++++
c+++++++ random effects        +++ 
c+++++++++++++++++++++++++++++++++

         acrate2=0.d0
         
         do ii=1,nrec

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

c++++++++++ generating a candidate
            ztz=sigmainv+1.d0/sigma2e(ii)
            zty=mu*sigmainv

            if(nfixed.eq.0)then
               tmp2=y(ii)
              else
               tmp1=0.d0
               do j=1,p 
                  tmp1=tmp1+x(ii,j)*beta(j)
               end do
               tmp2=y(ii)-tmp1
            end if
            zty=zty+tmp2/sigma2e(ii)
            ztz=1.d0/ztz
            tmp1=ztz*zty  

            thetac=rnorm(b(ii),sqrt(ztz))
            theta=b(ii)

c++++++++++ evaluating the likelihood
            loglikn=0.d0
            logliko=0.d0

            if(nfixed.eq.0)then
               tmp2=y(ii)
              else
               tmp1=0.d0
               do j=1,p 
                  tmp1=tmp1+x(ii,j)*beta(j)
               end do
               tmp2=y(ii)-tmp1
            end if
               
            logliko=dnrm(tmp2,theta,sqrt(sigma2e(ii)),1)
            loglikn=dnrm(tmp2,thetac,sqrt(sigma2e(ii)),1)

c++++++++++ evaluating the prior for the candidate

            call condupptprior(thetac,ii,maxm,mdzero,nrec,alpha,mu,
     &                         sigma,b,
     &                         whicho,whichn,logpriorn)

c++++++++++ evaluating the prior for the current value

            call condupptprior(theta,ii,maxm,mdzero,nrec,alpha,mu,
     &                         sigma,b,
     &                         whicho,whichn,logprioro)

c++++++++++ mh step
  
            ratio=loglikn-logliko+
     &            logpriorn-logprioro

            if(log(dble(runif())).lt.ratio)then
               acrate2=acrate2+1.d0
               b(ii)=thetac
            end if
         end do

         acrate(1)=acrate(1)+acrate2/dble(nrec)

         
c+++++++ update the log-likelihood for random effects

         call loglik_unippt(nrec,mdzero,maxm,alpha,mu,sigma,b,
     &                      whicho,whichn,logliko)

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ MH to update mu                                            +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(s0.gt.0.d0)then

c++++++++++ candidate
            logval=log(1.d0+dabs(mu))
            tmp1=exp(
     &         (tune1(1)/2.d0)+
     &         (tune1(2)/2.d0)*(logval-pilogestmu(1))
     &              )         
            muc=rnorm(mu,tmp1)

            logcgkn=dnrm(mu,muc,tmp1,1) 

            logval=log(1.d0+dabs(muc))
            tmp1=exp(
     &         (tune1(1)/2.d0)+
     &         (tune1(2)/2.d0)*(logval-pilogestmu(1))
     &              )         
            logcgko=dnrm(muc,mu,tmp1,1) 


c++++++++++ prior
            logpriorn=dnrm(muc, m0, sqrt(s0),1)
            logprioro=dnrm(mu , m0, sqrt(s0),1)

c++++++++++ evaluate log-likelihood
  
            call loglik_unippt(nrec,mdzero,maxm,alpha,muc,sigma,b,
     &                         whicho,whichn,loglikn)

c++++++++++ acceptance step

            ratio=loglikn+logpriorn-logliko-logprioro+
     &            logcgkn-logcgko

            acceptb=0

            if(log(dble(runif())).lt.ratio)then
               acceptb=1
            end if
            
            if(acceptb.eq.1)then
               acrate(2)=acrate(2)+1.d0
               aratemu(1)=aratemu(1)+1.d0
               mu=muc
               logliko=loglikn
            end if

            logval=log(1.d0+dabs(mu))
            pilogestmu(2)=pilogestmu(2)+logval
            
            if(logval.gt.pilogestmu(1))then
                aratemu(2)=aratemu(2)+1.d0
                if(acceptb.eq.1)then
                   aratemu(3)=aratemu(3)+1.d0
                end if
              else
                aratemu(4)=aratemu(4)+1.d0
                if(acceptb.eq.1)then
                   aratemu(5)=aratemu(5)+1.d0
                end if
            end if  

            pilogestmu(1)=pilogestmu(2)/(dble(iscan)+nscanp)

c++++++++++ do adapting

            skipm = skipm + 1
            if(skipm.eq.100)then
               counterm=counterm+1.d0
               aratemu(1)=aratemu(1)/dble(100)
               
c+++++++++++++ adapt a.
               if(aratemu(1).gt.0.44d0)then
                  tune1(1)=tune1(1)+
     &                     min(0.01d0,1.d0/sqrt(counterm))
                 else
                  tune1(1)=tune1(1)-
     &                     min(0.01d0,1.d0/sqrt(counterm))
               end if 

c+++++++++++++ adapt b.

               if(aratemu(3)*aratemu(4).lt.aratemu(5)*aratemu(2))
     &         then
                  tune1(2)=tune1(2)-
     &                     min(0.01d0,1.d0/sqrt(counterm))
                 else
                  tune1(2)=tune1(2)+
     &                     min(0.01d0,1.d0/sqrt(counterm))
               end if 

c+++++++++++++ prevent a from getting to extreme
               if(tune1(1).gt.maxa)then
                  tune1(1)=maxa
               end if   

               if(tune1(1).lt.-maxa)then
                  tune1(1)=-maxa
               end if   

c+++++++++++++ prevent b from getting to extreme
               if(tune1(2).gt.maxb)then
                  tune1(2)=maxb
               end if   

               if(tune1(2).lt.-maxb)then
                  tune1(2)=-maxb
               end if   

c+++++++++++++ set accetptance rate in batch to zero
               aratemu(1)=0.d0
               aratemu(2)=0.d0
               aratemu(3)=0.d0
               aratemu(4)=0.d0
               aratemu(5)=0.d0
               skipm=0
            end if    
         end if

c         call dblepr("mu",-1,mu,1)

c+++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Sampling log-sigma : Gamma prior
c+++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(tau1.gt.0.d0)then 

c++++++++++ candidate

            theta=log(sqrt(sigma))  
            logval=log(1.d0+dabs(theta))
            tmp1=exp(
     &         (tune2(1)/2.d0)+
     &         (tune2(2)/2.d0)*(logval-pilogestsig(1))
     &              )         
            thetac=rnorm(theta,tmp1)
            logcgkn=dnrm(theta,thetac,tmp1,1)

            logval=log(1.d0+dabs(thetac))
            tmp1=exp(
     &         (tune2(1)/2.d0)+
     &         (tune2(2)/2.d0)*(logval-pilogestsig(1))
     &              )         
            logcgko=dnrm(thetac,theta,tmp1,1)

c++++++++++ likelihood
            sdc=exp(thetac)

            call loglik_unippt(nrec,mdzero,maxm,alpha,mu,sdc**2,b,
     &                         whicho,whichn,loglikn)

c++++++++++ prior

            logpriorn=-tau1*thetac-tau2*exp(-2.d0*thetac)/2.d0
            logprioro=-tau1*theta-tau2*exp(-2.d0*theta)/2.d0

c++++++++++ acceptance step

            ratio=loglikn+logpriorn-logliko-logprioro+
     &            logcgkn-logcgko

            acceptb=0

            if(log(dble(runif())).lt.ratio)then
               acceptb=1
            end if
            
            if(acceptb.eq.1)then
               acrate(3)=acrate(3)+1.d0
               aratesig(1)=aratesig(1)+1.d0
               sigma=sdc**2
               sigmainv=1.d0/sigma
               logliko=loglikn
            end if

            theta=log(sqrt(sigma))  
            logval=log(1.d0+dabs(theta))
            pilogestsig(2)=pilogestsig(2)+logval
            
            if(logval.gt.pilogestsig(1))then
                aratesig(2)=aratesig(2)+1.d0
                if(acceptb.eq.1)then
                   aratesig(3)=aratesig(3)+1.d0
                end if
              else
                aratesig(4)=aratesig(4)+1.d0
                if(acceptb.eq.1)then
                   aratesig(5)=aratesig(5)+1.d0
                end if
            end if  

            pilogestsig(1)=pilogestsig(2)/(dble(iscan)+nscanp)

c++++++++++ do adapting

            skips = skips + 1
            if(skips.eq.100)then
               counters=counters+1.d0
               aratesig(1)=aratesig(1)/dble(100)
               
c+++++++++++++ adapt a.
               if(aratesig(1).gt.0.44d0)then
                  tune2(1)=tune2(1)+
     &                     min(0.01d0,1.d0/sqrt(counters))
                 else
                  tune2(1)=tune2(1)-
     &                     min(0.01d0,1.d0/sqrt(counters))
               end if 

c+++++++++++++ adapt b.

               if(aratesig(3)*aratesig(4).lt.aratesig(5)*aratesig(2))
     &         then
                  tune2(2)=tune2(2)-
     &                     min(0.01d0,1.d0/sqrt(counters))
                 else
                  tune2(2)=tune2(2)+
     &                     min(0.01d0,1.d0/sqrt(counters))
               end if 

c+++++++++++++ prevent a from getting to extreme
               if(tune2(1).gt.maxa)then
                  tune2(1)=maxa
               end if   

               if(tune2(1).lt.-maxa)then
                  tune2(1)=-maxa
               end if   

c+++++++++++++ prevent b from getting to extreme
               if(tune2(2).gt.maxb)then
                  tune2(2)=maxb
               end if   

               if(tune2(2).lt.-maxb)then
                  tune2(2)=-maxb
               end if   

c+++++++++++++ set accetptance rate in batch to zero
               aratesig(1)=0.d0
               aratesig(2)=0.d0
               aratesig(3)=0.d0
               aratesig(4)=0.d0
               aratesig(5)=0.d0               
               skips=0
            end if    
         end if

c         call dblepr("sigma",-1,sigma,1)

c++++++++++++++++++++++++++++++++++         
c+++++++ Precision parameter
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()
         
         if(aa0.gt.0.d0)then

c++++++++++ candidate

            theta=log(alpha)  
            logval=log(1.d0+dabs(theta))
            tmp1=exp(
     &         (tune3(1)/2.d0)+
     &         (tune3(2)/2.d0)*(logval-pilogesta(1))
     &              )         
            thetac=rnorm(theta,tmp1)
            alphac=exp(thetac) 

            logcgkn=dlnrm(alpha,thetac,tmp1,1)

            logval=log(1.d0+dabs(thetac))
            tmp1=exp(
     &         (tune3(1)/2.d0)+
     &         (tune3(2)/2.d0)*(logval-pilogesta(1))
     &              )         
            logcgko=dlnrm(alphac,theta,tmp1,1)

c++++++++++ evaluate log-prior for candidate value of the parameters

            call dgamma2(alphac,aa0,ab0,logpriorn)  

c++++++++++ evaluate log-prior for current value of parameters

            call dgamma2(alpha,aa0,ab0,logprioro)

c++++++++++ evaluate log-likelihood

            call loglik_unippt(nrec,mdzero,maxm,alphac,mu,sigma,b,
     &                        whicho,whichn,loglikn)

c++++++++++ acceptance step

            ratio=loglikn+logpriorn-logliko-logprioro+
     &            logcgkn-logcgko

            acceptb=0

            if(log(dble(runif())).lt.ratio)then
               acceptb=1
            end if
            
            if(acceptb.eq.1)then
               acrate(4)=acrate(4)+1.d0
               aratea(1)=aratea(1)+1.d0
               alpha=alphac
               logliko=loglikn
            end if

            theta=log(alpha)  
            logval=log(1.d0+dabs(theta))
            pilogesta(2)=pilogesta(2)+logval
            
            if(logval.gt.pilogesta(1))then
                aratea(2)=aratea(2)+1.d0
                if(acceptb.eq.1)then
                   aratea(3)=aratea(3)+1.d0
                end if
              else
                aratea(4)=aratea(4)+1.d0
                if(acceptb.eq.1)then
                   aratea(5)=aratea(5)+1.d0
                end if
            end if  

            pilogesta(1)=pilogesta(2)/(dble(iscan)+nscanp)

c++++++++++ do adapting

            skipa = skipa + 1
            if(skipa.eq.100)then
               countera=countera+1.d0
               aratea(1)=aratea(1)/dble(100)
               
c+++++++++++++ adapt a.
               if(aratea(1).gt.0.44d0)then
                  tune3(1)=tune3(1)+
     &                     min(0.01d0,1.d0/sqrt(countera))
                 else
                  tune3(1)=tune3(1)-
     &                     min(0.01d0,1.d0/sqrt(countera))
               end if 

c+++++++++++++ adapt b.

               if(aratea(3)*aratea(4).lt.aratea(5)*aratea(2))
     &         then
                  tune3(2)=tune3(2)-
     &                     min(0.01d0,1.d0/sqrt(countera))
                 else
                  tune3(2)=tune3(2)+
     &                     min(0.01d0,1.d0/sqrt(countera))
               end if 

c+++++++++++++ prevent a from getting to extreme
               if(tune3(1).gt.maxa)then
                  tune3(1)=maxa
               end if   

               if(tune3(1).lt.-maxa)then
                  tune3(1)=-maxa
               end if   

c+++++++++++++ prevent b from getting to extreme
               if(tune3(2).gt.maxb)then
                  tune3(2)=maxb
               end if   

               if(tune3(2).lt.-maxb)then
                  tune3(2)=-maxb
               end if   

c+++++++++++++ set accetptance rate in batch to zero
               aratea(1)=0.d0
               aratea(2)=0.d0
               aratea(3)=0.d0
               aratea(4)=0.d0
               aratea(5)=0.d0               
               skipa=0
            end if    

         end if 

c         call dblepr("alpha",-1,alpha,1)

c+++++++ save samples
         
         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1
c+++++++++++++ random effects

               do i=1,nrec
                  bsave(i)=bsave(i)+b(i)
                  randsave(isave,i)=b(i)
               end do

c+++++++++++++ predictive information

               call sampupptpred(maxm,mdzero,nrec,alpha,mu,
     &                           sigma,b,
     &                           whicho,whichn,theta)

               randsave(isave,nrec+1)=theta

c+++++++++++++ functional parameters

               call samplefuncupt(mdzero,maxm,nrec,alpha,b,mu,sigma,
     &                            theta,ztz)

               betar=mu+sqrt(sigma)*theta
               thetasave(isave,1)=betar

c+++++++++++++ regression coefficients

               if(nfixed.gt.0)then
                  do i=1,p
                     thetasave(isave,1+i)=beta(i)
                     betasave(i)=betasave(i)+beta(i)
                  end do
               end if   

c+++++++++++++ baseline mean

               thetasave(isave,1+nfixed+1)=mu

c+++++++++++++ baseline covariance

               thetasave(isave,1+nfixed+2)=sigma

c+++++++++++++ cluster information
               thetasave(isave,1+nfixed+3)=alpha

c+++++++++++++ cpo
               dbarc=0.d0
               do i=1,nrec
                  tmp1=0.d0
                  if(nfixed.gt.0)then
                     do j=1,p
                        tmp1=tmp1+x(i,j)*beta(j)
                     end do   
                  end if
                  tmp1=tmp1+b(i) 
                  tmp2=dnrm(y(i),tmp1,sqrt(sigma2e(i)),0)
                  cpo(i,1)=cpo(i,1)+1.0d0/tmp2  
                  cpo(i,2)=cpo(i,2)+tmp2                    
                  tmp2=dnrm(y(i),tmp1,sqrt(sigma2e(i)),1)
                  dbarc=dbarc+tmp2
               end do

c+++++++++++++ dic
               dbar=dbar-2.d0*dbarc
               
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
      
      do i=1,4
         acrate(i)=acrate(i)/dble(nscan)
      end do
      
      do i=1,nrec
         cpo(i,1)=dble(nsave)/cpo(i,1)
         cpo(i,2)=cpo(i,2)/dble(nsave)
      end do

      do i=1,p
         betasave(i)=betasave(i)/dble(nsave)
      end do

      do i=1,nrec
         bsave(i)=bsave(i)/dble(nsave)
      end do   

      dhat=0.d0
      lpml=0.d0
      do i=1,nrec
         tmp1=0.d0
         if(nfixed.gt.0)then
            do j=1,p
               tmp1=tmp1+x(i,j)*betasave(j)
            end do   
         end if
         tmp1=tmp1+bsave(i) 
         dhat=dhat+dnrm(y(i),tmp1,sqrt(sigma2e(i)),1)
         lpml=lpml+log(cpo(i,1))
      end do
      dhat=-2.d0*dhat

      dbar=dbar/dble(nsave)
      pd=dbar-dhat
      
      mc(1)=dbar
      mc(2)=dhat
      mc(3)=pd
      mc(4)=dbar+pd
      mc(5)=lpml

      mcmcad(1)=tune1(1)
      mcmcad(2)=tune1(2)
      mcmcad(3)=counterm
      mcmcad(4)=pilogestmu(1)
      mcmcad(5)=pilogestmu(2)
      mcmcad(6)=tune2(1)
      mcmcad(7)=tune2(2)
      mcmcad(8)=counters 
      mcmcad(9)=pilogestsig(1)
      mcmcad(10)=pilogestsig(2)

      mcmcad(11)=tune3(1)
      mcmcad(12)=tune3(2)
      mcmcad(13)=countera 
      mcmcad(14)=pilogesta(1)
      mcmcad(15)=pilogesta(2)

      mcmcad(16)=dble(nscan)
      
      return
      end
         

         
