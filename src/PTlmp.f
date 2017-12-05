c=======================================================================                      
      subroutine ptlmp(maxm,mdzero,ngrid,nrec,p,x,y,                    
     &                 a0b0,betapm,betapv,tau,m0,s0,                    
     &                 mcmc,nsave,propv,mcmcad,                         
     &                 seed,                                            
     &                 acrate,randsave,thetasave,cpo,f,                 
     &                 alpha,beta,mu,sigma2,v,                          
     &                 betac,iflag,vc,workm1,workm2,                    
     &                 workmh1,workv1,workv2,grid,whicho,whichn,xtx)    
c=======================================================================                      
c     # arguments = 40.
c
c     Subroutine `ptlmp' to run a Markov chain in the  
c     semiparametric linear regression model using a partially 
c     specified Mixture of Polya trees. 
c
c     Copyright: Alejandro Jara, 2006-2010.
c
c     Version 3.0: 
c
c     Last modification: 20-08-2007.
c
c     Changes and Bug fixes: 
c
c     Version 2.0 to Version 3.0:
c          - Allows for non median zero regression models.
c
c     Version 1.0 to Version 2.0:
c          - Uses vectors to keep the observations in each partition.
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
c        maxm        :  integer giving the number of binary partitions
c                       in the Finite Polya tree prior.
c        mdzero      :  integer indicating whether a median zero is 
c                       fitted or not.
c        ngrid       :  integer giving the size of the grid where
c                       the density estimate is evaluated.
c        nrec        :  integer giving the number of observations.
c        p           :  integer giving the number of fixed coefficients.
c        x           :  real matrix giving the design matrix for the 
c                       fixed effects, x(nrec,p). 
c        y           :  real vector giving the limits of the responses,
c                       y(nrec).
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
c        tau1, tau2  :  reals giving the hyperparameters of the prior 
c                       distribution for the inverse of the 
c                       variance of the normal baseline distribution,
c                       1/sigma ~ Gamma(tau1/2,tau2/2).
c        m0,s0       :  
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
c        mcmcad      :  real vector keeping elements for adaptive MH,
c                       mcmcad(13).
c        
c-----------------------------------------------------------------------
c
c---- Random numbers ---------------------------------------------------
c
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        
c-----------------------------------------------------------------------
c
c---- Output ----------------------------------------------------------- 
c
c        acrate      :  real giving the MH acceptance rate, acrate(4). 
c        cpo         :  real giving the cpo, cpo(nrec).
c        f           :  real vector giving the density estimate at the
c                       grid, f(ngrid).
c        randsave    :  real matrix containing the mcmc samples for
c                       the errors and prediction,
c                       randsave(nsave,nrec+1).
c        thetasave   :  real vector containing the mcmc sample for the
c                       regression parameters, betsave(nsave,p+2). 
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        alpha       :  real giving the current value of the precision
c                       parameter of the Polya tree process.
c        beta        :  real vector giving the current value of the 
c                       regression coefficients, beta(p).
c        mu          :  real giving the mean of the baseline 
c                       distribution.
c        sigma2      :  real giving the current value of the variance
c                       of the baseline distribution.
c        v           :  real vector giving the current value of the 
c                       errors, v(nrec).
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        alphac      :  real giving the candidate value of the precision
c                       parameter.
c        betac       :  real vector giving the candidate value of the 
c                       candidate for regression parameters, betac(p).
c        denom       :  real working variable.
c        dispcount   :  index. 
c        dlnrm       :  density of a lognormal normal distribution.
c        dnrm        :  density of a normal distribution.
c        grid        :  real vector giving the grid where the density
c                       estimate is evaluated, grid(ngrid).
c        i           :  index.
c        iflag       :  integer vector used to evaluate the prior
c                       distribution for the regression coefficients, 
c                       iflag(p).
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        logcgko     :  real working variable.
c        logcgkn     :  real working variable.
c        loglikec    :  real working variable.
c        loglikeo    :  real working variable.
c        logpriorc   :  real working variable.
c        logprioro   :  real working variable.
c        nscan       :  index.
c        ratio       :  real working variable.
c        rgamma      :  real gamma random number generator.
c        runif       :  real uniform random number generator.
c        sd          :  real working variable.
c        sdc         :  real working variable.
c        sisgma2c    :  real working variable.
c        sprint      :  integer function to print on screen.
c        skipcount   :  index. 
c        tmp1        :  real used to accumulate quantities. 
c        vc          :  real vector giving the candidate value of the 
c                       errors, vc(nrec).
c        vpred       :  real working variable.
c        whicho      :  integer vector giving the observation in each
c                       partitio, whicho(nrec).
c        whichn      :  integer vector giving the observation in each
c                       partitio, whichn(nrec).
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
c        xtx         :  real matrix used to update the fixed effects,
c                       xtx(p,p).
c
c=======================================================================                  
     
      implicit none 
      
c+++++Partially specified Polya Trees parameter      
      integer maxm,mdzero
     
c+++++Observed variables
      integer ngrid,nrec,p  
      double precision x(nrec,p),y(nrec)

c+++++Prior information
      double precision a0b0(2),aa0,ab0,betapm(p),betapv(p,p)
      double precision tau(2),tau1,tau2
      double precision m0,s0
      
c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay
      double precision propv(p,p)
      double precision mcmcad(18)

c+++++Random numbers
      integer seed(2),seed1,seed2

c+++++Stored output
      double precision acrate(4),randsave(nsave,nrec+1),
     1  thetasave(nsave,p+3)
      double precision cpo(nrec),f(ngrid)
      
c+++++Current values of the parameters
      double precision alpha,beta(p),mu,sigma2,v(nrec)
      
c+++++External Working space
      integer iflag(p)
      integer whicho(nrec),whichn(nrec)
      double precision betac(p)
      double precision grid(ngrid)
      double precision vc(nrec)
      double precision workm1(p,p),workm2(p,p),workmh1(p*(p+1)/2)
      double precision workv1(p),workv2(p)
      double precision xtx(p,p)

c+++++Internal Working space
      integer dispcount
      integer i
      integer iscan
      integer isave
      integer j
      integer nscan
      integer sprint
      integer skipcount
      double precision alphac
      double precision dlnrm
      double precision dnrm
      double precision logcgkn,logcgko
      double precision loglikec,loglikeo
      double precision logpriorc,logprioro
      double precision muc
      double precision ratio
      double precision rnorm
      double precision sd,sdc
      double precision theta,thetac
      double precision tmp1
      double precision vpred
      
      real runif

c+++++Adaptive MH
      integer acceptb
      double precision logval
      double precision maxa,maxb
      double precision nscanp

c+++++Adaptive MH for beta
      integer skipb
      double precision arateb
      double precision counterb
      double precision tune1

c+++++Adaptive MH for mu
      integer skipm
      double precision aratemu(5),counterm
      double precision pilogestmu(2)
      double precision tune2(2)

c+++++Adaptive MH for sigma2
      integer skips
      double precision aratesig(5),counters
      double precision pilogestsig(2)
      double precision tune3(2)

c+++++Adaptive MH for alpha
      integer skipa
      double precision aratea(5),countera
      double precision pilogesta(2)
      double precision tune4(2)

c+++++CPU time
      double precision sec00,sec0,sec1,sec
      
c++++ initialize variables

      aa0=a0b0(1)
      ab0=a0b0(2)

      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)

      tau1=tau(1)
      tau2=tau(2)
      
      sd=sqrt(sigma2)

c++++ set random number generator
      seed1=seed(1)
      seed2=seed(2)
      
      call setall(seed1,seed2)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*(nsave)

      maxa=100.d0
      maxb=100.d0

      arateb=0.d0
      skipb=0
      tune1=mcmcad(1)
      counterb=mcmcad(2)

      aratemu(1)=0.d0
      aratemu(2)=0.d0
      aratemu(3)=0.d0
      aratemu(4)=0.d0
      aratemu(5)=0.d0
      skipm=0
      tune2(1)=mcmcad(3)
      tune2(2)=mcmcad(4)
      counterm=mcmcad(5)
      pilogestmu(1)=mcmcad(6)
      pilogestmu(2)=mcmcad(7)
     
      aratesig(1)=0.d0
      aratesig(2)=0.d0
      aratesig(3)=0.d0
      aratesig(4)=0.d0
      aratesig(5)=0.d0      
      skips=0
      tune3(1)=mcmcad(8)
      tune3(2)=mcmcad(9)
      counters=mcmcad(10)
      pilogestsig(1)=mcmcad(11)
      pilogestsig(2)=mcmcad(12)

      aratea(1)=0.d0
      aratea(2)=0.d0
      aratea(3)=0.d0
      aratea(4)=0.d0
      aratea(5)=0.d0      
      skipa=0
      tune4(1)=mcmcad(13)
      tune4(2)=mcmcad(14)
      countera=mcmcad(15)
      pilogesta(1)=mcmcad(16)
      pilogesta(2)=mcmcad(17)

      nscanp=mcmcad(18)

      call cpu_time(sec0)
      sec00=0.d0

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ First computation of loglikelihood (to reduce CPU time)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      loglikeo=0.d0

      do i=1,nrec
         call rchkusr()   
         
         tmp1=0.d0
         do j=1,p
            tmp1=tmp1+x(i,j)*beta(j)
         end do
         v(i)=y(i)-tmp1
      end do
      
      call loglik_unippt(nrec,mdzero,maxm,alpha,mu,sigma2,v,
     &                   whicho,whichn,loglikeo)
      
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ Scanning the posterior distribution
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ MH to update the regression coefficients                   +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ sample candidates
         do i=1,p
            do j=1,p
               xtx(i,j)=tune1*propv(i,j)
            end do
         end do

         call rmvnorm(p,beta,xtx,workmh1,workv1,betac)
         
c+++++++ evaluate log-prior for candidate value of the parameters

         call dmvn(p,betac,betapm,betapv,logpriorc,workv1,workm1,
     &             workm2,workv2,iflag)  

c+++++++ evaluate log-prior for current value of parameters

         call dmvn(p,beta,betapm,betapv,logprioro,workv1,workm1,
     &             workm2,workv2,iflag)  

c+++++++ evaluate log-likelihood

         loglikec=0.d0

         do i=1,nrec
            call rchkusr()   
            tmp1=0.d0
            do j=1,p
               tmp1=tmp1+x(i,j)*betac(j)
            end do
            vc(i)=y(i)-tmp1
         end do

         call loglik_unippt(nrec,mdzero,maxm,alpha,mu,sigma2,vc,
     &                      whicho,whichn,loglikec)

c+++++++ acceptance step

         ratio=loglikec+logpriorc-loglikeo-logprioro

         if(log(dble(runif())).lt.ratio)then
            acrate(1)=acrate(1)+1.d0
            arateb=arateb+1.d0
            loglikeo=loglikec
            do i=1,p
               beta(i)=betac(i) 
            end do
            do i=1,nrec
               v(i)=vc(i)
            end do
         end if

c+++++++ do adapting

         skipb = skipb + 1
         if(skipb.eq.100)then
            counterb=counterb+1.d0
            arateb=arateb/dble(100)

            if(p.eq.1)then 
               if(arateb.gt.0.44d0)then
                  tune1=tune1+
     &                  min(0.01d0,1.d0/sqrt(counterb))
                 else
                  tune1=tune1-
     &                  min(0.01d0,1.d0/sqrt(counterb))
               end if
            end if 

            if(p.gt.1)then 
               if(arateb.gt.0.25d0)then
                  tune1=tune1+
     &                  min(0.01d0,1.d0/sqrt(counterb))
                 else
                  tune1=tune1-
     &                  min(0.01d0,1.d0/sqrt(counterb))
               end if
            end if 

c++++++++++ prevent from getting to extreme
            if(tune1.gt.maxa)then
               tune1=maxa
            end if   

            if(tune1.lt.0.0001d0)then
               tune1=0.0001d0
            end if   

            arateb=0.d0
            skipb=0 
         end if


c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ MH to update mu                                            +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(s0.gt.0.d0)then

c++++++++++ candidate
            logval=log(1.d0+dabs(mu))
            tmp1=exp(
     &         (tune2(1)/2.d0)+
     &         (tune2(2)/2.d0)*(logval-pilogestmu(1))
     &              )         
            muc=rnorm(mu,tmp1)

            logcgkn=dnrm(mu,muc,tmp1,1) 

            logval=log(1.d0+dabs(muc))
            tmp1=exp(
     &         (tune2(1)/2.d0)+
     &         (tune2(2)/2.d0)*(logval-pilogestmu(1))
     &              )         
            logcgko=dnrm(muc,mu,tmp1,1) 


c++++++++++ prior
            logpriorc=dnrm(muc, m0, sqrt(s0),1)
            logprioro=dnrm(mu , m0, sqrt(s0),1)

c++++++++++ evaluate log-likelihood

            loglikec=0.d0
 
            call loglik_unippt(nrec,mdzero,maxm,alpha,muc,sigma2,v,
     &                         whicho,whichn,loglikec)

c++++++++++ acceptance step

            ratio=loglikec+logpriorc-loglikeo-logprioro+
     &            logcgkn-logcgko

            acceptb=0

            if(log(dble(runif())).lt.ratio)then
               acceptb=1
            end if
            
            if(acceptb.eq.1)then
               acrate(2)=acrate(2)+1.d0
               aratemu(1)=aratemu(1)+1.d0
               mu=muc
               loglikeo=loglikec
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
                  tune2(1)=tune2(1)+
     &                     min(0.01d0,1.d0/sqrt(counterm))
                 else
                  tune2(1)=tune2(1)-
     &                     min(0.01d0,1.d0/sqrt(counterm))
               end if 

c+++++++++++++ adapt b.

               if(aratemu(3)*aratemu(4).lt.aratemu(5)*aratemu(2))
     &         then
                  tune2(2)=tune2(2)-
     &                     min(0.01d0,1.d0/sqrt(counterm))
                 else
                  tune2(2)=tune2(2)+
     &                     min(0.01d0,1.d0/sqrt(counterm))
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
               aratemu(1)=0.d0
               aratemu(2)=0.d0
               aratemu(3)=0.d0
               aratemu(4)=0.d0
               aratemu(5)=0.d0
               skipm=0
            end if    
         end if

c+++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Sampling log-sigma : Gamma prior
c+++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(tau1.gt.0.d0)then 

c++++++++++ candidate

            theta=log(sqrt(sigma2))  
            logval=log(1.d0+dabs(theta))
            tmp1=exp(
     &         (tune3(1)/2.d0)+
     &         (tune3(2)/2.d0)*(logval-pilogestsig(1))
     &              )         
            thetac=rnorm(theta,tmp1)
            logcgkn=dnrm(theta,thetac,tmp1,1)

            logval=log(1.d0+dabs(thetac))
            tmp1=exp(
     &         (tune3(1)/2.d0)+
     &         (tune3(2)/2.d0)*(logval-pilogestsig(1))
     &              )         
            logcgko=dnrm(thetac,theta,tmp1,1)

c++++++++++ likelihood
            sdc=exp(thetac)
 
            loglikec=0.d0

            call loglik_unippt(nrec,mdzero,maxm,alpha,mu,sdc**2,v,
     &                         whicho,whichn,loglikec)

c++++++++++ prior

            logpriorc=-tau1*thetac-tau2*exp(-2.d0*thetac)/2.d0
            logprioro=-tau1*theta-tau2*exp(-2.d0*theta)/2.d0

c++++++++++ acceptance step

            ratio=loglikec+logpriorc-loglikeo-logprioro+
     &            logcgkn-logcgko

            acceptb=0

            if(log(dble(runif())).lt.ratio)then
               acceptb=1
            end if
            
            if(acceptb.eq.1)then
               acrate(3)=acrate(3)+1.d0
               aratesig(1)=aratesig(1)+1.d0
               sd=sdc
               sigma2=sd**2
               loglikeo=loglikec
            end if

            theta=log(sqrt(sigma2))  
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
                  tune3(1)=tune3(1)+
     &                     min(0.01d0,1.d0/sqrt(counters))
                 else
                  tune3(1)=tune3(1)-
     &                     min(0.01d0,1.d0/sqrt(counters))
               end if 

c+++++++++++++ adapt b.

               if(aratesig(3)*aratesig(4).lt.aratesig(5)*aratesig(2))
     &         then
                  tune3(2)=tune3(2)-
     &                     min(0.01d0,1.d0/sqrt(counters))
                 else
                  tune3(2)=tune3(2)+
     &                     min(0.01d0,1.d0/sqrt(counters))
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
               aratesig(1)=0.d0
               aratesig(2)=0.d0
               aratesig(3)=0.d0
               aratesig(4)=0.d0
               aratesig(5)=0.d0               
               skips=0
            end if    
         end if

c++++++++++++++++++++++++++++++++++         
c+++++++ Precision parameter
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()
         
         if(aa0.gt.0.d0)then

c++++++++++ sample candidates

            theta=log(alpha)  
            logval=log(1.d0+dabs(theta))
            tmp1=exp(
     &         (tune4(1)/2.d0)+
     &         (tune4(2)/2.d0)*(logval-pilogesta(1))
     &              )         
            thetac=rnorm(theta,tmp1)
            alphac=exp(thetac) 

            logcgkn=dlnrm(alpha,thetac,tmp1,1)

            logval=log(1.d0+dabs(thetac))
            tmp1=exp(
     &         (tune4(1)/2.d0)+
     &         (tune4(2)/2.d0)*(logval-pilogesta(1))
     &              )         
            logcgko=dlnrm(alphac,theta,tmp1,1)

c++++++++++ evaluate log-prior for candidate value of the parameters

            call dgamma2(alphac,aa0,ab0,logpriorc)  

c++++++++++ evaluate log-prior for current value of parameters

            call dgamma2(alpha,aa0,ab0,logprioro)

c++++++++++ evaluate log-likelihood

            loglikec=0.d0

            call loglik_unippt(nrec,mdzero,maxm,alphac,mu,sigma2,v,
     &                         whicho,whichn,loglikec)

c++++++++++ acceptance step

            ratio=loglikec+logpriorc-loglikeo-logprioro+
     &            logcgkn-logcgko

            acceptb=0

            if(log(dble(runif())).lt.ratio)then
               acceptb=1
            end if
            
            if(acceptb.eq.1)then
               acrate(4)=acrate(4)+1.d0
               aratea(1)=aratea(1)+1.d0
               alpha=alphac
               loglikeo=loglikec
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
                  tune4(1)=tune4(1)+
     &                     min(0.01d0,1.d0/sqrt(countera))
                 else
                  tune4(1)=tune4(1)-
     &                     min(0.01d0,1.d0/sqrt(countera))
               end if 

c+++++++++++++ adapt b.

               if(aratea(3)*aratea(4).lt.aratea(5)*aratea(2))
     &         then
                  tune4(2)=tune4(2)-
     &                     min(0.01d0,1.d0/sqrt(countera))
                 else
                  tune4(2)=tune4(2)+
     &                     min(0.01d0,1.d0/sqrt(countera))
               end if 

c+++++++++++++ prevent a from getting to extreme
               if(tune4(1).gt.maxa)then
                  tune4(1)=maxa
               end if   

               if(tune4(1).lt.-maxa)then
                  tune4(1)=-maxa
               end if   

c+++++++++++++ prevent b from getting to extreme
               if(tune4(2).gt.maxb)then
                  tune4(2)=maxb
               end if   

               if(tune4(2).lt.-maxb)then
                  tune4(2)=-maxb
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

c+++++++++++++ precision parameter information
               thetasave(isave,p+3)=alpha

c+++++++++++++ errors and predictive information
               do j=1,nrec
                  randsave(isave,j)=v(j)
               end do

               call sampupptpred(maxm,mdzero,nrec,alpha,mu,
     &                           sigma2,v,
     &                           whicho,whichn,vpred)

               randsave(isave,nrec+1)=vpred
               
c+++++++++++++ cpo

               do i=1,nrec
                  loglikec=0.d0

                  vpred=v(i)
                  
                  call condupptprior(vpred,i,maxm,mdzero,nrec,alpha,
     &                               mu,sigma2,v,
     &                               whicho,whichn,loglikec)

                  cpo(i)=cpo(i)+1.d0/exp(loglikec)  
               end do

c+++++++++++++ density estimate

               do i=1,ngrid
               
                  loglikec=0.d0

                  call gridupptprior(grid(i),maxm,mdzero,nrec,
     &                               alpha,mu,sigma2,v,
     &                               whicho,whichn,loglikec)
                  f(i)=f(i)+exp(loglikec)  
               end do

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

      do i=1,4
         acrate(i)=acrate(i)/dble(nscan)        
      end do   
      
      do i=1,nrec
         cpo(i)=dble(nsave)/cpo(i)
      end do
      
      do i=1,ngrid
         f(i)=f(i)/dble(nsave)
      end do
      
      mcmcad(1)=tune1
      mcmcad(2)=counterb
      mcmcad(3)=tune2(1)
      mcmcad(4)=tune2(2)
      mcmcad(5)=counterm
      mcmcad(6)=pilogestmu(1)
      mcmcad(7)=pilogestmu(2)
      mcmcad(8)=tune3(1)
      mcmcad(9)=tune3(2)
      mcmcad(10)=counters 
      mcmcad(11)=pilogestsig(1)
      mcmcad(12)=pilogestsig(2)
      mcmcad(13)=tune4(1)
      mcmcad(14)=tune4(2)
      mcmcad(15)=countera 
      mcmcad(16)=pilogesta(1)
      mcmcad(17)=pilogesta(2)
      mcmcad(18)=dble(nscan)

      return
      end


