
c=======================================================================                  
      subroutine pbinary(link,nrec,p,sens,spec,x,yobs,
     &                   betapm,betapv, 
     &                   mcmc,nsave,propv,
     &                   acrate,thetasave,cpo,
     &                   beta,
     &                   betac,eta,etan,
     &                   iflag,
     &                   seed1,seed2,
     &                   workm1,workm2,
     &                   workmh1,workv1,workv2)
c=======================================================================
c
c     Version 1.0: 
c     Last modification: 01-07-2006.
c
c     Subroutine `pbinary' to run a Markov chain in a  
c     parametric binary regression model. 
c
c     Copyright: Alejandro Jara, 2006-2010.
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
c        link        :  integer giving the link to be considered,
c                       (1) logit, (2) probit, (3) cloglog, and
c                       (4) cauchy.
c        nrec        :  integer giving the number of observations.
c        p           :  integer giving the number of fixed coefficients.
c        sens        :  real vector of sensitivity, sens(nrec).
c        spec        :  real vector of specificity, spec(nrec).
c        x           :  real matrix giving the design matrix for the 
c                       fixed effects, x(nrec,p). 
c        yobs        :  integer vector giving the oberved binary 
c                       response, yobs(nrec).
c
c---- Prior information ------------------------------------------------
c
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
c        thetasave   :  real vector containing the mcmc sample for the
c                       regression parameters, betsave(nsave,p). 
c        cpo         :  real giving the cpo, cpo(nrec).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        beta        :  real vector giving the current value of the 
c                       regression coefficients, beta(p).
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        betac       :  real vector giving the current value of the 
c                       candidate for regression parameters, betac(p).
c        cdfcauchy   :  cdf of a cauchy distribution.
c        cdfnorm     :  cdf of a normal distribution.
c        dispcount   :  index. 
c        eta         :  real vector giving the linear predictor, 
c                       eta(nrec).
c        etan        :  real vector giving the linear predictor, 
c                       etan(nrec).
c        i           :  index. 
c        iflag       :  integer vector used to evaluate the prior
c                       distribution for the regression coefficients, 
c                       iflag(p).
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        logliko     :  real working variable.
c        loglikn     :  real working variable.
c        logprioro   :  real working variable.
c        logpriorn   :  real working variable.
c        nscan       :  index.
c        ok          :  integer indicator.
c        ratio       :  real working variable.
c        runif       :  real uniform random number generator.
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        skipcount   :  index. 
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
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
c
c=======================================================================                  

      implicit none

c+++++Constants
      double precision zero,one
      parameter(zero=0.d0)
      parameter(one =1.d0)

c+++++Observed variables
      integer link,nrec,p,yobs(nrec)
      double precision sens(nrec),spec(nrec)
      double precision x(nrec,p)

c+++++Prior information
      double precision betapm(p),betapv(p,p)

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay
      double precision propv(p,p)

c+++++Stored output
      double precision acrate
      double precision thetasave(nsave,p)
      double precision cpo(nrec)

c+++++Current values of the parameters
      double precision beta(p)
  
c+++++Working space
      integer dispcount
      integer i
      integer iflag(p)
      integer isave,iscan
      integer j
      integer nscan
      integer ok
      integer sprint  
      integer seed1,seed2
      integer skipcount
      
      double precision betac(p)
      double precision cdfcauchy,cdfnorm
      double precision eta(nrec),etan(nrec)
      double precision logliko,loglikn,logprioro,logpriorn
      double precision ratio
      double precision tmp1,tmp2
      double precision workm1(p,p),workm2(p,p),workmh1(p*(p+1)/2)
      double precision workv1(p),workv2(p)

      real runif
      
c+++++CPU time
      double precision sec00,sec0,sec1,sec
      
c++++ initialize variables

      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      
      tmp1=0.d0
      tmp2=0.d0

c++++ set random number generator

      call setall(seed1,seed2)


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
      
      call cpu_time(sec0)
      sec00=0.d0

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


c+++++++ evaluate log-likelihood for current and candidate value 
c+++++++ of parameters

        
         logliko=0.d0
         loglikn=0.d0
 
         do i=1,nrec

            tmp1=0.d0
            do j=1,p
               tmp1=tmp1+x(i,j)*betac(j)                   
            end do
            etan(i)=tmp1
         
            if(link.eq.1)then
               tmp1=exp(eta(i))/(1.d0+exp(eta(i)))
               tmp2=exp(etan(i))/(1.d0+exp(etan(i)))
            end if

            if(link.eq.2)then
               tmp1=cdfnorm(eta(i),0.d0,1.d0,1,0)
               tmp2=cdfnorm(etan(i),0.d0,1.d0,1,0)
            end if
            
            if(link.eq.3)then
               tmp1=1.d0-exp(-exp(eta(i)))
               tmp2=1.d0-exp(-exp(etan(i)))
            end if

            if(link.eq.4)then
               tmp1=cdfcauchy(eta(i),0.d0,1.d0,1,0)
               tmp2=cdfcauchy(etan(i),0.d0,1.d0,1,0)
            end if
            
            if(yobs(i).eq.1)then
               tmp1=tmp1*sens(i)+(1.d0-spec(i))*(1.d0-tmp1) 
               if(tmp1.lt.zero)go to 100
               if(tmp1.gt.one )go to 100
               tmp2=tmp2*sens(i)+(1.d0-spec(i))*(1.d0-tmp2) 
               if(tmp2.lt.zero)go to 100
               if(tmp2.gt.one )go to 100
               logliko=logliko+log(tmp1)
               loglikn=loglikn+log(tmp2)
            end if   
            if(yobs(i).eq.0)then
               tmp1=tmp1*sens(i)+(1.d0-spec(i))*(1.d0-tmp1) 
               if(tmp1.lt.zero)go to 100
               if(tmp1.gt.one )go to 100
               tmp2=tmp2*sens(i)+(1.d0-spec(i))*(1.d0-tmp2) 
               if(tmp2.lt.zero)go to 100
               if(tmp2.gt.one )go to 100
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

100      continue

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

c+++++++++++++ cpo, errors and predictive information
               do j=1,nrec

                  if(link.eq.1)then
                    tmp1=exp(eta(j))/(1.d0+exp(eta(j)))
                  end if

                  if(link.eq.2)then
                    tmp1=cdfnorm(eta(j),0.d0,1.d0,1,0)
                  end if
            
                  if(link.eq.3)then
                    tmp1=1.d0-exp(-exp(eta(j)))
                  end if

                  if(link.eq.4)then
                    tmp1=cdfcauchy(eta(j),0.d0,1.d0,1,0)
                  end if
         
                  tmp1=sens(j)*tmp1+(1.d0-spec(j))*(1.d0-tmp1)

                  if(yobs(j).eq.1)then
                    cpo(j)=cpo(j)+1.0d0/tmp1 
                  end if   
                  if(yobs(j).eq.0)then
                    tmp1=1.0d0-tmp1 
                    cpo(j)=cpo(j)+1.0d0/tmp1 
                  end if            
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
      
     
      acrate=acrate/dble(nscan)      
      
      do i=1,nrec
         cpo(i)=dble(nsave)/cpo(i)
      end do
      
      return
      end

