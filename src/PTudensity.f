
c=======================================================================                      
      subroutine ptdensityu(ngrid,nrec,y,
     &                      ab,murand,sigmarand,jfr,m0,s0,tau,
     &                      mcmcvec,nsave,tune1,tune2,tune3,
     &                      acrate,f,thetasave,cpo,
     &                      cpar,mu,sigma,
     &                      grid,seed,whicho,whichn)
c=======================================================================                      
c
c     Subroutine `ptdensityu' to run a Markov chain for univariate 
c     density estimation using a Mixture of Polya Tree prior. The
c     Polya Tree is centered in a N(mu,sigma2) distribution.
c
c     Copyright: Alejandro Jara, 2006-2010.
c
c     Version 3.0: 
c
c     Last modification: 25-01-2009.
c     
c     Changes and Bug fixes: 
c
c     Version 2.0 to Version 3.0:
c          - Centering parameters can be fixed.
c          - Proper prior can be used on the centering parameters.
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
c        ngrid       :  integer giving the size of the grid where
c                       the density estimate is evaluated.
c        nrec        :  integer giving the number of observations.
c        y           :  real vector giving the response variables,
c                       y(nrec).
c
c-----------------------------------------------------------------------
c
c---- Prior information ------------------------------------------------
c 
c        ca, cb      :  real giving the hyperparameters of the prior
c                       distribution for the precision parameter,
c                       c ~ Gamma(ca,cb). If ca<0 the precision 
c                       parameter is considered as a constant.
c        jfr         :  integer vector indicating whether Jeffery's
c                       prior is used for the centering parameters.
c        m0          :  real giving the mean of the normal prior
c                       for the centering mean.
c        s0          :  real giving the variance of the normal prior
c                       for the centering mean.
c        tau         :  real vector giving the hyperparameters of
c                       inverse gamma prior for the centering variance.
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
c        tune1       :  real giving the tuning parameter for MH of
c                       mean baseline.
c        tune2       :  real giving the tuning parameter for MH of
c                       variance baseline.
c        tune3       :  real giving the tuning parameter for MH of
c                       precision parameter.
c        
c-----------------------------------------------------------------------
c
c---- Output ----------------------------------------------------------- 
c
c        acrate      :  real vector giving the MH acceptance rate,
c                       acrate(3). 
c        f           :  real vector giving the density estimate at the
c                       grid, f(ngrid).
c        thetasave   :  real vector containing the mcmc sample for the
c                       parameters, betsave(nsave,3). 
c        cpo         :  real giving the cpo, cpo(nrec).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        cpar        :  real giving the current value of the precision
c                       parameter of the Polya Tree.
c        mu          :  real giving the current value of the 
c                       baseline mean.
c        sigma       :  real giving the he current value of the
c                       baseline standard deviation.
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        cdfnorm     :  cdf of a normal distribution.
c        cparc       :  real giving the value of the candidate
c                       for the precision parameter.
c        dispcount   :  index. 
c        dlnrm       :  density of a log-normal distribution.
c        dnrm        :  density of a normal distribution.
c        grid        :  real vector giving the grid where the density
c                       estimate is evaluated, grid(ngrid) .
c        i           :  index.
c        invcdfnorm  :  quantile function for a normal distribution.
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        je2         :  index.
c        k           :  index.   
c        k1          :  index.
c        k2          :  index.
c        logcgkn     :  real working variable.
c        logcgko     :  real working variable.
c        loglikn     :  real working variable.
c        logliko     :  real working variable.
c        logpriorn   :  real working variable.
c        logprioro   :  real working variable.
c        muc         :  real giving the value of the candidate
c                       for the baseline mean. 
c        nscan       :  index.
c        parti       :  index.
c        pprn        :  index.
c        quan        :  real working variable.
c        ratio       :  real working variable.
c        rbeta       :  real beta random number generator.
c        rtlnorm     :  real truncated log normal random number generator.
c        rnorm       :  real normal random number generator.
c        runif       :  real uniform random number generator.
c        s           :  real giving the sample variance. 
c        sec         :  cpu time working variable.
c        sec0        :  cpu time working variable.
c        sec00       :  cpu time working variable.
c        sec1        :  cpu time working variable.
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        sigmac      :  real giving the value of the candidate
c                       for the baseline standard deviation. 
c        sigma2      :  real giving the current value of
c                       the baseline variance. 
c        sigma2c     :  real giving the value of the candidate
c                       for the baseline variance. 
c        skipcount   :  index. 
c        sprint      :  integer function to print on screen.
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
c        whicho      :  integer vector giving the observation in each
c                       partition, whicho(nrec).
c        whichn      :  integer vector giving the observation in each
c                       partition, whichn(nrec).
c        ybar        :  real giving the sample mean.
c
c=======================================================================
      implicit none 

c+++++Data
      integer ngrid,nrec
      double precision y(nrec)

c+++++Prior information
      double precision ab(2),ca,cb
      integer murand,sigmarand,jfr(2)
      double precision m0,s0
      double precision tau(2)

c+++++MCMC parameters
      integer mcmcvec(3),nburn,nskip,nsave,ndisplay
      double precision tune1,tune2,tune3

c+++++Stored output
      double precision acrate(3),f(ngrid)
      double precision thetasave(nsave,3)
      double precision cpo(nrec)

c+++++Current values of the parameters
      double precision cpar,mu,sigma

c+++++Working space - CPU time
      double precision sec00,sec0,sec1,sec

c+++++Working space - Density
      double precision grid(ngrid)
      
c+++++Working space - Distributions
      double precision dnrm,dlnrm
      double precision invcdfnorm

c+++++Working space - General
      integer countero,countern
      integer i,j,je2,k,k1,k2,l
      integer nint,ok,parti,pprn,sprint
      integer whicho(nrec),whichn(nrec)
      double precision prob
      double precision quan
      double precision s,tmp1,tmp2,ybar

c+++++Working space - MCMC scans
      integer dispcount,isave,iscan,nscan,skipcount

c+++++Working space - MH steps
      double precision cparc
      double precision logcgkn,logcgko
      double precision loglikn,logliko
      double precision logpriorn,logprioro
      double precision muc,sigma2,sigmac,sigma2c
      double precision ratio

c+++++Working space - Random number generator
      integer seed(2),seed1,seed2
      double precision rnorm,rtlnorm
      real runif

c++++ initialize variables
      nburn=mcmcvec(1)
      nskip=mcmcvec(2)
      ndisplay=mcmcvec(3)
      
      ca=ab(1)
      cb=ab(2)
      
      sigma2=sigma**2
     
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
      
      call cpu_time(sec0)
      sec00=0.d0
  
      ybar=0.d0
      s=0.d0
      do i=1,nrec
         ybar=ybar+y(i) 
      end do
      ybar=ybar/dble(nrec)
      
      do i=1,nrec
         s=s+(y(i)-ybar)**2 
      end do
      s=s/dble(nrec)


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ First computation of loglikelihood (to reduce CPU time)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      logliko=0.d0

      do i=1,nrec

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++ first observation
         if(i.eq.1)then

            logliko=dnrm(y(1),mu, sigma, 1)

c+++++++ following observations
          else

            nint=2
            prob=1.d0/dble(nint)
            quan=invcdfnorm(prob,mu,sigma,1,0)

            countero=0
            
            if(y(i).le.quan) then
               parti=1 
               do l=1,i-1
                  if(y(l).le.quan)then
                     countero=countero+1
                     whicho(countero)=l
                  end if   
               end do
             else
               parti=2
               do l=1,i-1
                  if(y(l).gt.quan)then
                     countero=countero+1
                     whicho(countero)=l
                  end if   
               end do
            end if  

            logliko=logliko+
     &              log(2.d0*cpar+dble(2*countero))-
     &              log(2.d0*cpar+dble(i-1))

            if(countero.eq.0) go to 1

            ok=1
            j=2
            do while(ok.eq.1)
               nint=2**j
               je2=j**2
               prob=1.d0/dble(nint)
               
               k1=2*(parti-1)+1
               k2=2*(parti-1)+2
               quan=invcdfnorm(dble(k1)*prob,mu,sigma,1,0)
               
               if(y(i).le.quan)then
                 parti=k1
                 k=k1
                else
                 parti=k2
                 k=k2
               end if  

               countern=0
               
               if(k.eq.1)then
                  do l=1,countero
                     if(y(whicho(l)).le.quan)then
                        countern=countern+1
                        whichn(countern)=whicho(l)
                     end if   
                  end do
                else if(k.eq.nint)then
                  quan=invcdfnorm(dble(k-1)*prob,mu,sigma,1,0) 
                  do l=1,countero
                     if(y(whicho(l)).gt.quan)then
                        countern=countern+1
                        whichn(countern)=whicho(l)
                     end if   
                  end do
                else
                  tmp1=invcdfnorm(dble(k-1)*prob,mu,sigma,1,0)
                  tmp2=invcdfnorm(dble(k  )*prob,mu,sigma,1,0)

                  if(tmp1.ge.tmp2)then
                    call rexit("Error in the limits")
                  end if  
                
                  do l=1,countero
                     if(y(whicho(l)).gt.tmp1.and.
     &                  y(whicho(l)).le.tmp2)then
                        countern=countern+1
                        whichn(countern)=whicho(l)
                     end if   
                  end do
               end if
               
               logliko=logliko+
     &                  log(2.d0*cpar*dble(je2)+dble(2*countern))-
     &                  log(2.d0*cpar*dble(je2)+dble(  countero))

               if(countern.eq.0)then
                  ok=0
                else  
                  countero=countern
                  do l=1,countern
                     whicho(l)=whichn(l)
                  end do
                  j=j+1
               end if   
            end do

1           continue
            
            logliko=logliko+dnrm(y(i),mu,sigma,1)
            
         end if
      end do


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ Scanning the posterior distribution
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      do iscan=1,nscan
  
c+++++++ check if the user has requested an interrupt
         call rchkusr()
 
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Updating mu using a MH step                  +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
         if(murand.eq.1)then
            muc=rnorm(mu,tune1*sigma/sqrt(dble(nrec)))

            loglikn=0.d0

            do i=1,nrec

c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

c+++++++++++++ first observation
               if(i.eq.1)then

                  loglikn=dnrm(y(1),muc,sigma ,1)

c+++++++++++++ following observations
                 else

                  nint=2
                  prob=1.d0/dble(nint)
                  quan=invcdfnorm(prob,muc,sigma,1,0)

                  countero=0
                 
                  if(y(i).le.quan) then
                     parti=1 
                     do l=1,i-1
                        if(y(l).le.quan)then
                           countero=countero+1
                           whicho(countero)=l
                        end if   
                     end do
                   else
                     parti=2 
                     do l=1,i-1
                        if(y(l).gt.quan)then
                           countero=countero+1
                           whicho(countero)=l
                        end if   
                     end do
                  end if  

                  loglikn=loglikn+
     &                    log(2.d0*cpar+dble(2*countero))-
     &                    log(2.d0*cpar+dble(i-1))

                  if(countero.eq.0) go to 2

                  ok=1
                  j=2
                  do while(ok.eq.1)
                     nint=2**j
                     je2=j**2
                     prob=1.d0/dble(nint)

                     k1=2*(parti-1)+1
                     k2=2*(parti-1)+2
                     quan=invcdfnorm(dble(k1)*prob,muc,sigma,1,0)
               
                     if(y(i).le.quan)then
                       parti=k1
                       k=k1
                      else
                       parti=k2
                       k=k2
                     end if                      
                    
                     countern=0
                    
                     if(k.eq.1)then
                        do l=1,countero
                           if(y(whicho(l)).le.quan)then
                              countern=countern+1
                              whichn(countern)=whicho(l)
                           end if   
                        end do
                      else if(k.eq.nint)then
                        quan=invcdfnorm(dble(k-1)*prob,muc,sigma,1,0) 
                        do l=1,countero
                           if(y(whicho(l)).gt.quan)then
                              countern=countern+1
                              whichn(countern)=whicho(l)
                           end if   
                        end do
                      else
                        tmp1=invcdfnorm(dble(k-1)*prob,muc,sigma,1,0)
                        tmp2=invcdfnorm(dble(k  )*prob,muc,sigma,1,0)

                        if(tmp1.ge.tmp2)then
                          call rexit("Error in the limits")
                        end if  
                     
                        do l=1,countero
                           if(y(whicho(l)).gt.tmp1.and.
     &                        y(whicho(l)).le.tmp2)then
                              countern=countern+1
                              whichn(countern)=whicho(l)
                           end if   
                        end do
                     end if
                    
                     loglikn=loglikn+
     &                       log(2.d0*cpar*dble(je2)+dble(2*countern))-
     &                       log(2.d0*cpar*dble(je2)+dble(  countero))

                     if(countern.eq.0)then
                        ok=0
                      else  
                        countero=countern
                        do l=1,countern
                           whicho(l)=whichn(l)
                        end do
                        j=j+1
                     end if   
                  end do

2                 continue
                 
                  loglikn=loglikn+dnrm(y(i),muc,sigma,1)

              end if                 
            end do 

c++++++++++ acceptance step

            logpriorn=0.d0
            logprioro=0.d0

            if(jfr(1).eq.0)then
               logprioro=dnrm(mu, m0,sqrt(s0),1)
               logpriorn=dnrm(muc,m0,sqrt(s0),1)
            end if

            ratio=loglikn-logliko+logpriorn-logprioro

            if(log(dble(runif())).lt.ratio)then
               mu=muc
               logliko=loglikn
               acrate(1)=acrate(1)+1.d0
            end if

        end if

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Updating sigma using a MH step               +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(sigmarand.eq.1)then

            loglikn=0.d0

            sigma2c=rtlnorm(log(sigma2),tune2*0.1d0,0,0,.true.,.true.)
            sigmac=sqrt(sigma2c)

            logcgkn=dlnrm(sigma2 ,log(sigma2c),tune2*0.1d0,1) 
            logcgko=dlnrm(sigma2c,log(sigma2 ),tune2*0.1d0,1) 


            do i=1,nrec
         
c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

c+++++++++++++ first observation
               if(i.eq.1)then

                  loglikn=dnrm(y(1),mu,sigmac ,1)

c+++++++++++++ following observations
                else

                  nint=2
                  prob=1.d0/dble(nint)
                  quan=invcdfnorm(prob,mu,sigmac,1,0)

                  countero=0
                 
                  if(y(i).le.quan) then
                     parti=1
                     do l=1,i-1
                        if(y(l).le.quan)then
                           countero=countero+1
                           whicho(countero)=l
                        end if   
                     end do
                   else 
                     parti=2
                     do l=1,i-1
                        if(y(l).gt.quan)then
                           countero=countero+1
                           whicho(countero)=l
                        end if   
                     end do
                  end if  

                  loglikn=loglikn+
     &                    log(2.d0*cpar+dble(2*countero))-
     &                    log(2.d0*cpar+dble(i-1))

                  if(countero.eq.0) go to 3

                  ok=1
                  j=2
                  do while(ok.eq.1)
                     nint=2**j
                     je2=j**2
                     prob=1.d0/dble(nint)

                     k1=2*(parti-1)+1
                     k2=2*(parti-1)+2
                     quan=invcdfnorm(dble(k1)*prob,mu,sigmac,1,0)
               
                     if(y(i).le.quan)then
                       parti=k1
                       k=k1
                      else
                       parti=k2
                       k=k2
                     end if  
                    
                     countern=0
                    
                     if(k.eq.1)then
                        do l=1,countero
                           if(y(whicho(l)).le.quan)then
                              countern=countern+1
                              whichn(countern)=whicho(l)
                           end if   
                        end do
                      else if(k.eq.nint)then
                        quan=invcdfnorm(dble(k-1)*prob,mu,sigmac,1,0) 
                        do l=1,countero
                           if(y(whicho(l)).gt.quan)then
                              countern=countern+1
                              whichn(countern)=whicho(l)
                           end if   
                        end do
                      else
                        tmp1=invcdfnorm(dble(k-1)*prob,mu,sigmac,1,0)
                        tmp2=invcdfnorm(dble(k  )*prob,mu,sigmac,1,0)

                        if(tmp1.ge.tmp2)then
                           call rexit("Error in the limits")
                        end if  
                     
                        do l=1,countero
                           if(y(whicho(l)).gt.tmp1.and.
     &                        y(whicho(l)).le.tmp2)then
                              countern=countern+1
                              whichn(countern)=whicho(l)
                           end if   
                        end do
                     end if
                    
                     loglikn=loglikn+
     &                       log(2.d0*cpar*dble(je2)+dble(2*countern))-
     &                       log(2.d0*cpar*dble(je2)+dble(  countero))

                     if(countern.eq.0)then
                        ok=0
                      else  
                        countero=countern
                        do l=1,countern
                           whicho(l)=whichn(l)
                        end do
                        j=j+1
                     end if   
                  end do
 
3                 continue
                 
                  loglikn=loglikn+dnrm(y(i),mu,sigmac,1)

               end if
            end do 

c++++++++++ acceptance step

            logpriorn=-log(sigma2c)
            logprioro=-log(sigma2)

            if(jfr(2).eq.0)then
               logpriorn=-(0.5d0*tau(1)+1.d0)*log(sigma2c)-
     &                     0.5d0*tau(2)/sigma2c

               logprioro=-(0.5d0*tau(1)+1.d0)*log(sigma2)-
     &                     0.5d0*tau(2)/sigma2
            end if

            ratio=loglikn-logliko+logcgkn-logcgko+
     &            logpriorn-logprioro 


            if(log(dble(runif())).lt.ratio)then
               sigma=sigmac
               sigma2=sigma2c
               logliko=loglikn
               acrate(2)=acrate(2)+1.d0
            end if

         end if

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ MH to update the c parameter                 +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if(ca.gt.0.d0)then

c++++++++++ sample candidates

            cparc=rtlnorm(log(cpar),tune3*0.2,0,0,.true.,.true.)
            logcgkn=dlnrm(cpar ,log(cparc),tune3*0.2,1) 
            logcgko=dlnrm(cparc,log(cpar ),tune3*0.2,1) 

c++++++++++ evaluate log-prior for candidate value of the parameters

            call dgamma2(cparc,ca,cb,logpriorn)  

c++++++++++ evaluate log-prior for current value of parameters

            call dgamma2(cpar ,ca,cb,logprioro)

c++++++++++ evaluate log-likelihood

            do i=1,nrec

c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

c++++++++++ first observation
            if(i.eq.1)then

                 loglikn=dnrm(y(1),mu,sigma ,1)

c++++++++++ following observations
              else

                 nint=2
                 prob=1.d0/dble(nint)
                 quan=invcdfnorm(prob,mu,sigma,1,0)

                 countero=0
                 
                 if(y(i).le.quan) then
                    parti=1 
                    do l=1,i-1
                       if(y(l).le.quan)then
                          countero=countero+1
                          whicho(countero)=l
                       end if   
                    end do
                  else
                    parti=2
                    do l=1,i-1
                       if(y(l).gt.quan)then
                          countero=countero+1
                          whicho(countero)=l
                       end if   
                    end do
                 end if  

                 loglikn=loglikn+
     &                   log(2.d0*cparc+dble(2*countero))-
     &                   log(2.d0*cparc+dble(i-1))

                 if(countero.eq.0) go to 4

                 ok=1
                 j=2
                 do while(ok.eq.1)
                    nint=2**j
                    je2=j**2
                    prob=1.d0/dble(nint)

                    k1=2*(parti-1)+1
                    k2=2*(parti-1)+2
                    quan=invcdfnorm(dble(k1)*prob,mu,sigma,1,0)
               
                    if(y(i).le.quan)then
                      parti=k1
                      k=k1
                     else
                      parti=k2
                      k=k2
                    end if  
      
                    countern=0
                    
                    if(k.eq.1)then
                       do l=1,countero
                          if(y(whicho(l)).le.quan)then
                             countern=countern+1
                             whichn(countern)=whicho(l)
                          end if   
                       end do
                     else if(k.eq.nint)then
                       quan=invcdfnorm(dble(k-1)*prob,mu,sigma,1,0) 
                       do l=1,countero
                          if(y(whicho(l)).gt.quan)then
                             countern=countern+1
                             whichn(countern)=whicho(l)
                          end if   
                       end do
                     else
                       tmp1=invcdfnorm(dble(k-1)*prob,mu,sigma,1,0)
                       tmp2=invcdfnorm(dble(k  )*prob,mu,sigma,1,0)

                       if(tmp1.ge.tmp2)then
                         call rexit("Error in the limits")
                       end if  
                     
                       do l=1,countero
                          if(y(whicho(l)).gt.tmp1.and.
     &                       y(whicho(l)).le.tmp2)then
                             countern=countern+1
                             whichn(countern)=whicho(l)
                          end if   
                       end do
                    end if
                    
                    loglikn=loglikn+
     &                      log(2.d0*cparc*dble(je2)+dble(2*countern))-
     &                      log(2.d0*cparc*dble(je2)+dble(  countero))

                    if(countern.eq.0)then
                       ok=0
                     else  
                       countero=countern
                       do l=1,countern
                          whicho(l)=whichn(l)
                       end do
                       j=j+1
                    end if   
                 end do

4                continue
                 
                 loglikn=loglikn+dnrm(y(i),mu,sigma,1)

            end if
            end do 

c++++++++++ acceptance step
            ratio=loglikn+logpriorn-logliko-logprioro+
     &            logcgkn-logcgko

            if(log(dble(runif())).lt.ratio)then
               cpar=cparc
               acrate(3)=acrate(3)+1.d0
               logliko=loglikn
            end if            
            
         end if 

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Save samples                                 +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1
              
c+++++++++++++ mu
               thetasave(isave,1)=mu

c+++++++++++++ sigma
               thetasave(isave,2)=sigma

c+++++++++++++ c parameter
               thetasave(isave,3)=cpar

c+++++++++++++ cpo

               do i=1,nrec

c++++++++++++++++ check if the user has requested an interrupt
                  call rchkusr()

                  loglikn=0.d0 
               
                  nint=2
                  tmp1=1.d0/dble(nint)
                  quan=invcdfnorm(tmp1,mu,sigma,1,0)
                  
                  countero=0
                  
                  if(y(i).le.quan)then
                      parti=1
                      do l=1,nrec
                         if(y(l).le.quan.and.l.ne.i)then
                            countero=countero+1
                            whicho(countero)=l
                         end if   
                      end do
                    else 
                      parti=2
                      do l=1,nrec
                         if(y(l).gt.quan.and.l.ne.i)then
                            countero=countero+1
                            whicho(countero)=l
                         end if   
                      end do
                  end if  
        
                  loglikn=loglikn+
     &                 log(2.d0*cpar+dble(2*countero))-
     &                 log(2.d0*cpar+dble(nrec-1))

                  if(countero.eq.0) go to 5

                  ok=1
                  j=2
                  do while(ok.eq.1)
                     nint=2**j
                     je2=j**2
                     prob=1.d0/dble(nint)

                     k1=2*(parti-1)+1
                     k2=2*(parti-1)+2
                     quan=invcdfnorm(dble(k1)*prob,mu,sigma,1,0)
               
                     if(y(i).le.quan)then
                       parti=k1
                       k=k1
                      else
                       parti=k2
                       k=k2
                     end if  

                     countern=0
                    
                     if(k.eq.1)then
                        do l=1,countero
                           if(y(whicho(l)).le.quan.and.
     &                        whicho(l).ne.i)then
                              countern=countern+1
                              whichn(countern)=whicho(l)
                           end if   
                        end do
                      else if(k.eq.nint)then
                        quan=invcdfnorm(dble(k-1)*prob,mu,sigma,1,0) 
                        do l=1,countero
                           if(y(whicho(l)).gt.quan.and.
     &                        whicho(l).ne.i)then
                              countern=countern+1
                              whichn(countern)=whicho(l)
                           end if   
                        end do
                      else
                        tmp1=invcdfnorm(dble(k-1)*prob,mu,sigma,1,0)
                        tmp2=invcdfnorm(dble(k  )*prob,mu,sigma,1,0)

                        if(tmp1.ge.tmp2)then
                          call rexit("Error in the limits")
                        end if  
                     
                        do l=1,countero
                           if(whicho(l).ne.i)then
                           if(y(whicho(l)).gt.tmp1.and.
     &                        y(whicho(l)).le.tmp2)then
                              countern=countern+1
                              whichn(countern)=whicho(l)
                           end if
                           end if
                        end do
                     end if
                    
                     loglikn=loglikn+
     &                      log(2.d0*cpar*dble(je2)+dble(2*countern))-
     &                      log(2.d0*cpar*dble(je2)+dble(  countero))

                     if(countern.eq.0)then
                        ok=0
                      else  
                        countero=countern
                        do l=1,countern
                           whicho(l)=whichn(l)
                        end do
                        j=j+1
                     end if   
                  end do

5                 continue
                 
                  loglikn=loglikn+dnrm(grid(i),mu,sigma,1)

                  cpo(i)=cpo(i)+1.d0/exp(loglikn)
               end do

c+++++++++++++ density 

               do i=1,ngrid
               
c++++++++++++++++ check if the user has requested an interrupt
                  call rchkusr()

                  loglikn=0.d0 
               
                  nint=2
                  tmp1=1.d0/dble(nint)
                  quan=invcdfnorm(tmp1,mu,sigma,1,0)
                  
                  countero=0
                  
                  if(grid(i).le.quan)then
                      parti=1
                      do l=1,nrec
                         if(y(l).le.quan)then
                            countero=countero+1
                            whicho(countero)=l
                         end if   
                      end do
                    else
                      parti=2
                      do l=1,nrec
                         if(y(l).gt.quan)then
                            countero=countero+1
                            whicho(countero)=l
                         end if   
                      end do
                  end if  
        
                  loglikn=loglikn+
     &                 log(2.d0*cpar+dble(2*countero))-
     &                 log(2.d0*cpar+dble(  nrec))

                  if(countero.eq.0) go to 6

                  ok=1
                  j=2
                  do while(ok.eq.1)
                     nint=2**j
                     je2=j**2
                     prob=1.d0/dble(nint)

                     k1=2*(parti-1)+1
                     k2=2*(parti-1)+2
                     quan=invcdfnorm(dble(k1)*prob,mu,sigma,1,0)
               
                     if(grid(i).le.quan)then
                       parti=k1
                       k=k1
                      else
                       parti=k2
                       k=k2
                     end if  

                     countern=0
                    
                     if(k.eq.1)then
                        do l=1,countero
                           if(y(whicho(l)).le.quan)then
                              countern=countern+1
                              whichn(countern)=whicho(l)
                           end if   
                        end do
                      else if(k.eq.nint)then
                        quan=invcdfnorm(dble(k-1)*prob,mu,sigma,1,0) 
                        do l=1,countero
                           if(y(whicho(l)).gt.quan)then
                              countern=countern+1
                              whichn(countern)=whicho(l)
                           end if   
                        end do
                      else
                        tmp1=invcdfnorm(dble(k-1)*prob,mu,sigma,1,0)
                        tmp2=invcdfnorm(dble(k  )*prob,mu,sigma,1,0)

                        if(tmp1.ge.tmp2)then
                          call rexit("Error in the limits")
                        end if  
                     
                        do l=1,countero
                           if(y(whicho(l)).gt.tmp1.and.
     &                        y(whicho(l)).le.tmp2)then
                              countern=countern+1
                              whichn(countern)=whicho(l)
                           end if   
                        end do
                     end if
                    
                     loglikn=loglikn+
     &                      log(2.d0*cpar*dble(je2)+dble(2*countern))-
     &                      log(2.d0*cpar*dble(je2)+dble(  countero))

                     if(countern.eq.0)then
                        ok=0
                      else  
                        countero=countern
                        do l=1,countern
                           whicho(l)=whichn(l)
                        end do
                        j=j+1
                     end if   
                  end do

6                 continue
                 
                  loglikn=loglikn+dnrm(grid(i),mu,sigma,1)
                 
                  f(i)=f(i)+exp(loglikn)

               end do

c+++++++++++++ print
               skipcount = 0
               if(dispcount.ge.ndisplay)then
                  call cpu_time(sec1)
                  sec00=sec00+(sec1-sec0)
                  sec=sec00
                  sec0=sec1
                  pprn=sprint(isave,nsave,sec)
                  dispcount=0
               end if   
            end if         
         end if
      end do

      do i=1,3       
         acrate(i)=acrate(i)/dble(nscan)      
      end do   
     
      do i=1,nrec
         cpo(i)=dble(nsave)/cpo(i)
      end do

      do i=1,ngrid
         f(i)=f(i)/dble(nsave)
      end do


      return
      end

