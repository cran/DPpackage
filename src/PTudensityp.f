
c=======================================================================                      
      subroutine ptdensityup(ngrid,nrec,y,
     &                       ab,murand,sigmarand,jfr,m0,s0,tau,
     &                       nlevel,ninter,
     &                       mcmcvec,nsave,tune1,tune2,tune3,
     &                       acrate,f,thetasave,cpo,
     &                       cpar,mu,sigma,
     &                       grid,intpn,intpo,
     &                       accums,assign,counter,endp,
     &                       prob,rvecs,
     &                       seed)
c=======================================================================                      
c
c     Subroutine `fptdensityu' to run a Markov chain for univariate 
c     density estimation using a Mixture of Finite Polya Tree prior. The
c     Polya Tree is centered in a N(mu,sigma2) distribution.
c
c     This subroutine is based on the mpt FORTRAN program of 
c     Tim Hanson.
c
c     Copyright: Alejandro Jara and Tim Hanson, 2006-2010.
c
c     Version 3.0: 
c
c     Last modification: 6-1-2010.
c
c     Changes and Bug fixes: 
c
c     Version 2.0 to Version 3.0:
c          - Centering parameters can be fixed.
c          - Proper prior can be used on the centering parameters.
c
c     Version 1.0 to Version 2.0:
c          - Fixed bug in computation of MH ratio for precision parameter.
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
c     The authors's contact information:
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
c      Tim Hanson
c      Division of Biostatistics
c      University of Minnesota
c      School of Public Health
c      A460 Mayo Building, 
c      MMC 303
c      420 Delaware St SE
c      Minneapolis, MN 55455
c      Voice: 612-626-7075  URL  : http://www.biostat.umn.edu/~hanson/
c      Fax  : 612-626-0660  Email: hanson@biostat.umn.edu
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
c        ninter      :  integer giving the number of final intervals 
c                       in the Finite Polya tree prior.
c        nlevel      :  integer giving the number of binary partitions
c                       in the Finite Polya tree prior.
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
c        mu          :  real vector giving the current value of the 
c                       baseline mean.
c        sigma       :  real giving the he current value of the
c                       baseline standard deviation.
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        accums      :  real matrix giving the accumulated log 
c                       probabilities used in the computation of
c                       final interval probabilities,
c                       accums(nlevel,ninter).
c        assign      :  integer matrix giving the possition of each
c                       observation in each partition,
c                       assign(nrec,nlevel).
c        cdfnorm     :  cdf of a normal distribution.
c        counter     :  integer matrix giving the number of subjects
c                       in each binary partition, counter(nlevel,ninter).
c        cparc       :  real giving the value of the candidate
c                       for the precision parameter.
c        dbet        :  density of a beta distribution.
c        dispcount   :  index. 
c        dlnrm       :  density of a log-normal distribution.
c        dnrm        :  density of a normal distribution.
c        endp        :  real vector giving the end of the intervals,
c                       endp(ninter-1).
c        grid        :  real vector giving the grid where the density
c                       estimate is evaluated, grid(ngrid) .
c        i           :  index.
c        intlp       :  index.
c        intpn       :  integer vector giving the interval possition
c                       for the current value of parameters for
c                       the observations, intpn(nrec).
c        intpo       :  integer vector giving the interval possition
c                       for the candidate value of parameters for
c                       the observations, intpo(nrec).
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
c        nint        :  index
c        npoints     :  index.
c        nscan       :  index.
c        pprn        :  index.
c        prob        :  real vector giving the probability of 
c                       the intervals, prob(ninter).
c        quan        :  real working variable.
c        ratio       :  real working variable.
c        rbeta       :  real beta random number generator.
c        rtlnorm     :  real truncated log normal random number generator.
c        rnorm       :  real normal random number generator.
c        runif       :  real uniform random number generator.
c        rvecs       :  real matrix giving the random vectors for the
c                       Polya tree,  rvecs(nlevel,ninter).
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
c        tmp3        :  real used to accumulate quantities.
c        uni         :  real working variable.
c        ybar        :  real giving the sample mean.
c
c=======================================================================

      implicit none 

c+++++Data
      integer ngrid,nrec
      double precision y(nrec)

c+++++Prior information
      integer nlevel,ninter
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
      double precision dbet,dnrm,dlnrm
      double precision invcdfnorm

c+++++Working space - General
      integer i,j,je2,k,k1,k2
      integer nint,npoints,pprn,sprint
      double precision quan
      double precision s,tmp1,tmp2,tmp3,ybar

c+++++Working space - MCMC scans
      integer dispcount,isave,iscan,nscan,skipcount

c+++++Working space - MH steps
      integer intlp
      integer intpn(nrec),intpo(nrec)
      double precision cparc
      double precision logcgkn,logcgko
      double precision loglikn,logliko
      double precision logpriorn,logprioro
      double precision muc,sigma2,sigmac,sigma2c
      double precision ratio

c+++++Working space - Polya tree parameters
      integer assign(nrec,nlevel)
      integer counter(nlevel,ninter)
      double precision accums(nlevel,ninter)
      double precision endp(ninter-1)
      double precision prob(ninter)
      double precision rvecs(nlevel,ninter)

c+++++Working space - Random number generator
      integer seed(2),seed1,seed2
      double precision rbeta,rnorm,rtlnorm
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
      
      npoints=ninter-1
  
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
  
      do iscan=1,nscan
  
c+++++++ check if the user has requested an interrupt
         call rchkusr()
 
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Generating from the full conditionals for Polya Tree        ++ 
c+++++++ parmeters                                                   ++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
         do i=1,nlevel
            nint=2**i
            do j=1,nint
               counter(i,j)=0
               accums(i,j)=0.d0
               rvecs(i,j)=0.d0
            end do
         end do   
         
         do i=1,nrec
            nint=2
            tmp1=1.d0/dble(nint)
            quan=invcdfnorm(tmp1,mu,sigma,1,0)
            if(y(i).le.quan)then
                assign(i,1)=1
                counter(1,1)=counter(1,1)+1
              else
                assign(i,1)=2
                counter(1,2)=counter(1,2)+1
            end if  
        
            do j=2,nlevel
               nint=2**j
               tmp1=1.d0/dble(nint)            
               k=assign(i,j-1)
               k1=2*(k-1)+1
               k2=2*(k-1)+2
               
               quan=invcdfnorm(dble(k1)*tmp1,mu,sigma,1,0) 
               
               if(y(i).le.quan)then
                  assign(i,j)=k1
                  counter(j,k1)=counter(j,k1)+1
                 else
                  assign(i,j)=k2
                  counter(j,k2)=counter(j,k2)+1
               end if  
            end do
            
            intpo(i)=assign(i,nlevel)
         end do
 
c+++++++ generating (X0,X1),(X00,X01,X10,X11),.... 
 
         tmp1=cpar+dble(counter(1,1))
         tmp2=cpar+dble(counter(1,2))
         tmp3=rbeta(tmp1,tmp2)
         rvecs(1,1)=tmp3
         rvecs(1,2)=1.d0-tmp3
         accums(1,1)=log(tmp3)
         accums(1,2)=log(1.d0-tmp3)

         do i=1,nlevel-1
            nint=2**i
            je2=(i+1)**2
            do j=1,nint
               k1=2*(j-1)+1
               k2=2*(j-1)+2            
               tmp1=cpar*dble(je2)+dble(counter(i+1,k1))
               tmp2=cpar*dble(je2)+dble(counter(i+1,k2))
               tmp3=rbeta(tmp1,tmp2)
               rvecs(i+1,k1)=tmp3
               rvecs(i+1,k2)=1.d0-tmp3
               accums(i+1,k1)=log(tmp3)+accums(i,j)
               accums(i+1,k2)=log(1.d0-tmp3)+accums(i,j)               
             end do
          end do   
 
          do i=1,ninter
             prob(i)=exp(accums(nlevel,i))
          end do
          
c++++++++ end points
  
          tmp1=1.d0/dble(ninter)  
          do i=1,npoints
             endp(i)=invcdfnorm(dble(i)*tmp1,mu,sigma,1,0)
          end do

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Updating mu using a MH step                  +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(murand.eq.1)then

            muc=rnorm(mu,tune1*sigma/sqrt(dble(nrec)))

            logliko=0.d0
            loglikn=0.d0

            do i=1,nrec
         
c+++++++++++++ possition according to the candidate

               nint=2
               tmp1=1.d0/dble(nint)
               quan=invcdfnorm(tmp1,muc,sigma,1,0)
               if(y(i).le.quan)then
                  intlp=1
                 else
                  intlp=2
               end if  
        
               do j=2,nlevel
                  nint=2**j
                  tmp1=1.d0/dble(nint)            
                  k=intlp
                  k1=2*(k-1)+1
                  k2=2*(k-1)+2
               
                  quan=invcdfnorm(dble(k1)*tmp1,muc,sigma,1,0) 
               
                  if(y(i).le.quan)then
                     intlp=k1
                    else
                     intlp=k2
                  end if  
               end do
            
               intpn(i)=intlp

c+++++++++++++ likelihood current

               tmp1=prob(intpo(i))*dble(ninter)*dnrm(y(i),mu,sigma,0)

               logliko=logliko+log(tmp1)

c+++++++++++++ likelihood candidate

               tmp2=prob(intpn(i))*dble(ninter)*dnrm(y(i),muc,sigma,0)
            
               loglikn=loglikn+log(tmp2)

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
               do i=1,nrec
                  intpo(i)=intpn(i)
               end do
               acrate(1)=acrate(1)+1.d0
            end if
          
         end if

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Updating sigma using a MH step               +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(sigmarand.eq.1)then

            logliko=0.d0
            loglikn=0.d0

            sigma2c=rtlnorm(log(sigma2),tune2*0.1d0,0,0,.true.,.true.)
            sigmac=sqrt(sigma2c)

            logcgkn=dlnrm(sigma2 ,log(sigma2c),tune2*0.1d0,1) 
            logcgko=dlnrm(sigma2c,log(sigma2 ),tune2*0.1d0,1) 


            do i=1,nrec
c+++++++++++++ possition according to the candidate

               nint=2
               tmp1=1.d0/dble(nint)
               quan=invcdfnorm(tmp1,mu,sigmac,1,0)
               if(y(i).le.quan)then
                  intlp=1
                 else
                  intlp=2
               end if  
        
               do j=2,nlevel
                  nint=2**j
                  tmp1=1.d0/dble(nint)            
                  k=intlp
                  k1=2*(k-1)+1
                  k2=2*(k-1)+2
               
                  quan=invcdfnorm(dble(k1)*tmp1,mu,sigmac,1,0) 
                 
                  if(y(i).le.quan)then
                     intlp=k1
                    else
                     intlp=k2
                  end if  
               end do
            
               intpn(i)=intlp

c+++++++++++++ likelihood current

               tmp1=prob(intpo(i))*dble(ninter)*dnrm(y(i),mu,sigma,0)

               logliko=logliko+log(tmp1)

c+++++++++++++ likelihood candidate

               tmp2=prob(intpn(i))*dble(ninter)*dnrm(y(i),mu,sigmac,0)
            
               loglikn=loglikn+log(tmp2)

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
               do i=1,nrec
                  intpo(i)=intpn(i)
               end do
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

c+++++++++++ evaluate log-likelihood

            logliko=0.d0
            loglikn=0.d0

            tmp1=cpar
            tmp2=cpar
            logliko=logliko+dbet(rvecs(1,1),tmp1,tmp2,1)

            tmp1=cparc
            tmp2=cparc
            loglikn=loglikn+dbet(rvecs(1,1),tmp1,tmp2,1)
            
            do i=1,nlevel-1
               nint=2**i
               je2=(i+1)**2
               do j=1,nint
                  k1=2*(j-1)+1
                  k2=2*(j-1)+2            
                  tmp1=cpar*dble(je2)
                  tmp2=cpar*dble(je2)
                  logliko=logliko+dbet(rvecs(i+1,k1),tmp1,tmp2,1)
                  
                  tmp1=cparc*dble(je2)
                  tmp2=cparc*dble(je2)
                  loglikn=loglikn+dbet(rvecs(i+1,k1),tmp1,tmp2,1)
                end do
             end do   

c++++++++++ acceptance step
            ratio=loglikn+logpriorn-logliko-logprioro+
     &            logcgkn-logcgko

            if(log(dble(runif())).lt.ratio)then
               cpar=cparc
               acrate(3)=acrate(3)+1.d0
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
                  tmp1=prob(intpo(i))*dble(ninter)*
     &                 dnrm(y(i),mu,sigma,0)
                  cpo(i)=cpo(i)+1.0d0/tmp1 
               end do

c+++++++++++++ density 

               do i=1,ngrid
                  nint=2
                  tmp1=1.d0/dble(nint)
                  quan=invcdfnorm(tmp1,mu,sigma,1,0)
                  if(grid(i).le.quan)then
                     intlp=1
                    else
                     intlp=2
                  end if  
        
                  do j=2,nlevel
                     nint=2**j
                     tmp1=1.d0/dble(nint)            
                     k=intlp
                     k1=2*(k-1)+1
                     k2=2*(k-1)+2
                     quan=invcdfnorm(dble(k1)*tmp1,mu,sigma,1,0) 
                     if(grid(i).le.quan)then
                        intlp=k1
                      else
                        intlp=k2
                     end if  
                  end do
                  f(i)=f(i)+prob(intlp)*dble(ninter)*
     &                 dnrm(grid(i),mu,sigma,0)
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

  
  
