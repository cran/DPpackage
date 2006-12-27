
c=======================================================================                      
      subroutine fptraschpoi(datastr,imiss,ngrid,nmissing,nsubject,p,y, #7
     &                       roffset,                                   #1
     &                       ninter,nlevel,                             #2
     &                       a0b0,b0,prec,sb,tau1,tau2,m,s,             #8
     &                       mcmc,nsave,tune3,tune4,tune5,              #5
     &                       acrate,cpo,f,faccum,randsave,thetasave,    #6
     &                       alpha,b,beta,mu,sigma,                     #5
     &                       accums,assignb,betac,counter,endp,         #5
     &                       iflag,intpn,intpo,prob,rvecs,seed,         #6
     &                       work1,work2,work3,                         #3 
     &                       workmh1,workv1,workv2,workv3,              #4
     &                       xtx,xty,grid)                              #3   
c=======================================================================                      
c
c     Subroutine `fptraschpoi' to run a Markov chain in the  
c     semiparametric Rasch Poisson Count model using a Polya tree prior
c     for the random effect distribution. 
c
c     Copyright: Alejandro Jara Vallejos, 2006-2007
c
c     Version 2.0: 
c     Last modification: 01-02-2007.
c
c     Changes and Bug fixes: 
c
c     Version 1.0 to Version 2.0:
c          - Correction in computation of MH ratio for random effects.
c          - Add offset.
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
c        datastr     :  integer matrix giving the position of missing 
c                       data points, datastr(nmissing,2)
c        ngrid       :  integer giving the size of the grid where
c                       the cdf estimate is evaluated.
c        imiss       :  integer indicating whether missing data are
c                       present (1) or absent (0). 
c        nmissing    :  integer giving the number of missing 
c                       observations.
c        nsubject    :  integer giving the number of subjects.
c        p           :  integer giving the number of items.
c        roffset     :  real matrix giving the offsets,
c                       roffset(nsubject,p).
c        y           :  integer matrix giving the response variable,
c                       y(nsubject,p).
c
c-----------------------------------------------------------------------
c
c---- Prior information ------------------------------------------------
c 
c        aa0, ab0    :  real giving the hyperparameters of the prior
c                       distribution for the precision parameter,
c                       alpha ~ Gamma(aa0,ab0). If aa0<0 the precision 
c                       parameter is considered as a constant.
c        b0          :  real vector giving the prior mean of difficulty
c                       parameters, b0(p-1).
c        prec        :  real matrix giving the prior precision matrix
c                       for the difficulty parameters, prec(p-1,p-1).
c        sb          :  real vector giving the product of the prior 
c                       precision and prior mean for the difficulty
c                       parameters, sb(p-1).
c        m           :  real vector giving the prior mean for the 
c                       baseline mean, m.
c        ninter      :  integer giving the number of final intervals 
c                       in the Finite Polya tree prior.
c        nlevel      :  integer giving the number of binary partitions
c                       in the Finite Polya tree prior.
c        s           :  real vector giving the prior variance for the 
c                       baseline mean, s.
c        tau1, tau2  :  reals giving the hyperparameters of the prior 
c                       distribution for the inverse of the random
c                       variance, 1/sigma2 ~ Gamma(tau1/2,tau2/2).
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
c        tune3       :  real giving the tuning parameter for MH of
c                       mean baseline.
c        tune4       :  real giving the tuning parameter for MH of
c                       variance baseline.
c        tune5       :  real giving the tuning parameter for MH of
c                       precision parameter.
c        
c-----------------------------------------------------------------------
c
c---- Output -----------------------------------------------------------
c
c        acrate      :  real vector giving the MH acceptance rate, 
c                       acrate(5). 
c        cpo         :  real matrix giving the cpo, cpo(nsubject,p).
c        f           :  real vector giving the density estimate at the
c                       grid, f(ngrid).
c        faccum      :  real vector giving the cdf estimate at the
c                       grid, faccum(ngrid).
c        randsave    :  real matrix containing the mcmc samples for
c                       the random effects and prediction,
c                       randsave(nsave,nsubject+1)
c        thetasave   :  real matrix containing the mcmc samples for
c                       the fixed effects,and mean and covariance of
c                       the baseline distribution, 
c                       thetsave(nsave,p+2).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        alpha       :  real giving the current value of the precision
c                       parameter of the Polya tree.
c        b           :  real vector giving the current value of the 
c                       random effects, b(nsubject).
c        beta        :  real vector giving the current value of the 
c                       difficulty parameters, beta(p-1).
c        mu          :  real giving the mean of the normal 
c                       base line distribution for the random effects,
c                       mu.
c        sigma       :  real giving the current value of the
c                       standard deviation for normal base line 
c                       distribution for the random effects,
c                       sigma.
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        accums      :  real matrix giving the accumulated log 
c                       probabilities used in the computation of
c                       final interval probabilities,
c                       accums(nlevel,ninter).
c        acrate2     :  real used to calculate the acceptance rate. 
c        alphac      :  real giving the current value of the candidate
c                       for the precision parameter.
c        assignb     :  integer matrix giving the possition of each
c                       observation in each partition,
c                       assignb(nrec,nlevel).
c        betac       :  real vector giving the current value of the 
c                       candidate for difficulty parameters, 
c                       betac(p-1).
c        counter     :  integer matrix giving the number of subjects
c                       in each binary partition, counter(nlevel,ninter).
c        dbet        :  density of a beta distribution.
c        dispcount   :  index. 
c        dlnrm       :  density of a log-normal distribution.
c        dnrm        :  density of a normal distribution.
c        endp        :  real vector giving the end of the intervals,
c                       endp(ninter-1).
c        grid        :  real vector giving the grid where the density
c                       estimate is evaluated, grid(ngrid) .
c        i           :  index. 
c        ii          :  index. 
c        imax        :  index.
c        imin        :  index.
c        iflag       :  integer vector used to invert the of the lhs
c                       least square solution for the difficulties,
c                       iflag(p-1).
c        intlp       :  index.
c        intpn       :  integer vector giving the interval possition
c                       for the current value of the random effects,
c                       intpn(nsubject).
c        intpo       :  integer vector giving the interval possition
c                       for the candidate value of the random effects,
c                       intpo(nsubject).
c        invcdfnorm  :  quantile function for a normal distribution.
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        je2         :  index.
c        k           :  index.   
c        l           :  index.
c        k1          :  index.   
c        k2          :  index.  
c        logcgkn     :  real working variable.
c        logcgko     :  real working variable.
c        loglikn     :  real working variable.
c        logliko     :  real working variable.
c        logpriorn   :  real working variable.
c        logprioro   :  real working variable.
c        nint        :  index
c        npoints     :  index.
c        nscan       :  index.
c        prob        :  real vector giving the probability of 
c                       the intervals, prob(ninter).
c        quan        :  real working variable.
c        ratio       :  real working variable.
c        rbeta       :  real beta random number generator.
c        rtlnorm     :  real truncated log normal random number generator.
c        runif       :  real uniform random number generator.
c        rvecs       :  real matrix giving the random vectors for the
c                       Polya tree,  rvecs(nlevel,ninter).
c        sec         :  cpu time working variable.
c        sec0        :  cpu time working variable.
c        sec00       :  cpu time working variable.
c        sec1        :  cpu time working variable.
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        skipcount   :  index. 
c        sprint      :  integer function to print on screen.
c        thetac      :  real giving the current value of the candidate
c                       for the quantile parameter.
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
c        tmp3        :  real used to accumulate quantities.
c        work1       :  real matrix used to update the difficulties,
c                       work1(p-1,p-1).
c        work2       :  real matrix used to update the difficulties,
c                       work2(p-1,p-1).
c        work3       :  real matrix used to update the difficulties,
c                       work3(p-1,p-1).
c        workmh1     :  real vector used to update the difficulties,
c                       workmh1((p-1)*p/2).
c        workv1      :  real vector used to update the difficulties,
c                       workv1(p-1).
c        workv2      :  real vector used to update the difficulties,
c                       workv2(p-1).
c        workv3      :  real vector used to update the difficulties,
c                       workv3(p-1).
c        xtx         :  real matrix givind the product X^tX, 
c                       xtx(p-1,p-1).
c        xty         :  real vector used to save the product 
c                       Xt(Y-Zb), xty(p-1).
c
c=======================================================================

      implicit none 

c+++++Data
      integer imiss,ngrid,nmissing,nsubject,p
      integer datastr(nmissing,2),y(nsubject,p)
      real*8 roffset(nsubject,p)
      
c+++++Prior 
      integer ninter,nlevel
      real*8 aa0,ab0,a0b0(2),b0(p-1),prec(p-1,p-1)
      real*8 sb(p-1)
      real*8 tau1,tau2,m,s

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay
      real*8 tune3,tune4,tune5

c+++++Output
      real*8 acrate(5)
      real*8 cpo(nsubject,p),f(ngrid),faccum(ngrid)
      real*8 randsave(nsave,nsubject+1)
      real*8 thetasave(nsave,p+2)

c+++++Current values of the parameters
      real*8 alpha,beta(p-1),b(nsubject)
      real*8 mu,sigma

c+++++Working space - General
      integer i,ii,intlp,j,je2,k,k1,k2
      integer nint,npoints,sprint
      real*8 dbet,dlnrm
      real*8 quan
      real*8 tmp1,tmp2,tmp3

c+++++Working space - Random effects
      integer imax,imin
      real*8 grid(ngrid),thetac
      real*8 zty,ztz,ztzinv

c+++++Working space - Random effects Distribution
      integer assignb(nsubject,nlevel)
      integer counter(nlevel,ninter)
      integer intpn(nsubject),intpo(nsubject)
      real*8 accums(nlevel,ninter)
      real*8 endp(ninter-1)
      real*8 prob(ninter)
      real*8 rvecs(nlevel,ninter)
      real*8 dnrm,dpoiss,cdfnorm,invcdfnorm

c+++++Working space - RNG
      integer rpois,seed(2),seed1,seed2
      real*8 rbeta,rgamma,rnorm,rtlnorm
      real runif

c+++++Working space - MCMC
      integer iscan,isave,nscan
      integer skipcount,dispcount

c+++++Working space - Difficulty parameters
      integer iflag(p-1)
      real*8 betac(p-1)
      real*8 detlog
      real*8 xtx(p-1,p-1),xty(p-1)
      real*8 workmh1((p-1)*p/2)
      real*8 work1(p-1,p-1),work2(p-1,p-1),work3(p-1,p-1)
      real*8 workv1(p-1),workv2(p-1),workv3(p-1)

c+++++Working space - GLM part
      integer yij
      real*8 acrate2
      real*8 eta,etac,gprime,gprimec,mean,meanc,offset,ytilde,ytildec

c+++++Working space - MH 
      integer nu
      real*8 alphac,sigmac,muc
      real*8 logcgkn,logcgko
      real*8 loglikn,logliko
      real*8 logpriorn,logprioro
      real*8 ratio,ssb

c+++++CPU time
      real*8 sec00,sec0,sec1,sec

c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      aa0=a0b0(1)
      ab0=a0b0(2)
      
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

      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++++
c+++++++ missing data
c++++++++++++++++++++++++++++++++++
         if(imiss.eq.1)then
            do ii=1,nmissing
               i=datastr(ii,1)
               j=datastr(ii,2)
               if(j.eq.1)then
                 eta=b(i)+roffset(i,j)
                else
                 eta=b(i)-beta(j-1)+roffset(i,j)
               end if
               
               mean=exp(eta)
               y(i,j)=rpois(mean)
            end do
         end if


c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++++
c+++++++ difficulty parameters
c++++++++++++++++++++++++++++++++++

         do i=1,p-1
            do j=1,p-1
               xtx(i,j)=prec(i,j)
               work1(i,j)=0.d0
               work2(i,j)=0.d0
               work3(i,j)=0.d0
            end do
            xty(i)=sb(i)
            workv1(i)=0.d0
            workv2(i)=0.d0
            workv3(i)=0.d0
            iflag(i)=0
         end do

         logliko=0.d0         
         
         do i=1,nsubject
            do j=1,p-1
               yij=y(i,j+1)            
               eta=b(i)-beta(j)+roffset(i,j+1) 
               offset=b(i)+roffset(i,j+1) 
               mean=exp(eta)
               gprime=exp(-eta)
               ytilde=eta+(dble(yij)-mean)*gprime-offset
               
               xtx(j,j)=xtx(j,j)+1.d0/gprime
               xty(j)=xty(j)-ytilde/gprime

               logliko=logliko+dpoiss(dble(yij),mean,1)               
            end do
         end do

         do i=1,p-1
            do j=1,p-1
               work1(i,j)=xtx(i,j)          
            end do
         end do

         call invdet(work1,p-1,work2,detlog,iflag,workv1)
            
         do i=1,p-1
            tmp1=0.d0
            do j=1,p-1
               tmp1=tmp1+work2(i,j)*xty(j) 
            end do
            workv2(i)=tmp1
         end do

         call rmvnorm(p-1,workv2,work2,workmh1,workv3,betac)


c+++++++ evaluating the candidate generating kernel

         call dmvn(p-1,betac,workv2,work2,logcgko,
     &             workv1,work1,work3,workv3,iflag)                 

  
c+++++++ prior ratio

         logprioro=0.d0
         logpriorn=0.d0
         
         do i=1,p-1
            do j=1,p-1
               logpriorn=logpriorn+(betac(i)-b0(i))* 
     &                    prec(i,j)      *
     &                   (betac(j)-b0(j))

               logprioro=logprioro+(beta(i) -b0(i))* 
     &                    prec(i,j)      *
     &                   (beta(j) -b0(j))
            end do
         end do
         
         logpriorn=-0.5d0*logpriorn
         logprioro=-0.5d0*logprioro
            
c+++++++ candidate generating kernel contribution

         do i=1,p-1
            do j=1,p-1
               xtx(i,j)=prec(i,j)
               work1(i,j)=0.d0
               work2(i,j)=0.d0
               work3(i,j)=0.d0
            end do
            xty(i)=sb(i)
            workv1(i)=0.d0
            workv2(i)=0.d0
            workv3(i)=0.d0
            iflag(i)=0
         end do

         loglikn=0.d0         

         do i=1,nsubject
            do j=1,p-1
               yij=y(i,j+1)
               etac=b(i)-betac(j)+roffset(i,j+1)  
               offset=b(i)+roffset(i,j+1) 
               meanc=exp(etac)
               gprimec=exp(-etac)
               ytildec=etac+(dble(yij)-meanc)*gprimec-offset
               
               xtx(j,j)=xtx(j,j)+1.d0/gprimec
               xty(j)=xty(j)-ytildec/gprimec
               
               loglikn=loglikn+dpoiss(dble(yij),meanc,1)
            end do
         end do

         do i=1,p-1
            do j=1,p-1
               work1(i,j)=xtx(i,j)          
            end do
         end do
         
         call invdet(work1,p-1,work2,detlog,iflag,workv1)

         do i=1,p-1
            tmp1=0.d0
            do j=1,p-1
               tmp1=tmp1+work2(i,j)*xty(j) 
            end do
            workv2(i)=tmp1
         end do

c+++++++ evaluating the candidate generating kernel
            
         call dmvn(p-1,beta,workv2,work2,logcgkn,
     &             workv1,work1,work3,workv3,iflag)                 
 

c+++++++ mh step
           
         ratio=dexp(loglikn-logliko+logcgkn-logcgko+
     &              logpriorn-logprioro)

         if(dble(runif()).lt.ratio)then
            acrate(1)=acrate(1)+1.d0
            do i=1,p-1
               beta(i)=betac(i) 
            end do
         end if
       

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Generating from the full conditionals for RE´s distribution ++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         do i=1,nlevel
            nint=2**i
            do j=1,nint
               counter(i,j)=0
               accums(i,j)=0.d0
               rvecs(i,j)=0.d0
            end do
         end do   

         do i=1,nsubject
            nint=2
            tmp1=1.d0/dble(nint)
            quan=invcdfnorm(tmp1,mu,sigma,1,0)
            if(b(i).le.quan)then
                assignb(i,1)=1
                counter(1,1)=counter(1,1)+1
              else
                assignb(i,1)=2
                counter(1,2)=counter(1,2)+1
            end if  
        
            do j=2,nlevel
               nint=2**j
               tmp1=1.d0/dble(nint)            
               k=assignb(i,j-1)
               k1=2*(k-1)+1
               k2=2*(k-1)+2
               
               quan=invcdfnorm(dble(k1)*tmp1,mu,sigma,1,0)
               
               if(b(i).le.quan)then
                  assignb(i,j)=k1
                  counter(j,k1)=counter(j,k1)+1
                 else
                  assignb(i,j)=k2
                  counter(j,k2)=counter(j,k2)+1
               end if  
            end do
         end do


c+++++++ generating (Y0,Y1),(Y00,Y01,Y10,Y11),.... 

         tmp1=alpha+counter(1,1)
         tmp2=alpha+counter(1,2)
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
               tmp1=alpha*dble(je2)+dble(counter(i+1,k1))
               tmp2=alpha*dble(je2)+dble(counter(i+1,k2))
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

c+++++++ end points

         tmp1=1.d0/dble(ninter)  
         do i=1,npoints
            endp(i)=invcdfnorm(dble(i)*tmp1,mu,sigma,1,0)             
         end do

c+++++++ check if the user has requested an interrupt
         call rchkusr()


c++++++++++++++++++++++++++++++++++         
c+++++++ random effects 
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         acrate2=0.d0

         do i=1,nsubject
            ztz=0.d0
            ztzinv=0.d0
            zty=0.d0
            
            logliko=0.d0                     
           
            do j=1,p
               if(j.eq.1)then
                  eta=b(i)+roffset(i,j)
                  offset=roffset(i,j)
                 else
                  eta=b(i)-beta(j-1)+roffset(i,j)
                  offset=roffset(i,j)-beta(j-1)
               end if
               
               yij=y(i,j)
               mean=exp(eta)
               gprime=exp(-eta)
               ytilde=eta+(dble(yij)-mean)*gprime-offset    

               ztz=ztz+1.d0/gprime
               zty=zty+ytilde/gprime
               
               logliko=logliko+dpoiss(dble(yij),mean,1)               
            end do

            ztz=ztz+(1.d0/sigma**2)
            ztzinv=1.d0/ztz

            zty=zty+(1.d0/sigma**2)*mu

            tmp2=ztzinv*zty
 
            thetac=rnorm(b(i),sqrt(ztzinv))

            logcgko=dnrm(thetac,b(i),sqrt(ztzinv),1)

c++++++++++ prior ratio

            nint=2
            tmp1=1.d0/dble(nint)
            quan=invcdfnorm(tmp1,mu,sigma,1,0)
            if(thetac.le.quan)then
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
            
               if(thetac.le.quan)then
                 intlp=k1
               else
                 intlp=k2
               end if  
            end do

            tmp1=prob(intlp)*dble(ninter)*dnrm(thetac,mu,sigma,0) 
            logpriorn=log(tmp1)

            intlp=assignb(i,nlevel)
            tmp2=prob(intlp)*dble(ninter)*dnrm(b(i),mu,sigma,0) 
            logprioro=log(tmp2)

c++++++++++ candidate generating kernel contribution

            loglikn=0.d0                     

            ztz=0.d0
            ztzinv=0.d0
            zty=0.d0
             
            do j=1,p
               if(j.eq.1)then
                  etac=thetac+roffset(i,j)
                  offset=roffset(i,j)
                 else
                  etac=thetac-beta(j-1)+roffset(i,j)
                  offset=roffset(i,j)-beta(j-1)
               end if
               
               yij=y(i,j)
               meanc=exp(etac)
               gprimec=exp(-etac)
               ytildec=etac+(dble(yij)-meanc)*gprimec-offset    

               ztz=ztz+1.d0/gprimec
               zty=zty+ytildec/gprimec
               
               loglikn=loglikn+dpoiss(dble(yij),meanc,1)               
            end do

            ztz=ztz+(1.d0/sigma**2)
            ztzinv=1.d0/ztz

            zty=zty+(1.d0/sigma**2)*mu

            tmp2=ztzinv*zty

            logcgkn=dnrm(b(i),thetac,sqrt(ztzinv),1)

c++++++++++ mh step
           
            ratio=dexp(loglikn-logliko+
     &                 logcgkn-logcgko+            
     &                 logpriorn-logprioro)

            if(dble(runif()).lt.ratio)then
               acrate2=acrate2+1.d0
               b(i)=thetac
            end if
         end do

         acrate(2)=acrate(2)+acrate2/dble(nsubject)


c++++++++++++++++++++++++++++++++++++++         
c+++++++ Base line distribution-mu
c++++++++++++++++++++++++++++++++++++++

         if(s.gt.0.d0)then

         muc=rnorm(mu,tune3*sigma/sqrt(dble(nsubject)))
         
c+++++++ evaluate log-prior for candidate value of the parameters

         logpriorn=dnrm(muc,m,sqrt(s),1)  

c+++++++ evaluate log-prior for current value of parameters

         logprioro=dnrm(mu,m,sqrt(s),1)  

c+++++++ evaluate log-likelihood for current and candidate value 
c+++++++ of parameters

         logliko=0.d0
         loglikn=0.d0

         do i=1,nsubject

c++++++++++ possition lp current

            nint=2
            tmp1=1.d0/dble(nint)
            quan=invcdfnorm(tmp1,mu,sigma,1,0) 
                        
            if(b(i).le.quan)then
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
               
               if(b(i).le.quan)then
                  intlp=k1
                 else
                  intlp=k2
               end if  
            end do
            
            intpo(i)=intlp
            
c++++++++++ possition lp candidate

            nint=2
            tmp1=1.d0/dble(nint)
            quan=invcdfnorm(tmp1,muc,sigma,1,0) 
            
            if(b(i).le.quan)then
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
               
               if(b(i).le.quan)then
                  intlp=k1
                 else
                  intlp=k2
               end if  
            end do
            
            intpn(i)=intlp

c++++++++++ likelihood current

            tmp1=prob(intpo(i))*dble(ninter)*dnrm(b(i),mu,sigma,0) 
            logliko=logliko+log(tmp1)

c++++++++++ likelihood candidate

            tmp2=prob(intpn(i))*dble(ninter)*dnrm(b(i),muc,sigma,0) 
            loglikn=loglikn+log(tmp2)

         end do


c+++++++ acceptance/rejection step

         ratio=dexp(loglikn+logpriorn-logliko-logprioro)

         if(dble(runif()).lt.ratio)then
            mu=muc
            do i=1,nsubject
               intpo(i)=intpn(i)
            end do
            acrate(3)=acrate(3)+1.d0
         end if
         
         end if


c++++++++++++++++++++++++++++++++++++++         
c+++++++ Base line distribution-sigma
c++++++++++++++++++++++++++++++++++++++

         if(tau1.gt.0.d0)then

         nu=(dble(nsubject))*tune4
         ssb=(sigma**2)*dble(nu)

         sigmac=sqrt(1.d0/
     &          rgamma(0.5d0*(dble(nu)),0.5d0*ssb))

         logcgko=
     &      dble(nu)*log(sigma)-(0.5*dble(nu)+1.d0)*log(sigmac**2)
     &     -ssb/(2.d0*(sigmac**2))
         

         ssb=(sigmac**2)*dble(nu)
         logcgkn=
     &      dble(nu)*log(sigmac)-(0.5*dble(nu)+1.d0)*log(sigma**2)
     &     -ssb/(2.d0*(sigma**2))


c+++++++ evaluate log-prior for candidate value of the parameters

         call dgamma2(1.d0/(sigmac**2),0.5d0*tau1,0.5d0*tau2,logpriorn)  

c+++++++ evaluate log-prior for current value of parameters

         call dgamma2(1.d0/(sigma**2),0.5d0*tau1,0.5d0*tau2,logprioro)  


c+++++++ evaluate log-likelihood for current and candidate value 
c+++++++ of parameters

         logliko=0.d0
         loglikn=0.d0

         do i=1,nsubject

c++++++++++ possition lp candidate

            nint=2
            tmp1=1.d0/dble(nint)
            quan=invcdfnorm(tmp1,mu,sigmac,1,0) 
            
            if(b(i).le.quan)then
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
               
               if(b(i).le.quan)then
                  intlp=k1
                 else
                  intlp=k2
               end if  
            end do
            
            intpn(i)=intlp

c++++++++++ likelihood current

            tmp1=prob(intpo(i))*dble(ninter)*dnrm(b(i),mu,sigma,0) 
            logliko=logliko+log(tmp1)

c++++++++++ likelihood candidate

            tmp2=prob(intpn(i))*dble(ninter)*dnrm(b(i),mu,sigmac,0) 
            loglikn=loglikn+log(tmp2)

         end do


c+++++++ acceptance/rejection step

         ratio=dexp(loglikn+logpriorn-logliko-logprioro+
     &              logcgkn-logcgko)

         if(dble(runif()).lt.ratio)then
            sigma=sigmac
            do i=1,nsubject
               intpo(i)=intpn(i)
            end do
            acrate(4)=acrate(4)+1.d0
         end if
         
         end if


c++++++++++++++++++++++++++++++++++         
c+++++++ Precision parameter
c++++++++++++++++++++++++++++++++++

         if(aa0.gt.0.d0)then
c++++++++++ sample candidates

            alphac=rtlnorm(log(alpha),tune5*0.1d0,0,0,.true.,.true.)
            logcgkn=dlnrm(alpha ,log(alphac),tune5*0.1d0,1) 
            logcgko=dlnrm(alphac,log(alpha ),tune5*0.1d0,1) 

c++++++++++ evaluate log-prior for candidate value of the parameters

            call dgamma2(alphac,aa0,ab0,logpriorn)  

c++++++++++ evaluate log-prior for current value of parameters

            call dgamma2(alpha,aa0,ab0,logprioro)

c+++++++++++ evaluate log-likelihood

            tmp1=alpha
            tmp2=alpha
            logliko=dbet(rvecs(1,1),tmp1,tmp2,1)
            
            tmp1=alphac
            tmp2=alphac
            loglikn=dbet(rvecs(1,1),tmp1,tmp2,1)

            do i=1,nlevel-1
               nint=2**i
               je2=(i+1)**2
               do j=1,nint
                  k1=2*(j-1)+1
                  k2=2*(j-1)+2            
                  tmp1=alpha*dble(je2)
                  tmp2=alpha*dble(je2)
                  logliko=logliko+dbet(rvecs(i+1,k1),tmp1,tmp2,1)

                  tmp1=alphac*je2
                  tmp2=alphac*je2
                  loglikn=loglikn+dbet(rvecs(i+1,k1),tmp1,tmp2,1)
               end do
            end do   

c++++++++++ acceptance step
            ratio=dexp(loglikn+logpriorn-logliko-logprioro+
     &                 logcgkn-logcgko)

            if(dble(runif()).lt.ratio)then
               alpha=alphac
               acrate(5)=acrate(5)+1.d0
            end if            
            
         end if 


c++++++++++++++++++++++++++++++++++         
c+++++++ save samples
c++++++++++++++++++++++++++++++++++         
         
         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1

c+++++++++++++ difficulty parameters

               do i=1,p-1
                  thetasave(isave,i)=beta(i)
               end do

c+++++++++++++ baseline mean

               thetasave(isave,p)=mu

c+++++++++++++ baseline stansard deviation

               thetasave(isave,p+1)=sigma

c+++++++++++++ cluster information

               thetasave(isave,p+2)=alpha

c+++++++++++++ random effects

               do i=1,nsubject
                  randsave(isave,i)=b(i)
               end do

c+++++++++++++ predictive information

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
            
                  if(intlp.eq.1)then
                     tmp1=0.d0
                     tmp2=prob(1)*
     &                    ( dble(ninter)*
     &                      cdfnorm(grid(i),mu,sigma,1,0)
     &                    )
                    else
                     imin=1
                     imax=intlp-1
                     tmp1=0.d0
                     do j=imin,imax
                        tmp1=tmp1+prob(j)
                     end do
                     tmp2=prob(intlp)*
     &                   ( dble(ninter)*
     &                     cdfnorm(grid(i),mu,sigma,1,0)-
     &                     dble(intlp-1)
     &                   )
                     tmp2=tmp1+tmp2                    
                  end if  
                  faccum(i)=faccum(i)+tmp2

                  f(i)=f(i)+prob(intlp)*dble(ninter)*
     &                 dnrm(grid(i),mu,sigma,0)

               end do

               tmp3=dble(runif())
               tmp1=prob(1)
               j=1
               do while(tmp3.gt.tmp1.and.j.lt.ninter)
                  j=j+1
                  tmp1=tmp1+prob(j)
               end do
               tmp2=(tmp3-tmp1+dble(j)*prob(j))/(dble(ninter)*prob(j))
               thetac=invcdfnorm(tmp2,mu,sigma,1,0) 
               
               randsave(isave,nsubject+1)=thetac

c+++++++++++++ cpo

               do i=1,nsubject
                  do j=1,p
                     yij=y(i,j)
                     if(j.eq.1)then
                       eta=b(i)+roffset(i,j)
                      else
                       eta=b(i)-beta(j-1)+roffset(i,j)
                     end if  
                     mean=exp(eta)
                     tmp1=dpoiss(dble(yij),mean,1)
                     cpo(i,j)=cpo(i,j)+1.0d0/exp(tmp1)   
                  end do
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
      
      do i=1,5
         acrate(i)=acrate(i)/dble(nscan)    
      end do   
      
      do i=1,nsubject
         do j=1,p
            cpo(i,j)=dble(nsave)/cpo(i,j)
         end do   
      end do

      do i=1,ngrid
         f(i)=f(i)/dble(nsave)
         faccum(i)=faccum(i)/dble(nsave)
      end do

      return
      end

