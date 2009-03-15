
c=======================================================================                      
      subroutine fptrasch(datastrm,imiss,ngrid,nmissing,nsubject,p,y,
     &                    ninter,nlevel,  
     &                    a0b0,b0,m,prec,s,tau1,tau2,
     &                    mcmc,nsave,tune3,tune4,tune5,
     &                    acrate,cpo,f,faccum,randsave,thetasave,
     &                    alpha,b,beta,mu,sigma,
     &                    accums,assignb,betac,propvf,bc,propvr,
     &                    counter,endp,
     &                    intpn,intpo,prob,rvecs,seed,workmh1,workv1,
     &                    grid)
c=======================================================================                      
c
c     Version 1.0: 
c     Last modification: 30-09-2006.
c
c     Subroutine `fptrasch' to run a Markov chain in the  
c     semiparametric Rasch model using a Polya tree prior
c     for the random effect distribution. 
c
c     Copyright: Alejandro Jara, 2006-2009.
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
c      Facultad de Ciencias Físicas y Matemáticas
c      Universidad de Concepción
c      Avenida Esteban Iturra S/N
c      Barrio Universitario
c      Concepción
c      Chile
c      Voice: +56-41-2203163  URL  : http://www2.udec.cl/~ajarav
c      Fax  : +56-41-2251529  Email: ajarav@udec.cl
c
c---- Data -------------------------------------------------------------
c 
c        datastrm    :  integer matrix giving the position of missing 
c                       data points, datastrm(nmissing,2)
c        imiss       :  integer indicating whether missing data are
c                       present (1) or absent (0). 
c        ngrid       :  integer giving the size of the grid where
c                       the cdf estimate is evaluated.
c        nmissing    :  integer giving the number of missing 
c                       observations.
c        nsubject    :  integer giving the number of subjects.
c        p           :  integer giving the number of items.
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
c        m           :  real vector giving the prior mean for the 
c                       baseline mean, m.
c        ninter      :  integer giving the number of final intervals 
c                       in the Finite Polya tree prior.
c        nlevel      :  integer giving the number of binary partitions
c                       in the Finite Polya tree prior.
c        prec        :  real matrix giving the prior precision matrix
c                       for the difficulty parameters, prec(p-1,p-1).
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
c                       assignb(nsubject,nlevel).
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
c        workmh1     :  real vector used to update the fixed effects,
c                       workmh1((p-1)*p/2).
c        workv1      :  real vector used to update the fixed effects,
c                       workv1(p-1).
c
c=======================================================================

      implicit none 

c+++++Data
      integer imiss,ngrid,nmissing,nsubject,p
      integer datastrm(nmissing,2),y(nsubject,p)

c+++++Prior 
      integer ninter,nlevel
      real*8 aa0,ab0,a0b0(2),b0(p-1),m,prec(p-1,p-1)
      real*8 s
      real*8 tau1,tau2

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
      real*8 eta,dbet,dbin,dlnrm,dnrm
      real*8 mean,quan
      real*8 tmp1,tmp2,tmp3

c+++++Working space - Random effects
      integer imax,imin
      real*8 bc(nsubject),grid(ngrid),thetac,propvr
      real*8 targetr

c+++++Working space - Random effects Distribution
      integer assignb(nsubject,nlevel)
      integer counter(nlevel,ninter)
      integer intpn(nsubject),intpo(nsubject)
      real*8 accums(nlevel,ninter)
      real*8 endp(ninter-1)
      real*8 prob(ninter)
      real*8 rvecs(nlevel,ninter)
      real*8 cdfnorm,invcdfnorm

c+++++Working space - RNG
      integer evali,seed(2),seed1,seed2
      real*8 rbeta,rnorm,rtlnorm
      real runif

c+++++Working space - MCMC
      integer iscan,isave,nscan
      integer skipcount,dispcount

c+++++Working space - Difficulty parameters
      real*8 betac(p-1)
      real*8 targetd,propvf(p-1,p-1)
      real*8 workmh1((p-1)*p/2),workv1(p-1)
      
c+++++Working space - GLM part
      integer yij
      real*8 acrate2
      real*8 logp

c+++++Working space - MH 
      real*8 alphac,muc,sigmac,logcgkc,logcgko,logliko,loglikc,ratio
      real*8 logpriorc,logprioro 

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
      
      do i=1,nsubject
         bc(i)=b(i)
      end do
      

      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++++
c+++++++ missing data
c++++++++++++++++++++++++++++++++++
         if(imiss.eq.1)then
            do ii=1,nmissing
               i=datastrm(ii,1)
               j=datastrm(ii,2)
               if(j.eq.1)then
                 eta=b(i)
                else
                 eta=b(i)-beta(j-1)
               end if
               
               mean=exp(eta)/(1.d0+exp(eta))
               
               call rbinom(1,mean,evali)
               y(i,j)=evali
            end do
         end if


c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++++
c+++++++ difficulty parameters
c++++++++++++++++++++++++++++++++++

         call rmvnorm(p-1,beta,propvf,workmh1,workv1,betac)
         
         logp=targetd(y,nsubject,p,betac,b,b0,prec)        
         logp=logp-targetd(y,nsubject,p,beta,b,b0,prec)        

c+++++++ mh step

         if(log(dble(runif())).lt.logp)then
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

         tmp1=alpha+dble(counter(1,1))
         tmp2=alpha+dble(counter(1,2))
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

            thetac=rnorm(b(i),sqrt(propvr))
            
            bc(i)=thetac
         
            logp=targetr(i,nsubject,p,bc,beta,nlevel,
     &                   ninter,prob,mu,sigma,y)        
     
            logp=logp-targetr(i,nsubject,p,b,beta,nlevel,
     &                        ninter,prob,mu,sigma,y)        

c++++++++++ mh step

            if(log(dble(runif())).lt.logp)then
               acrate2=acrate2+1.d0
               b(i)=thetac
              else
               bc(i)=b(i)
            end if
         end do

         acrate(2)=acrate(2)+acrate2/dble(nsubject)


c++++++++++++++++++++++++++++++++++++++         
c+++++++ Base line distribution-mu
c++++++++++++++++++++++++++++++++++****

         if(s.gt.0.d0)then

         muc=rnorm(mu,tune3*0.025d0)
         
         logcgkc=dnrm(mu ,muc,tune3*0.025d0,1) 
         
         logcgko=dnrm(muc,mu ,tune3*0.025d0,1)

c+++++++ evaluate log-prior for candidate value of the parameters

         logpriorc=dnrm(muc,m,sqrt(s),1)  

c+++++++ evaluate log-prior for current value of parameters

         logprioro=dnrm(mu,m,sqrt(s),1)  

c+++++++ evaluate log-likelihood for current and candidate value 
c+++++++ of parameters

         logliko=0.d0
         loglikc=0.d0

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
            loglikc=loglikc+log(tmp2)

         end do


c+++++++ aceptance/rejection step

         ratio=dexp(loglikc+logpriorc-logliko-logprioro+
     &              logcgkc-logcgko)

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

         sigmac=rtlnorm(log(sigma),tune4*0.025d0,0,0,.true.,.true.)

         logcgkc=dlnrm(sigma ,log(sigmac),tune4*0.025d0,1) 
         
         logcgko=dlnrm(sigmac,log(sigma) ,tune4*0.025d0,1)

c+++++++ evaluate log-prior for candidate value of the parameters

         call dgamma2(1.d0/(sigmac**2),0.5d0*tau1,0.5d0*tau2,logpriorc)

c+++++++ evaluate log-prior for current value of parameters

         call dgamma2(1.d0/(sigma**2),0.5d0*tau1,0.5d0*tau2,logprioro)


c+++++++ evaluate log-likelihood for current and candidate value 
c+++++++ of parameters

         logliko=0.d0
         loglikc=0.d0

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
            loglikc=loglikc+log(tmp2)

         end do


c+++++++ aceptance/rejection step

         ratio=dexp(loglikc+logpriorc-logliko-logprioro+
     &              logcgkc-logcgko)

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
            logcgkc=dlnrm(alpha ,log(alphac),tune5*0.1d0,1) 
            logcgko=dlnrm(alphac,log(alpha ),tune5*0.1d0,1) 

c++++++++++ evaluate log-prior for candidate value of the parameters

            call dgamma2(alphac,aa0,ab0,logpriorc)  

c++++++++++ evaluate log-prior for current value of parameters

            call dgamma2(alpha,aa0,ab0,logprioro)

c+++++++++++ evaluate log-likelihood

            tmp1=alpha
            tmp2=alpha
            logliko=dbet(rvecs(1,1),tmp1,tmp2,1)
            
            tmp1=alphac
            tmp2=alphac
            loglikc=dbet(rvecs(1,1),tmp1,tmp2,1)

            do i=1,nlevel-1
               nint=2**i
               je2=(i+1)**2
               do j=1,nint
                  k1=2*(j-1)+1
                  k2=2*(j-1)+2            
                  tmp1=alpha*dble(je2)
                  tmp2=alpha*dble(je2)
                  logliko=logliko+dbet(rvecs(i+1,k1),tmp1,tmp2,1)

                  tmp1=alphac*dble(je2)
                  tmp2=alphac*dble(je2)
                  loglikc=loglikc+dbet(rvecs(i+1,k1),tmp1,tmp2,1)
               end do
            end do   

c++++++++++ acceptance step
            ratio=dexp(loglikc+logpriorc-logliko-logprioro+
     &                 logcgkc-logcgko)

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

c+++++++++++++ baseline standard mean

               thetasave(isave,p)=mu

c+++++++++++++ baseline standard deviation

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
                       eta=b(i)
                      else
                       eta=b(i)-beta(j-1)
                     end if  
                     mean=exp(eta)/(1.d0+exp(eta))
                     tmp1=dbin(dble(yij),1.d0,mean,1)
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


c=======================================================================
      double precision function targetd(y,nsubject,p,beta,b,b0,prec)
c=======================================================================
c     calculates the logarithm of the full conditional distribution
c     of the difficulty parameters in a Rasch model.
      implicit none
      integer i,j,nsubject,p
      integer y(nsubject,p)
      real*8 b(nsubject),beta(p-1),b0(p-1),prec(p-1,p-1)
      real*8 eta,mean,loglik,logpr
      real*8 dbin,tmp1

      loglik=0.d0
      logpr=0.d0
      
      do i=1,nsubject
         do j=1,p-1
            eta=b(i)-beta(j) 

            mean=exp(eta)/(1.d0+exp(eta))
            loglik=loglik+dbin(dble(y(i,j+1)),1.d0,mean,1)
         end do
      end do


      tmp1=0.d0
        
      do i=1,p-1
         do j=1,p-1
            tmp1=tmp1+(beta(i)-b0(i))* 
     &                 prec(i,j)      *
     &                (beta(j)-b0(j))

         end do
      end do
      logpr=-0.5d0*tmp1
      
      targetd=loglik+logpr      
      
      return
      end


c=======================================================================
      double precision function targetr(index,nsubject,p,b,beta,nlevel,
     &                                  ninter,prob,mu,sigma,y)        
c=======================================================================
c     calculates the logarithm of the full conditional distribution
c     of the random effects in a Rasch model.
      implicit none
      integer index,j,nint,ninter,nlevel,nsubject,p
      integer y(nsubject,p),intlp,k1,k2,k
      real*8 b(nsubject),beta(p-1)
      real*8 eta,mean,loglik,logpr
      real*8 dbin,tmp1
      real*8 dnrm,invcdfnorm,prob(ninter),mu,sigma,quan

      loglik=0.d0
      logpr=0.d0
      
      do j=1,p
         if(j.eq.1)then
            eta=b(index)
           else
            eta=b(index)-beta(j-1)
         end if
         mean=exp(eta)/(1.d0+exp(eta))
         loglik=loglik+dbin(dble(y(index,j)),1.d0,mean,1)
      end do


      nint=2
      tmp1=1.d0/dble(nint)
      quan=invcdfnorm(tmp1,mu,sigma,1,0)
      if(b(index).le.quan)then
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
      
         if(b(index).le.quan)then
           intlp=k1
         else
           intlp=k2
         end if  
      end do

      tmp1=prob(intlp)*dble(ninter)*dnrm(b(index),mu,sigma,0) 
      tmp1=log(tmp1)

      logpr=tmp1
      
      targetr=loglik+logpr      
      
      return
      end

