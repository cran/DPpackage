
c=======================================================================                      
      subroutine ptraschpoi(datastr,imiss,nmissing,nsubject,p,y, 
     &                      roffset,
     &                      ngrid,grid,                                   
     &                      maxm,ntprob,ntsets,                             
     &                      a0b0,b0,prec,sb,tau1,tau2,m,s,             
     &                      mcmc,nsave,algo,              
     &                      acrate,cpo,cpov,
     &                      densave,cdfsave,randsave,thetasave,bf,    
     &                      alpha,b,beta,mu,sigma2,probmat, 
     &                      seed,         
     &                      ptcount,kvec,prob,   
     &                      betac,workvp,workmhp,iflagp,              
     &                      xtx,xty)
c=======================================================================                    
c
c     Subroutine `ptraschpoi' to run a Markov chain in the  
c     semiparametric Rasch Poisson Count model using a Polya tree prior
c     for the random effect distribution. 
c
c     Copyright: Alejandro Jara, 2006-2010.
c
c     Version 3.0: 
c     Last modification: 04-09-2009.
c
c     Changes and Bug fixes: 
c
c     Version 2.0 to Version 3.0:
c          - Difficulty parameters can be sampled using MH or slice sampler.
c          - Random effects and PT parameters are sampled using a slice sampler.
c          - Optimize PT part of the model.
c          - Add computation of functionals.
c          - Samples of the density and CDF are returned now.  
c          - BF is computed using the Savage-Dickey ratio.
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
c        grid        :  real vector giving the grid, grid(ngrid).
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
c        ntsets      :  integer giving the number of final intervals 
c                       in the Finite Polya tree prior.
c        maxm        :  integer giving the number of binary partitions
c                       in the Finite Polya tree prior.
c        ntprob      :  integer giving the number of total conditional
c                       probabilities.
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
c        algo        :  integer indicating the algorithm used to sample
c                       difficulty parameters, algo=1 MH-WLS, 
c                       algo=2 Slice sampling in each coordinate.
c        
c-----------------------------------------------------------------------
c
c---- Output -----------------------------------------------------------
c
c        acrate      :  real giving the MH acceptance rate for  
c                       difficulty parameters. 
c        cpo         :  real matrix giving the cpo, cpo(nsubject,p).
c        densave     :  real matrix giving the debsity samples,
c                       densave(nsave,ngrid).
c        cdfsave     :  real matrix giving the CDF samples,
c                       cdfsave(nsave,ngrid).
c        randsave    :  real matrix containing the mcmc samples for
c                       the random effects and prediction,
c                       randsave(nsave,nsubject+1)
c        thetasave   :  real matrix containing the mcmc samples for
c                       the fixed effects,and mean and covariance of
c                       the baseline distribution, 
c                       thetsave(nsave,p+4).
c        bf          :  real giving the Bayes factor.
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
c        sigma2      :  real giving the current value of the
c                       variance for normal base line 
c                       distribution for the random effects,
c                       sigma.
c        probmat     :  real matrix giving the PT conditional 
c                       probabilities, probmat(ntsets,maxm).
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        acrate2     :  real used to calculate the acceptance rate. 
c        alphac      :  real giving the current value of the candidate
c                       for the precision parameter.
c        betac       :  real vector giving the current value of the 
c                       candidate for difficulty parameters, 
c                       betac(p-1).
c        dbet        :  density of a beta distribution.
c        dispcount   :  index. 
c        dlnrm       :  density of a log-normal distribution.
c        dnrm        :  density of a normal distribution.
c        i           :  index. 
c        ii          :  index. 
c        imax        :  index.
c        imin        :  index.
c        iflagp      :  integer vector used to invert the of the lhs
c                       least square solution for the difficulties,
c                       iflagp(p-1).
c        invcdfnorm  :  quantile function for a normal distribution.
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        je2         :  index.
c        k           :  index.   
c        l           :  index.
c        k1          :  index.   
c        k2          :  index.  
c        kvec        :  integer working vector, kvec(maxm).
c        logcgkn     :  real working variable.
c        logcgko     :  real working variable.
c        loglikn     :  real working variable.
c        logliko     :  real working variable.
c        logpriorn   :  real working variable.
c        logprioro   :  real working variable.
c        nscan       :  index.
c        ptcount     :  integer working vector, ptcount(ntprob).
c        prob        :  real working vector, prob(ntsets).
c        ratio       :  real working variable.
c        rbeta       :  real beta random number generator.
c        rtnorm      :  real truncated normal random number generator.
c        runif       :  real uniform random number generator.
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
c        workmhp     :  real vector used to update the difficulties,
c                       workmhp((p-1)*p/2).
c        workvp      :  real vector used to update the difficulties,
c                       workvp(p-1).
c        xtx         :  real matrix givind the product X^tX, 
c                       xtx(p-1,p-1).
c        xty         :  real vector used to save the product 
c                       Xt(Y-Zb), xty(p-1).
c
c=======================================================================

      implicit none 

c+++++Data
      integer imiss,nmissing,nsubject,p
      integer datastr(nmissing,2),y(nsubject,p)
      double precision roffset(nsubject,p)

c+++++Density and CDF
      integer ngrid
      double precision grid(ngrid)
      
c+++++Prior 
      integer maxm,ntprob,ntsets
      double precision aa0,ab0,a0b0(2),b0(p-1),prec(p-1,p-1)
      double precision sb(p-1)
      double precision tau1,tau2,m,s

c+++++MCMC parameters
      integer algo
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Output
      double precision acrate
      double precision cpo(nsubject,p)
      double precision cpov(nsubject)
      double precision densave(nsave,ngrid),cdfsave(nsave,ngrid)
      double precision randsave(nsave,nsubject+1)
      double precision thetasave(nsave,p+4)
      double precision bf

c+++++Current values of the parameters
      double precision alpha,beta(p-1),b(nsubject)
      double precision mu,sigma2
      double precision probmat(ntsets,maxm)

c+++++Working space - Polya tree
      integer ptcount(ntprob)
      integer kvec(maxm)
      double precision prob(ntsets)

c+++++Working space - Difficulty parameters
      integer iflagp(p-1)
      double precision betac(p-1)
      double precision xtx(p-1,p-1),xty(p-1)
      double precision workmhp((p-1)*p/2)
      double precision workvp(p-1)

c+++++RN
      integer seed(2),seed1,seed2

c+++++Working space - General
      integer i,i1,ii,j,jj,j1,j2,k1,evali
      integer n1,n2,sprint
      double precision dbet
      double precision tmp1,tmp2,tmp3,tmp4

c+++++Working space - Random effects
      double precision thetac

c+++++Working space - Ditributions and RNG
      integer rpois
      double precision rbeta,rtnorm
      real runif
      double precision dnrm,dpoiss,cdfnorm,invcdfnorm

c+++++Working space - MCMC
      integer iscan,isave,nscan
      integer skipcount,dispcount

c+++++Working space - GLM part
      integer yij
      double precision eta,gprime,mean,offset,ytilde

c+++++Polya tree
      integer kphi
      double precision mureal
      double precision mu2real
      double precision sigma2real
      double precision liminf,limsup

c+++++Working space - MH 
      double precision logcgkn,logcgko
      double precision loglikn,logliko
      double precision logpriorn,logprioro
      double precision ratio
      double precision logdenom
      double precision lognumer
      double precision numerBF
      double precision denomBF

c+++++Working space slice sampling
      double precision rexpo,re,uwork
      double precision logy,xx0,xx1,llim,rlim
      double precision grlim,gllim,gxx1

c+++++CPU time
      double precision sec00,sec0,sec1,sec

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

c++++ set probabilities

      numerBF=0.d0
      denomBF=0.d0
 
      do i=1,ntsets
         kphi=i
         do j1=1,maxm
            j2=maxm-j1+1
            kvec(j2)=kphi  
            kphi=int(real(kphi+1)/2.0)
         end do

         logprioro=0.d0
         do j1=1,maxm
            logprioro=logprioro+log(probmat(kvec(j1),j1))    
         end do
         prob(i)=exp(logprioro)
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

      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++++
c+++++++ missing data
c++++++++++++++++++++++++++++++++++
         if(imiss.eq.1)then

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

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

c++++++++++++++++++++++++++++++++++
c+++++++ difficulty parameters
c++++++++++++++++++++++++++++++++++

         if(algo.eq.1)then

         do i=1,p-1
            do j=1,p-1
               xtx(i,j)=prec(i,j)
            end do
            xty(i)=sb(i)
         end do

         logliko=0.d0         
         
         do i=1,nsubject

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

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

         call inverse(xtx,p-1,iflagp)      
            
         do i=1,p-1
            tmp1=0.d0
            do j=1,p-1
               tmp1=tmp1+xtx(i,j)*xty(j) 
            end do
            workvp(i)=tmp1
         end do

         call rmvnorm(p-1,workvp,xtx,workmhp,xty,betac)

c+++++++ evaluating the candidate generating kernel

         call dmvnd(p-1,betac,workvp,xtx,logcgko,iflagp)
  
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
            end do
            xty(i)=sb(i)
         end do

         loglikn=0.d0         

         do i=1,nsubject

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            do j=1,p-1
               yij=y(i,j+1)
               eta=b(i)-betac(j)+roffset(i,j+1)  
               offset=b(i)+roffset(i,j+1) 
               mean=exp(eta)
               gprime=exp(-eta)
               ytilde=eta+(dble(yij)-mean)*gprime-offset
               
               xtx(j,j)=xtx(j,j)+1.d0/gprime
               xty(j)=xty(j)-ytilde/gprime
               
               loglikn=loglikn+dpoiss(dble(yij),mean,1)
            end do
         end do

         call inverse(xtx,p-1,iflagp)      
            
         do i=1,p-1
            tmp1=0.d0
            do j=1,p-1
               tmp1=tmp1+xtx(i,j)*xty(j) 
            end do
            workvp(i)=tmp1
         end do

c+++++++ evaluating the candidate generating kernel

         call dmvnd(p-1,beta,workvp,xtx,logcgkn,iflagp)
            
c+++++++ mh step
           
         ratio=loglikn-logliko+logcgkn-logcgko+
     &              logpriorn-logprioro

         if(log(runif()).lt.ratio)then
            acrate=acrate+1.d0
            do i=1,p-1
               beta(i)=betac(i) 
            end do
         end if

c+++++++ SLICE sampler

         else

         do j=1,p-1

            evali=1
            xx0=beta(j)
            call targetrashpd(j,nsubject,p,y,b,beta,
     &                       roffset,prec,b0,tmp1)
            re=rexpo(1.d0)
            logy=tmp1-re

            uwork=dble(runif())  
            llim=xx0-uwork
            rlim=llim+1.d0

            evali=evali+1
            beta(j)=llim
            call targetrashpd(j,nsubject,p,y,b,beta,
     &                       roffset,prec,b0,gllim)

            evali=evali+1
            beta(j)=rlim
            call targetrashpd(j,nsubject,p,y,b,beta,
     &                       roffset,prec,b0,grlim)


            do while(gllim.gt.logy)
               llim=llim-1.d0

               evali=evali+1
               beta(j)=llim
               call targetrashpd(j,nsubject,p,y,b,beta,
     &                          roffset,prec,b0,gllim)
            end do 

            do while(grlim.gt.logy)
               rlim=rlim+1.d0

               evali=evali+1
               beta(j)=rlim
               call targetrashpd(j,nsubject,p,y,b,beta,
     &                          roffset,prec,b0,grlim)
            end do 

            xx1=llim+(rlim-llim)*dble(runif())

            evali=evali+1
            beta(j)=xx1
            call targetrashpd(j,nsubject,p,y,b,beta,
     &                       roffset,prec,b0,gxx1)

            do while(gxx1.lt.logy)
               if(xx1.gt.xx0)rlim=xx1
               if(xx1.lt.xx0)llim=xx1

               xx1=llim+(rlim-llim)*dble(runif())

               evali=evali+1
               beta(j)=xx1
               call targetrashpd(j,nsubject,p,y,b,beta,
     &                          roffset,prec,b0,gxx1)

            end do

            beta(j)=xx1

         end do
         
         end if

c++++++++++++++++++++++++++++++++++         
c+++++++ random effects 
c++++++++++++++++++++++++++++++++++

         do i=1,ntprob
            ptcount(i)=0
         end do

         do i=1,nsubject

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            evali=1
            xx0=b(i)
            call targetreslice(i,nsubject,p,y,roffset,maxm,ntsets,
     &                         prob,mu,sigma2,xx0,beta,tmp1)
            re=rexpo(1.d0)
            logy=tmp1-re

            uwork=dble(runif())  
            llim=xx0-uwork
            rlim=llim+1.d0

            evali=evali+1
            call targetreslice(i,nsubject,p,y,roffset,maxm,ntsets,
     &                         prob,mu,sigma2,llim,beta,gllim)

            evali=evali+1
            call targetreslice(i,nsubject,p,y,roffset,maxm,ntsets,
     &                         prob,mu,sigma2,rlim,beta,grlim)


            do while(gllim.gt.logy)
               llim=llim-1.d0

               evali=evali+1
               call targetreslice(i,nsubject,p,y,roffset,maxm,ntsets,
     &                            prob,mu,sigma2,llim,beta,gllim)
            end do 

            do while(grlim.gt.logy)
               rlim=rlim+1.d0

               evali=evali+1
               call targetreslice(i,nsubject,p,y,roffset,maxm,ntsets,
     &                            prob,mu,sigma2,rlim,beta,grlim)
            end do 

            xx1=llim+(rlim-llim)*dble(runif())

            evali=evali+1
            call targetreslice(i,nsubject,p,y,roffset,maxm,ntsets,
     &                         prob,mu,sigma2,xx1,beta,gxx1)

            do while(gxx1.lt.logy)
               if(xx1.gt.xx0)rlim=xx1
               if(xx1.lt.xx0)llim=xx1

               xx1=llim+(rlim-llim)*dble(runif())

               evali=evali+1
               call targetreslice(i,nsubject,p,y,roffset,maxm,ntsets,
     &                            prob,mu,sigma2,xx1,beta,gxx1)

            end do

            b(i)=xx1

            tmp1=(b(i)-mu)/sqrt(sigma2) 
            if(tmp1.gt.4.0d0)then
               tmp2=0.999968d0
              else if(tmp1.lt.-4.0d0)then
               tmp2=0.000032d0
              else 
               tmp2=cdfnorm(tmp1,0.d0,1.d0,1,0)
            end if

            kphi=int(real(2**maxm)*tmp2+1)

            do j1=1,maxm
               j2=maxm-j1+1
               kvec(j2)=kphi  
               kphi=int(real(kphi+1)/2.0)
            end do
     
            j=0
            do j1=1,maxm
               j2=j+kvec(j1)
               j=j+2**j1
               ptcount(j2)=ptcount(j2)+1
            end do

         end do


c+++++++++++++++++++++++++++++++++++++         
c+++++++ PT probabilities 
c+++++++++++++++++++++++++++++++++++++

         i1=0
         do i=1,maxm
            j1=2**(i-1)
            do j=1,j1
               k1=i1+j

               ii=(k1-1)*2+1
               jj=(k1-1)*2+2
  
               n1=ptcount(ii)
               n2=ptcount(jj)

               tmp2=alpha*dble(i**2)

               tmp1=rbeta(tmp2+dble(n1),tmp2+dble(n2))

               ii=(j-1)*2+1
               jj=(j-1)*2+2
               probmat(ii,i)=tmp1
               probmat(jj,i)=1.d0-tmp1
            end do
            i1=i1+j1
         end do

         do i=1,ntsets
            kphi=i
            do j1=1,maxm
               j2=maxm-j1+1
               kvec(j2)=kphi  
               kphi=int(real(kphi+1)/2.0)
            end do
            logprioro=0.d0
            do j1=1,maxm
               logprioro=logprioro+log(probmat(kvec(j1),j1))    
            end do
            prob(i)=exp(logprioro)
         end do

c++++++++++++++++++++++++++++++++++++++         
c+++++++ Base line distribution mu
c++++++++++++++++++++++++++++++++++++++

         if(s.gt.0.d0)then
           
            evali=1
            xx0=mu
            call llikaslicept(maxm,ntsets,prob,
     &                        nsubject,b,xx0,sigma2,tmp1)
            re=rexpo(1.d0)
            logy=dnrm(xx0,m,sqrt(s),1)
            logy=logy+tmp1-re
c            call intpr("evaluation #",-1,evali,1)

            uwork=dble(runif())  
            llim=xx0-uwork
            rlim=llim+1.d0

            evali=evali+1
            call llikaslicept(maxm,ntsets,prob,
     &                        nsubject,b,llim,sigma2,tmp1)
            gllim=dnrm(llim,m,sqrt(s),1)
            gllim=gllim+tmp1
c            call intpr("evaluation #",-1,evali,1)

            evali=evali+1
            call llikaslicept(maxm,ntsets,prob,
     &                        nsubject,b,rlim,sigma2,tmp1)
            grlim=dnrm(rlim,m,sqrt(s),1)
            grlim=grlim+tmp1
c            call intpr("evaluation #",-1,evali,1)

c            call dblepr("xx0",-1,xx0,1)
c            call dblepr("llim",-1,llim,1)
c            call dblepr("rlim",-1,rlim,1)
c            call dblepr("logy",-1,logy,1)
c            call dblepr("gllim",-1,gllim,1)
c            call dblepr("grlim",-1,grlim,1)

            do while(gllim.gt.logy)
               llim=llim-1.d0

               evali=evali+1
               call llikaslicept(maxm,ntsets,prob,
     &                           nsubject,b,llim,sigma2,tmp1)
               gllim=dnrm(llim,m,sqrt(s),1)
               gllim=gllim+tmp1
c               call intpr("evaluation l#",-1,evali,1)
            end do 

            do while(grlim.gt.logy)
               rlim=rlim+1.d0

               evali=evali+1
               call llikaslicept(maxm,ntsets,prob,
     &                           nsubject,b,rlim,sigma2,tmp1)
               grlim=dnrm(rlim,m,sqrt(s),1)
               grlim=grlim+tmp1
c               call intpr("evaluation r#",-1,evali,1)
            end do 

            xx1=llim+(rlim-llim)*dble(runif())

            evali=evali+1
            call llikaslicept(maxm,ntsets,prob,
     &                        nsubject,b,xx1,sigma2,tmp1)
            gxx1=dnrm(xx1,m,sqrt(s),1)
            gxx1=gxx1+tmp1
c            call intpr("evaluation c#",-1,evali,1)

c            call dblepr("xx1",-1,xx1,1)
c            call dblepr("gxx1",-1,gxx1,1)

            do while(gxx1.lt.logy)
               if(xx1.gt.xx0)rlim=xx1
               if(xx1.lt.xx0)llim=xx1

               xx1=llim+(rlim-llim)*dble(runif())

               evali=evali+1
               call llikaslicept(maxm,ntsets,prob,
     &                           nsubject,b,xx1,sigma2,tmp1)
               gxx1=dnrm(xx1,m,sqrt(s),1)
               gxx1=gxx1+tmp1

c               call intpr("evaluation c#",-1,evali,1)
c               call dblepr("xx1",-1,xx1,1)
c               call dblepr("gxx1",-1,gxx1,1)
            end do

            mu=xx1

c            call dblepr("mu",-1,mu,1)

         end if

c++++++++++++++++++++++++++++++++++++++         
c+++++++ Base line distribution-sigma2
c++++++++++++++++++++++++++++++++++++++

         if(tau1.gt.0.d0)then

            evali=evali+1
            xx0=log(sigma2)
            re=rexpo(1.d0)
            logy=-re
            logy=logy-(0.5d0*tau1+1.d0)*xx0-
     &                 0.5d0*tau2/exp(xx0)
            call llikaslicept(maxm,ntsets,prob,
     &                        nsubject,b,mu,exp(xx0),tmp1)
            logy=logy+tmp1

           
            uwork=dble(runif())  
            llim=xx0-uwork
            rlim=llim+1.d0

            evali=evali+1
            gllim=-(0.5d0*tau1+1.d0)*llim-
     &              0.5d0*tau2/exp(llim)
            call llikaslicept(maxm,ntsets,prob,
     &                        nsubject,b,mu,exp(llim),tmp1)
            gllim=gllim+tmp1


            evali=evali+1
            grlim=-(0.5d0*tau1+1.d0)*rlim-
     &              0.5d0*tau2/exp(rlim)
            call llikaslicept(maxm,ntsets,prob,
     &                        nsubject,b,mu,exp(rlim),tmp1)
            grlim=grlim+tmp1


            do while(gllim.gt.logy)
               llim=llim-1.d0

               evali=evali+1
               gllim=-(0.5d0*tau1+1.d0)*llim-
     &                 0.5d0*tau2/exp(llim)
               call llikaslicept(maxm,ntsets,prob,
     &                        nsubject,b,mu,exp(llim),tmp1)
               gllim=gllim+tmp1
            end do 

            do while(grlim.gt.logy)
               rlim=rlim+1.d0

               evali=evali+1
               grlim=-(0.5d0*tau1+1.d0)*rlim-
     &                 0.5d0*tau2/exp(rlim)
               call llikaslicept(maxm,ntsets,prob,
     &                        nsubject,b,mu,exp(rlim),tmp1)
               grlim=grlim+tmp1
            end do 

            xx1=llim+(rlim-llim)*dble(runif())

            evali=evali+1
            gxx1=-(0.5d0*tau1+1.d0)*xx1-
     &             0.5d0*tau2/exp(xx1)
            call llikaslicept(maxm,ntsets,prob,
     &                        nsubject,b,mu,exp(xx1),tmp1)
            gxx1=gxx1+tmp1


            do while(gxx1.lt.logy)
               if(xx1.gt.xx0)rlim=xx1
               if(xx1.lt.xx0)llim=xx1

               xx1=llim+(rlim-llim)*dble(runif())

               evali=evali+1
               gxx1=-(0.5d0*tau1+1.d0)*xx1-
     &                0.5d0*tau2/exp(xx1)
               call llikaslicept(maxm,ntsets,prob,
     &                           nsubject,b,mu,exp(xx1),tmp1)
               gxx1=gxx1+tmp1
            end do

            sigma2=exp(xx1)

c            call dblepr("sigma2",-1,sigma2,1)

         end if


c++++++++++++++++++++++++++++++++++         
c+++++++ Precision parameter
c++++++++++++++++++++++++++++++++++

         if(aa0.gt.0.d0)then

            evali=1
            xx0=alpha
            re=rexpo(1.d0)
            logy=-re

            call targetcslice(aa0,ab0,xx0,maxm,ntprob,
     &                        ntsets,ptcount,probmat,tmp1)
            logy=logy+tmp1

           
            uwork=dble(runif())*0.25d0  
            llim=xx0-uwork
            rlim=llim+0.25d0

            if(llim.lt.0.01d0)llim=0.01d0 

            evali=evali+1
            call targetcslice(aa0,ab0,llim,maxm,ntprob,
     &                        ntsets,ptcount,probmat,gllim)

            evali=evali+1
            call targetcslice(aa0,ab0,rlim,maxm,ntprob,
     &                        ntsets,ptcount,probmat,grlim)

            do while(gllim.gt.logy)
               llim=llim-0.25d0

               if(llim.lt.0.01d0)then
                  llim=0.01d0
                  gllim=logy-1.d0
                 else   
                  evali=evali+1
                  call targetcslice(aa0,ab0,llim,maxm,ntprob,
     &                           ntsets,ptcount,probmat,gllim)
               end if
            end do 

            do while(grlim.gt.logy)
               rlim=rlim+0.25d0

               evali=evali+1
               call targetcslice(aa0,ab0,rlim,maxm,ntprob,
     &                           ntsets,ptcount,probmat,grlim)
            end do 

c            call intpr("evaluation #",-1,evali,1)
c            call dblepr("xx0",-1,xx0,1)
c            call dblepr("llim",-1,llim,1)
c            call dblepr("rlim",-1,rlim,1)
c            call dblepr("logy",-1,logy,1)
c            call dblepr("gllim",-1,gllim,1)
c            call dblepr("grlim",-1,grlim,1)

            xx1=llim+(rlim-llim)*dble(runif())

            evali=evali+1
            call targetcslice(aa0,ab0,xx1,maxm,ntprob,
     &                        ntsets,ptcount,probmat,gxx1)

c            call intpr("evaluation c#",-1,evali,1)
c            call dblepr("xx1",-1,xx1,1)
c            call dblepr("gxx1",-1,gxx1,1)

            do while(gxx1.lt.logy)
               if(xx1.gt.xx0)rlim=xx1
               if(xx1.lt.xx0)llim=xx1
               if(llim.lt.0.01d0)llim=0.01d0 

               xx1=llim+(rlim-llim)*dble(runif())

               evali=evali+1
               call targetcslice(aa0,ab0,xx1,maxm,ntprob,
     &                           ntsets,ptcount,probmat,gxx1)

c               call intpr("evaluation c#",-1,evali,1)
c               call dblepr("xx1",-1,xx1,1)
c               call dblepr("gxx1",-1,gxx1,1)
            end do

            alpha=xx1

c            call dblepr("alpha",-1,alpha,1)

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

c+++++++++++++ mean and variances of random effects
c+++++++++++++ and random sample

               call simdisc(prob,ntsets,ntsets,evali)

               tmp1=1.d0/dble(ntsets)
               mureal=0.d0
               mu2real=0.d0
               liminf=-999.d9
               limsup=invcdfnorm(tmp1,0.d0,1.d0,1,0)

               do i=1,ntsets
                  logprioro=dble(maxm)*log(2.d0)
                  logprioro=logprioro+log(prob(i))

                  if(i.eq.1)then
                     tmp1=dnrm(limsup,0.d0,1.d0,0)
                     mureal=mureal-exp(logprioro)*tmp1

                     tmp3=1.d0/dble(ntsets)
                     tmp4=limsup*tmp1
                     mu2real=mu2real+exp(logprioro)*(tmp3-tmp4)

                     if(evali.eq.i)then
                        thetac=rtnorm(mu,sqrt(sigma2),0.d0,
     &                         limsup,.true.,.false.) 
                     end if                     

                  else if(i.eq.ntsets)then
                     tmp2=dnrm(liminf,0.d0,1.d0,0)
                     mureal=mureal+exp(logprioro)*tmp2

                     tmp3=1.d0/dble(ntsets)
                     tmp4=liminf*tmp2
                     mu2real=mu2real+exp(logprioro)*(tmp3+tmp4)

                     if(evali.eq.i)then
                        thetac=rtnorm(mu,sqrt(sigma2),liminf,
     &                         0.d0,.false.,.true.) 
                     end if                     

                  else
                     tmp1=dnrm(limsup,0.d0,1.d0,0)
                     tmp2=dnrm(liminf,0.d0,1.d0,0)
                     mureal=mureal-exp(logprioro)*(tmp1-tmp2)

                     tmp3=1.d0/dble(ntsets)
                     tmp4=limsup*tmp1-liminf*tmp2
                     mu2real=mu2real+exp(logprioro)*(tmp3-tmp4)

                     if(evali.eq.i)then
                        thetac=rtnorm(mu,sqrt(sigma2),liminf,
     &                         limsup,.false.,.false.) 
                     end if                     
                  end if

c                  call dblepr("liminf",-1,liminf,1)
c                  call dblepr("limsup",-1,limsup,1)

                  liminf=limsup
                  tmp1=dble(i+1)/dble(ntsets)
                  limsup=invcdfnorm(tmp1,0.d0,1.d0,1,0)

               end do

               randsave(isave,nsubject+1)=thetac

               sigma2real=sigma2*(mu2real-mureal*mureal)
               mureal=mu+sqrt(sigma2)*mureal

c               call dblepr("mu",-1,mureal,1)
c               call dblepr("sigma2",-1,sigma2real,1)

               thetasave(isave,p)=mureal
               thetasave(isave,p+1)=sigma2real

c+++++++++++++ baseline mean

               thetasave(isave,p+2)=mu

c+++++++++++++ baseline stansard deviation

               thetasave(isave,p+3)=sigma2

c+++++++++++++ cluster information

               thetasave(isave,p+4)=alpha

c+++++++++++++ random effects

               do i=1,nsubject
                  randsave(isave,i)=b(i)
               end do

c+++++++++++++ density

               do i=1,ngrid
                  tmp1=(grid(i)-mu)/sqrt(sigma2) 
                  if(tmp1.gt.4.0d0)then
                     tmp2=0.999968d0
                    else if(tmp1.lt.-4.0d0)then
                     tmp2=0.000032d0
                    else 
                     tmp2=cdfnorm(tmp1,0.d0,1.d0,1,0)
                  end if

                  kphi=int(real(2**maxm)*tmp2+1)
                  logprioro=dnrm(grid(i),mu,sqrt(sigma2),1)     
                  logprioro=logprioro+dble(maxm)*log(2.d0)
                  logprioro=logprioro+log(prob(kphi))
                  densave(isave,i)=exp(logprioro)

                  logprioro=log(prob(kphi))
                  if(kphi.eq.1)then
                     tmp2=exp(logprioro)*
     &                    dble(ntsets)*
     &                    cdfnorm(grid(i),mu,sqrt(sigma2),1,0)
                    else
                     tmp1=0.d0
                     do j=1,kphi-1
                        tmp1=tmp1+prob(j)    
                     end do

                     tmp2=exp(logprioro)*
     &                    (dble(ntsets)*
     &                     cdfnorm(grid(i),mu,sqrt(sigma2),1,0)-
     &                     dble(kphi-1))
 
                     tmp2=tmp1+tmp2                    
                  end if  

                  cdfsave(isave,i)=tmp2

               end do


c+++++++++++++ cpo
               do i=1,nsubject
                  tmp2=0.d0
                  do j=1,p
                     yij=y(i,j)
                     if(j.eq.1)then
                       eta=b(i)+roffset(i,j)
                      else
                       eta=b(i)-beta(j-1)+roffset(i,j)
                     end if  
                     mean=exp(eta)
                     tmp1=dpoiss(dble(yij),mean,1)
                     cpo(i,j)=cpo(i,j)+1.d0/exp(tmp1)
                     tmp2=tmp2+tmp1    
                  end do
                  cpov(i)=cpov(i)+1.d0/exp(tmp2)
               end do

c+++++++++++++ Savage-Dickey ratio
               i1=0
               logdenom=0.d0
               lognumer=0.d0

               do i=1,maxm
                  j1=2**(i-1)
                  do j=1,j1
                     k1=i1+j

                     ii=(k1-1)*2+1
                     jj=(k1-1)*2+2
  
                     n1=ptcount(ii)
                     n2=ptcount(jj)

                     tmp2=alpha*dble(i**2)

                     if(i.eq.1)then
                        logdenom=0.d0
                        lognumer=0.d0
                       else 
                        logdenom=logdenom+
     &                   dbet(0.5d0,tmp2+dble(n1),tmp2+dble(n2),1)

                        lognumer=lognumer+
     &                   dbet(0.5d0,tmp2,tmp2,1)
                     end if
                  end do
                  i1=i1+j1
               end do
               numerBF=numerBF+exp(lognumer)
               denomBF=denomBF+exp(logdenom)

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
      
      do i=1,nsubject
         do j=1,p
            cpo(i,j)=dble(nsave)/cpo(i,j)
         end do 
         cpov(i)=dble(nsave)/cpov(i)  
      end do

      numerBF=numerBF/dble(nsave)
      denomBF=denomBF/dble(nsave)
      bf=numerBF/denomBF 

      return
      end

c=================================================================================
      subroutine llikaslicept(maxm,ntsets,prob,nsubject,b,mu,sigma2,out)
c=================================================================================
c     A.J., 2009
c=================================================================================
      implicit none
      integer maxm,ntsets,nsubject
      double precision prob(ntsets),b(nsubject)
      double precision mu,sigma2

      integer i,kphi
      double precision dnrm,cdfnorm,tmp1,tmp2

      double precision out

      out=0.d0
      do i=1,nsubject
c+++++++ check if the user has requested an interrupt
         call rchkusr()

         tmp1=(b(i)-mu)/sqrt(sigma2) 
         if(tmp1.gt.4.0d0)then
            tmp2=0.999968d0
          else if(tmp1.lt.-4.0d0)then
            tmp2=0.000032d0
          else 
            tmp2=cdfnorm(tmp1,0.d0,1.d0,1,0)
         end if
         kphi=int(real(2**maxm)*tmp2+1)
         out=out+dnrm(b(i),mu,sqrt(sigma2),1)     
         out=out+log(prob(kphi))
      end do
      return
      end


c=================================================================================
      subroutine targetreslice(i,nsubject,p,y,roffset,maxm,ntsets,prob,
     &                         mu,sigma2,bi,beta,out)
c=================================================================================
c     A.J., 2009
c=================================================================================
      implicit none
      integer i
      integer nsubject,p,y(nsubject,p)
      integer maxm,ntsets
      double precision roffset(nsubject,p) 
      double precision mu,sigma2
      double precision bi,beta(p-1)
      double precision prob(ntsets)
 
      integer j,yij,kphi
      double precision cdfnorm,dnrm,dpoiss,eta,mean,tmp1,tmp2

      double precision out

      out=0.d0
      do j=1,p
         if(j.eq.1)then
            eta=bi+roffset(i,j)
           else
            eta=bi-beta(j-1)+roffset(i,j)
         end if
               
         yij=y(i,j)
         mean=exp(eta)
         out=out+dpoiss(dble(yij),mean,1)               
      end do

      tmp1=(bi-mu)/sqrt(sigma2) 
      if(tmp1.gt.4.0d0)then
         tmp2=0.999968d0
        else if(tmp1.lt.-4.0d0)then
         tmp2=0.000032d0
        else 
         tmp2=cdfnorm(tmp1,0.d0,1.d0,1,0)
      end if
      kphi=int(real(2**maxm)*tmp2+1)
      out=out+dnrm(bi,mu,sqrt(sigma2),1)     
      out=out+log(prob(kphi))
 
      return
      end

c=================================================================================
      subroutine targetcslice(aa0,ab0,alpha,maxm,ntprob,ntsets,
     &                        ptcount,probmat,out)
c=================================================================================
c     A.J., 2009
c=================================================================================
      implicit none
      integer maxm,ntprob,ntsets
      integer ptcount(ntprob)
      double precision aa0,ab0,alpha
      double precision probmat(ntsets,maxm)

      integer i,i1,j,j1,k1,ii,jj,n1,n2
      double precision tmp1,tmp1c,tmp2,tmp4,lbetaf

      double precision out

      out=(aa0-1.d0)*log(alpha)-
     &     ab0*alpha

      i1=0
      do i=1,maxm
         j1=2**(i-1)
         do j=1,j1
            k1=i1+j

            ii=(k1-1)*2+1
            jj=(k1-1)*2+2
  
            n1=ptcount(ii)
            n2=ptcount(jj)

            ii=(j-1)*2+1
            jj=(j-1)*2+2
            tmp1=probmat(ii,i)
            tmp1c=probmat(jj,i)

            tmp2=alpha*2.d0**i
            tmp4=(tmp2+real(n1)-1.d0)*log(tmp1)
            tmp4=tmp4+(tmp2+real(n2)-1.d0)*log(tmp1c)
            tmp4=tmp4-lbetaf(tmp2+dble(n1),tmp2+dble(n2))

            out=out+tmp4
         end do
         i1=i1+j1
      end do

c      call dblepr("out",-1,out,1)
      return
      end
  

c=================================================================================
      subroutine targetrashpd(j,nsubject,p,y,b,beta,
     &                       roffset,prec,b0,out)
c=================================================================================
c     A.J., 2009
c=================================================================================
      implicit none
      integer j,nsubject,p
      integer y(nsubject,p)
      double precision roffset(nsubject,p)
      double precision b(nsubject)
      double precision beta(p-1) 
      double precision prec(p-1,p-1)
      double precision b0(p-1) 

      integer yij,i,jj 
      double precision dpoiss,eta,mean,tmp1

      double precision out

      out=0.d0         
         
      do i=1,nsubject

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         yij=y(i,j+1)            
         eta=b(i)-beta(j)+roffset(i,j+1) 
         mean=exp(eta)
         out=out+dpoiss(dble(yij),mean,1)               
      end do

      tmp1=0.d0
      do i=1,p-1
         do jj=1,p-1
            tmp1=tmp1+(beta(i)-b0(i))* 
     &               prec(i,jj)      *
     &              (beta(jj)-b0(jj))

         end do
      end do
      out=out-0.5d0*tmp1
   
      return
      end

