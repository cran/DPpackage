c=======================================================================                      
      subroutine ptglmmgam
     &                (datastr,maxni,nrec,nsubject,nfixed,p,q,subject,
     &                 x,y,z,
     &                 a0b0,mu0,prec1,prec2,sb,tinv,
     &                 mcmc,nsave,
     &                 acrate,cpo,randsave,thetasave,
     &                 curr,b,beta,betar,mu,sigma,ortho,
     &                 betac,
     &                 iflagp,workmp1,workmp2,workmhp1,workvp1,
     &                 xtx,xty,
     &                 iflagr,parti,whicho,whichn,bz,bzc,limw,linf,lsup,
     &                 propvr,sigmainv,theta,thetac,
     &                 workmhr,workmr,workmr1,workmr2,workvr,ybar,
     &                 sigmac,sigmainvc,workmhr2,
     &                 massi,pattern,betasave,bsave)
c=======================================================================                      
c     # of arguments = 64.
c
c     Subroutine `ptglmmgam' to run a Markov chain in a semiparametric 
c     gamma mixed effect model, using a Mixture of Multivariate Polya 
c     trees prior for the distribution of the random effects.
c
c     Copyright: Alejandro Jara and Timothy Hanson, 2007-2010.
c
c     The parametrization considered here is:
c     log p(Y)= v* ((-1/mu)*Y - log(mu)) + c(Y,v),
c     such that:
c              - E(Y)  =mu
c              - Var(Y)= (1/v)*(mu^2)
c
c     Note that the commonly used parametrization is phi=(1/v). 
c     A Gamma(tau1/2,tau2/2) is specified on v
c
c     Version 1.0:
c
c     Last modification: 01-07-2008.
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
c     The authors' contact information:
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
c        datastr     :  integer matrix giving the number of measurements
c                       and the location in y of the observations for 
c                       each subject, datastr(nsubject,maxni+1)
c        maxni       :  integer giving the maximum number of 
c                       measurements for subject.
c        nrec        :  integer giving the number of observations.
c        nsubject    :  integer giving the number of subjects.
c        nfixed      :  integer giving the number of fixed effects,
c                       if nfixed is 0 then p=1.
c        p           :  integer giving the number of fixed coefficients.
c        q           :  integer giving the number of random effects.
c        subject     :  integer vector giving the subject for each.
c                       observation, subject(nsubject).
c        x           :  real matrix giving the design matrix for the 
c                       fixed effects (also offset), x(nrec,p+1). 
c        y           :  real vector giving the response variable,
c                       y(nrec).
c        z           :  real matrix giving the design matrix for the 
c                       random effects, z(nrec,q). 
c
c-----------------------------------------------------------------------
c
c---- Prior information ------------------------------------------------
c 
c        aa0, ab0    :  real giving the hyperparameters of the prior
c                       distribution for the precision parameter,
c                       alpha ~ Gamma(aa0,ab0). If aa0<0 the precision 
c                       parameter is considered as a constant.
c        b0          :  real vector giving the prior mean
c                       for the fixed effects, b0(p)
c        m           :  integer giving the number of binary partitions
c                       in each margin of the Multivariate
c                       Polya tree prior.
c        mu0         :  real vector giving the prior mean
c                       for the baseline mean, mu0(q).
c        nu0         :  integer giving the degrees of freedom for the
c                       inverted-Wishart prior distribution for the
c                       covariance matrix of the random effects
c                       (This is for the base line).
c        prec1       :  real matrix giving the prior precision matrix
c                       for the fixed effects, prec1(p,p).
c        prec2       :  real matrix giving the prior precision matrix
c                       for the baseline mean, prec2(q,q).
c        sb          :  real vector giving the product of the prior 
c                       precision and prior mean for the fixed effects
c                       sb(p,2).
c        tinv        :  real matrix giving the scale matrix for the
c                       inverted-Wishart prior distribution for the
c                       covariance matrix of the random effects, 
c                       sigma ~ Inv-Wishart(nu0,tinv^{-1}), such that 
c                       E(sigma)=(1/(nu0-q-1)) * tinv 
c                       (This is for the base line distribution).
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
c        samplef     :  integer indicating whether the functionals
c                       must be sampled (1) or not (0).          
c        
c-----------------------------------------------------------------------
c
c---- Output -----------------------------------------------------------
c
c        acrate      :  real vector giving the MH acceptance rate, 
c                       acrate(6). 
c        cpo         :  real giving the cpo. 
c        randsave    :  real matrix containing the mcmc samples for
c                       the random effects and prediction,
c                       randsave(nsave,q*(nsubject+1)).
c        thetasave   :  real matrix containing the mcmc samples for
c                       the averaged random effects, fixed effects, 
c                       error variance, and mean and covariance of
c                       the baseline distribution, 
c                       thetsave(nsave,q+nfixed+q+2*nuniq(Sigma)+2+q*q).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        cpar        :  real giving the current value of the precision
c                       parameter of the PT.
c        disp        :  real giving the value of the dispersion 
c                       parameter.
c        b           :  real matrix giving the current value of the 
c                       random effects, b(nsubject,q).
c        beta        :  real vector giving the current value of the 
c                       fixed effects, beta(p).
c        betar       :  real vector giving the current value of the 
c                       averaged random effects, betar(q).
c        mu          :  real vector giving the mean of the normal 
c                       base line distribution for the random effects,
c                       mu(q).
c        sigma       :  real matrix giving the current value of the
c                       covariance matrix for normal base line 
c                       distribution for the random effects,
c                       sigma(q,q).
c        ortho       :  real matrix giving the current value of the
c                       orthogonal matrix, ortho(q,q).
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        acrate2     :  real working varible
c        betac       :  real vector giving the candidate value of the 
c                       fixed effects, betac(p).
c        bz          :  real matrix giving the current value of the 
c                       standarized random effects, bz(nsubject,q).
c        bzc         :  real matrix giving the candidate value of the 
c                       standarized random effects, bz(nsubject,q).
c        cparc       :  real giving the value of the candidate
c                       for the precision parameter.
c        detlogl     :  real used to save the log-determinant in a
c                       matrix inversion process.
c        detloglc    :  real used to save the log-determinant in a
c                       matrix inversion process.
c        dispc       :  real giving the value of the candidate
c                       for the dispersion parameter.
c        dispcount   :  index. 
c        dlnrm       :  density of a log-normal distribution.
c        i           :  index. 
c        ii          :  index. 
c        iflagp      :  integer vector used to invert the of the lhs
c                       least square solution for the fixed effects,
c                       iflagp(p).
c        iflagr      :  integer vector used to invert the of the lhs
c                       least square solution for the random effects,
c                       iflagr(q).
c        ihmssf      :  integer function to evaluate the position of an
c                       element in a matrix based on a half-stored 
c                       version.
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        k           :  index. 
c        l           :  index.
c        limw        :  real vector giving the limits of partitions, 
c                       limw(q).
c        linf        :  real vector giving the limits of partitions, 
c                       linf(q).
c        logcgko     :  real working variable.
c        logcgkn     :  real working variable.
c        loglikec    :  real working variable.
c        loglikeo    :  real working variable.
c        logpriorc   :  real working variable.
c        logprioro   :  real working variable.
c        lsup        :  real vector giving the limits of partitions, 
c                       lsup(q).
c        massi       :  integer vector giving the number of RE
c                       in each element of the partition, massi(2**q).
c        narea       :  integer giving the total number of areas per 
c                       partition, narea=2**q.
c        ni          :  integer indicator used in updating the state. 
c        nscan       :  integer indicating the total number of MCMC
c                       scans.
c        nu          :  real working variable. 
c        parti       :  integer vector giving the partition,
c                       parti(q). 
c        pattern     :  integer vector giving the pattern of an observation,
c                       pattern(q). 
c        propvr      :  real matrix used to update the random effects,
c                       propvr(q,q).
c        ratio       :  real working variable.
c        rgamma      :  real gamma random number generator.
c        rtlnorm     :  real truncated log normal random number generator.
c        runif       :  real uniform random number generator.
c        sec         :  cpu time working variable.
c        sec0        :  cpu time working variable.
c        sec00       :  cpu time working variable.
c        sec1        :  cpu time working variable.
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        sigmac      :  real matrix giving the candidate value of the
c                       baseline covariance matrix, sigmac(q,q).
c        sigmainv    :  real matrix giving the inverse of the current
c                       value of the baseline covariance matrix, 
c                       sigmainv(q,q).
c        sigmainvc   :  real matrix giving the inverse of the candidate
c                       value of the baseline covariance matrix, 
c                       sigmainvc(q,q).
c        skipcount   :  index. 
c        sprint      :  integer function to print on screen.
c        sse         :  real used to save the SS of the errors.
c        theta       :  real vector used to save randomnly generated
c                       random effects, theta(q).
c        thetac      :  real vector used to save randomnly generated
c                       random effects, thetac(q).
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
c        whicho      :  integer vector giving the rand. eff. in each
c                       partition, whicho(nsubject).
c        whichn      :  integer vector giving the rand. eff. in each
c                       partition, whichn(nsubject).
c        workmhp1    :  real vector used to update the fixed effects,
c                       workmhp1(p*(p+1)/2)
c        workmhr     :  real vector used to update the random effects,
c                       workmhr(q*(q+1)/2)
c        workmhr2    :  real vector used to update the baseline cov,
c                       workmhr2(q*(q+1)/2)
c        workmp1     :  real matrix used to update the fixed effects,
c                       workmp1(p,p).
c        workmp2     :  real matrix used to update the fixed effects,
c                       workmp2(p,p).
c        workmr      :  real matrix used to update the random effects,
c                       workmr(q,q).
c        workmr1     :  real matrix used to update the random effects,
c                       workmr1(q,q).
c        workmr2     :  real matrix used to update the random effects,
c                       workmr2(q,q).
c        workvp1     :  real vector used to update the fixed effects,
c                       workvp1(p)
c        workvr      :  real vector used to update the random effects,
c                       workvr(p)
c        xtx         :  real matrix givind the product X^tX, xtx(p,p).
c        xty         :  real vector used to save the product 
c                       Xt(Y-Zb), xty(p).
c        ybar        :  real vector used to save the mean of random  
c                       effects and the probabilities in each 
c                       partition area, ybar(2**q).
c
c=======================================================================                  
      implicit none 

c+++++Data
      integer maxni,nrec,nsubject,nfixed,p,q,subject(nrec)
      integer datastr(nsubject,maxni+1)
      double precision x(nrec,p+1),z(nrec,q),xtx(p,p),y(nrec)
      
c+++++Prior 
      integer fixed,m,murand,sigmarand,typepr
      double precision aa0,ab0,a0b0(10),mu0(q),nu0,prec1(p,p),prec2(q,q)
      double precision sb(p,2)
      double precision tau1,tau2,tinv(q,q)      

c+++++MCMC parameters
      integer mcmc(12),nburn,nskip,nsave,ndisplay,nbase,samplef
      double precision tune1,tune2,tune3,tune4,tune5

c+++++Output
      double precision acrate(7),cpo(nrec,2)
      double precision randsave(nsave,q*(nsubject+1))
      double precision thetasave(nsave,q+nfixed+q+(q*(q+1))+2+q*q)

c+++++Current values of the parameters
      double precision cpar,curr(7),disp,beta(p),b(nsubject,q)
      double precision betar(q)
      double precision mu(q),sigma(q,q)
      double precision ortho(q,q)

c+++++Working space - External
      integer iflagp(p) 
      double precision betac(p)
      double precision workmp1(p,p),workmp2(p,p)
      double precision workmhp1(p*(p+1)/2)
      double precision workvp1(p)
      double precision xty(p)

      integer iflagr(q) 
      integer parti(q)
      integer whicho(nsubject),whichn(nsubject)      
      double precision bz(nsubject,q),bzc(nsubject,q)
      double precision limw(q),linf(q),lsup(q)
      double precision propvr(q,q)
      double precision sigmainv(q,q)
      double precision theta(q),thetac(q)
      double precision workmhr(q*(q+1)/2)
      double precision workmr(q,q)
      double precision workmr1(q,q),workmr2(q,q)
      double precision workvr(q)
      double precision ybar(2**q)
      
      double precision sigmac(q,q),sigmainvc(q,q)
      double precision workmhr2(q*(q+1)/2)

      integer massi(2**q)
      integer pattern(q)

c++++ model's performance
      double precision mc(5)
      double precision betasave(p+1),bsave(nsubject,q)

c+++++Working space - RNG
      integer seed1,seed2
      
c+++++Working space - Internal
      integer baseskip
      integer dispcount
      integer i,j,k,l
      integer ihmssf
      integer iscan,isave
      integer narea,ni,nscan,nu
      integer skipcount
      integer sprint
      double precision acrate2,cparc
      double precision detlogl,detloglc,dispnew,dlnrm
      double precision eta,gprime
      double precision mean
      double precision logcgkn,logcgko
      double precision loglikn,logliko
      double precision logpriorn,logprioro
      double precision offset
      double precision ratio,rnorm,rtlnorm,runif
      double precision slogmu,slogy,symu
      double precision targetp
      double precision tmp1,tmp2
      double precision trigamm
      double precision yij,ytilde

c++++ model's performance
      double precision dbarc,dbar,dhat,pd,lpml

c+++++CPU time
      double precision sec00,sec0,sec1,sec

c+++++Adaptive MH for mu
      integer countermh
      double precision mumh(100),sigmamh(100,100)

c+++++Adaptive MH for sigma
      integer nadaptive
      parameter(nadaptive=2000)
      integer adaptives,sigmaskip
      double precision aratesigma

c+++++Adaptive MH for c
      integer adaptivec,cskip
      double precision aratec

c+++++Adaptive MH for partition
      integer adaptivep,pskip
      double precision aratep

c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      nbase=mcmc(4)      
      m=mcmc(5)
      seed1=mcmc(6)
      seed2=mcmc(7)
      typepr=mcmc(8)
      murand=mcmc(9)
      sigmarand=mcmc(10)
      fixed=mcmc(11)
      samplef=mcmc(12)

      aa0=a0b0(1)
      ab0=a0b0(2)
      nu0=a0b0(3)
      tune1=a0b0(4)
      tune2=a0b0(5)
      tune3=a0b0(6)
      tune4=a0b0(7)
      tune5=a0b0(8)
      tau1=a0b0(9)
      tau2=a0b0(10)

      cpar=curr(1)
      disp=curr(2)

      narea=2**q
      
c++++ set random number generator
      call setall(seed1,seed2)

c++++ transforming random effects and calculate log-likelihood
c++++ for the baseline covariance matrix

      logliko=0.d0

      call rhaar2(workmr,ortho,q,workmr1)

      do i=1,q
         mumh(i)=0.d0
         do j=1,q
            sigmamh(i,j)=0.d0
            workmr(i,j)=sigma(i,j)
         end do
      end do
      call inversedet(workmr,q,iflagr,detlogl)

      do i=1,q
         do j=1,q
            workmr(i,j)=0.d0
            sigmainv(i,j)=0.d0
         end do
      end do
      call cholesky(q,sigma,workmhr)
      do i=1,q
         do j=1,i
            workmr(i,j)=workmhr(ihmssf(i,j,q))
         end do
      end do

      do i=1,q
         do j=1,q
            tmp1=0.d0
            do k=1,q
               tmp1=tmp1+workmr(i,k)*workmr1(k,j)
            end do 
            sigmainv(i,j)=tmp1  
         end do
      end do

      call inverse(sigmainv,q,iflagr)      

      call loglikpt_mucan(m,q,nsubject,parti,
     &                    whicho,whichn,b,bz,cpar,detlogl,
     &                    linf,lsup,mu,sigmainv,
     &                    theta,fixed,logliko)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      dbar=0.d0
      isave=0
      skipcount=0
      dispcount=0
      baseskip=0
      countermh=0
      
      adaptives=0
      aratesigma=0.d0
      sigmaskip=0
      if(sigmarand.eq.1.and.tune2.lt.0.d0)then
         adaptives=1
         tune2=10.0
         nburn=nburn+nadaptive
      end if  

      adaptivec=0
      aratec=0.d0
      cskip=0
      if(aa0.gt.0.d0.and.tune3.lt.0.d0)then
         adaptivec=1
         tune3=1.0
         nburn=nburn+nadaptive
      end if  

      adaptivep=0
      aratep=0.d0
      pskip=0
      if(typepr.eq.1.and.tune4.lt.0.d0)then
         adaptivep=1
         tune4=1.d0
         nburn=nburn+nadaptive
      end if  

      
      nscan=nburn+(nskip+1)*(nsave)

      call cpu_time(sec0)
      sec00=0.d0
      
      do iscan=1,nscan


c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++
c+++++++ fixed effects        +++
c++++++++++++++++++++++++++++++++

         if(nfixed.gt.0)then
            do i=1,p
               do j=1,p
                  xtx(i,j)=prec1(i,j)
               end do
               xty(i)=sb(i,1)
            end do

            logliko=0.d0

            do i=1,nrec
               eta=0.d0
               offset=0.d0
               mean=0.d0
               gprime=0.d0
               yij=y(i)
               
               do j=1,p
                  eta=eta+x(i,j)*beta(j)
               end do
               
               do j=1,q
                  eta=eta+z(i,j)*b(subject(i),j) 
                  offset=offset+z(i,j)*b(subject(i),j) 
               end do

               eta=eta+x(i,p+1)
               offset=offset+x(i,p+1)

               mean=exp(eta)
               gprime=exp(-eta)

               ytilde=eta+(yij-mean)*gprime-offset
               
               do j=1,p
                  do k=1,p
                     xtx(j,k)=xtx(j,k)+x(i,j)*x(i,k)
                  end do
                  xty(j)=xty(j)+x(i,j)*ytilde
               end do
               
               call dgamma(y(i),mean,disp,tmp1) 
               logliko=logliko+tmp1
            end do

            call inverse(xtx,p,iflagp)      

            do i=1,p
               tmp1=0.d0
               do j=1,p
                  tmp1=tmp1+xtx(i,j)*xty(j) 
               end do
               workvp1(i)=tmp1
            end do

            call rmvnorm(p,workvp1,xtx,workmhp1,xty,betac)

c++++++++++ evaluating the candidate generating kernel

            call dmvn2(p,betac,workvp1,xtx,logcgko,
     &                 xty,workmp1,workmp2,iflagp)                 
            
c++++++++++ evaluating the likelihood

            do i=1,p
               do j=1,p
                  xtx(i,j)=prec1(i,j)
               end do
               xty(i)=sb(i,1)
            end do

            loglikn=0.d0

            do i=1,nrec
               eta=0.d0
               offset=0.d0
               mean=0.d0
               gprime=0.d0
               yij=y(i)
            
               do j=1,p
                  eta=eta+x(i,j)*betac(j)
               end do
               
               do j=1,q
                  eta=eta+z(i,j)*b(subject(i),j) 
                  offset=offset+z(i,j)*b(subject(i),j) 
               end do

               eta=eta+x(i,p+1)
               offset=offset+x(i,p+1)
               
               mean=exp(eta)
               gprime=exp(-eta)
               
               ytilde=eta+(yij-mean)*gprime-offset
               
               do j=1,p
                  do k=1,p
                     xtx(j,k)=xtx(j,k)+x(i,j)*x(i,k)
                  end do
                  xty(j)=xty(j)+x(i,j)*ytilde
               end do

               call dgamma(y(i),mean,disp,tmp1) 
               loglikn=loglikn+tmp1
            end do

            call inverse(xtx,p,iflagp)      

            do i=1,p
               tmp1=0.d0
               do j=1,p
                  tmp1=tmp1+xtx(i,j)*xty(j) 
               end do
               workvp1(i)=tmp1
            end do

c++++++++++ evaluating the candidate generating kernel

            call dmvn2(p,beta,workvp1,xtx,logcgkn,
     &                 xty,workmp1,workmp2,iflagp)                 

c++++++++++ prior ratio
            logprioro=0.d0
            logpriorn=0.d0
         
            do i=1,p
               do j=1,p
                  logpriorn=logpriorn+(betac(i)-sb(i,2))* 
     &                      prec1(i,j)       *
     &                      (betac(j)-sb(j,2))

                  logprioro=logprioro+(beta(i)-sb(i,2))* 
     &                      prec1(i,j)      *
     &                      (beta(j)-sb(j,2))

               end do
            end do
      
            logpriorn=-0.5d0*logpriorn
            logprioro=-0.5d0*logprioro

c++++++++++ mh step

            ratio=loglikn-logliko+logcgkn-logcgko+
     &            logpriorn-logprioro

            if(log(dble(runif())).lt.ratio)then
               acrate(1)=acrate(1)+1.d0
               do i=1,p
                  beta(i)=betac(i) 
               end do
            end if
         end if            


c+++++++++++++++++++++++++++++++++
c+++++++ random effects        +++ 
c+++++++++++++++++++++++++++++++++

         acrate2=0.d0

         do i=1,q
            ybar(i)=0.d0
            do j=1,q
               workmr(i,j)=sigma(i,j)
            end do
            workvr(i)=0.d0
            iflagr(i)=0
         end do

         call inverse(workmr,q,iflagr)      

         do i=1,nsubject

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            do j=1,q
               theta(j)=b(i,j)
            end do   

c++++++++++ generating a candidate
            do j=1,q
               tmp1=0.d0
               do k=1,q
                  propvr(j,k)=workmr(j,k)
               end do
               workvr(j)=0.d0
               iflagr(j)=0
            end do

            ni=datastr(i,1) 

            logliko=0.d0
            
            do j=1,ni
               eta=0.d0 
               mean=0.d0
               gprime=0.d0
               offset=0.d0

               yij=y(datastr(i,j+1))
            
               do k=1,p
                  eta=eta+x(datastr(i,j+1),k)*beta(k)
                  offset=offset+x(datastr(i,j+1),k)*beta(k)
               end do
            
               do k=1,q
                  eta=eta+z(datastr(i,j+1),k)*theta(k)
               end do

               eta=eta+x(datastr(i,j+1),p+1)
               offset=offset+x(datastr(i,j+1),p+1)
               
               mean=exp(eta)
               gprime=exp(-eta)
               
               ytilde=eta+(yij-mean)*gprime-offset
            
               do k=1,q
                  do l=1,q
                     propvr(k,l)=propvr(k,l)+
     &                        z(datastr(i,j+1),k)*
     &                        z(datastr(i,j+1),l)
                  end do
               end do

               call dgamma(yij,mean,disp,tmp1) 
               logliko=logliko+tmp1
            end do

            call inverse(propvr,q,iflagr)      

            call rmvnorm(q,theta,propvr,workmhr,workvr,thetac)

c++++++++++ evaluating the candidate generating kernel

            call dmvn2(q,thetac,theta,propvr,logcgkn,
     &                 workvr,workmr1,workmr2,iflagr)

c++++++++++ evaluating the likelihood

            do j=1,q
               tmp1=0.d0
               do k=1,q
                  propvr(j,k)=workmr(j,k)
               end do
               workvr(j)=0.d0
               iflagr(j)=0
            end do

            ni=datastr(i,1) 

            loglikn=0.d0
            
            do j=1,ni
               eta=0.d0 
               mean=0.d0
               gprime=0.d0
               offset=0.d0

               yij=y(datastr(i,j+1))            
            
               do k=1,p
                  eta=eta+x(datastr(i,j+1),k)*beta(k)
                  offset=offset+x(datastr(i,j+1),k)*beta(k)
               end do
            
               do k=1,q
                  eta=eta+z(datastr(i,j+1),k)*thetac(k)
               end do

               eta=eta+x(datastr(i,j+1),p+1)
               offset=offset+x(datastr(i,j+1),p+1)
               
               mean=exp(eta)
               gprime=exp(-eta)
               
               ytilde=eta+(yij-mean)*gprime-offset
            
               do k=1,q
                  do l=1,q
                     propvr(k,l)=propvr(k,l)+
     &                        z(datastr(i,j+1),k)*
     &                        z(datastr(i,j+1),l)
                  end do
               end do

               call dgamma(yij,mean,disp,tmp1) 
               loglikn=loglikn+tmp1
            end do

            call inverse(propvr,q,iflagr)      

c++++++++++ evaluating the candidate generating kernel

            call dmvn2(q,theta,thetac,propvr,logcgkn,
     &                 workvr,workmr1,workmr2,iflagr)


c++++++++++ evaluating the prior

            logprioro=0.d0
            logpriorn=0.d0

            do j=1,q
               limw(j)=bz(i,j)
            end do

            call condptprior(limw,i,nsubject,q,bz,cpar,m,detlogl,
     &                       linf,lsup,parti,whicho,whichn,
     &                       fixed,logprioro)

            do j=1,q
               tmp1=0.d0
               do k=1,q
                  tmp1=tmp1+sigmainv(j,k)*(thetac(k)-mu(k))   
               end do
               limw(j)=tmp1
            end do

            call condptprior(limw,i,nsubject,q,bz,cpar,m,detlogl,
     &                       linf,lsup,parti,whicho,whichn,
     &                       fixed,logpriorn)

c++++++++++ mh step
  
            ratio=loglikn-logliko+
     &            logcgkn-logcgko+
     &            logpriorn-logprioro

            if(log(dble(runif())).lt.ratio)then
               acrate2=acrate2+1.d0
               do j=1,q
                  b(i,j)=thetac(j)
                  bz(i,j)=limw(j)
               end do
            end if

         end do

         acrate(2)=acrate(2)+acrate2/dble(nsubject)


c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Updating Dispersion parameter using a MH step      +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         slogy=0.d0
         symu=0.d0
         slogmu=0.d0

         do i=1,nrec
             eta=0.d0
             mean=0.d0
             do j=1,p
                eta=eta+x(i,j)*beta(j)
             end do
               
             do j=1,q
                eta=eta+z(i,j)*b(subject(i),j) 
             end do
             
             eta=eta+x(i,p+1)                    
             
             mean=exp(eta)
         
             slogy=slogy+log(y(i))
             symu=symu+y(i)/mean
             slogmu=slogmu+eta
         end do

c+++++++ generating a candidate

         tmp1=1.d0/((disp**2)*dble(nrec)*(trigamm(disp)-1.d0/disp))
         tmp2=sqrt(tune5*tmp1)

         dispnew=rtlnorm(log(disp),tmp2,0,0,.true.,.true.)

c+++++++ evaluating the candidate generating kernel

         logcgko=dlnrm(dispnew,log(disp),tmp2,1) 

         tmp1=1.d0/
     &        ((dispnew**2)*dble(nrec)*(trigamm(dispnew)-1.d0/dispnew))
         tmp2=sqrt(tune5*tmp1)

         logcgkn=dlnrm(disp,log(dispnew),tmp2,1) 

c+++++++ evaluating the full conditional

         loglikn=targetp(nrec,tau1,tau2,dispnew,slogy,symu,slogmu)
         logliko=targetp(nrec,tau1,tau2,disp   ,slogy,symu,slogmu)  
         
         ratio=loglikn-logliko+logcgkn-logcgko

         if(log(dble(runif())).lt.ratio)then
            acrate(3)=acrate(3)+1.d0
            disp=dispnew
         end if


         baseskip = baseskip + 1
         if(baseskip.ge.nbase)then
         countermh=countermh+1

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Updating mu using a MH step                        +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++ update the log-likelihood for random effects

         if(murand.eq.1.or.sigmarand.eq.1.or.aa0.gt.0.d0.or.
     &      typepr.eq.1)then            
            call loglikpt_updatet(m,q,nsubject,parti,
     &                         whicho,whichn,bz,cpar,detlogl,
     &                         linf,lsup,
     &                         fixed,logliko)
         end if

         if(murand.eq.1)then

c++++++++++ generating a candidate

            if(tune1.gt.0.d0)then 

               do i=1,q
                  do j=1,q
                     propvr(i,j)=(tune1)*sigma(i,j)/dble(nsubject)
                  end do
               end do
               call rmvnorm(q,mu,propvr,workmhr,workvr,theta)

             else  
               if(iscan.le.nburn/2)then
                  do i=1,q
                     do j=1,q
                        propvr(i,j)=(0.01d0)*sigma(i,j)/dble(nsubject)
                     end do
                  end do
                  call rmvnorm(q,mu,propvr,workmhr,workvr,theta)
                else
                  ratio=dble(runif())
                  if(ratio.le.0.25)then
                     do i=1,q
                        do j=1,q
                           propvr(i,j)=(5.4264d0/dble(q))*sigmamh(i,j)
                        end do
                     end do
                     call rmvnorm(q,mu,propvr,workmhr,workvr,theta)
                   else if(ratio.le.0.5)then
                     do i=1,q
                        do j=1,q
                          propvr(i,j)=sigma(i,j)/dble(nsubject)
                        end do
                     end do
                     call rmvnorm(q,mu,propvr,workmhr,workvr,theta)
                   else
                     do i=1,q
                        do j=1,q
                          propvr(i,j)=0.01d0*sigma(i,j)/dble(nsubject)
                        end do
                     end do
                     call rmvnorm(q,mu,propvr,workmhr,workvr,theta)
                  end if  
               end if  
            end if

c++++++++++ evaluating priors 

            logpriorn=0.d0
            logprioro=0.d0
         
            do i=1,q
               do j=1,q
                  logpriorn=logpriorn+(theta(i)-mu0(i))* 
     &                                 prec2(i,j)    *
     &                                (theta(j)-mu0(j))

                  logprioro=logprioro+(mu(i)-mu0(i))* 
     &                                 prec2(i,j)    *
     &                                (mu(j)-mu0(j))

               end do
            end do
         
            logpriorn=-0.5d0*logpriorn
            logprioro=-0.5d0*logprioro

c++++++++++ evaluating likelihood for muc

            call loglikpt_mucan(m,q,nsubject,parti,
     &                          whicho,whichn,b,bzc,cpar,detlogl,
     &                          linf,lsup,theta,sigmainv,
     &                          thetac,fixed,loglikn)


c++++++++++ acceptance step

            ratio=loglikn-logliko+logpriorn-logprioro

            if(log(dble(runif())).lt.ratio)then
               do i=1,q
                  mu(i)=theta(i)
               end do
               do i=1,nsubject
                  do j=1,q
                     bz(i,j)=bzc(i,j)
                  end do   
               end do
               logliko=loglikn
               acrate(4)=acrate(4)+1.d0
            end if
            
c++++++++++ addapting the parameters for the MH algorithm
            if(countermh.eq.1)then
               do i=1,q
                  mumh(i)=mu(i)
                  do j=1,q
                     sigmamh(i,j)=sigma(i,j)/dble(nsubject)
                  end do
               end do
             else
               do i=1,q
                  do j=1,q
                     sigmamh(i,j)=sigmamh(i,j)+(1.d0/dble(countermh))*
     &               ( 
     &                (mu(i)-mumh(i))*(mu(j)-mumh(j)) - 
     &                (dble(countermh)/dble(countermh-1))*sigmamh(i,j)
     &               ) 
                  end do
               end do
            
               do i=1,q
                  mumh(i)=mumh(i)+(1.d0/dble(countermh))*(mu(i)-mumh(i))
               end do
            end if
            
         end if
         
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Updating sigma using a MH step                     +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(sigmarand.eq.1)then
         
c++++++++++ Addaptive MH

            if(adaptives.eq.1)then  
               sigmaskip = sigmaskip + 1
               if(sigmaskip.eq.100)then
                  aratesigma=aratesigma/dble(100)

                  if(iscan.le.nadaptive)then  
                     if(q.eq.1)then
                        if(aratesigma.lt.0.44)then
                           tune2=exp(log(tune2)+(0.44-aratesigma))
                         else
                           tune2=exp(log(tune2)-(aratesigma-0.44))
                        end if  
                       else
                        if(aratesigma.lt.0.234)then
                           tune2=exp(log(tune2)+(0.234-aratesigma))
                         else
                           tune2=exp(log(tune2)-(aratesigma-0.234))
                        end if  
                     end if  

                   else 
                     if(q.eq.1)then
                        if(aratesigma.lt.0.44)then
                           tune2=exp(log(tune2)+
     &                 min(0.01d0,1.d0/sqrt(dble(iscan-nadaptive))))
                          else
                           tune2=exp(log(tune2)-
     &                 min(0.01d0,1.d0/sqrt(dble(iscan-nadaptive))))
                        end if  
                       else
                        if(aratesigma.lt.0.234)then
                           tune2=exp(log(tune2)+
     &                 min(0.01d0,1.d0/sqrt(dble(iscan-nadaptive))))
                          else
                           tune2=exp(log(tune2)-
     &                 min(0.01d0,1.d0/sqrt(dble(iscan-nadaptive))))
                        end if  
                     end if  
                  end if
                  
                  nu=(dble(nsubject))*tune2
                  if(nu.le.(q+1))tune2=dble(q+2)/dble(nsubject)
                  
                  sigmaskip=0
                  aratesigma=0.d0
               end if
            end if

c++++++++++ START: Simple MH
c++++++++++ generating the candidate value

            nu=(dble(nsubject))*tune2
         
            do i=1,q
               do j=1,q
                  sigmac(i,j)=dble(nu-q-1)*sigma(i,j)
               end do
            end do

            call riwishart(q,nu,sigmac,workmr2,workmr,workvr,
     &                     workmhr,workmhr2,iflagr)


c++++++++++ evaluating the candidate generating kernel

            do i=1,q
               do j=1,q
                  propvr(i,j)=dble(nu-q-1)*sigma(i,j)
               end do
            end do

            call diwishart(q,nu,sigmac,propvr,workmr1,workmr2,workvr,
     &                     iflagr,logcgko)        

            do i=1,q
               do j=1,q
                  propvr(i,j)=dble(nu-q-1)*sigmac(i,j)
               end do
            end do

            call diwishart(q,nu,sigma,propvr,workmr1,workmr2,workvr,
     &                     iflagr,logcgkn)        
   
c++++++++++ ENDS: Simple MH

c++++++++++ evaluating the prior

            call diwishart(q,int(nu0),sigmac,tinv,workmr1,workmr2,
     &                     workvr,iflagr,logpriorn)        

            call diwishart(q,int(nu0),sigma,tinv,workmr1,workmr2,
     &                     workvr,iflagr,logprioro)        

c++++++++++ evaluating likelihood for sigmac

            call rhaar2(workmr,ortho,q,workmr1)

            call loglikpt_covarcan2(m,q,nsubject,iflagr,parti,
     &                              whicho,whichn,b,bzc,cpar,detloglc,
     &                              linf,lsup,mu,sigmac,sigmainvc,
     &                              workmr1,workvr,workmhr,workmr,
     &                              loglikn,fixed)

c++++++++++ acceptance step
         
            ratio=loglikn-logliko+logcgkn-logcgko+
     &            logpriorn-logprioro

            if(log(dble(runif())).lt.ratio)then
               do i=1,q
                  do j=1,q
                     sigma(i,j)=sigmac(i,j)
                     sigmainv(i,j)=sigmainvc(i,j)
                  end do
               end do
               do i=1,nsubject
                  do j=1,q
                     bz(i,j)=bzc(i,j)
                  end do   
               end do 
               detlogl=detloglc
               logliko=loglikn
               acrate(5)=acrate(5)+1.d0
               
               if(adaptives.eq.1)aratesigma=aratesigma+1.d0
               
            end if
         end if

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ MH to update the c parameter                 +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(aa0.gt.0.d0)then

c++++++++++ Addaptive MH

            if(adaptivec.eq.1)then  
               cskip = cskip + 1
               if(cskip.eq.100)then
                  aratec=aratec/dble(100)
                  if(iscan.le.nadaptive)then  
                     if(aratec.lt.0.44)then
                        tune3=exp(log(tune3)+(0.44-aratec))
                      else
                        tune3=exp(log(tune3)-(aratec-0.44))
                     end if  
                   else 
                     if(aratec.gt.0.44)then
                        tune3=exp(log(tune3)+
     &                        min(0.01d0,1.d0/sqrt(dble(iscan))))
                       else
                        tune3=exp(log(tune3)-
     &                        min(0.01d0,1.d0/sqrt(dble(iscan))))
                     end if 
                  end if    
                  cskip=0
                  aratec=0.d0
               end if
            end if

c++++++++++ sample candidates

            cparc=rtlnorm(log(cpar),tune3*1.0,0,0,.true.,.true.)
            logcgkn=dlnrm(cpar ,log(cparc),tune3*1.0,1) 
            logcgko=dlnrm(cparc,log(cpar ),tune3*1.0,1) 

c++++++++++ evaluate log-prior for candidate value of the parameters

            call dgamma2(cparc,aa0,ab0,logpriorn)  

c++++++++++ evaluate log-prior for current value of parameters

            call dgamma2(cpar ,aa0,ab0,logprioro)


c++++++++++ evaluate log-likelihood for candidate value of 
c++++++++++ the parameters

            call loglikpt_cparcan(m,q,nsubject,iflagr,parti,
     &                            whicho,whichn,bz,cparc,detlogl,
     &                            linf,lsup,
     &                            theta,fixed,loglikn)

c++++++++++ acceptance step
            ratio=loglikn+logpriorn-logliko-logprioro+
     &            logcgkn-logcgko

            if(log(dble(runif())).lt.ratio)then
               cpar=cparc
               acrate(6)=acrate(6)+1.d0
               logliko=loglikn

               if(adaptivec.eq.1)aratec=aratec+1.d0
               
            end if            
         end if

         baseskip=0
         end if

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Updating the partition                       +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(typepr.eq.1)then

c++++++++++ Addaptive MH

            if(adaptivep.eq.1)then  
               pskip = pskip + 1
               if(pskip.eq.100)then
                  aratep=aratep/dble(100)
                  if(iscan.le.nadaptive)then  
                     if(aratep.lt.0.234)then
                        tune4=exp(log(tune4)+(0.234-aratep))
                      else
                        tune4=exp(log(tune4)-(aratep-0.234))
                     end if  
                   else 
                     if(aratep.gt.0.234)then
                        tune4=exp(log(tune4)+
     &                        min(0.01d0,1.d0/sqrt(dble(iscan))))
                       else
                        tune4=exp(log(tune4)-
     &                        min(0.01d0,1.d0/sqrt(dble(iscan))))
                     end if 
                  end if    
                  pskip=0
                  aratep=0.d0
               end if
            end if

c            call rhaar(q,workmr1,propvr)

            do i=1,q
               do j=1,q   
                   workmr2(i,j)=rnorm(ortho(i,j),tune4*0.05d0)
               end do
            end do
            call rhaar2(workmr,workmr2,q,propvr)

            call loglikpt_covarcan2(m,q,nsubject,iflagr,parti,
     &                              whicho,whichn,b,bzc,cpar,detloglc,
     &                              linf,lsup,mu,sigma,sigmainvc,
     &                              propvr,workvr,workmhr,workmr,
     &                              loglikn,fixed)

c++++++++++ acceptance step
         
            ratio=loglikn-logliko

            if(log(dble(runif())).lt.ratio)then
               acrate(7)=acrate(7)+1.d0
               do i=1,q
                  do j=1,q
                     ortho(i,j)=workmr2(i,j)
                     sigmainv(i,j)=sigmainvc(i,j)
                  end do
               end do

               do i=1,nsubject
                  do j=1,q
                     bz(i,j)=bzc(i,j)
                  end do   
               end do 
               detlogl=detloglc
               logliko=loglikn

               if(adaptivep.eq.1)aratep=aratep+1.d0

            end if
         end if 

c++++++++++++++++++++++++++++++++++         
c+++++++ save samples
c++++++++++++++++++++++++++++++++++         

         curr(1)=cpar
         curr(2)=disp

         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1

c+++++++++++++ random effects

               k=0
               do i=1,nsubject
                  do j=1,q
                     bsave(i,j)=bsave(i,j)+b(i,j)                  
                     k=k+1
                     randsave(isave,k)=b(i,j)
                  end do   
               end do


c+++++++++++++ predictive information

               call sampredptun(narea,q,nsubject,parti,m,ybar,
     &               massi,pattern,iflagr,whichn,whicho,bz,
     &               cpar,limw,linf,lsup,thetac,fixed)

               do i=1,q
                  do j=1,q
                     propvr(i,j)=0.d0
                     workmr(i,j)=0.d0
                  end do
               end do
               call cholesky(q,sigma,workmhr)
               do i=1,q
                  do j=1,i
                     propvr(i,j)=workmhr(ihmssf(i,j,q))
                  end do
               end do

               call rhaar2(workmr,ortho,q,workmr1)

               do i=1,q
                  do j=1,q
                     tmp1=0.d0
                     do k=1,q
                        tmp1=tmp1+propvr(i,k)*workmr1(k,j) 
                     end do 
                     workmr(i,j)=tmp1
                  end do
               end do

               k=nsubject*q
               do i=1,q
                  tmp1=0.d0
                  do j=1,q
                     tmp1=tmp1+workmr(i,j)*thetac(j)   
                  end do
                  theta(i)=tmp1+mu(i)
                  k=k+1
                  randsave(isave,k)=theta(i)
               end do

c+++++++++++++ functional parameter

               if(samplef.eq.1)then

               call samplefuncpt(fixed,m,q,nsubject,cpar,bz,theta,
     &                           iflagr,pattern,thetac,workmr1)  

               do i=1,q
                  tmp1=0.d0
                  do j=1,q
                     tmp1=tmp1+workmr(i,j)*thetac(j)   
                  end do
                  betar(i)=tmp1+mu(i)
               end do

               do i=1,q
                  do j=1,q
                     tmp1=0.d0  
                     do k=1,q
                        tmp1=tmp1+workmr(i,k)*workmr1(k,j)
                     end do
                     propvr(i,j)=tmp1
                  end do
               end do

               do i=1,q
                  do j=1,q
                     tmp1=0.d0  
                     do k=1,q
                        tmp1=tmp1+propvr(i,k)*workmr(j,k)
                     end do
                     workmr1(i,j)=tmp1
                  end do
               end do
               

c+++++++++++++ regression coefficients

               do i=1,q
                  thetasave(isave,i)=betar(i)
               end do

               end if  

               if(nfixed.gt.0)then
                  do i=1,p
                     thetasave(isave,q+i)=beta(i)
                     betasave(i)=betasave(i)+beta(i)
                  end do
               end if   

c+++++++++++++ dispersion parameter

               thetasave(isave,q+nfixed+1)=1.d0/disp
               betasave(p+1)=betasave(p+1)+disp

c+++++++++++++ baseline mean

               do i=1,q
                  thetasave(isave,q+nfixed+1+i)=mu(i)
               end do

c+++++++++++++ baseline covariance

               k=0
               do i=1,q
                  do j=i,q
                     k=k+1
                     thetasave(isave,q+nfixed+1+q+k)=sigma(i,j)
                  end do
               end do

c+++++++++++++ precision parameter
               k=(q*(q+1)/2)  
               thetasave(isave,q+nfixed+1+q+k+1)=cpar

c+++++++++++++ type of partition
               call rhaar2(workmr,ortho,q,workmr2)

               k=(q*(q+1)/2)+1
               do i=1,q
                  do j=1,q
                     k=k+1
                     thetasave(isave,q+nfixed+1+q+k)=workmr2(i,j) 
                  end do
               end do   

c+++++++++++++ random effects variance
               l=0
               do i=1,q
                  do j=i,q
                     l=l+1
                     thetasave(isave,q+nfixed+1+q+k+l)=workmr1(i,j)
                  end do
               end do   

c+++++++++++++ cpo
               dbarc=0.d0
               do i=1,nrec
                  yij=y(i)
                  eta=0.d0
                  do j=1,p
                     eta=eta+x(i,j)*beta(j)
                  end do
                  do j=1,q
                     eta=eta+z(i,j)*b(subject(i),j)
                  end do
                  eta=eta+x(i,p+1)                  
                  
                  mean=exp(eta)
                  
                  call dgamma(yij,mean,disp,tmp1) 
                  cpo(i,1)=cpo(i,1)+1.0d0/exp(tmp1)
                  cpo(i,2)=cpo(i,2)+exp(tmp1)

                  dbarc=dbarc+tmp1                  
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
      
      do i=1,3
         acrate(i)=acrate(i)/dble(nscan)
      end do
      do i=4,6
         acrate(i)=acrate(i)*dble(nbase)/dble(nscan)
      end do
      acrate(7)=acrate(7)/dble(nscan)
      
      do i=1,nrec
         cpo(i,1)=dble(nsave)/cpo(i,1)
         cpo(i,2)=cpo(i,2)/dble(nsave)                           
      end do

      do i=1,p+1
         betasave(i)=betasave(i)/dble(nsave)
      end do

      do i=1,nsubject
         do j=1,q
            bsave(i,j)=bsave(i,j)/dble(nsave)
         end do
      end do   

      dhat=0.d0
      lpml=0.d0
      do i=1,nrec
         yij=y(i)
         eta=0.d0
         do j=1,p
            eta=eta+x(i,j)*betasave(j)
         end do
         do j=1,q
            eta=eta+z(i,j)*bsave(subject(i),j)
         end do
         eta=eta+x(i,p+1)                  
         
         mean=exp(eta)
         call dgamma(yij,mean,betasave(p+1),tmp1) 

         dhat=dhat+tmp1   
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
            
      do i=1,5
         curr(i+2)=mc(i)
      end do
      
      return
      end


