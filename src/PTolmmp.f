c=======================================================================                      
      subroutine ptolmmp
     &                (datastr,maxni,nrec,nsubject,ncateg,p,q,nfixed,   #8 
     &                 subject,                                         #1  
     &                 x,xtx,yr,y,z,                                    #5
     &                 a0b0,mu0,prec1,prec2,sb,tinv,                    #6
     &                 mcmc,nsave,                                      #2
     &                 acrate,cpo,randsave,thetasave,                   #4
     &                 curr,b,beta,mu,sigma,                            #5 
     &                 iflagp,res,workmp1,workmhp1,workvp1,             #5
     &                 xty,                                             #1
     &                 iflagr,parti,whicho,whichn,bz,bzc,limw,linf,lsup,#9
     &                 propvr,sigmainv,theta,thetac,                    #4
     &                 workmhr,workmr,workmr1,workmr2,workvr,ybar,ztz,  #7
     &                 zty,                                             #1
     &                 sigmac,sigmainvc,workmhr2,                       #3 
     &                 massi,pattern,betasave,bsave)                    #4
c=======================================================================                      
c     # of arguments = 65.
c
c     Subroutine `ptolmmp' to run a Markov chain in a semiparametric 
c     probit mixed effect model for ordinal data, using a Mixture of 
c     Multivariate Polya trees prior for the distribution of the random 
c     effects.
c
c     Copyright: Alejandro Jara, 2007
c
c     Version 1.0:
c
c     Last modification: 30-04-2007.
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
c     Alejandro Jara
c     Biostatistical Centre
c     Katholieke Universiteit Leuven
c     U.Z. Sint-Rafa�l
c     Kapucijnenvoer 35
c     B-3000 Leuven
c     Voice: +32 (0)16 336892 
c     Fax  : +32 (0)16 337015 
c     URL  : http://student.kuleuven.be/~s0166452/
c     Email: Alejandro.JaraVallejos@med.kuleuven.be
c
c---- Data -------------------------------------------------------------
c 
c        datastr     :  integer matrix giving the number of measurements
c                       and the location in y of the observations for 
c                       each subject, datastr(nsubject,maxni+1)
c        maxni       :  integer giving the maximum number of 
c                       measurements for subject.
c        ncateg      :  integer giving the number of levels in the
c                       ordinal response. 
c        nrec        :  integer giving the number of observations.
c        nsubject    :  integer giving the number of subjects.
c        nfixed      :  integer giving the number of fixed effects,
c                       if nfixed is 0 then p=1.
c        p           :  integer giving the number of fixed coefficients.
c        q           :  integer giving the number of random effects.
c        subject     :  integer vector giving the subject for each.
c                       observation, subject(nsubject).
c        typep       :  type of partition to be considered.
c                       1=chol,2=square root 1, 3=square root 2.
c        x           :  real matrix giving the design matrix for the 
c                       fixed effects, x(nrec,p). 
c        xtx         :  real matrix givind the product X^tX, xtx(p,p).
c        y           :  real vector giving the latent data,
c                       y(nrec).
c        yr          :  integer vector giving the response variable,
c                       yr(nrec).
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
c                       precision and prior mean for the fixed effects,
c                       sb(p).
c        tinv        :  real matrix giving the scale matrix for the
c                       inverted-Wishart prior distribution for the
c                       covariance matrix of the random effects, 
c                       sigma ~ Inv-Wishart(nu0,tinv^{-1}), such that 
c                       E(sigma)=(1/(nu0-q-1)) * tinv 
c                       (This is for the base line distribution)
c
c-----------------------------------------------------------------------
c
c---- MCMC parameters --------------------------------------------------
c
c        nburn       :  integer giving the number of burn-in scans.
c        ndisplay    :  integer giving the number of saved scans to be
c                       displayed on screen.
c        nskip       :  integer giving the thinning interval.
c        nsave       :  integer giving the number of scans to be saved.\
c        nbase       :  integer giving the the number of scans where 
c                       the baseline distribution and the precision
c                       parameter are sampled.
c        
c-----------------------------------------------------------------------
c
c---- Output -----------------------------------------------------------
c
c        acrate      :  real vector giving the MH acceptance rate, 
c                       acrate(4). 
c        cpo         :  real giving the cpo. 
c        randsave    :  real matrix containing the mcmc samples for
c                       the random effects and prediction,
c                       randsave(nsave,q*(nsubject+1)).
c        thetasave   :  real matrix containing the mcmc samples for
c                       the averaged random effects, fixed effects, 
c                       error variance, and mean and covariance of
c                       the baseline distribution, 
c                       thetsave(nsave,q+nfixed+q+2*nuniq(Sigma)+1+
c                       ncateg-2).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        cpar        :  real giving the current value of the precision
c                       parameter of the PT.
c        b           :  real matrix giving the current value of the 
c                       random effects, b(nsubject,q).
c        beta        :  real vector giving the current value of the 
c                       fixed effects, beta(p).
c        betar       :  real vector giving the current value of the 
c                       averaged random effects, betar(q) (contained
c                       now in the vector curr).
c        mu          :  real vector giving the mean of the normal 
c                       base line distribution for the random effects,
c                       mu(q).
c        sigma       :  real matrix giving the current value of the
c                       covariance matrix for normal base line 
c                       distribution for the random effects,
c                       sigma(q,q).
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        acrate2     :  real working varible
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
c        dispcount   :  index. 
c        dlnrm       :  density of a log-normal distribution.
c        dnrm        :  density of a normal distribution.
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
c        res         :  real vector used to save the residual effects,
c                       res(nrec).
c        rtlnorm     :  real truncated log normal random number generator.
c        rtnorm      :  real truncated normal random number generator.
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
c        tmp3        :  real used to accumulate quantities.
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
c        xty         :  real vector used to save the product 
c                       Xt(Y-Zb), xty(p).
c        ybar        :  real vector used to save the mean of random  
c                       effects and the probabilities in each 
c                       partition area, ybar(2**q).
c        zty         :  real vector used to save the product 
c                       Zt(Y-Xbeta), zty(q).
c        ztz         :  real matrix used to save the product 
c                       ZtSigma^1Z, ztz(q,q).
c
c=======================================================================                  
      implicit none 

c+++++Data
      integer maxni,ncateg,nrec,nsubject,nfixed,p,q,subject(nrec)
      integer datastr(nsubject,maxni+1),yr(nrec)
      real*8 y(nrec),x(nrec,p),z(nrec,q),xtx(p,p)	
      
c+++++Prior 
      integer fixed,m,murand,sigmarand,typep
      real*8 aa0,ab0,a0b0(6),mu0(q),nu0,prec1(p,p),prec2(q,q)
      real*8 sb(p)
      real*8 tinv(q,q)      

c+++++MCMC parameters
      integer mcmc(11),nburn,nskip,nsave,ndisplay,nbase
      real*8 tune1,tune2,tune3

c+++++Output
      real*8 acrate(4),cpo(nrec,2)
      real*8 randsave(nsave,q*(nsubject+1))
      real*8 thetasave(nsave,q+nfixed+q+(q*(q+1))+1+ncateg-2)

c+++++Current values of the parameters
      real*8 curr(ncateg+q+5),cpar,beta(p),b(nsubject,q)
      real*8 mu(q),sigma(q,q)

c++++ model�s performance
      real*8 mc(5)
      real*8 betasave(p+ncateg-1),bsave(nsubject,q)

c+++++Working space - External
      integer iflagp(p) 
      real*8 res(nrec)
      real*8 workmp1(p,p)
      real*8 workmhp1(p*(p+1)/2)
      real*8 workvp1(p)
      real*8 xty(p)

      integer iflagr(q) 
      integer parti(q)
      integer whicho(nsubject),whichn(nsubject)      
      real*8 bz(nsubject,q),bzc(nsubject,q)
      real*8 limw(q),linf(q),lsup(q)
      real*8 propvr(q,q)
      real*8 sigmainv(q,q)
      real*8 theta(q),thetac(q)
      real*8 workmhr(q*(q+1)/2)
      real*8 workmr(q,q)
      real*8 workmr1(q,q),workmr2(q,q)
      real*8 workvr(q)
      real*8 ybar(2**q)
      real*8 ztz(q,q),zty(q)
      
      real*8 sigmac(q,q),sigmainvc(q,q)
      real*8 workmhr2(q*(q+1)/2)

      integer massi(2**q)
      integer pattern(q)

c+++++Working space - RNG
      integer seed1,seed2

c++++ model�s performance
      real*8 dbarc,dbar,dhat,pd,lpml
      
c+++++Working space - Internal
      integer baseskip
      integer dispcount
      integer i,j,k,l
      integer ihmssf
      integer iscan,isave
      integer narea,ni,nscan,nu
      integer skipcount
      integer sprint
      real*8 acrate2,cdfnorm,cparc
      real*8 detlogl,detloglc,dlnrm,dnrm
      real*8 logcgkn,logcgko
      real*8 loglikn,logliko
      real*8 logpriorn,logprioro
      real*8 ratio,rtnorm,rtlnorm,runif,tmp1,tmp2,tmp3
      logical ainf,asup

c+++++CPU time
      real*8 sec00,sec0,sec1,sec

c+++++Adaptive MH for mu
      integer countermh
      real*8 mumh(100),sigmamh(100,100)

c+++++Adaptive MH for sigma
      integer nadaptive
      parameter(nadaptive=2000)
      integer adaptives,sigmaskip
      real*8 aratesigma

c+++++Adaptive MH for c
      integer adaptivec,cskip
      real*8 aratec

c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      nbase=mcmc(4)      
      m=mcmc(5)
      seed1=mcmc(6)
      seed2=mcmc(7)
      typep=mcmc(8)
      murand=mcmc(9)
      sigmarand=mcmc(10)
      fixed=mcmc(11)

      aa0=a0b0(1)
      ab0=a0b0(2)
      nu0=a0b0(3)
      tune1=a0b0(4)
      tune2=a0b0(5)
      tune3=a0b0(6)
      
      narea=2**q

      cpar=curr(ncateg)
      
c++++ set random number generator
      call setall(seed1,seed2)

c++++ transforming random effects and calculate log-likelihood
c++++ for the baseline covariance matrix

      logliko=0.d0

      do i=1,q
         mumh(i)=0.d0      
         do j=1,q
            sigmamh(i,j)=0.d0         
            workmr(i,j)=sigma(i,j)
         end do
      end do

      call invdet(workmr,q,sigmainv,detlogl,iflagr,workvr)

      if(typep.eq.1)then
         do i=1,q
            do j=1,q
               workmr(i,j)=0.d0
               sigmainv(i,j)=0.d0
            end do
         end do
         call cholesky(q,sigma,workmhr)
         do i=1,q
            do j=1,i
               sigmainv(i,j)=workmhr(ihmssf(i,j,q))
            end do
         end do
         call inverse(sigmainv,q,iflagr)      
         
       else if(typep.eq.2)then
         do i=1,q
            do j=1,q
               workmr(i,j)=sigma(i,j)
               propvr(i,j)=0.d0
               sigmainv(i,j)=0.d0
            end do
         end do
         call eigenv(q,q,workmr,workvr,zty,ztz)

         do i=1,q
            propvr(i,i)=sqrt(workvr(i))
         end do
         
         do i=1,q
            do j=1,q
               tmp1=0.d0
               do k=1,q
                  tmp1=tmp1+ztz(i,k)*propvr(k,j)
               end do
               sigmainv(i,j)=tmp1
            end do
         end do
         call inverse(sigmainv,q,iflagr)      
      
       else
         do i=1,q
            do j=1,q
               workmr(i,j)=sigma(i,j)
               propvr(i,j)=0.d0
               sigmainv(i,j)=0.d0
            end do
         end do
         call eigenv(q,q,workmr,workvr,zty,ztz)

         do i=1,q
            propvr(i,i)=sqrt(workvr(i))
         end do

         do i=1,q
            do j=1,q
               tmp1=0.d0
               do k=1,q
                  tmp1=tmp1+ztz(i,k)*propvr(k,j)
               end do
               workmr(i,j)=tmp1
            end do
         end do

         do i=1,q
            do j=1,q
               tmp1=0.d0
               do k=1,q
                  tmp1=tmp1+workmr(i,k)*ztz(j,k)
               end do
               sigmainv(i,j)=tmp1
            end do
         end do
         call inverse(sigmainv,q,iflagr)      
      end if 

      call loglikpt_mucan(m,q,nsubject,parti,
     &                    whicho,whichn,b,bz,cpar,detlogl,
     &                    linf,lsup,mu,sigmainv,
     &                    theta,workmhr,workmr,workvr,
     &                    fixed,logliko)
     

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
      
      nscan=nburn+(nskip+1)*(nsave)

      call cpu_time(sec0)
      sec00=0.d0
      
      do iscan=1,nscan


c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++++
c+++++++ latent variable
c++++++++++++++++++++++++++++++++++

         do i=1,nrec
            tmp1=0.d0
            if(nfixed.gt.0)then
              do j=1,p
                 tmp1=tmp1+x(i,j)*beta(j)
              end do
            end if
            
            do j=1,q
               tmp1=tmp1+z(i,j)*b(subject(i),j) 
            end do
            
            if(yr(i).eq.1)then
              ainf=.true.
              asup=.false.
              y(i)=rtnorm(tmp1,1.d0,0.d0,curr(1),ainf,asup) 
            end if
            
            if(yr(i).eq.ncateg)then
              ainf=.false.
              asup=.true.
              y(i)=rtnorm(tmp1,1.d0,curr(ncateg-1),0.d0,ainf,asup) 
            end if

            do j=2,ncateg-1
               if(yr(i).eq.j)then
                  ainf=.false.
                  asup=.false.
                  y(i)=rtnorm(tmp1,1.d0,curr(j-1),curr(j),
     &                        ainf,asup) 
               end if
            end do
         end do


c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++
c+++++++ fixed effects        +++
c++++++++++++++++++++++++++++++++

         if(nfixed.eq.0)go to 1
            do i=1,p
               xty(i)=sb(i)
            end do

            do i=1,nrec
               tmp1=0.d0
               do j=1,q
                  tmp1=tmp1+z(i,j)*b(subject(i),j) 
               end do
               tmp1=y(i)-tmp1
             
               do j=1,p
                  xty(j)=xty(j)+x(i,j)*tmp1
               end do
            end do

            do i=1,p
               do j=1,p
                  workmp1(i,j)=xtx(i,j)+prec1(i,j)          
               end do
            end do

            call inverse(workmp1,p,iflagp)      

            do i=1,p
               tmp1=0.d0
               do j=1,p
                  tmp1=tmp1+workmp1(i,j)*xty(j) 
               end do
               workvp1(i)=tmp1
            end do

            call rmvnorm(p,workvp1,workmp1,workmhp1,xty,beta)

1        continue            


         if(nfixed.eq.0)then
             do i=1,nrec
                res(i)=y(i) 
             end do
           else
             do i=1,nrec
                tmp1=0.d0
                do j=1,p
                   tmp1=tmp1+x(i,j)*beta(j)    
                end do
                res(i)=y(i)-tmp1
             end do
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
                  propvr(j,k)=0.d0 
                  ztz(j,k)=workmr(j,k)
                  tmp1=tmp1+workmr(j,k)*mu(k)
               end do
               zty(j)=tmp1
               workvr(j)=0.d0
               iflagr(j)=0
            end do

            ni=datastr(i,1) 
            logliko=0.d0

            do j=1,q
               do k=1,q
                  tmp1=0.d0
                  do l=1,ni
                     tmp1=tmp1+z(datastr(i,l+1),j)*
     &                         z(datastr(i,l+1),k)
                  end do
                  ztz(j,k)=ztz(j,k)+tmp1
               end do
            end do

            do j=1,q
               tmp1=0.d0
               do k=1,ni
                  tmp1=tmp1+z(datastr(i,k+1),j)*
     &                      res(datastr(i,k+1))
               end do
               zty(j)=zty(j)+tmp1
            end do                  


            call inverse(ztz,q,iflagr)      

c            do j=1,q
c               limw(j)=0.d0
c               tmp1=0.d0
c               do k=1,q
c                  tmp1=tmp1+ztz(j,k)*zty(k) 
c               end do
c               limw(j)=tmp1
c            end do
c            call rmvnorm(q,limw,ztz,workmhr,workvr,thetac)

            call rmvnorm(q,theta,ztz,workmhr,workvr,thetac)

c++++++++++ evaluating the candidate generating kernel

c            call dmvn2(q,thetac,limw,ztz,logcgko,
c     &                 workvr,workmr1,workmr2,iflagr)
c
c            call dmvn2(q,theta,limw,ztz,logcgkn,
c     &                 workvr,workmr1,workmr2,iflagr)

c++++++++++ evaluating the likelihood
            loglikn=0.d0
            logliko=0.d0

            do j=1,ni

               tmp1=0.d0
               do k=1,p
                  tmp1=tmp1+x(datastr(i,j+1),k)*beta(k)   
               end do
               
               do k=1,q
                  tmp1=tmp1+z(datastr(i,j+1),k)*theta(k)   
               end do
               
               logliko=logliko+dnrm(y(datastr(i,j+1)),
     &                              tmp1,1.d0,1)

               tmp1=0.d0
               do k=1,p
                  tmp1=tmp1+x(datastr(i,j+1),k)*beta(k)   
               end do

               do k=1,q
                  tmp1=tmp1+z(datastr(i,j+1),k)*thetac(k)   
               end do
               
               loglikn=loglikn+dnrm(y(datastr(i,j+1)),
     &                              tmp1,1.d0,1)

            end do

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
c     &            logcgkn-logcgko+
     &            logpriorn-logprioro

            if(log(dble(runif())).lt.ratio)then
               acrate2=acrate2+1.d0
               do j=1,q
                  b(i,j)=thetac(j)
                  bz(i,j)=limw(j)
               end do
            end if

         end do

         acrate(1)=acrate(1)+acrate2/dble(nsubject)


         baseskip = baseskip + 1
         if(baseskip.ge.nbase)then
         countermh=countermh+1

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Updating mu using a MH step                        +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do i=1,q
            curr(ncateg+i)=0.d0
         end do

c         do i=1,nsubject
c            do j=1,q
c               curr(ncateg+j)=curr(ncateg+j)+b(i,j)
c            end do
c         end do
c
c         do i=1,q
c            curr(ncateg+i)=curr(ncateg+i)/dble(nsubject)
c         end do


c+++++++ update the log-likelihood for random effects
            
         call loglikpt_updatet(m,q,nsubject,parti,
     &                         whicho,whichn,bz,cpar,detlogl,
     &                         linf,lsup,
     &                         fixed,logliko)

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
     &                              prec2(i,j)    *
     &                             (theta(j)-mu0(j))

                  logprioro=logprioro+(mu(i)-mu0(i))* 
     &                              prec2(i,j)    *
     &                             (mu(j)-mu0(j))

               end do
            end do
         
            logpriorn=-0.5d0*logpriorn
            logprioro=-0.5d0*logprioro

c++++++++++ evaluating likelihood for muc

            call loglikpt_mucan(m,q,nsubject,parti,
     &                          whicho,whichn,b,bzc,cpar,detlogl,
     &                          linf,lsup,theta,sigmainv,
     &                          thetac,workmhr,workmr,workvr,
     &                          fixed,loglikn)


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
               acrate(2)=acrate(2)+1.d0
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
     &                      min(0.01,1.d0/sqrt(dble(iscan-nadaptive))))
                          else
                           tune2=exp(log(tune2)-
     &                      min(0.01,1.d0/sqrt(dble(iscan-nadaptive))))
                        end if  
                       else
                        if(aratesigma.lt.0.234)then
                           tune2=exp(log(tune2)+
     &                      min(0.01,1.d0/sqrt(dble(iscan-nadaptive))))
                          else
                           tune2=exp(log(tune2)-
     &                      min(0.01,1.d0/sqrt(dble(iscan-nadaptive))))
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

            call riwishart(q,nu,sigmac,sigmainvc,workmr,workvr,
     &                     workmhr,workmhr2,iflagr)

c++++++++++ evaluating the candidate generating kernel

            do i=1,q
               do j=1,q
                  ztz(i,j)=dble(nu-q-1)*sigma(i,j)
               end do
            end do

            call diwishart(q,nu,sigmac,ztz,workmr1,workmr2,workvr,
     &                     iflagr,logcgko)        

            do i=1,q
               do j=1,q
                  ztz(i,j)=dble(nu-q-1)*sigmac(i,j)
               end do
            end do

            call diwishart(q,nu,sigma,ztz,workmr1,workmr2,workvr,
     &                     iflagr,logcgkn)        

c++++++++++ ENDS: Simple MH

c++++++++++ evaluating the prior

            call diwishart(q,int(nu0),sigmac,tinv,workmr1,workmr2,
     &                     workvr,iflagr,logpriorn)        

            call diwishart(q,int(nu0),sigma,tinv,workmr1,workmr2,
     &                     workvr,iflagr,logprioro)        

c++++++++++ evaluating likelihood for sigmac

            call loglikpt_covarcan(m,q,nsubject,iflagr,parti,
     &                             whicho,whichn,b,bzc,cpar,detloglc,
     &                             linf,lsup,mu,sigmac,sigmainvc,
     &                             theta,workmhr,workmr,workvr,
     &                             loglikn,typep,workmr1,workmr2,
     &                             fixed)
     
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
               acrate(3)=acrate(3)+1.d0
               
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
     &                        min(0.01,1.d0/sqrt(dble(iscan))))                     
                       else
                        tune3=exp(log(tune3)-
     &                        min(0.01,1.d0/sqrt(dble(iscan))))                     
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
     &                            theta,workmhr,workmr,workvr,
     &                            fixed,loglikn)

c++++++++++ acceptance step
            ratio=loglikn+logpriorn-logliko-logprioro+
     &            logcgkn-logcgko

            if(log(dble(runif())).lt.ratio)then
               cpar=cparc
               acrate(4)=acrate(4)+1.d0
               logliko=loglikn

               if(adaptivec.eq.1)aratec=aratec+1.d0
               
            end if            
         end if

         baseskip=0
         end if

c++++++++++++++++++++++++++++++++++         
c+++++++ cutoff points
c++++++++++++++++++++++++++++++++++
         
         do i=2,ncateg-1
            tmp1=-100000.d0
            tmp2=+100000.d0
            
            do j=1,nrec
               if(yr(j).eq.(i+1))then
                  if(y(j).lt.tmp2)tmp2=y(j)
               end if
               if(yr(j).eq.i)then
                  if(y(j).gt.tmp1)tmp1=y(j)
               end if
            end do
            
            if(tmp2.lt.tmp1)then
              call rexit("Error in the limtis")
            end if
            
            curr(i)=tmp1+dble(runif())*(tmp2-tmp1)
         end do

c++++++++++++++++++++++++++++++++++         
c+++++++ save samples
c++++++++++++++++++++++++++++++++++         

         curr(ncateg)=cpar
         
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
     &               cpar,limw,linf,lsup,workmr,workmhr,
     &               workvr,thetac,fixed)

               if(typep.eq.1)then
                  do i=1,q
                     do j=1,q
                        workmr(i,j)=0.d0
                     end do
                  end do
                  call cholesky(q,sigma,workmhr)
                  do i=1,q
                     do j=1,i
                        workmr(i,j)=workmhr(ihmssf(i,j,q))
                     end do
                  end do

                else if(typep.eq.2)then
                  do i=1,q
                     do j=1,q
                        workmr(i,j)=0.d0
                        propvr(i,j)=0.d0
                     end do
                  end do
                  call eigenv(q,q,sigma,workvr,zty,ztz)

                  do i=1,q
                     propvr(i,i)=sqrt(workvr(i))
                  end do
         
                  do i=1,q
                     do j=1,q
                        tmp1=0.d0
                        do k=1,q
                           tmp1=tmp1+ztz(i,k)*propvr(k,j)
                        end do
                        workmr(i,j)=tmp1
                     end do
                  end do
                else
                  do i=1,q
                     do j=1,q
                        workmr(i,j)=0.d0
                        propvr(i,j)=0.d0
                        workmr1(i,j)=0.d0
                     end do
                  end do
                  call eigenv(q,q,sigma,workvr,zty,ztz)
                  do i=1,q
                     propvr(i,i)=sqrt(workvr(i))
                  end do

                  do i=1,q
                     do j=1,q
                        tmp1=0.d0
                        do k=1,q
                           tmp1=tmp1+ztz(i,k)*propvr(k,j)
                        end do
                        workmr1(i,j)=tmp1
                     end do
                  end do

                  do i=1,q
                     do j=1,q
                        tmp1=0.d0
                        do k=1,q
                           tmp1=tmp1+workmr1(i,k)*ztz(j,k)
                        end do
                        workmr(i,j)=tmp1
                     end do
                  end do
               end if 

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

               call samplefuncpt(fixed,m,q,nsubject,cpar,bz,theta,
     &                           iflagr,pattern,thetac,ztz)  

               do i=1,q
                  tmp1=0.d0
                  do j=1,q
                     tmp1=tmp1+workmr(i,j)*thetac(j)   
                  end do
                  curr(ncateg+i)=tmp1+mu(i)
               end do

               do i=1,q
                  do j=1,q
                     tmp1=0.d0  
                     do k=1,q
                        tmp1=tmp1+workmr(i,k)*ztz(k,j)                  
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
                  thetasave(isave,i)=curr(ncateg+i)
               end do

               if(nfixed.gt.0)then
                  do i=1,p
                     thetasave(isave,q+i)=beta(i)
                     betasave(i)=betasave(i)+beta(i)                     
                  end do
               end if   

c+++++++++++++ baseline mean

               do i=1,q
                  thetasave(isave,q+nfixed+i)=mu(i)
               end do

c+++++++++++++ baseline covariance

               k=0
               do i=1,q
                  do j=i,q
                     k=k+1
                     thetasave(isave,q+nfixed+q+k)=sigma(i,j)
                  end do
               end do

c+++++++++++++ cutoff points
               k=(q*(q+1)/2)    
               do i=2,ncateg-1
                  thetasave(isave,q+nfixed+q+k+i-1)=curr(i)
               end do
               
               do i=1,ncateg-1
                  betasave(p+i)=betasave(p+i)+curr(i)                                 
               end do

c+++++++++++++ precision parameter
               k=(q*(q+1)/2)+ncateg-2    
               thetasave(isave,q+nfixed+q+k+1)=cpar

c+++++++++++++ random effects variance
               l=0
               do i=1,q
                  do j=i,q
                     l=l+1
                     thetasave(isave,q+nfixed+q+k+1+l)=workmr1(i,j)
                  end do
               end do   

c+++++++++++++ cpo
               dbarc=0.d0
               do i=1,nrec
                  tmp1=0.d0
                  tmp2=0.d0
                  tmp3=0.d0
                  
                  if(nfixed.gt.0)then
                     do j=1,p
                        tmp1=tmp1+x(i,j)*beta(j)
                     end do
                  end if   
		  do j=1,q
		     tmp1=tmp1+z(i,j)*b(subject(i),j)
                  end do
                  
                  if(yr(i).eq.1)then
                     tmp1=curr(1)-tmp1
                     tmp3=cdfnorm(tmp1,0.d0,1.d0,1,0)
                  end if
            
                  if(yr(i).eq.ncateg)then
                     tmp1=curr(ncateg-1)-tmp1
                     tmp3=cdfnorm(tmp1,0.d0,1.d0,0,0)                  
                  end if

                  do j=2,ncateg-1
                     if(yr(i).eq.j)then
                        tmp2=curr(j)-tmp1
                        tmp3=curr(j-1)-tmp1
                        tmp3=cdfnorm(tmp2,0.d0,1.d0,1,0)-
     &                       cdfnorm(tmp3,0.d0,1.d0,1,0)                  
                     end if
                  end do                  
                  
                  cpo(i,1)=cpo(i,1)+1.0d0/tmp3
                  cpo(i,2)=cpo(i,2)+tmp3                  
                  
                  dbarc=dbarc+log(tmp3)
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
      
      do i=1,1
         acrate(i)=acrate(i)/dble(nscan)
      end do
      do i=2,4
         acrate(i)=acrate(i)*dble(nbase)/dble(nscan)
      end do
      
      do i=1,nrec
         cpo(i,1)=dble(nsave)/cpo(i,1)
         cpo(i,2)=cpo(i,2)/dble(nsave)                                             
      end do

      do i=1,p+ncateg-1
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
         tmp1=0.d0
         tmp2=0.d0
         tmp3=0.d0
         
         if(nfixed.gt.0)then
            do j=1,p
               tmp1=tmp1+x(i,j)*betasave(j)
            end do
         end if   
	  do j=1,q
	     tmp1=tmp1+z(i,j)*bsave(subject(i),j)
         end do
         
         if(yr(i).eq.1)then
            tmp1=betasave(p+1)-tmp1
            tmp3=cdfnorm(tmp1,0.d0,1.d0,1,0)
         end if
      
         if(yr(i).eq.ncateg)then
            tmp1=betasave(p+ncateg-1)-tmp1
            tmp3=cdfnorm(tmp1,0.d0,1.d0,0,0)                  
         end if

         do j=2,ncateg-1
            if(yr(i).eq.j)then
               tmp2=betasave(p+j)-tmp1
               tmp3=betasave(p+j-1)-tmp1
               tmp3=cdfnorm(tmp2,0.d0,1.d0,1,0)-
     &              cdfnorm(tmp3,0.d0,1.d0,1,0)                  
            end if
         end do                  
         dhat=dhat+log(tmp3)
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
         curr(ncateg+q+i)=mc(i)      
      end do
      
      return
      end
         