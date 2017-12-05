c=======================================================================                      
      subroutine dpmglmmpois(
     &                  maxni,nrec,nsubject,nfixed,p,q,subject,datastr, 
     &                  y,x,z,roffset,                                  
     &                  a0b0,prec,sb,nu,tinv1,smu,psiinv,tinv2,         
     &                  mcmc,nsave,                                       
     &                  ncluster,ss,alpha,beta,b,betar,mu,              
     &                  sigma,sigmainv,mub,sigmab,sigmabinv,mc,         
     &                  acrate,cpo,randsave,thetasave,musave,clustsave, 
     &                  iflagp,betac,workmhp,workvp,xtx,xty,            
     &                  iflagr,theta,thetac,workmhr,workmhr2,workmr,    
     &                  workvr,ztz,zty,cstrt,ccluster,prob,quadf,       
     &                  seed,betasave,bsave)
c=======================================================================                      
c     # of arguments = 63.
c
c     Subroutine `dpmglmmpois' to run a Markov chain in a semiparametric 
c     poisson mixed effect model, using a Dirichlet process mixture of 
c     normals prior for the distribution of the random effects.
c
c     Copyright: Alejandro Jara, 2007-2010.
c
c     Version 1.0:
c
c     Last modification: 20-04-2007.
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
c        roffset     :  real vector giving the real offset for each
c                       observation, roffset(nrec).
c        subject     :  integer vector giving the subject for each.
c                       observation, subject(nsubject).
c        x           :  real matrix giving the design matrix for the 
c                       fixed effects, x(nrec,p). 
c        y           :  integer vector giving the response variable,
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
c        nu01        :  integer giving the degrees of freedom for the
c                       inverted-Wishart prior distribution for the
c                       covariance matrix of normal kernel.
c        nu02        :  integer giving the degrees of freedom for the
c                       inverted-Wishart prior distribution for the
c                       centering distribution.
c        prec        :  real matrix giving the prior precision matrix
c                       for the fixed effects, prec(p,p).
c        psiinv      :  real matrix giving the prior precision matrix
c                       for the baseline mean, psiinv(q,q).
c        sb          :  real vector giving the product of the prior 
c                       precision and prior mean for the fixed effects,
c                       sb(p).
c        smu         :  real vector giving the product of the prior 
c                       precision and prior mean for the baseline mean,
c                       smu(q).
c        tinv1       :  real matrix giving the scale matrix for the
c                       inverted-Wishart prior distribution for the
c                       covariance matrix of the normal kernel.
c        tinv2       :  real matrix giving the scale matrix for the
c                       inverted-Wishart prior distribution for the
c                       covariance matrix of centering distribution.
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
c        
c-----------------------------------------------------------------------
c
c---- Output -----------------------------------------------------------
c
c        cpo         :  real giving the cpo. 
c        clustsave   :  integer matrix containing the cardinality of 
c                       each cluster in each of the mcmc samples,
c                       clustsave(nsave,nsubject). 
c        randsave    :  real matrix containing the mcmc samples for
c                       the random effects and prediction,
c                       randsave(nsave,q*(nsubject+1)).
c        musave      :  real matrix containing the mcmc samples for the
c                       cluster locations, musave(nsave,q*nsubject).
c        thetasave   :  real matrix containing the mcmc samples for
c                       the averaged random effects, fixed effects, 
c                       error variance, the normal kernel variance, 
c                       and mean and covariance ofthe baseline 
c                       distribution, the number of clusters, and the 
c                       precision parameter, 
c                       thetsave(nsave,q+nfixed+q+(q*(q+1))+2).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        alpha       :  real giving the current value of the precision
c                       parameter of the Dirichlet process.
c        b           :  real matrix giving the current value of the 
c                       random effects, b(nsubject,q).
c        beta        :  real vector giving the current value of the 
c                       fixed effects, beta(p).
c        betar       :  real vector giving the current value of the 
c                       averaged random effects, betar(q).
c        mu         :   real matrix giving the cluster locations
c                       mu(nsubject,q).
c        mub         :  real vector giving the mean of the normal 
c                       base line distribution for the random effects,
c                       mub(q).
c        ncluster    :  integer giving the number of clusters in the
c                       random effects.
c        sigma       :  real matrix giving the current value of the
c                       covariance matrix for normal kernel,
c                       sigma(q,q).
c        sigmab      :  real matrix giving the current value of the
c                       covariance matrix for normal base line 
c                       distribution for the random effects,
c                       sigma(q,q).
c        sigmainv    :  real matrix used to save the inverse of the
c                       kernel covariance matrix,
c                       sigmainv(q,q).
c        sigmabinv   :  real matrix used to save the inverse of the
c                       centering covariance matrix,
c                       sigmainv(q,q).
c        ss          :  integer vector giving the cluster label for 
c                       each subject, ss(nsubject).
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        ccluster    :  integer vector indicating the number of
c                       subjects in each cluster, ccluster(nsubject).
c        cstrt       :  integer matrix used to save the cluster
c                       structure, cstrt(nsubject,nsubject).
c        dispcount   :  index. 
c        dnrm        :  density of a normal distribution.
c        evali       :  integer indicator used in updating the state.
c        i           :  index. 
c        ii          :  index. 
c        iflagp      :  integer vector used to invert the of the lhs
c                       least square solution for the fixed effects,
c                       iflagp(p).
c        iflagr      :  integer vector used to invert the of the lhs
c                       least square solution for the random effects,
c                       iflagr(q).
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        k           :  index. 
c        l           :  index.
c        m           :  index.
c        prob        :  real vector used to update the cluster 
c                       structure, prob(nsubject+1).
c        quadf       :  real matrix used to save the bilinear product
c                       of random effects, quadf(q,q).
c        ni          :  integer indicator used in updating the state. 
c        ns          :  integer indicator used in updating the state. 
c        nscan       :  integer indicating the total number of MCMC
c                       scans.
c        sec         :  cpu time working variable.
c        sec0        :  cpu time working variable.
c        sec00       :  cpu time working variable.
c        sec1        :  cpu time working variable.
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        since       :  index.
c        skipcount   :  index. 
c        sse         :  real used to save the SS of the errors.
c        theta       :  real vector used to save randomnly generated
c                       random effects, theta(q).
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
c        tmp3        :  real used to accumulate quantities.
c        workmr      :  real matrix used to update the random effects,
c                       workmr(q,q).
c        workmhp     :  real vector used to update the fixed effects,
c                       workmhp(p*(p+1)/2).
c        workmhr     :  real vector used to update the random effects,
c                       workmhr(q*(q+1)/2).
c        workmhr2    :  real vector used to update the random effects,
c                       workmhr2(q*(q+1)/2).
c        workvp      :  real vector used to update the fixed effects,
c                       workvp(p).
c        workvr      :  real vector used to update the random effects,
c                       workvr(p).
c        xtx         :  real matrix givind the product X^tX, xtx(p,p).
c        xty         :  real vector used to save the product 
c                       Xt(Y-Zb), xty(p).
c        y           :  real vector giving the latent data,
c                       y(nrec).
c        zty         :  real vector used to save the product 
c                       Zt(Y-Xbeta), zty(q).
c        ztz         :  real matrix used to save the product 
c                       ZtSigma^1Z, ztz(q,q).
c
c=======================================================================                  
      implicit none 

c+++++Data
      integer maxni,nrec,nsubject,nfixed,p,q,subject(nrec)
      integer datastr(nsubject,maxni+1),y(nrec)
      double precision roffset(nrec),x(nrec,p),z(nrec,q)
      
c+++++Prior 
      integer nu01,nu02,nu(2) 
      double precision aa0,ab0,a0b0(2) 
      double precision prec(p,p),sb(p,2)
      double precision tinv1(q,q)
      double precision smu(q),psiinv(q,q)
      double precision tinv2(q,q)
      
c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Output
      integer clustsave(nsave,nsubject) 
      double precision acrate(2),cpo(nrec,2)
      double precision randsave(nsave,q*(nsubject+1))
      double precision musave(nsave,q*nsubject)
      double precision thetasave(nsave,q+nfixed+q+(q*(q+1))+2)

c+++++Current values of the parameters
      integer ncluster,ss(nsubject)
      double precision alpha
      double precision beta(p)
      double precision b(nsubject,q)
      double precision betar(q)
      double precision mu(nsubject,q)
      double precision sigma(q,q),sigmainv(q,q)
      double precision mub(q)
      double precision sigmab(q,q),sigmabinv(q,q)

c+++++Seeds
      integer seed(2),seed1,seed2

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++External working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++fixed effects
      integer iflagp(p) 
      double precision betac(p)
      double precision workmhp(p*(p+1)/2)
      double precision workvp(p)
      double precision xtx(p,p)
      double precision xty(p)

c+++++random effects
      integer iflagr(q) 
      double precision ztz(q,q),zty(q)
      double precision theta(q),thetac(q)
      double precision workmr(q,q)      
      double precision workmhr(q*(q+1)/2),workmhr2(q*(q+1)/2)
      double precision workvr(q)      

c+++++DPM
      integer cstrt(nsubject,nsubject)
      integer ccluster(nsubject)
      double precision prob(nsubject+1)

c+++++Kernel
      double precision quadf(q,q)

c++++ model's performance
      double precision mc(5)
      double precision betasave(p),bsave(nsubject,q)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Internal working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++General
      integer evali,ii,i,j,k,l,m,ni,ns 
      integer ok
      integer since,sprint 
      integer yij      
      double precision acrate2
      double precision eta,gprime
      double precision mean
      double precision logcgkn,logcgko
      double precision loglikn,logliko
      double precision logpriorn,logprioro
      double precision offset
      double precision ratio
      double precision tmp1,tmp2,tmp3
      double precision ytilde

c+++++MCMC
      integer dispcount,isave,iscan,nscan,skipcount 

c+++++RNG and distributions
      double precision dpoiss
      real runif

c+++++DP (functional parameter)
      double precision eps,rbeta,weight
      parameter(eps=0.01)

c++++ model's performance
      double precision dbarc,dbar,dhat,pd,lpml

c+++++CPU time
      double precision sec00,sec0,sec1,sec

c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)

      aa0=a0b0(1)
      ab0=a0b0(2)

      nu01=nu(1)
      nu02=nu(2)
      
c++++ set random number generator
      seed1=seed(1)
      seed2=seed(2)
      call setall(seed1,seed2)


c++++ set configurations
      do i=1,nsubject
         ccluster(ss(i))=ccluster(ss(i))+1
         cstrt(ss(i),ccluster(ss(i)))=i
      end do
 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      dbar=0.d0
      isave=0
      skipcount=0
      dispcount=0
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
                  xtx(i,j)=prec(i,j)
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

               eta=eta+roffset(i)
               offset=offset+roffset(i)

               mean=exp(eta)
               gprime=exp(-eta)

               ytilde=eta+(dble(yij)-mean)*gprime-offset
               
               do j=1,p
                  do k=1,p
                     xtx(j,k)=xtx(j,k)+x(i,j)*x(i,k)*mean
                  end do
                  xty(j)=xty(j)+x(i,j)*ytilde*mean
               end do
               
               logliko=logliko+dpoiss(dble(yij),mean,1)
            end do

            call inverse(xtx,p,iflagp)      

            do i=1,p
               tmp1=0.d0
               do j=1,p
                  tmp1=tmp1+xtx(i,j)*xty(j) 
               end do
               workvp(i)=tmp1
            end do

            call rmvnorm(p,workvp,xtx,workmhp,xty,betac)

c++++++++++ evaluating the candidate generating kernel

            call dmvnd(p,betac,workvp,xtx,logcgko,iflagp)

c++++++++++ evaluating the likelihood

            do i=1,p
               do j=1,p
                  xtx(i,j)=prec(i,j)
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

               eta=eta+roffset(i)
               offset=offset+roffset(i)
               
               mean=exp(eta)
               gprime=exp(-eta)
               
               ytilde=eta+(dble(yij)-mean)*gprime-offset
               
               do j=1,p
                  do k=1,p
                     xtx(j,k)=xtx(j,k)+x(i,j)*x(i,k)*mean
                  end do
                  xty(j)=xty(j)+x(i,j)*ytilde*mean
               end do

               loglikn=loglikn+dpoiss(dble(yij),mean,1)
            end do

            call inverse(xtx,p,iflagp)      

            do i=1,p
               tmp1=0.d0
               do j=1,p
                  tmp1=tmp1+xtx(i,j)*xty(j) 
               end do
               workvp(i)=tmp1
            end do

c++++++++++ evaluating the candidate generating kernel

            call dmvnd(p,beta,workvp,xtx,logcgkn,iflagp)

c++++++++++ prior ratio
            logprioro=0.d0
            logpriorn=0.d0
         
            do i=1,p
               do j=1,p
                  logpriorn=logpriorn+(betac(i)-sb(i,2))* 
     &                      prec(i,j)       *
     &                      (betac(j)-sb(j,2))

                  logprioro=logprioro+(beta(i)-sb(i,2))* 
     &                      prec(i,j)      *
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
 
         do ii=1,nsubject

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            do j=1,q
               theta(j)=b(ii,j)
            end do   

c++++++++++ generating a candidate

            do j=1,q
               tmp1=0.d0
               do k=1,q
                  ztz(j,k)=sigmainv(j,k)
                  tmp1=tmp1+mu(ss(ii),k)*sigmainv(j,k)
               end do
               zty(j)=tmp1
            end do

            ni=datastr(ii,1) 

            logliko=0.d0
            
            do j=1,ni
               eta=0.d0 
               mean=0.d0
               gprime=0.d0
               offset=0.d0

               yij=y(datastr(ii,j+1))
            
               do k=1,p
                  eta=eta+x(datastr(ii,j+1),k)*beta(k)
                  offset=offset+x(datastr(ii,j+1),k)*beta(k)
               end do
            
               do k=1,q
                  eta=eta+z(datastr(ii,j+1),k)*theta(k)
               end do

               eta=eta+roffset(datastr(ii,j+1))
               offset=offset+roffset(datastr(ii,j+1))
               
               mean=exp(eta)
               gprime=exp(-eta)
               
               ytilde=eta+(dble(yij)-mean)*gprime-offset
            
               do k=1,q
                  do l=1,q
                     ztz(k,l)=ztz(k,l)+
     &                        z(datastr(ii,j+1),k)*
     &                        z(datastr(ii,j+1),l)*mean
                  end do
                  zty(k)=zty(k)+z(datastr(ii,j+1),k)*ytilde*mean
               end do

               logliko=logliko+dpoiss(dble(yij),mean,1)

            end do

            call inverse(ztz,q,iflagr)      

            do j=1,q
               tmp1=0.d0
               do k=1,q
                  tmp1=tmp1+ztz(j,k)*zty(k) 
               end do
               workvr(j)=tmp1
            end do

            call rmvnorm(q,workvr,ztz,workmhr,zty,thetac)
            
c++++++++++ evaluating the candidate generating kernel

            call dmvnd(q,thetac,workvr,ztz,logcgko,iflagr)

c++++++++++ evaluating the likelihood

            do j=1,q
               tmp1=0.d0
               do k=1,q
                  ztz(j,k)=sigmainv(j,k)
                  tmp1=tmp1+mu(ss(ii),k)*sigmainv(j,k)
               end do
               zty(j)=tmp1
            end do

            ni=datastr(ii,1) 

            loglikn=0.d0
            
            do j=1,ni
               eta=0.d0 
               mean=0.d0
               gprime=0.d0
               offset=0.d0

               yij=y(datastr(ii,j+1))            
            
               do k=1,p
                  eta=eta+x(datastr(ii,j+1),k)*beta(k)
                  offset=offset+x(datastr(ii,j+1),k)*beta(k)
               end do
            
               do k=1,q
                  eta=eta+z(datastr(ii,j+1),k)*thetac(k)
               end do

               eta=eta+roffset(datastr(ii,j+1))
               offset=offset+roffset(datastr(ii,j+1))
               
               mean=exp(eta)
               gprime=exp(-eta)
               
               ytilde=eta+(dble(yij)-mean)*gprime-offset
            
               do k=1,q
                  do l=1,q
                     ztz(k,l)=ztz(k,l)+
     &                        z(datastr(ii,j+1),k)*
     &                        z(datastr(ii,j+1),l)*mean
                  end do
                  zty(k)=zty(k)+z(datastr(ii,j+1),k)*ytilde*mean
               end do

               loglikn=loglikn+dpoiss(dble(yij),mean,1)

            end do

            call inverse(ztz,q,iflagr)      

            do j=1,q
               tmp1=0.d0
               do k=1,q
                  tmp1=tmp1+ztz(j,k)*zty(k) 
               end do
               workvr(j)=tmp1
            end do
            
c++++++++++ evaluating the candidate generating kernel

            call dmvnd(q,theta,workvr,ztz,logcgkn,iflagr)

c++++++++++ evaluating the prior

            logprioro=0.d0
            logpriorn=0.d0

            do j=1,q
               do k=1,q
                  logpriorn=logpriorn+(thetac(j)-mu(ss(ii),j))* 
     &                      sigmainv(j,k)      *
     &                      (thetac(k)-mu(ss(ii),k))

                  logprioro=logprioro+(theta(j)-mu(ss(ii),j))* 
     &                      sigmainv(j,k)      *
     &                      (theta(k)-mu(ss(ii),k))

               end do
            end do
      
            logpriorn=-0.5d0*logpriorn
            logprioro=-0.5d0*logprioro

c++++++++++ mh step
  
            ratio=loglikn-logliko+
     &            logcgkn-logcgko+
     &            logpriorn-logprioro

            if(log(dble(runif())).lt.ratio)then
               acrate2=acrate2+1.d0
               do j=1,q
                  b(ii,j)=thetac(j)
               end do
            end if

         end do

         acrate(2)=acrate(2)+acrate2/dble(nsubject)


c+++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ a) Polya Urn based on a collapsed state
c+++++++++++++++++++++++++++++++++++++++++++++++++

         do i=1,nsubject
         
            ns=ccluster(ss(i))

c++++++++++ subject in cluster with more than 1 observations
             
            if(ns.gt.1)then
          
               j=1
               ok=0
               do while(ok.eq.0.and.j.le.ns)
                  if(cstrt(ss(i),j).eq.i)ok=j
                  j=j+1
               end do
   
               do j=ok,ns-1
                  cstrt(ss(i),j)=cstrt(ss(i),j+1)
               end do
          
               ccluster(ss(i))=ccluster(ss(i))-1 

               do j=1,ncluster
                  
                  do k=1,q
                     tmp1=0.d0
                     do l=1,q
                        ztz(k,l)=sigmabinv(k,l)+
     &                           dble(ccluster(j))*sigmainv(k,l)
                        tmp1=tmp1+sigmabinv(k,l)*mub(l)
                     end do
                     zty(k)=tmp1
                  end do
                  
                  call inverse(ztz,q,iflagr) 
                  
                  do k=1,ccluster(j)
                     do l=1,q
                        tmp1=0.d0 
                        do m=1,q
                           tmp1=tmp1+sigmainv(l,m)*b(cstrt(j,k),m)   
                        end do
                        zty(l)=zty(l)+tmp1
                     end do
                  end do 
                  
                  do k=1,q
                     tmp1=0.d0
                     do l=1,q
                        tmp1=tmp1+ztz(k,l)*zty(l)                     
                     end do
                     theta(k)=tmp1
                  end do
                  
                  do k=1,q
                     workvr(k)=b(i,k)
                     do l=1,q
                        ztz(k,l)=ztz(k,l)+sigma(k,l)
                     end do
                  end do
                  
                  call dmvnd(q,workvr,theta,ztz,tmp1,iflagr)        

                  prob(j)=dble(ccluster(j))*exp(tmp1)
               end do
               
               do k=1,q
                  workvr(k)=b(i,k)
                  theta(k)=mub(k)
                  do l=1,q
                     ztz(k,l)=sigma(k,l)+sigmab(k,l)
                  end do
               end do
               
               call dmvnd(q,workvr,theta,ztz,tmp1,iflagr)        
                    
               prob(ncluster+1)=alpha*exp(tmp1)

               call simdisc(prob,nsubject+1,ncluster+1,evali)

               ss(i)=evali
               
               ccluster(evali)=ccluster(evali)+1
               
               cstrt(evali,ccluster(evali))=i
               
               if(evali.gt.ncluster)then
                  ncluster=ncluster+1
               end if
            end if


c++++++++++ subject in cluster with only 1 observation
             
            if(ns.eq.1)then
                
               since=ss(i)
                
               if(since.lt.ncluster)then
                   call relabeldpm(i,since,nsubject,q,ncluster,
     &                             ccluster,ss,cstrt)                   
               end if

               ccluster(ncluster)=ccluster(ncluster)-1 
               ncluster=ncluster-1

               do j=1,ncluster
                  do k=1,q
                     tmp1=0.d0
                     do l=1,q
                        ztz(k,l)=sigmabinv(k,l)+
     &                           dble(ccluster(j))*sigmainv(k,l)
                        tmp1=tmp1+sigmabinv(k,l)*mub(l)
                     end do
                     zty(k)=tmp1
                  end do
                  
                  call inverse(ztz,q,iflagr) 
                  
                  do k=1,ccluster(j)
                     do l=1,q
                        tmp1=0.d0 
                        do m=1,q
                           tmp1=tmp1+sigmainv(l,m)*b(cstrt(j,k),m)   
                        end do
                        zty(l)=zty(l)+tmp1
                     end do
                  end do 
                  
                  do k=1,q
                     tmp1=0.d0
                     do l=1,q
                        tmp1=tmp1+ztz(k,l)*zty(l)                     
                     end do
                     theta(k)=tmp1
                  end do
                  
                  do k=1,q
                     workvr(k)=b(i,k)
                     do l=1,q
                        ztz(k,l)=ztz(k,l)+sigma(k,l)
                     end do
                  end do
                  
                  call dmvnd(q,workvr,theta,ztz,tmp1,iflagr)
                  
                  prob(j)=dble(ccluster(j))*exp(tmp1)
               end do

               do k=1,q
                  workvr(k)=b(i,k)
                  theta(k)=mub(k)
                  do l=1,q
                     ztz(k,l)=sigma(k,l)+sigmab(k,l)
                  end do
               end do
               
               call dmvnd(q,workvr,theta,ztz,tmp1,iflagr) 
               
               prob(ncluster+1)=alpha*exp(tmp1)
               
               call simdisc(prob,nsubject+1,ncluster+1,evali)
               
               ss(i)=evali
               
               ccluster(evali)=ccluster(evali)+1
               
               cstrt(evali,ccluster(evali))=i
               
               if(evali.gt.ncluster)then
                  ncluster=ncluster+1
               end if
            end if

         end do

c++++++++++++++++++++++++++++++
c+++++++ b) Resampling step
c++++++++++++++++++++++++++++++

         do ii=1,ncluster

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            ns=ccluster(ii)
            
            do k=1,q
               tmp1=0.d0
               do l=1,q
                  ztz(k,l)=sigmabinv(k,l)+
     &                     dble(ns)*sigmainv(k,l)
                  tmp1=tmp1+sigmabinv(k,l)*mub(l)
               end do
               zty(k)=tmp1
            end do
            
            call inverse(ztz,q,iflagr) 
            
            do k=1,ns
               do l=1,q
                  tmp1=0.d0 
                  do m=1,q
                     tmp1=tmp1+sigmainv(l,m)*b(cstrt(ii,k),m)   
                  end do
                  zty(l)=zty(l)+tmp1
               end do
            end do 

            do k=1,q
               tmp1=0.d0
               do l=1,q
                  tmp1=tmp1+ztz(k,l)*zty(l)                     
               end do
               workvr(k)=tmp1
            end do
            
            call rmvnorm(q,workvr,ztz,workmhr,zty,theta)

            do k=1,q
               mu(ii,k)=theta(k)
            end do
         end do


c++++++++++++++++++++++++++++++
c+++++++ Kernel variance
c++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do i=1,q
            do j=1,q
               quadf(i,j)=0.d0
            end do
         end do
         
         do i=1,nsubject
            do j=1,q
               do k=1,q
                  quadf(j,k)=quadf(j,k)+               
     &                (b(i,j)-mu(ss(i),j))*(b(i,k)-mu(ss(i),k))
               end do
            end do            
         end do

         do i=1,q
            do j=1,q
               quadf(i,j)=quadf(i,j)+tinv1(i,j)
            end do
         end do

         call riwishart(q,nu01+nsubject,quadf,ztz,workmr,workvr,
     &                  workmhr,workmhr2,iflagr)
         do i=1,q
            do j=1,q
               sigma(i,j)=quadf(i,j)
               sigmainv(i,j)=ztz(i,j)
            end do
         end do


c++++++++++++++++++++++++++++++++++         
c+++++++ Base line distribution
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do i=1,q
            zty(i)=smu(i)
            do j=1,q
               ztz(i,j)=(sigmabinv(i,j)*dble(ncluster))+psiinv(i,j)
            end do
         end do

         call inverse(ztz,q,iflagr) 

         do i=1,ncluster
            do j=1,q
               tmp1=0.d0
               do k=1,q
                  tmp1=tmp1+sigmabinv(j,k)*mu(i,k)
               end do
               zty(j)=zty(j)+tmp1
            end do
         end do
     
         do i=1,q
            tmp1=0.d0
            do j=1,q
               tmp1=tmp1+ztz(i,j)*zty(j)
            end do
            workvr(i)=tmp1
         end do

         call rmvnorm(q,workvr,ztz,workmhr,zty,theta)


c+++++++ check if the user has requested an interrupt
         call rchkusr()
     
         do i=1,q
            mub(i)=theta(i)
            do j=1,q
               quadf(i,j)=0.d0
            end do
         end do
         
         do i=1,ncluster
            do j=1,q
               do k=1,q
                  quadf(j,k)=quadf(j,k)+               
     &                 (mu(i,j)-mub(j))*(mu(i,k)-mub(k))
               end do
            end do
         end do

         do i=1,q
            do j=1,q
               quadf(i,j)=quadf(i,j)+tinv2(i,j)
            end do
         end do

         call riwishart(q,nu02+ncluster,quadf,ztz,workmr,workvr,
     &                  workmhr,workmhr2,iflagr)

         do i=1,q
            do j=1,q
               sigmab(i,j)=quadf(i,j)
               sigmabinv(i,j)=ztz(i,j)
            end do
         end do

c++++++++++++++++++++++++++++++++++         
c+++++++ Precision parameter
c++++++++++++++++++++++++++++++++++
         if(aa0.gt.0.d0)then
            call samalph(alpha,aa0,ab0,ncluster,nsubject)
         end if 


c++++++++++++++++++++++++++++++++++         
c+++++++ save samples
c++++++++++++++++++++++++++++++++++         

         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1

c+++++++++++++ random effects
               k=0
               do i=1,ncluster
                  do j=1,q
                     k=k+1
                     musave(isave,k)=mu(i,j)
                  end do
                  clustsave(isave,i)=ccluster(i)
               end do

               k=0
               do i=1,nsubject
                  do j=1,q
                     bsave(i,j)=bsave(i,j)+b(i,j)                  
                     k=k+1
                     randsave(isave,k)=b(i,j)
                  end do   
               end do

c+++++++++++++ predictive information

               do i=1,ncluster
                  prob(i)=dble(ccluster(i))/(alpha+dble(nsubject))
               end do
               prob(ncluster+1)=alpha/(alpha+dble(nsubject))

               call simdisc(prob,nsubject+1,ncluster+1,evali)
               
               if(evali.le.ncluster)then
                  do j=1,q
                     theta(j)=mu(evali,j)
                  end do
                else
                  call rmvnorm(q,mub,sigmab,workmhr,workvr,theta)
               end if
               call rmvnorm(q,theta,sigma,workmhr,workvr,zty)
               
               do i=1,q
                  k=k+1
                  randsave(isave,k)=zty(i) 
               end do

c+++++++++++++ functional parameters
               
               tmp1=rbeta(1.d0,alpha+dble(nsubject))
               do i=1,q
                  betar(i)=tmp1*theta(i)
               end do
               tmp2=tmp1
               weight=(1.d0-tmp1)
               
               do while((1.d0-tmp2).gt.eps)
                  tmp3=rbeta(1.d0,alpha+dble(nsubject))
                  tmp1=weight*tmp3
                  weight=weight*(1.d0-tmp3)

                  do i=1,ncluster
                     prob(i)=dble(ccluster(i))/(alpha+dble(nsubject))
                  end do
                  prob(ncluster+1)=alpha/(alpha+dble(nsubject))

                  call simdisc(prob,nsubject+1,ncluster+1,evali)
               
                  if(evali.le.ncluster)then
                     do j=1,q
                        theta(j)=mu(evali,j)
                     end do
                   else
                     call rmvnorm(q,mub,sigmab,workmhr,workvr,theta)
                  end if

                  do i=1,q
                     betar(i)=betar(i)+tmp1*theta(i)
                  end do
                  tmp2=tmp2+tmp1
               end do

               do i=1,ncluster
                  prob(i)=dble(ccluster(i))/(alpha+dble(nsubject))
               end do
               prob(ncluster+1)=alpha/(alpha+dble(nsubject))

               call simdisc(prob,nsubject+1,ncluster+1,evali)
               
               if(evali.le.ncluster)then
                  do j=1,q
                     theta(j)=mu(evali,j)
                  end do
                else
                  call rmvnorm(q,mub,sigmab,workmhr,workvr,theta)
               end if
               
               tmp1=weight

               do i=1,q
                  betar(i)=betar(i)+tmp1*theta(i)
               end do

c+++++++++++++ regression coefficients

               do i=1,q
                  thetasave(isave,i)=betar(i)
               end do

               if(nfixed.gt.0)then
                  do i=1,p
                     thetasave(isave,q+i)=beta(i)
                     betasave(i)=betasave(i)+beta(i)
                  end do
               end if   

c+++++++++++++ kernel variance
               k=0
               do i=1,q
                  do j=i,q
                     k=k+1
                     thetasave(isave,q+nfixed+k)=sigma(i,j)
                  end do
               end do

c+++++++++++++ baseline mean
               k=(q*(q+1)/2) 
               do i=1,q
                  thetasave(isave,q+nfixed+k+i)=mub(i)
               end do

c+++++++++++++ baseline covariance
               l=(q*(q+1)/2)
               k=0
               do i=1,q
                  do j=i,q
                     k=k+1
                     thetasave(isave,q+nfixed+l+q+k)=sigmab(i,j)
                  end do
               end do

c+++++++++++++ cluster information
               k=(q*(q+1)/2)  
               thetasave(isave,q+nfixed+k+q+k+1)=ncluster
               thetasave(isave,q+nfixed+k+q+k+2)=alpha

c+++++++++++++ cpo
               dbarc=0.d0
               do i=1,nrec
                  yij=y(i)
                  eta=0.d0
                  if(nfixed.gt.0)then
                     do j=1,p
                        eta=eta+x(i,j)*beta(j)
                     end do
                  end if   
                  do j=1,q
                     eta=eta+z(i,j)*b(subject(i),j)
                  end do
                  eta=eta+roffset(i)               
                  mean=exp(eta)
                  
                  tmp1=dpoiss(dble(yij),mean,0)
                  cpo(i,1)=cpo(i,1)+1.0d0/tmp1
                  cpo(i,2)=cpo(i,2)+tmp1
                  
                  tmp1=dpoiss(dble(yij),mean,1)
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
      
      do i=1,2
         acrate(i)=acrate(i)/dble(nscan)
      end do
      
      do i=1,nrec
         cpo(i,1)=dble(nsave)/cpo(i,1)
         cpo(i,2)=cpo(i,2)/dble(nsave)
      end do

      do i=1,p
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
         eta=eta+roffset(i)                  
         
         mean=exp(eta)
         dhat=dhat+dpoiss(dble(yij),mean,1)
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
            
      
      return
      end
         
