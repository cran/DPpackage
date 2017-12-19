
c=======================================================================                      
      subroutine sprasch(datastr,imiss,ngrid,nmissing,nsubject,p,y,
     &                   roffset, 
     &                   a0b0,b0,prec,psiinv,sb,smu,tau1,tau2,
     &                   mcmc,nsave,
     &                   acrate,cpo,cpov,
     &                   cdfsave,randsave,thetasave,
     &                   alpha,b,bclus,beta,mu,ncluster,sigma2,
     &                   sigmainv,ss,
     &                   betac,ccluster,cstrt,fsavet,iflagp,prob,
     &                   seed,
     &                   workmhp,workvp,
     &                   xtx,xty,grid,workcpo)
c=======================================================================                      
c
c     Subroutine `sprasch' to run a Markov chain in the  
c     semiparametric Rasch model. In this routine, inference 
c     is based on  the Polya urn representation of Dirichlet process.
c     The algorithm 8 with m=1 of Neal (2000) is used to sample the 
c     configurations.
c
c     Copyright: Alejandro Jara, 2006 - 2010.
c
c     Version 3.0: 
c
c     Last modification: 05-09-2009.
c
c     Changes and Bug fixes: 
c
c     Version 2.0 to Version 3.0:
c          - Add computation of functionals.
c          - Improve efficiency of DP part.
c          - CDF samples are returned now.
c
c     Version 1.0 to Version 2.0:
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
c        imiss       :  integer indicating whether missing data are
c                       present (1) or absent (0). 
c        ngrid       :  integer giving the size of the grid where
c                       the cdf estimate is evaluated.
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
c        psiinv      :  real giving the prior precision
c                       for the baseline mean, psiinv.
c        sb          :  real vector giving the product of the prior 
c                       precision and prior mean for the difficulty
c                       parameters, sb(p-1).
c        smu         :  real vector giving the product of the prior 
c                       precision and prior mean for the baseline mean,
c                       smu.
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
c        
c-----------------------------------------------------------------------
c
c---- Output -----------------------------------------------------------
c
c        acrate      :  real vector giving the MH acceptance rate, 
c                       acrate(2). 
c        cpo         :  real matrix giving the cpo, cpo(nsubject,p).
c        cdfsave     :  real matrix containing the mcmc samples for the 
c                       cdf, cdfsave(nsave,ngrid).
c        randsave    :  real matrix containing the mcmc samples for
c                       the random effects and prediction,
c                       randsave(nsave,nsubject+1)
c        thetasave   :  real matrix containing the mcmc samples for
c                       the fixed effects,and mean and covariance of
c                       the baseline distribution, 
c                       thetsave(nsave,p+5).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        alpha       :  real giving the current value of the precision
c                       parameter of the Dirichlet process.
c        b           :  real vector giving the current value of the 
c                       random effects, b(nsubject).
c        bclus       :  real vector giving the current value of the 
c                       different values of random effects, 
c                       bclus(nsubject).
c        beta        :  real vector giving the current value of the 
c                       difficulty parameters, beta(p-1).
c        mu          :  real giving the mean of the normal 
c                       base line distribution for the random effects,
c                       mu.
c        ncluster    :  integer giving the number of clusters in the
c                       random effects.
c        sigma       :  real giving the current value of the
c                       standard deviation for normal base line 
c                       distribution for the random effects,
c                       sigma.
c        sigmainv    :  real used to save the inverse of the base line 
c                       variance for the random effects,
c                       sigmainv.
c        ss          :  integer vector giving the cluster label for 
c                       each subject, ss(nsubject).
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        acrate2     :  real used to calculate the acceptance rate. 
c        betac       :  real vector giving the current value of the 
c                       candidate for difficulty parameters, 
c                       betac(p-1).
c        ccluster    :  integer vector indicating the number of
c                       subjects in each cluster, ccluster(nsubject).
c        detlog      :  real used to save the log-determinant in a
c                       matrix inversion process.
c        dispcount   :  index. 
c        evali       :  integer indicator used in updating the state.
c        fsavet      :  real vector used to store the cdf, 
c                       fsavet(ngrid).
c        grid        :  real vector giving the grid where the density
c                       estimate is evaluated, grid(ngrid) .
c        i           :  index. 
c        ii          :  index. 
c        iflag       :  integer vector used to invert the of the lhs
c                       least square solution for the difficulties,
c                       iflag(p-1).
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        k           :  index. 
c        l           :  index.
c        prob        :  real vector used to update the cluster 
c                       structure, prob(nsubject+100).
c        ns          :  integer indicator used in updating the state. 
c        nscan       :  integer indicating the total number of MCMC
c                       scans.
c        runif       :  uniform random number generator.
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        since       :  index.
c        skipcount   :  index. 
c        theta       :  real used to save randomnly generated
c                       random effects, theta.
c        thetac      :  real used to save randomnly generated
c                       random effects, thetac.
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
c        workmhp     :  real vector used to update the difficulties,
c                       workmhp((p-1)*p/2).
c        workvp      :  real vector used to update the difficulties,
c                       workvp(p-1).
c        xtx         :  real matrix givind the product X^tX, 
c                       xtx(p-1,p-1).
c        xty         :  real vector used to save the product 
c                       Xt(Y-Zb), xty(p-1).
c        zty         :  real used to save the product 
c                       Zt(Y-Xbeta), zty.
c        ztz         :  real used to save the product 
c                       ZtSigma^1Z, ztz.
c        ztzinv      :  real used to save the inverted 
c                       ztz, ztzinv.
c=======================================================================                  
      implicit none 

c+++++Data
      integer imiss,ngrid,nmissing,nsubject,p
      integer datastr(nmissing,2),y(nsubject,p)
      double precision roffset(nsubject,p)

c+++++Prior 
      double precision aa0,ab0,a0b0(2),b0(p-1),prec(p-1,p-1),psiinv
      double precision sb(p-1),smu
      double precision tau1,tau2

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Output
      double precision acrate(2)
      double precision cpo(nsubject,p)
      double precision cpov(nsubject)
      double precision cdfsave(nsave,ngrid)
      double precision randsave(nsave,nsubject+1)
      double precision thetasave(nsave,p+5)

c+++++Current values of the parameters
      integer ncluster,ss(nsubject)
      double precision alpha,beta(p-1),b(nsubject)
      double precision bclus(nsubject+100)
      double precision mu,sigma2,sigmainv

c+++++Working space - Loops
      integer ii,i,j,k

c+++++Working space - Random effects
      double precision fsavet(ngrid)
      double precision dnrm,grid(ngrid),thetac
      double precision zty,ztz

c+++++Working space - RNG
      integer seed(2),seed1,seed2
      double precision rgamma,rnorm
      real runif

c+++++Working space - MCMC
      integer iscan,isave,nscan
      integer sprint,skipcount,dispcount
      double precision dbin
      double precision logcgkn,logcgko
      double precision loglikn,logliko
      double precision logpriorn,logprioro
      double precision ratio
      
c+++++Working space - Configurations
      integer ccluster(nsubject),evali,isample,ok 
      integer cstrt(nsubject,nsubject)
      integer ns
      integer since
      double precision prob(nsubject+2)
      double precision tmp1,tmp2,tmp3
      double precision theta

c+++++Working space - Difficulty parameters
      integer iflagp(p-1)
      double precision betac(p-1)
      double precision xtx(p-1,p-1),xty(p-1)
      double precision workmhp((p-1)*p/2)
      double precision workvp(p-1)

c+++++Working space - CPO
      double precision workcpo(nsubject,p)

c+++++Working space - GLM part
      integer yij
      double precision acrate2
      double precision eta,gprime,mean,offset,ytilde

c+++++CPU time
      double precision sec00,sec0,sec1,sec

c++++ DP (functional parameter)
      double precision eps,rbeta,weight
      double precision mureal,sigma2real,mu2real
      parameter(eps=0.01)

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
     
c++++ cluster structure
      do i=1,nsubject
         ccluster(ss(i))=ccluster(ss(i))+1
         cstrt(ss(i),ccluster(ss(i)))=i
      end do
      sigmainv=1.d0/sigma2 
      
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
            do ii=1,nmissing
               i=datastr(ii,1)
               j=datastr(ii,2)
               if(j.eq.1)then
                 eta=b(i)+roffset(i,j)
                else
                 eta=b(i)-beta(j-1)+roffset(i,j)
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

         do i=1,p-1
            do j=1,p-1
               xtx(i,j)=prec(i,j)
            end do
            xty(i)=sb(i)
         end do

         logliko=0.d0         
         
         do i=1,nsubject
            do j=1,p-1
               yij=y(i,j+1)
               eta=b(i)-beta(j)+roffset(i,j+1) 
               offset=b(i)+roffset(i,j+1) 

               mean=exp(eta)/(1.d0+exp(eta))
               logliko=logliko+dbin(dble(yij),1.d0,mean,1)

               tmp1=mean*(1.0d0-mean)
               gprime=1.d0/tmp1

               ytilde=eta+(dble(yij)-mean)*gprime-offset
               xtx(j,j)=xtx(j,j)+1.d0/gprime
               xty(j)=xty(j)-ytilde/gprime

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
            do j=1,p-1
               yij=y(i,j+1)
               eta=b(i)-betac(j)+roffset(i,j+1)  
               offset=b(i)+roffset(i,j+1) 
               mean=exp(eta)/(1.d0+exp(eta))
               loglikn=loglikn+dbin(dble(yij),1.d0,mean,1)

               tmp1=mean*(1.0d0-mean)
               gprime=1.d0/tmp1

               ytilde=eta+(dble(yij)-mean)*gprime-offset
               xtx(j,j)=xtx(j,j)+1.d0/gprime
               xty(j)=xty(j)-ytilde/gprime
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
     &         logpriorn-logprioro

         if(log(runif()).lt.ratio)then
            acrate(1)=acrate(1)+1.d0
            do i=1,p-1
               beta(i)=betac(i) 
            end do
         end if
            
c+++++++ check if the user has requested an interrupt
         call rchkusr()


c++++++++++++++++++++++++++++++++++         
c+++++++ random effects 
c++++++++++++++++++++++++++++++++++

c++++++++++++++++++++++++++++++
c+++++++ a) Polya Urn 
c++++++++++++++++++++++++++++++

         do i=1,nsubject
         
            ns=ccluster(ss(i))

            if(ns.gt.1)then
               isample=1
               ccluster(ss(i))=ccluster(ss(i))-1 

               j=1
               ok=0
               do while(ok.eq.0.and.j.le.ns)
                  if(cstrt(ss(i),j).eq.i)ok=j
                  j=j+1
               end do

               do j=ok,ns-1
                  cstrt(ss(i),j)=cstrt(ss(i),j+1)
               end do

             else
               isample=2
               since=ss(i)
    
               if(since.lt.ncluster)then

                  call relabelrasch(i,since,nsubject,ncluster,
     &                              ccluster,ss,cstrt,bclus)

               end if
               ccluster(ncluster)=ccluster(ncluster)-1 
               ncluster=ncluster-1
            end if 

            do j=1,ncluster
                  
               tmp1=0.d0
                  
               yij=y(i,1)
               eta=bclus(j)+roffset(i,1) 
               mean=exp(eta)/(1.d0+exp(eta))
               tmp1=tmp1+dbin(dble(yij),1.d0,mean,1)

               do k=1,p-1
                  yij=y(i,k+1)
                  eta=bclus(j)-beta(k)+roffset(i,k+1) 
                  mean=exp(eta)/(1.d0+exp(eta))
                  tmp1=tmp1+dbin(dble(yij),1.d0,mean,1)
               end do

               prob(j)=exp(log(dble(ccluster(j)))+
     &                        tmp1)
            end do
               
            if(isample.eq.1)then
               theta=rnorm(mu,sqrt(sigma2))
               bclus(ncluster+1)=theta
            end if 
            
            theta=bclus(ncluster+1)
            tmp1=0.d0
               
            yij=y(i,1)
            eta=theta+roffset(i,1) 
            mean=exp(eta)/(1.d0+exp(eta))
            tmp1=tmp1+dbin(dble(yij),1.d0,mean,1)

            do k=1,p-1
               yij=y(i,k+1)
               eta=theta-beta(k)+roffset(i,k+1) 
               mean=exp(eta)/(1.d0+exp(eta))
               tmp1=tmp1+dbin(dble(yij),1.d0,mean,1)
            end do

            prob(ncluster+1)=alpha*exp(tmp1)
               
            call simdisc(prob,nsubject+100,ncluster+1,evali)

            ss(i)=evali
            ccluster(evali)=ccluster(evali)+1
            cstrt(evali,ccluster(evali))=i

            if(evali.gt.ncluster)then
               ncluster=ncluster+1
            end if

         end do


c++++++++++++++++++++++++++++++
c+++++++ b) Resampling step
c++++++++++++++++++++++++++++++

         acrate2=0.d0
         
         do i=1,ncluster

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            ztz=1.d0/sigma2
            zty=mu/sigma2

            logliko=0.d0                     

            ns=ccluster(i)
            do j=1,ns

c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

               ii=cstrt(i,j)

               do k=1,p
                  if(k.eq.1)then
                     eta=bclus(i)+roffset(ii,k)
                     offset=roffset(ii,k)
                    else
                     eta=bclus(i)-beta(k-1)+roffset(ii,k)
                     offset=roffset(ii,k)-beta(k-1)
                  end if
                    
                  yij=y(ii,k)

                  mean=exp(eta)/(1.d0+exp(eta))
                  logliko=logliko+dbin(dble(yij),1.d0,mean,1)

                  tmp1=mean*(1.0d0-mean)
                  gprime=1.d0/tmp1

                  ytilde=eta+(dble(yij)-mean)*gprime-offset
                  ztz=ztz+1.d0/gprime
                  zty=zty+ytilde/gprime

               end do
            end do

            ztz=1.d0/ztz
            zty=ztz*zty

            thetac=rnorm(zty,sqrt(ztz))

c++++++++++ evaluating the candidate generating kernel

            logcgko=dnrm(thetac,zty,sqrt(ztz),1)

c++++++++++ prior ratio

            logprioro=0.d0
            logpriorn=0.d0
         
            logpriorn=(thetac-mu)* 
     &            sigmainv*
     &           (thetac-mu)

            logprioro=(bclus(i)-mu)* 
     &            sigmainv*
     &           (bclus(i)-mu)
            
            logpriorn=-0.5d0*logpriorn
            logprioro=-0.5d0*logprioro


c++++++++++ candidate generating kernel contribution

            ztz=1.d0/sigma2
            zty=mu/sigma2

            loglikn=0.d0                     
            
            ns=ccluster(i)
            do j=1,ns

c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

               ii=cstrt(i,j)

               do k=1,p
                  if(k.eq.1)then
                     eta=thetac+roffset(ii,k)
                     offset=roffset(ii,k)
                    else
                     eta=thetac-beta(k-1)+roffset(ii,k)
                     offset=roffset(ii,k)-beta(k-1)
                  end if
                    
                  yij=y(ii,k)

                  mean=exp(eta)/(1.d0+exp(eta))
                  loglikn=loglikn+dbin(dble(yij),1.d0,mean,1)

                  tmp1=mean*(1.0d0-mean)
                  gprime=1.d0/tmp1

                  ytilde=eta+(dble(yij)-mean)*gprime-offset
                  ztz=ztz+1.d0/gprime
                  zty=zty+ytilde/gprime
               end do
            end do

            ztz=1.d0/ztz
            zty=ztz*zty

c++++++++++ evaluating the candidate generating kernel

            theta=bclus(i)
            logcgkn=dnrm(theta,zty,sqrt(ztz),1)                 

c++++++++++ mh step
           
            ratio=loglikn-logliko+logcgkn-logcgko+
     &            logpriorn-logprioro

            if(log(runif()).lt.ratio)then
               acrate2=acrate2+1.d0
               bclus(i)=thetac
            end if
         end do

         acrate(2)=acrate(2)+acrate2/dble(ncluster)

         do i=1,nsubject
            b(i)=bclus(ss(i))
         end do


c++++++++++++++++++++++++++++++++++         
c+++++++ Base line distribution
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(psiinv.gt.0.d0)then
            tmp2=smu
            tmp1=1.d0/(sigmainv*dble(ncluster)+psiinv)
            do i=1,ncluster
               tmp2=tmp2+sigmainv*bclus(i)
            end do
            mean=tmp1*tmp2
            mu=rnorm(mean,sqrt(tmp1))
         end if 

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(tau1.gt.0.d0)then
            tmp1=0.d0
            do i=1,ncluster
               tmp1=tmp1+(bclus(i)-mu)*(bclus(i)-mu)                   
            end do

            sigma2=1.d0/
     &             rgamma(0.5d0*(dble(ncluster)+tau1),0.5d0*(tmp1+tau2))
     
            sigmainv=1.d0/sigma2
         end if   
         
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

c+++++++++++++ difficulty parameters

               do i=1,p-1
                  thetasave(isave,i)=beta(i)
               end do

c+++++++++++++ Partially sampling the DP for mean and variance

               do i=1,ngrid
                  fsavet(i)=0.d0
               end do

               mureal=0.d0 
               mu2real=0.d0
               sigma2real=0.d0

               do i=1,ncluster
                   prob(i)=real(ccluster(i))/(alpha+real(nsubject))
               end do
               prob(ncluster+1)=alpha/(alpha+dble(nsubject))
               call simdisc(prob,nsubject+100,ncluster+1,evali)

               if(evali.le.ncluster)then
                  theta=bclus(evali)
               end if
               if(evali.eq.ncluster+1)then 
                  theta=rnorm(mu,sqrt(sigma2))
               end if

               tmp1=rbeta(1.d0,alpha+dble(nsubject))
               tmp2=tmp1
               weight=(1.d0-tmp1)

               mureal=tmp1*theta 
               mu2real=tmp1*theta*theta
               do j=1,ngrid
                  if(theta.le.grid(j))then
                     fsavet(j)=fsavet(j)+tmp1
                  end if
               end do

               do while((1.0-tmp2).gt.eps)
                  call rchkusr()

                  tmp3=rbeta(1.d0,alpha+real(nsubject))
                  tmp1=weight*tmp3
                  weight=weight*(1.0-tmp3)

                  call simdisc(prob,nsubject+100,ncluster+1,evali)

                  if(evali.le.ncluster)then
                     theta=bclus(evali)
                  end if
                  if(evali.eq.ncluster+1)then 
                     theta=rnorm(mu,sqrt(sigma2))
                  end if

                  mureal=mureal+tmp1*theta 
                  mu2real=mu2real+tmp1*theta*theta

                  do j=1,ngrid
                     if(theta.le.grid(j))then
                        fsavet(j)=fsavet(j)+tmp1
                     end if
                  end do
                  tmp2=tmp2+tmp1
               end do

               tmp3=1.d0
               tmp1=weight*tmp3
               call simdisc(prob,nsubject+100,ncluster+1,evali)

               if(evali.le.ncluster)then
                  theta=bclus(evali)
               end if
               if(evali.eq.ncluster+1)then 
                  theta=rnorm(mu,sqrt(sigma2))
               end if

               mureal=mureal+tmp1*theta 
               mu2real=mu2real+tmp1*theta*theta

               do j=1,ngrid
                  if(theta.le.grid(j))then
                     fsavet(j)=fsavet(j)+tmp1
                  end if
               end do

               sigma2real=mu2real-mureal*mureal

               thetasave(isave,p)=mureal
               thetasave(isave,p+1)=sigma2real

               do i=1,ngrid
                  cdfsave(isave,i)=fsavet(i)
               end do

c+++++++++++++ baseline mean

               thetasave(isave,p+2)=mu

c+++++++++++++ baseline stansard deviation

               thetasave(isave,p+3)=sigma2

c+++++++++++++ cluster information

               thetasave(isave,p+4)=ncluster
               thetasave(isave,p+5)=alpha

c+++++++++++++ random effects

               do i=1,nsubject
                  randsave(isave,i)=b(i)
               end do

c+++++++++++++ predictive information

               do i=1,ncluster
                  prob(i)=dble(ccluster(i))/(alpha+dble(nsubject))
               end do
               prob(ncluster+1)=alpha/(alpha+dble(nsubject))

               call simdisc(prob,nsubject+100,ncluster+1,evali)
               
               if(evali.le.ncluster)then
                  theta=bclus(evali)
               end if
               if(evali.gt.ncluster)then
                  theta=rnorm(mu,sqrt(sigma2))
               end if
               randsave(isave,nsubject+1)=theta

c+++++++++++++ cpo

               tmp3=0.d0
               do i=1,nsubject
                  tmp2=0.d0
                  do j=1,p
                     yij=y(i,j)
                     if(j.eq.1)then
                       eta=b(i)+roffset(i,j)
                      else
                       eta=b(i)-beta(j-1)+roffset(i,j)
                     end if  

                     mean=exp(eta)/(1.d0+exp(eta))
                     tmp1=dbin(dble(yij),1.d0,mean,1)
                     cpo(i,j)=cpo(i,j)+1.0d0/exp(tmp1)   
                     tmp2=tmp2+tmp1
                     tmp3=tmp3+log(dble(isave)/cpo(i,j))
                  end do
                  cpov(i)=cpov(i)+1.d0/exp(tmp2)
               end do

c               call dblepr("LPML",-1,tmp3,1)

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
      
      acrate(1)=acrate(1)/dble(nscan)    
      acrate(2)=acrate(2)/dble(nscan)    
      
      do i=1,nsubject
         do j=1,p
            cpo(i,j)=dble(nsave)/cpo(i,j)
         end do   
         cpov(i)=dble(nsave)/cpov(i)  
      end do
             
      return
      end


