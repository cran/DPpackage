
c=======================================================================                      
      subroutine spraschpoi(datastr,imiss,ngrid,nmissing,nsubject,p,y,
     &                      roffset, 
     &                      a0b0,b0,prec,psiinv,sb,smu,tau1,tau2,
     &                      mcmc,nsave,
     &                      acrate,cpo,faccum,randsave,thetasave,
     &                      alpha,b,bclus,beta,mu,ncluster,sigma,
     &                      sigmainv,ss,
     &                      betac,ccluster,fsavet,iflag,prob,seed,
     &                      work1,work2,work3,
     &                      workmh1,workv1,workv2,workv3,
     &                      xtx,xty,grid)
c=======================================================================                      
c
c     Subroutine `spraschpoi' to run a Markov chain in the  
c     semiparametric Rasch Poissin Count model. In this routine, inference 
c     is based on  the Polya urn representation of Dirichlet process.
c     The algorithm 8 with m=1 of Neal (2000) is used to sample the 
c     configurations.
c
c     Copyright: Alejandro Jara Vallejos, 2006 - 2007
c
c     Version 3.0: 
c
c     Last modification: 01-02-2007.
c
c     Changes and Bug fixes: 
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
c        faccum      :  real vector giving the cdf estimate at the
c                       grid, faccum(ngrid).
c        randsave    :  real matrix containing the mcmc samples for
c                       the random effects and prediction,
c                       randsave(nsave,nsubject+1)
c        thetasave   :  real matrix containing the mcmc samples for
c                       the fixed effects,and mean and covariance of
c                       the baseline distribution, 
c                       thetsave(nsave,p+3).
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
c        fsavet      :  integer vector used to store the number of
c                       observations to calculate the predictive
c                       cdf distribution, fsavet(ngrid).
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
c                       structure, prob(nsubject+2).
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
      real*8 roffset(nsubject,p)

c+++++Prior 
      real*8 aa0,ab0,a0b0(2),b0(p-1),prec(p-1,p-1),psiinv
      real*8 sb(p-1),smu
      real*8 tau1,tau2

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Output
      real*8 acrate(2)
      real*8 cpo(nsubject,p),faccum(ngrid)
      real*8 randsave(nsave,nsubject+1)
      real*8 thetasave(nsave,p+3)

c+++++Current values of the parameters
      integer ncluster,ss(nsubject)
      real*8 alpha,beta(p-1),b(nsubject)
      real*8 bclus(nsubject)
      real*8 mu,sigma,sigmainv

c+++++Working space - Loops
      integer ii,i,j,k

c+++++Working space - Random effects
      integer fsavet(ngrid)
      real*8 cdfnorm,dnrm,grid(ngrid),thetac
      real*8 zty,ztz,ztzinv

c+++++Working space - RNG
      integer rpois,seed(2),seed1,seed2
      real*8 rgamma,rnorm
      real runif

c+++++Working space - MCMC
      integer iscan,isave,nscan
      integer sprint,skipcount,dispcount
      real*8 dpoiss
      real*8 logcgkn,logcgko
      real*8 loglikn,logliko
      real*8 logpriorn,logprioro
      real*8 ratio
      
c+++++Working space - Configurations
      integer ccluster(nsubject),evali 
      integer ns
      integer since
      real*8 prob(nsubject+2)
      real*8 tmp1,tmp2
      real*8 theta

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
     
c++++ cluster structure
      do i=1,nsubject
         ccluster(ss(i))=ccluster(ss(i))+1
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
         
         call invdet(xtx,p-1,work2,detlog,iflag,workv1)
            
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
         
         call invdet(xtx,p-1,work2,detlog,iflag,workv1)

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


c++++++++++++++++++++++++++++++++++         
c+++++++ random effects 
c++++++++++++++++++++++++++++++++++

c++++++++++++++++++++++++++++++
c+++++++ a) Polya Urn 
c++++++++++++++++++++++++++++++

         do i=1,nsubject
         
            ns=ccluster(ss(i))

c++++++++++ observation in cluster with more than 1 element
             
            if(ns.gt.1)then
 
               ccluster(ss(i))=ccluster(ss(i))-1 
               
               do j=1,ncluster
                  
                  tmp1=0.d0
                  
                  yij=y(i,1)
                  eta=bclus(j)+roffset(i,1) 
                  mean=exp(eta)
                  tmp1=tmp1+dpoiss(dble(yij),mean,1)

                  do k=1,p-1
                     yij=y(i,k+1)
                     eta=bclus(j)-beta(k)+roffset(i,k+1) 
                     mean=exp(eta)
                     tmp1=tmp1+dpoiss(dble(yij),mean,1)
                  end do

                  prob(j)=exp(log(dble(ccluster(j)))+
     &                        tmp1)
               end do
               
               theta=rnorm(mu,sigma)

               tmp1=0.d0
               
               yij=y(i,1)
               eta=theta+roffset(i,1) 
               mean=exp(eta)
               tmp1=tmp1+dpoiss(dble(yij),mean,1)

               do k=1,p-1
                  yij=y(i,k+1)
                  eta=theta-beta(k)+roffset(i,k+1) 
                  mean=exp(eta)
                  tmp1=tmp1+dpoiss(dble(yij),mean,1)
               end do

               prob(ncluster+1)=exp(log(alpha)+tmp1)
               
               call simdisc(prob,nsubject+2,ncluster+1,evali)

               if(evali.le.ncluster)then
                  ss(i)=evali
                  ccluster(evali)=ccluster(evali)+1
               end if   
               if(evali.gt.ncluster)then
                  ncluster=ncluster+1
                  ss(i)=ncluster
                  ccluster(ncluster)=1
                  bclus(ncluster)=theta
               end if               
            end if

c++++++++++ subject in cluster with only 1 observation

            if(ns.eq.1)then
                
               since=ss(i)
                
               if(since.lt.ncluster)then
                   call relabelraschpoi(i,since,nsubject,ncluster,
     &                          ccluster,ss,bclus,theta)                   
	       end if

               ccluster(ncluster)=ccluster(ncluster)-1 
               ncluster=ncluster-1

               do j=1,ncluster
                  
                  tmp1=0.d0
                  yij=y(i,1)
                  eta=bclus(j)+roffset(i,1) 
                  mean=exp(eta)
                  tmp1=tmp1+dpoiss(dble(yij),mean,1)

                  do k=1,p-1
                     yij=y(i,k+1)
                     eta=bclus(j)-beta(k)+roffset(i,k+1) 
                     mean=exp(eta)
                     tmp1=tmp1+dpoiss(dble(yij),mean,1)
                  end do

                  prob(j)=exp(log(dble(ccluster(j)))+
     &                        tmp1)
               end do
               
               theta=b(i)

               tmp1=0.d0
               
               yij=y(i,1)
               eta=theta+roffset(i,1) 
               mean=exp(eta)
               tmp1=tmp1+dpoiss(dble(yij),mean,1)

               do k=1,p-1
                  yij=y(i,k+1)
                  eta=theta-beta(k)+roffset(i,k+1) 
                  mean=exp(eta)
                  tmp1=tmp1+dpoiss(dble(yij),mean,1)
               end do

               prob(ncluster+1)=exp(log(alpha)+tmp1)
               
               call simdisc(prob,nsubject+2,ncluster+1,evali)

               if(evali.le.ncluster)then
                  ss(i)=evali
                  ccluster(evali)=ccluster(evali)+1
               end if   
               if(evali.gt.ncluster)then
                  ncluster=ncluster+1
                  ss(i)=ncluster
                  ccluster(ncluster)=1
                  bclus(ncluster)=theta
               end if      
	    
	    end if

         end do


c++++++++++++++++++++++++++++++
c+++++++ b) Resampling step
c++++++++++++++++++++++++++++++

         acrate2=0.d0
         
         do ii=1,ncluster

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            ztz=0.d0
            ztzinv=0.d0
            zty=0.d0

            logliko=0.d0                     
            
            do i=1,nsubject
               if(ii.eq.ss(i))then
                 do j=1,p
                    if(j.eq.1)then
                       eta=bclus(ii)+roffset(i,j)
                       offset=roffset(i,j)
                      else
                       eta=bclus(ii)-beta(j-1)+roffset(i,j)
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
               end if
            end do

            ztz=ztz+sigmainv
            ztzinv=1.d0/ztz

            zty=zty+sigmainv*mu

            tmp2=ztzinv*zty
 
            thetac=rnorm(tmp2,sqrt(ztzinv))

c++++++++++ evaluating the candidate generating kernel

            logcgko=dnrm(thetac,tmp2,sqrt(ztzinv),1)

c++++++++++ prior ratio

            logprioro=0.d0
            logpriorn=0.d0
         
            logpriorn=(thetac-mu)* 
     &            sigmainv*
     &           (thetac-mu)

            logprioro=(bclus(ii)-mu)* 
     &            sigmainv*
     &           (bclus(ii)-mu)
            
            logpriorn=-0.5d0*logpriorn
            logprioro=-0.5d0*logprioro


c++++++++++ candidate generating kernel contribution

            ztz=0.d0
            ztzinv=0.d0
            zty=0.d0

            loglikn=0.d0                     
            
            do i=1,nsubject
               if(ii.eq.ss(i))then
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
               end if
            end do
   
            ztz=ztz+sigmainv
            ztzinv=1.d0/ztz

            zty=zty+sigmainv*mu

            tmp2=ztzinv*zty

c++++++++++ evaluating the candidate generating kernel

            theta=bclus(ii)
            logcgkn=dnrm(theta,tmp2,sqrt(ztzinv),1)                 

c++++++++++ mh step
           
            ratio=dexp(loglikn-logliko+logcgkn-logcgko+
     &                 logpriorn-logprioro)

            if(dble(runif()).lt.ratio)then
               acrate2=acrate2+1.d0
               bclus(ii)=thetac
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

            tmp2=1.d0/
     &             rgamma(0.5d0*(dble(ncluster)+tau1),0.5d0*(tmp1+tau2))
     
            sigma=sqrt(tmp2)
            sigmainv=1.d0/tmp2
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

c+++++++++++++ baseline mean

               thetasave(isave,p)=mu

c+++++++++++++ baseline stansard deviation

               thetasave(isave,p+1)=sigma

c+++++++++++++ cluster information

               thetasave(isave,p+2)=ncluster
               thetasave(isave,p+3)=alpha

c+++++++++++++ random effects

               do i=1,nsubject
                  randsave(isave,i)=b(i)
               end do


c+++++++++++++ predictive information

               do i=1,ngrid
                  fsavet(i)=0
               end do

               do i=1,ncluster
                  prob(i)=dble(ccluster(i))/(alpha+dble(nsubject))
                  
                  j=1
                  do while(bclus(i).le.grid(j).and.j.le.ngrid)
                     fsavet(j)=fsavet(j)+ccluster(i)
                     j = j + 1
                  end do

               end do

               do j=1,ngrid
                  faccum(j)=faccum(j)+
     &             ( alpha*cdfnorm(grid(j),mu,sigma,1,0) +
     &               dble(fsavet(j)))/
     &              (alpha+dble(nsubject))  

               end do

               prob(ncluster+1)=alpha/(alpha+dble(nsubject))

               call simdisc(prob,nsubject+2,ncluster+1,evali)
               
               if(evali.le.ncluster)then
                  theta=bclus(evali)
               end if
               if(evali.gt.ncluster)then
                  theta=rnorm(mu,sigma)
               end if
               
               randsave(isave,nsubject+1)=theta

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
      
      acrate(1)=acrate(1)/dble(nscan)    
      acrate(2)=acrate(2)/dble(nscan)    
      
      do i=1,nsubject
         do j=1,p
            cpo(i,j)=dble(nsave)/cpo(i,j)
         end do   
      end do

      do i=1,ngrid
         faccum(i)=faccum(i)/dble(nsave)
      end do
             
      return
      end


c=======================================================================      
      subroutine relabelraschpoi(ind,since,nsubject,ncluster,ccluster,
     &                           ss,bclus,theta)
c=======================================================================
c     relabel the clusters after the elimination of one of them
c     A.J.V., 2006
      implicit none
      integer i,j,ind,since,nsubject,ncluster,ccluster(nsubject)
      integer ss(nsubject)
      real*8 bclus(nsubject),theta

      theta=bclus(since)
      
      do i=since+1,ncluster
         do j=1,nsubject
            if(ss(j).eq.i)then
               ss(j)=i-1 
            end if
         end do
         bclus(i-1)=bclus(i)
         ccluster(i-1)=ccluster(i)
      end do
      
      ss(ind)=ncluster
      
      bclus(ncluster)=theta
      ccluster(ncluster)=1
      
      return
      end  
      