
c=======================================================================                      
      subroutine dpmmeta(nrec,nfixed,p,                                 
     &                   y,x,sigma2e,                                   
     &                   a0b0,prec,sb,tau,smu,psiinv,                   
     &                   mcmc,nsave,                                     
     &                   ncluster,ss,alpha,beta,b,mu,                   
     &                   sigma,mub,sigmab,mc,                           
     &                   cpo,randsave,thetasave,musave,clustsave,       
     &                   iflagp,workmhp,workmp,workvp,xty,              
     &                   cstrt,ccluster,prob,                           
     &                   seed,betasave,bsave)                           
c=======================================================================                      
c     # of arguments = 40.
c
c     Subroutine `dpmmeta' to run a Markov chain in the semiparametric 
c     meta-analytic linear mixed model using a Dirichlet Process Mixture
c     of normals prior for the distributions of the random effecs. 
c
c     Copyright: Alejandro Jara, 2007-2010.
c
c     Version 1.0:
c
c     Last modification: 16-05-2007.
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
c        nrec        :  integer giving the number of observations.
c        nfixed      :  integer giving the number of fixed effects,
c                       if nfixed is 0 then p=1.
c        p           :  integer giving the number of fixed coefficients.
c        x           :  real matrix giving the design matrix for the 
c                       fixed effects, x(nrec,p). 
c        xtx         :  real matrix givind the product X^tX, xtx(p,p).
c        y           :  real vector giving the response variable,
c                       y(nrec).
c        sigma2e     :  real vector giving the value of the residual
c                       variances, sigma2e(nrec).
c
c-----------------------------------------------------------------------
c
c---- Prior information ------------------------------------------------
c 
c        aa0, ab0    :  real giving the hyperparameters of the prior
c                       distribution for the precision parameter,
c                       alpha ~ Gamma(aa0,ab0). If aa0<0 the precision 
c                       parameter is considered as a constant.
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
c        tau01,tau02 :  reals giving the hyperparameters of the prior 
c                       distribution for the inverse of kernel
c                       variance, 1/sigma ~ Gamma(tau01/2,tau02/2).
c        tau11,tau12 :  reals giving the hyperparameters of the prior 
c                       distribution for the inverse of centering
c                       variance, 1/sigmab ~ Gamma(tau11/2,tau12/2).
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
c                       clustsave(nsave,nrec). 
c        randsave    :  real matrix containing the mcmc samples for
c                       the random effects and prediction,
c                       randsave(nsave,(nrec+1)).
c        musave      :  real matrix containing the mcmc samples for the
c                       cluster locations, musave(nsave,nrec).
c        thetasave   :  real matrix containing the mcmc samples for
c                       the averaged random effects, fixed effects, 
c                       error variance, the normal kernel variance, 
c                       and mean and covariance ofthe baseline 
c                       distribution, the number of clusters, and the 
c                       precision parameter, 
c                       thetsave(nsave,nfixed+6).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        alpha       :  real giving the current value of the precision
c                       parameter of the Dirichlet process.
c        b           :  real matrix giving the current value of the 
c                       random effects, b(nrec).
c        beta        :  real vector giving the current value of the 
c                       fixed effects, beta(p).
c        mu         :   real vector giving the cluster locations
c                       mu(nrec).
c        mub         :  real giving the mean of the normal 
c                       base line distribution for the random effects,
c                       mub.
c        ncluster    :  integer giving the number of clusters in the
c                       random effects.
c        sigma       :  real giving the current value of the
c                       variance for normal kernel.
c        sigmab      :  real giving the current value of the
c                       variance for normal base line 
c                       distribution for the random effects.
c        ss          :  integer vector giving the cluster label for 
c                       each subject, ss(nrec).
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        ccluster    :  integer vector indicating the number of
c                       subjects in each cluster, ccluster(nrec).
c        cstrt       :  integer matrix used to save the cluster
c                       structure, cstrt(nrec,nrec).
c        dispcount   :  index. 
c        dnrm        :  density of a normal distribution.
c        evali       :  integer indicator used in updating the state.
c        i           :  index. 
c        ii          :  index. 
c        iflagp      :  integer vector used to invert the of the lhs
c                       least square solution for the fixed effects,
c                       iflagp(p).
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        k           :  index. 
c        l           :  index.
c        m           :  index.
c        prob        :  real vector used to update the cluster 
c                       structure, prob(nrec+1).
c        ns          :  integer indicator used in updating the state. 
c        nscan       :  integer indicating the total number of MCMC
c                       scans.
c        rgamma      :  real gamma random number generator.
c        sec         :  cpu time working variable.
c        sec0        :  cpu time working variable.
c        sec00       :  cpu time working variable.
c        sec1        :  cpu time working variable.
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        since       :  index.
c        skipcount   :  index. 
c        theta       :  real used to save randomnly generated
c                       random effects, theta.
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
c        tmp3        :  real used to accumulate quantities.
c        workmp      :  real matrix used to update the fixed effects,
c                       workmp(p,p).
c        workmhp     :  real vector used to update the fixed effects,
c                       workmhp(p*(p+1)/2).
c        workvp      :  real vector used to update the fixed effects,
c                       workvp(p).
c        xty         :  real vector used to save the product 
c                       Xt(Y-Zb), xty(p).
c=======================================================================                  
      implicit none 

c+++++Data
      integer nrec,nfixed,p
      double precision y(nrec),x(nrec,p),sigma2e(nrec)

c+++++Prior 
      integer murand,sigmarand
      double precision aa0,ab0,a0b0(2) 
      double precision prec(p,p),sb(p)
      double precision tau01,tau02,tau11,tau12,tau(4)
      double precision smu,psiinv

c+++++MCMC parameters
      integer mcmc(5),nburn,nskip,nsave,ndisplay

c+++++Current values of the parameters
      integer ncluster,ss(nrec)
      double precision alpha
      double precision beta(p)
      double precision b(nrec)
      double precision mu(nrec)
      double precision sigma
      double precision mub
      double precision sigmab

c+++++Output
      integer clustsave(nsave,nrec) 
      double precision cpo(nrec,2)
      double precision randsave(nsave,nrec+1)
      double precision musave(nsave,nrec)
      double precision thetasave(nsave,nfixed+6)

c+++++Seeds
      integer seed(2),seed1,seed2

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++External working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++fixed effects
      integer iflagp(p)
      double precision workmhp(p*(p+1)/2),workmp(p,p),workvp(p)
      double precision xty(p)

c+++++DPM
      integer cstrt(nrec,nrec)
      integer ccluster(nrec)
      double precision prob(nrec+1)

c++++ model's performance
      double precision mc(5)
      double precision betasave(p),bsave(nrec)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Internal working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++General
      integer evali,ii,i,j,k,ns 
      integer ok
      integer since,sprint 
      double precision betar
      double precision theta,tmp1,tmp2,tmp3
      double precision ztz,zty
      double precision sigmainv,sigmabinv

c+++++MCMC
      integer dispcount,isave,iscan,nscan,skipcount 

c+++++RNG and distributions
      double precision dnrm,rnorm,rgamma

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
      murand=mcmc(4)
      sigmarand=mcmc(5)

      tau01=tau(1)
      tau02=tau(2)
      tau11=tau(3)
      tau12=tau(4)

      aa0=a0b0(1)
      ab0=a0b0(2)
      
c++++ set random number generator
      seed1=seed(1)
      seed2=seed(2)
      call setall(seed1,seed2)

c++++ set configurations
      do i=1,nrec
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


      sigmainv=1.d0/sigma
      sigmabinv=1.d0/sigmab

      call cpu_time(sec0)
      sec00=0.d0
      
      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++++
c+++++++ fixed effects
c++++++++++++++++++++++++++++++++++

         if(nfixed.gt.0)then
            do i=1,p
               xty(i)=sb(i)
               workvp(i)=0.d0
               do j=1,p
                  workmp(i,j)=prec(i,j)
               end do
            end do

            do i=1,nrec
               tmp1=y(i)-b(i) 

               do j=1,p
                  xty(j)=xty(j)+x(i,j)*(tmp1/sigma2e(i))
               end do
               
               do j=1,p
                  do k=1,p 
                     workmp(j,k)=workmp(j,k)+x(i,j)*x(i,k)/sigma2e(i)
                  end do   
               end do
            end do

            call inverse(workmp,p,iflagp) 

            do i=1,p
               tmp1=0.d0
               do j=1,p
                  tmp1=tmp1+workmp(i,j)*xty(j) 
               end do
               workvp(i)=tmp1
            end do

            call rmvnorm(p,workvp,workmp,workmhp,xty,beta)
         end if

c++++++++++++++++++++++++++++++++++         
c+++++++ random effects 
c++++++++++++++++++++++++++++++++++

         do ii=1,nrec

c++++++++++ check if the user has requested an interrupt
            call rchkusr()
            ztz=sigmainv+1.d0/sigma2e(ii)
            zty=mu(ss(ii))*sigmainv

            if(nfixed.eq.0)then
               tmp2=y(ii)
              else
               tmp1=0.d0
               do j=1,p 
                  tmp1=tmp1+x(ii,j)*beta(j)
               end do
               tmp2=y(ii)-tmp1
            end if
            zty=zty+tmp2/sigma2e(ii)
            ztz=1.d0/ztz
            tmp1=ztz*zty  
            theta=rnorm(tmp1,sqrt(ztz))
            b(ii)=theta
         end do
         
c+++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ a) Polya Urn based on a collapsed state
c+++++++++++++++++++++++++++++++++++++++++++++++++

         do i=1,nrec
         
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
                  ztz=sigmabinv+dble(ccluster(j))*sigmainv
                  zty=sigmabinv*mub
                  ztz=1.d0/ztz
                  do k=1,ccluster(j)
                     zty=zty+sigmainv*b(cstrt(j,k))   
                  end do 
                  tmp1=ztz*zty
                  ztz=ztz+sigma
                  tmp2=dnrm(b(i),tmp1,sqrt(ztz),1)
                  prob(j)=exp(log(dble(ccluster(j)))+tmp2)
               end do
               
               tmp1=mub
               ztz=sigma+sigmab
               tmp2=dnrm(b(i),tmp1,sqrt(ztz),1)
               prob(ncluster+1)=exp(log(alpha)+tmp2)

               call simdisc(prob,nrec+1,ncluster+1,evali)

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
                   call relabeldpm(i,since,nrec,1,ncluster,
     &                             ccluster,ss,cstrt)                   
               end if

               ccluster(ncluster)=ccluster(ncluster)-1 
               ncluster=ncluster-1

               do j=1,ncluster
                  ztz=sigmabinv+dble(ccluster(j))*sigmainv
                  zty=sigmabinv*mub
                  ztz=1.d0/ztz
                  do k=1,ccluster(j)
                     zty=zty+sigmainv*b(cstrt(j,k))   
                  end do 
                  tmp1=ztz*zty
                  ztz=ztz+sigma
                  tmp2=dnrm(b(i),tmp1,sqrt(ztz),1)
                  prob(j)=exp(log(dble(ccluster(j)))+tmp2)
               end do
               
               tmp1=mub
               ztz=sigma+sigmab
               tmp2=dnrm(b(i),tmp1,sqrt(ztz),1)
               prob(ncluster+1)=exp(log(alpha)+tmp2)

               call simdisc(prob,nrec+1,ncluster+1,evali)

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
            
            ztz=sigmabinv+dble(ns)*sigmainv
            zty=sigmabinv*mub
            ztz=1.d0/ztz
            
            do k=1,ns
               zty=zty+sigmainv*b(cstrt(ii,k))    
            end do 

            tmp1=ztz*zty

            theta=rnorm(tmp1,sqrt(ztz))            
            
            mu(ii)=theta
         end do

c++++++++++++++++++++++++++++++
c+++++++ Kernel variance
c++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         tmp1=0.d0
         do i=1,nrec
            tmp1=tmp1+(b(i)-mu(ss(i)))*(b(i)-mu(ss(i)))
         end do


         sigma=1.d0/
     &         rgamma(0.5d0*(dble(nrec)+tau01),0.5d0*(tmp1+tau02))

         sigmainv=1.d0/sigma 

c++++++++++++++++++++++++++++++++++         
c+++++++ Base line distribution
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(murand.eq.1)then
            zty=smu
            ztz=(sigmabinv*dble(ncluster))+psiinv
            ztz=1.d0/ztz

            do i=1,ncluster
               zty=zty+sigmabinv*mu(i) 
            end do
            tmp1=ztz*zty

            theta=rnorm(tmp1,sqrt(ztz))            
            mub=theta
         end if   

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(sigmarand.eq.1)then
            tmp1=0.d0
            do i=1,ncluster
               tmp1=tmp1+(mu(i)-mub)*(mu(i)-mub)
            end do

            sigmab=1.d0/
     &          rgamma(0.5d0*(dble(ncluster)+tau11),0.5d0*(tmp1+tau12))

            sigmabinv=1.d0/sigmab 
         end if   

c++++++++++++++++++++++++++++++++++         
c+++++++ Precision parameter
c++++++++++++++++++++++++++++++++++
         if(aa0.gt.0.d0)then
            call samalph(alpha,aa0,ab0,ncluster,nrec)
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
               do i=1,ncluster
                  musave(isave,i)=mu(i)
                  clustsave(isave,i)=ccluster(i)
               end do

               do i=1,nrec
                  bsave(i)=bsave(i)+b(i)                  
                  randsave(isave,i)=b(i)
               end do

c+++++++++++++ predictive information

               do i=1,ncluster
                  prob(i)=dble(ccluster(i))/(alpha+dble(nrec))
               end do
               prob(ncluster+1)=alpha/(alpha+dble(nrec))

               call simdisc(prob,nrec+1,ncluster+1,evali)
               
               if(evali.le.ncluster)then
                  theta=mu(evali)
                else
                  theta=rnorm(mub,sqrt(sigmab))  
               end if
               tmp1=rnorm(theta,sqrt(sigma))  
               randsave(isave,nrec+1)=tmp1

c+++++++++++++ functional parameters
               
               tmp1=rbeta(1.d0,alpha+dble(nrec))
               betar=tmp1*theta
               tmp2=tmp1
               weight=(1.d0-tmp1)
               
               do while((1.d0-tmp2).gt.eps)
                  tmp3=rbeta(1.d0,alpha+dble(nrec))
                  tmp1=weight*tmp3
                  weight=weight*(1.d0-tmp3)

                  do i=1,ncluster
                     prob(i)=dble(ccluster(i))/(alpha+dble(nrec))
                  end do
                  prob(ncluster+1)=alpha/(alpha+dble(nrec))

                  call simdisc(prob,nrec+1,ncluster+1,evali)
               
                  if(evali.le.ncluster)then
                     theta=mu(evali)
                   else
                     theta=rnorm(mub,sqrt(sigmab))  
                  end if
                  betar=betar+tmp1*theta
                  tmp2=tmp2+tmp1
               end do

               do i=1,ncluster
                  prob(i)=dble(ccluster(i))/(alpha+dble(nrec))
               end do
               prob(ncluster+1)=alpha/(alpha+dble(nrec))

               call simdisc(prob,nrec+1,ncluster+1,evali)
               
               if(evali.le.ncluster)then
                  theta=mu(evali)
                else
                  theta=rnorm(mub,sqrt(sigmab))
               end if
               
               tmp1=weight
               betar=betar+tmp1*theta

c+++++++++++++ regression coefficients

               thetasave(isave,1)=betar

               if(nfixed.gt.0)then
                  do i=1,p
                     thetasave(isave,1+i)=beta(i)
                     betasave(i)=betasave(i)+beta(i)
                  end do
               end if   

c+++++++++++++ kernel variance
               thetasave(isave,1+nfixed+1)=sigma

c+++++++++++++ baseline mean
               thetasave(isave,1+nfixed+2)=mub

c+++++++++++++ baseline covariance
               thetasave(isave,1+nfixed+3)=sigmab

c+++++++++++++ cluster information
               thetasave(isave,1+nfixed+4)=ncluster
               thetasave(isave,1+nfixed+5)=alpha

c+++++++++++++ cpo
               dbarc=0.d0
               do i=1,nrec
                  tmp1=0.d0
                  if(nfixed.gt.0)then
                     do j=1,p
                        tmp1=tmp1+x(i,j)*beta(j)
                     end do   
                  end if
                  tmp1=tmp1+b(i) 
                  tmp2=dnrm(y(i),tmp1,sqrt(sigma2e(i)),0)
                  cpo(i,1)=cpo(i,1)+1.0d0/tmp2  
                  cpo(i,2)=cpo(i,2)+tmp2                    
                  tmp2=dnrm(y(i),tmp1,sqrt(sigma2e(i)),1)
                  dbarc=dbarc+tmp2
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
      
      do i=1,nrec
         cpo(i,1)=dble(nsave)/cpo(i,1)
         cpo(i,2)=cpo(i,2)/dble(nsave)
      end do

      do i=1,p
         betasave(i)=betasave(i)/dble(nsave)
      end do

      do i=1,nrec
         bsave(i)=bsave(i)/dble(nsave)
      end do   

      dhat=0.d0
      lpml=0.d0
      do i=1,nrec
         tmp1=0.d0
         if(nfixed.gt.0)then
            do j=1,p
               tmp1=tmp1+x(i,j)*betasave(j)
            end do   
         end if
         tmp1=tmp1+bsave(i) 
         dhat=dhat+dnrm(y(i),tmp1,sqrt(sigma2e(i)),1)
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
