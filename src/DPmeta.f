
c=======================================================================                      
      subroutine dpmeta(nrec,nfixed,p,                                  
     &                  x,y,sigma2e,                                    
     &                  a0b0,prec,psiinv,sb,smu,tau,                    
     &                  mcmc,nsave,randsave,thetasave,cpo,alpha,b,      
     &                  bclus,beta,betar,mu,ncluster,sigma,ss,mc,       
     &                  cstrt,ccluster,iflag,prob,                      
     &                  res,seed,work1,workmh1,workv1,xty,              
     &                  betasave,bsave)                                 
c=======================================================================                   
c     # of arguments = 39.
c
c     Subroutine `dpmeta' to run a Markov chain in the semiparametric 
c     meta-analytic linear mixed model using a Dirichlet Process prior 
c     for the distributions of the random effecs. 
c     In this routine, inference is 
c     based on the Polya urn representation of Dirichlet process.
c
c     Copyright: Alejandro Jara, 2007-2010.
c
c     Version 1.0: 
c
c     Last modification: 16-04-2007.
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
c        z           :  real matrix giving the design matrix for the 
c                       random effects, z(nrec,q). 
c        sigma2e     :  real vector giving the value of the residual
c                       variances, sigma2e(nrec).
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
c        psiinv      :  real giving the prior precision matrix
c                       for the baseline mean.
c        sb          :  real vector giving the product of the prior 
c                       precision and prior mean for the fixed effects,
c                       sb(p).
c        smu         :  real giving the product of the prior 
c                       precision and prior mean for the baseline mean,
c                       smu.
c        tau1, tau2  :  reals giving the hyperparameters of the prior 
c                       distribution for the inverse of the baseline
c                       variance, 1/sigma ~ Gamma(tau1/2,tau2/2).
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
c        randsave    :  real matrix containing the mcmc samples for
c                       the random effects and prediction,
c                       randsave(nsave,nrec+1).
c        thetasave   :  real matrix containing the mcmc samples for
c                       the averaged random effects, fixed effects, 
c                       error variance, and mean and covariance of
c                       the baseline distribution, 
c                       thetsave(nsave,1+nfixed+1+1+1+2).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        alpha       :  real giving the current value of the precision
c                       parameter of the Dirichlet process.
c        b           :  real vector giving the current value of the 
c                       random effects, b(nrec).
c        bclus       :  real vector giving the current value of the 
c                       different values of random effects, 
c                       bclus(nrec).
c        beta        :  real vector giving the current value of the 
c                       fixed effects, beta(p).
c        betar       :  real giving the current value of the 
c                       averaged random effects.
c        mu          :  real giving the mean of the normal 
c                       base line distribution for the random effects.
c        ncluster    :  integer giving the number of clusters in the
c                       random effects.
c        sigma       :  real giving the current value of the
c                       variance for normal base line 
c                       distribution for the random effects.
c        ss          :  integer vector giving the cluster label for 
c                       each subject, ss(nrec).
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        ccluster    :  integer vector indicating the number of
c                       subjects in each cluster, ccluster(nrec).
c        cstrt       :  integer matrix used to save the cluster
c                       structure, cstrt(nrec,nrec).
c        detlog      :  real used to save the log-determinant in a
c                       matrix inversion process.
c        dispcount   :  index. 
c        evali       :  integer indicator used in updating the state.
c        i           :  index. 
c        ii          :  index. 
c        iflag       :  integer vector used to invert the of the lhs
c                       least square solution for the fixed effects,
c                       iflag(p).
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        prob        :  real vector used to update the cluster 
c                       structure, prob(nrec+1).
c        ns          :  integer indicator used in updating the state. 
c        nscan       :  integer indicating the total number of MCMC
c                       scans.
c        res         :  real vector used to save the residual effects,
c                       res(nrec).
c        rgamma      :  real gamma random number generator.
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        since       :  index.
c        skipcount   :  index. 
c        sse         :  real used to save the SS of the errors.
c        theta       :  real used to save randomnly generated
c                       random effects.
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
c        tmp3        :  real used to accumulate quantities.
c        tpi         :  real parameter used to evaluate the normal 
c                       density.
c        work1       :  real matrix used to update the fixed effects,
c                       work1(p,p).
c        workmh1     :  real vector used to update the fixed effects,
c                       workmh1(p*(p+1)/2)
c        workv1      :  real vector used to update the fixed effects,
c                       workv1(p)
c        xty         :  real vector used to save the product 
c                       Xt(Y-Zb), xty(p).
c=======================================================================                  
      implicit none 

c+++++Data
      integer nrec,nfixed,p
      double precision y(nrec),x(nrec,p),sigma2e(nrec)
     
c+++++Prior 
      integer murand,sigmarand
      double precision aa0,ab0,a0b0(2),prec(p,p),psiinv
      double precision sb(p),smu
      double precision tau1,tau2,tau(2)

c+++++MCMC parameters
      integer mcmc(5),nburn,nskip,nsave,ndisplay

c+++++Output
      double precision cpo(nrec,2)
      double precision randsave(nsave,nrec+1)
      double precision thetasave(nsave,nfixed+5)

c+++++Current values of the parameters
      integer ncluster,ss(nrec)
      double precision alpha,beta(p),b(nrec)
      double precision betar,bclus(nrec)
      double precision mu,sigma

c+++++Seeds
      integer seed(2),seed1,seed2

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++External working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++fixed effects
      integer iflag(p)
      double precision work1(p,p)
      double precision workmh1(p*(p+1)/2)
      double precision workv1(p)
      double precision xty(p)

c+++++DP
      integer cstrt(nrec,nrec)
      integer ccluster(nrec)
      double precision prob(nrec+1)

c+++++Residuals
      double precision res(nrec)

c++++ model's performance
      double precision mc(5)
      double precision betasave(p),bsave(nrec)
      
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Internal working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++General
      integer evali,ii,i,j,k,ns,ok 
      integer since,sprint 
      double precision sigmainv,sse,theta
      double precision tmp1,tmp2,tmp3
      double precision ztz,zty
      
c+++++MCMC
      integer dispcount,isave,iscan,nscan,skipcount 

c+++++RNG and distributions
      double precision dnrm,rgamma,rnorm

c+++++DP
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
      
      tau1=tau(1)
      tau2=tau(2)
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

      sigmainv=1.d0/sigma

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

c++++++++++++++++++++++++++++++++++
c+++++++ fixed effects
c++++++++++++++++++++++++++++++++++

         if(nfixed.gt.0)then
            do i=1,p
               xty(i)=sb(i)
               workv1(i)=0.d0
               do j=1,p
                  work1(i,j)=prec(i,j)
               end do
            end do

            do i=1,nrec
               tmp1=y(i)-b(i) 

               do j=1,p
                  xty(j)=xty(j)+x(i,j)*(tmp1/sigma2e(i))
               end do
               
               do j=1,p
                  do k=1,p 
                     work1(j,k)=work1(j,k)+x(i,j)*x(i,k)/sigma2e(i)
                  end do   
               end do
            end do

            call inverse(work1,p,iflag) 

            do i=1,p
               tmp1=0.d0
               do j=1,p
                  tmp1=tmp1+work1(i,j)*xty(j) 
               end do
               workv1(i)=tmp1
            end do

            call rmvnorm(p,workv1,work1,workmh1,xty,beta)
         end if

         
c++++++++++++++++++++++++++++++++++         
c+++++++ random effects 
c++++++++++++++++++++++++++++++++++

c++++++++++++++++++++++++++++++
c+++++++ a) Polya Urn 
c++++++++++++++++++++++++++++++

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
                  tmp1=dnrm(res(i),bclus(j),sqrt(sigma2e(i)),1)
                  prob(j)=exp(log(dble(ccluster(j)))+tmp1)
               end do
               tmp1=dnrm(res(i),mu,sqrt(sigma2e(i)+sigma),1)
               prob(ncluster+1)=exp(log(alpha)+tmp1)

               call simdisc(prob,nrec+1,ncluster+1,evali)

               ss(i)=evali
               
               ccluster(evali)=ccluster(evali)+1

               cstrt(evali,ccluster(evali))=i
               
               if(evali.gt.ncluster)then
                  ncluster=ncluster+1
                  ztz=1.d0/sigma2e(i)
                  zty=res(i)/sigma2e(i)
                  ztz=1.d0/(ztz+sigmainv)
                  zty=zty+sigmainv*mu   
                  tmp1=ztz*zty  
                  theta=rnorm(tmp1,sqrt(ztz))
                  bclus(evali)=theta
               end if
            end if


c++++++++++ subject in cluster with only 1 observation
             
            if(ns.eq.1)then
                
               since=ss(i)
                
               if(since.lt.ncluster)then
                   call relabelmeta(i,since,nrec,ncluster,
     &                             cstrt,ccluster,ss,bclus)
               end if

               ccluster(ncluster)=ccluster(ncluster)-1 
               ncluster=ncluster-1

               do j=1,ncluster
                  tmp1=dnrm(res(i),bclus(j),sqrt(sigma2e(i)),1)
                  prob(j)=exp(log(dble(ccluster(j)))+tmp1)
               end do
               tmp1=dnrm(res(i),mu,sqrt(sigma2e(i)+sigma),1)
               prob(ncluster+1)=exp(log(alpha)+tmp1)

               call simdisc(prob,nrec+1,ncluster+1,evali)
               
               ss(i)=evali
               
               ccluster(evali)=ccluster(evali)+1

               cstrt(evali,ccluster(evali))=i
               
               if(evali.gt.ncluster)then
                  ncluster=ncluster+1
                  ztz=1.d0/sigma2e(i)
                  zty=res(i)/sigma2e(i)
                  ztz=1.d0/(ztz+sigmainv)
                  zty=zty+sigmainv*mu   
                  tmp1=ztz*zty  
                  theta=rnorm(tmp1,sqrt(ztz))
                  bclus(evali)=theta
               end if
            end if
         end do
         
c++++++++++++++++++++++++++++++
c+++++++ b) Resampling step
c++++++++++++++++++++++++++++++

         do ii=1,ncluster

c++++++++++ check if the user has requested an interrupt
            call rchkusr()
    
            ztz=0.d0
            zty=0.d0

            ns=ccluster(ii)

            do i=1,ns
               ztz=ztz+1.d0/sigma2e(cstrt(ii,i))
               zty=zty+res(cstrt(ii,i))/sigma2e(cstrt(ii,i))
            end do
            
            ztz=1.d0/(ztz+sigmainv)
            zty=zty+sigmainv*mu   

            tmp1=ztz*zty  
            
            theta=rnorm(tmp1,sqrt(ztz))
            bclus(ii)=theta
            
            do i=1,ns
               b(cstrt(ii,i))=theta
            end do
         end do

c++++++++++++++++++++++++++++++++++         
c+++++++ Base line distribution
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(murand.eq.1)then

            tmp1=1.d0/(sigmainv*dble(ncluster))+psiinv
            tmp2=smu
            do i=1,ncluster
               tmp2=tmp2+sigmainv*bclus(i)       
            end do
            tmp3=tmp2*tmp1
            mu=rnorm(tmp3,sqrt(tmp1))     

         end if

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(sigmarand.eq.1)then
         
            sse=0.d0
         
            do i=1,ncluster
               sse=sse+(bclus(i)-mu)*(bclus(i)-mu) 
            end do
            
            sigma=1.d0/
     &           rgamma(0.5d0*(dble(ncluster)+tau1),0.5d0*(sse+tau2))
     
            sigmainv=1.d0/sigma
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
                  theta=bclus(evali)
               end if
               if(evali.gt.ncluster)then
                  theta=rnorm(mu,sqrt(sigma))  
               end if
               randsave(isave,nrec+1)=theta

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
                     theta=bclus(evali)                 
                  end if
                  if(evali.gt.ncluster)then
                     theta=rnorm(mu,sqrt(sigma))  
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
                  theta=bclus(evali)
               end if
               if(evali.gt.ncluster)then
                  theta=rnorm(mu,sqrt(sigma))  
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

c+++++++++++++ baseline mean

               thetasave(isave,1+nfixed+1)=mu

c+++++++++++++ baseline covariance

               thetasave(isave,1+nfixed+2)=sigma

c+++++++++++++ cluster information
               thetasave(isave,1+nfixed+3)=ncluster
               thetasave(isave,1+nfixed+4)=alpha

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
         

