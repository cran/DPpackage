
c=======================================================================                      
      subroutine sprasch(datastr,imiss,ngrid,nmissing,nsubject,p,y,
     &                   a0b0,b0,prec,psiinv,smu,tau1,tau2,
     &                   mcmc,nsave,
     &                   acrate,cpo,faccum,randsave,thetasave,
     &                   alpha,b,bclus,beta,mu,ncluster,sigma,
     &                   sigmainv,ss,     
     &                   betac,ccluster,fsavet,propvf,propvr,
     &                   prob,seed1,seed2,workmh1,workv1,grid)
c=======================================================================                      
c
c     Subroutine `sprasch' to run a Markov chain in the  
c     semiparametric Rasch model using a Dirichlet process prior
c     for the random effect distribution. In this routine, inference 
c     is based on  the Polya urn representation of Dirichlet process.
c     The algorithm 8 with m=1 of Neal (2000) is used to sample the 
c     configurations. 
c
c     Copyright: Alejandro Jara, 2006-2009
c
c     Version 1.0: 
c     Last modification: 09-04-2007.
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
c        dispcount   :  index. 
c        evali       :  integer indicator used in updating the state.
c        fsavet      :  integer vector used to store the number of
c                       observations to calculate the predictive
c                       cdf distribution, fsavet(ngrid).
c        grid        :  real vector giving the grid where the density
c                       estimate is evaluated, grid(ngrid) .
c        i           :  index. 
c        ii          :  index. 
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        k           :  index. 
c        l           :  index.
c        prob        :  real vector used to update the cluster 
c                       structure, prob(nsubject+2).
c        propvf      :  real matrix giving the variance of the normal
c                       proposal for the difficulty parameters,
c                       propvf(p-1,p-1).
c        propvr      :  real giving the variance of the normal
c                       proposal for the random effects, propvr. 
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
c        workmh1     :  real vector used to update the difficulties,
c                       workmh1((p-1)*p/2).
c        workv1      :  real vector used to update the difficulties,
c                       workv1(p-1).
c
c=======================================================================                  

      implicit none 

c+++++Data
      integer imiss,ngrid,nmissing,nsubject,p
      integer datastr(nmissing,2),y(nsubject,p)

c+++++Prior 
      real*8 aa0,ab0,a0b0(2),b0(p-1),prec(p-1,p-1),psiinv
      real*8 smu
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

c+++++Working space - General
      integer i,ii,j,k
      integer sprint
      real*8 eta,cdfnorm,dbin
      real*8 mean
      real*8 tmp1,tmp2

c+++++Working space - Random effects
      integer fsavet(ngrid)
      real*8 grid(ngrid),thetac,propvr

c+++++Working space - Configurations
      integer ccluster(nsubject),evali 
      integer ns
      integer since
      real*8 prob(nsubject+2)
      real*8 theta

c+++++Working space - RNG
      integer seed1,seed2
      real*8 rgamma,rnorm
      real runif

c+++++Working space - MCMC
      integer iscan,isave,nscan
      integer skipcount,dispcount

c+++++Working space - Difficulty parameters
      real*8 betac(p-1)
      real*8 propvf(p-1,p-1)
      real*8 workmh1((p-1)*p/2),workv1(p-1)
      
c+++++Working space - GLM part
      integer yij
      real*8 acrate2
      real*8 logp

c+++++CPU time
      real*8 sec00,sec0,sec1,sec

c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      aa0=a0b0(1)
      ab0=a0b0(2)
      
c++++ set random number generator
      
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

         logp=0.d0 

         call rmvnorm(p-1,beta,propvf,workmh1,workv1,betac)

         do i=1,nsubject
            do j=1,p-1
               eta=b(i)-betac(j) 
               mean=exp(eta)/(1.d0+exp(eta))
               logp=logp+dbin(dble(y(i,j+1)),1.d0,mean,1)
            end do
         
            do j=1,p-1
               eta=b(i)-beta(j) 
               mean=exp(eta)/(1.d0+exp(eta))
               logp=logp-dbin(dble(y(i,j+1)),1.d0,mean,1)
            end do
         end do

         tmp1=0.d0
         tmp2=0.d0
         
         do i=1,p-1
            do j=1,p-1
               tmp1=tmp1+(betac(i)-b0(i))* 
     &              prec(i,j)      *
     &              (betac(j)-b0(j))

               tmp2=tmp2+(beta(i)-b0(i))* 
     &              prec(i,j)      *
     &              (beta(j)-b0(j))

            end do
         end do
      
         logp=logp-0.5d0*tmp1+0.5d0*tmp2 
         
c+++++++ mh step

         if(log(dble(runif())).lt.logp)then
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
                  eta=bclus(j)

                  mean=exp(eta)/(1.d0+exp(eta))
                  tmp1=tmp1+dbin(dble(yij),1.d0,mean,1)

                  do k=1,p-1
                     yij=y(i,k+1)
                     eta=bclus(j)-beta(k)
                     mean=exp(eta)/(1.d0+exp(eta))                     

                     tmp1=tmp1+dbin(dble(yij),1.d0,mean,1)
                  end do

                  prob(j)=exp(log(dble(ccluster(j)))+
     &                        tmp1)
               end do
               
               theta=rnorm(mu,sigma)

               tmp1=0.d0
               
               yij=y(i,1)
               eta=theta
               mean=exp(eta)/(1.d0+exp(eta))
               
               tmp1=tmp1+dbin(dble(yij),1.d0,mean,1)

               do k=1,p-1
                  yij=y(i,k+1)
                  eta=theta-beta(k)
                  mean=exp(eta)/(1.d0+exp(eta))
               
                  tmp1=tmp1+dbin(dble(yij),1.d0,mean,1)
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
                   call relabelrasch(i,since,nsubject,ncluster,
     &                          ccluster,ss,bclus)                   
               end if

               ccluster(ncluster)=ccluster(ncluster)-1 
               ncluster=ncluster-1


               do j=1,ncluster
                  
                  tmp1=0.d0
                  
                  yij=y(i,1)
                  eta=bclus(j)

                  mean=exp(eta)/(1.d0+exp(eta))                     

                  tmp1=tmp1+dbin(dble(yij),1.d0,mean,1)
                  
                  do k=1,p-1
                     yij=y(i,k+1)
                     eta=bclus(j)-beta(k)

                     mean=exp(eta)/(1.d0+exp(eta))                     

                     tmp1=tmp1+dbin(dble(yij),1.d0,mean,1)

                  end do

                  prob(j)=exp(log(dble(ccluster(j)))+
     &                        tmp1)
               end do
               
               theta=b(i)

               tmp1=0.d0
               
               yij=y(i,1)
               eta=theta
               mean=exp(eta)/(1.d0+exp(eta))                     

               tmp1=tmp1+dbin(dble(yij),1.d0,mean,1)
               
               do k=1,p-1
                  yij=y(i,k+1)
                  eta=theta-beta(k)

                  mean=exp(eta)/(1.d0+exp(eta))                     

                  tmp1=tmp1+dbin(dble(yij),1.d0,mean,1)
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
         
            thetac=rnorm(bclus(ii),sqrt(propvr))
            logp=0.d0

            do i=1,nsubject
               if(ii.eq.ss(i))then
                 do j=1,p
                    if(j.eq.1)then
                       eta=bclus(ii)
                      else
                       eta=bclus(ii)-beta(j-1)
                    end if
                    yij=y(i,j)
                    mean=exp(eta)/(1.d0+exp(eta))
                    logp=logp-dbin(dble(yij),1.d0,mean,1)

                    if(j.eq.1)then
                       eta=thetac
                      else
                       eta=thetac-beta(j-1)
                    end if
                    yij=y(i,j)
                    mean=exp(eta)/(1.d0+exp(eta))
                    logp=logp+dbin(dble(yij),1.d0,mean,1)
                 end do
               end if
            end do
            
c++++++++++ prior ratio

            tmp1=(thetac-mu)* 
     &            sigmainv*
     &           (thetac-mu)

            tmp2=(bclus(ii)-mu)* 
     &            sigmainv*
     &           (bclus(ii)-mu)
            
            logp=logp-0.5d0*tmp1+0.5d0*tmp2 
            
c++++++++++ mh step

            if(log(dble(runif())).lt.logp)then
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
      
      do i=1,2
         acrate(i)=acrate(i)/dble(nscan)    
      end do   
      
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
      subroutine relabelrasch(ind,since,nsubject,ncluster,ccluster,
     &                        ss,bclus)
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
