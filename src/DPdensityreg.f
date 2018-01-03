c=======================================================================                      
      subroutine dpdenregr(nrec,nx,nvar,nmissi,nmiss,z,missp,npred,
     &                     xpred,ngrid,grid,  
     &                     a0b0,k0,nuvec,s2inv, 
     &                     s2invm2, 
     &                     psiinv2,tau,mcmc,nsave,
     &                     cpo,thetasave,denspm,denspl,densph, 
     &                     meanfpm,meanfpl,meanfph, 
     &                     alpha,m1,muclus,ncluster,
     &                     psi1,psiinv1,s1,sigmaclus,ss,ccluster,cstrt,
     &                     iflag,num,denom,fs,fm, 
     &                     muwork,prob,seed,sigmawork,sigworkinv,theta,
     &                     workm1,workm2,workm3,workmh1,workmh2,workv1,
     &                     workv2,workv3,ywork,
     &                     iflagx,workvx,workmx,worksam,
     &                     numcpo,denomcpo)
c=======================================================================                      
c     # of arguments = 63.
c
c     Subroutine `dpdenregr' to run a Markov chain in the DP mixture of  
c     normals model for conditional density estimation. 
c
c     In this routine, inference is based on the 
c     Polya urn representation of the Dirichlet process. The algorithm
c     8 of Neal (2000) is used with m=1. 
c
c     Copyright: Alejandro Jara, 2008-2010.
c
c     Version 1.0: 
c
c     Last modification: 25-06-2008.
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
c        nmissi      :  integer indicating whether the data contain
c                       (1) or not (0) missing data. 
c        nmiss       :  integer giving the number of missing data 
c                       points. 
c        missp       :  integer matrix giving the location of the
c                       missing data points, missp(nmiss,2).
c        nrec        :  integer giving the number of observations.
c        nvar        :  integer giving the number of variables.
c        nx          :  integer giving the number of predictors.
c        z           :  real matrix giving the response variables,
c                       z(nrec,nvar).
c
c---- Predictions ------------------------------------------------------
c 
c        cband       :  integer value indicating whether the 
c                       credible bands need to be computed or not.
c        ngrid       :  integer giving the size of the grid where
c                       the densities will be evaluated.
c        grid        :  real vector giving the grid of values for
c                       the response, grid(ngrid).
c        npred       :  integer giving the number of predictions.
c        tband       :  integer indicating the type of credible 
c                       band that need to be computed.
c        xpred       :  real matrix giving the value of the
c                       predictors for prediction, xpred(npred,nx).
c
c-----------------------------------------------------------------------
c
c---- Prior information ------------------------------------------------
c 
c        aa0, ab0    :  reals giving the hyperparameters of the prior
c                       distribution for the precision parameter,
c                       alpha ~ Gamma(aa0,ab0). If aa0<0 the precision 
c                       parameter is considered as a constant.
c        nu1         :  integer giving the degrees of freedom for the
c                       inverted-Wishart component of the baseline
c                       distribution.
c        nu2         :  integer giving the degrees of freedom for the
c                       inverted-Wishart prior distribution for the
c                       covariance matrix of the normal component of
c                       the baseline distribution.
c        m1rand      :  integer indicating wether m1 should be held 
c                       fix, 0, or random, 1.
c        s2inv       :  real matrix giving the precision of the normal
c                       prior distribution on the mean of the normal 
c                       component of the baseline distribution,
c                       s2inv(nvar,nvar).
c        s2invm2     :  real vector giving the the product of precision 
c                       matrix and the prior mean of the normal
c                       prior distribution on the mean of the normal 
c                       component of the baseline distribution,
c                       s2ivm2(nvar).
c        tau1, tau2  :  reals giving the hyperparameters of the prior 
c                       distribution for k0, k0 ~ Gamma(tau1/2,tau2/2).
c        psiinv2     :  real matrix giving the inverse of the scale 
c                       matrix for the inverted-Wishart prior on the
c                       variance matrix of the normal component of 
c                       the baseline distribution, psiinv2(nvar,nvar).
c
c        NOTE        :  the inverted-Wishart here is parametrized,
c                       sigma ~ Inv-Wishart(nu0,tinv^{-1}), such that 
c                       E(sigma)=(1/(nu0-q-1)) * tinv
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
c        cpo         :  real matrix giving the cpo's. 
c        thetasave   :  real matrix containing the mcmc samples for
c                       the parameters, 
c        meansave    :  real matrix containing the mcmc samples for
c                       the conditional means. 
c                       thetasave(nsave,nvar+nvar*(nvar+1)/2+3).
c        denspm      :  real matrix giving the posterior means for the
c                       conditional densities, denspm(npred,ngrid).
c        denspl      :  real matrix giving the lower limit for the
c                       95%HPD for the conditional densities, 
c                       denspl(npred,ngrid).
c        densph      :  real matrix giving the upper limit for the
c                       95%HPD for the conditional densities, 
c                       densph(npred,ngrid).
c        meanfpm     :  real vector giving the posterior means for the
c                       conditional means, meanfpm(npred).
c        meanfpl     :  real vector giving the lower limit for the
c                       95%HPD for the conditional means, 
c                       meanfpl(npred).
c        meanfph     :  real vector giving the upper limit for the
c                       95%HPD for the conditional means, 
c                       meanfph(npred).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        alpha       :  real giving the current value of the precision
c                       parameter of the Dirichlet process.
c        k0          :  real giving the precision parameter for the 
c                       normal baseline
c        m1          :  real vector giving the mean of the normal 
c                       component of the baseline distribution, m1(nvar)
c        muclus      :  real matrix giving the current value of the 
c                       means, muclus(nrec+2,nvar).
c        ncluster    :  integer giving the number of clusters in the
c                       data.
c        psi1        :  real matrix giving the scale matrix for the
c                       inverted-Wishart component of the baseline
c                       distribution, psi1(nvar,nvar).
c        psiinv1     :  real matrix giving the inverse of the scale 
c                       matrix for the inverted-Wishart component of 
c                       the baseline distribution, psiinv1(nvar,nvar).
c        sigmaclus   :  real matrix giving the current value of the
c                       variances, sigmaclus(nrec+2,nvar*(nvar+1)/2) .
c        ss          :  integer vector giving the cluster label for 
c                       each record, ss(nrec).
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        ccluster    :  integer vector indicating the number of
c                       subjects in each cluster, ccluster(nrec).
c        cstrt       :  integer matrix used to save the cluster
c                       structure, cstrt(nsave,nsave).
c        count       :  index.           
c        detlog      :  real used to save the log-determinant in an
c                       matrix inversion process.
c        denom       :  real matrix used to evaluate the conditional
c                       densities, denom(npred).
c        dispcount   :  index. 
c        dnrm        :  density of a normal distribution.
c        evali       :  integer indicator used in updating the state.
c        fs          :  real vector used to evaluate the conditional
c                       densities, fs(ngrid).
c        fm          :  real vector used to evaluate the conditional
c                       means, fs(npred).
c        i           :  index. 
c        ii          :  index. 
c        ihmssf      :  integer function to determine the position of a
c                       half-stored matrix.
c        iflag       :  integer vector used to evaluate the mvn density,
c                       iflag(nvar).
c        iflagx      :  integer vector used to evaluate the mvn density,
c                       iflagx(nx).
c        isave       :  index. 
c        iscan       :  index.
c        j           :  index. 
c        k           :  index. 
c        l           :  index. 
c        l1          :  index. 
c        l2          :  index. 
c        muwork      :  real vector used to save the mean,
c                       one observation, muwork(nvar).
c        num         :  real matrix used to evaluate the conditional 
c                       densities, num(npred,ngrid).  
c        ns          :  integer indicator used in updating the state. 
c        nscan       :  integer indicating the total number of MCMC
c                       scans.
c        nuniqs      :  integer giving the dimension of the half-stored
c                       covariance matrix.
c        nuwork      :  index.
c        prob        :  real vector used to update the cluster 
c                       structure, prob(nrec+2).
c        rgamma      :  gamma random number generator
c        s1          :  real matrix giving the covariance matrix of 
c                       the normal component of the baseline 
c                       distribution, s1(nvar,nvar).
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        seed3       :  seed for random number generation.
c        sigmawork   :  real matrix used to save the variance of
c                       one observation, sigmawork(nvar,nvar).
c        sigworkinv  :  real matrix used to save the inverse of the
c                       variance of one observation, 
c                       sigworkinv(nvar,nvar).
c        since       :  index.
c        skipcount   :  index. 
c        theta       :  real vector used to save randomnly generated
c                       mean vector, theta(nvar).
c        tmp1        :  real working variable. 
c        tmp2        :  real working variable.
c        workm1      :  real matrix used to update the cluster 
c                       structure, workm1(nvar,nvar).
c        workm2      :  real matrix used to update the cluster 
c                       structure, workm2(nvar,nvar).
c        workm3      :  real matrix used to update the cluster 
c                       structure, workm3(nvar,nvar).
c        workmx      :  real matrix used to evaluate the
c                       conditional densities, workmn(nx,nx).
c        workmh1     :  real vector used to update the cluster
c                       structure, workmh1(nvar*(nvar+1)/2).
c        workmh2     :  real vector used to update the cluster
c                       structure, workmh2(nvar*(nvar+1)/2).
c        worksam     :  real matrix used to store samples,
c                       worksam(ngrid).
c        workv1      :  real vector used to update the cluster
c                       structure, workv1(nvar).
c        workv2      :  real vector used to update the cluster
c                       structure, workv2(nvar).
c        workv3      :  real vector used to update the cluster
c                       structure, workv3(nvar).
c        workvx      :  real vector used to evaluate the
c                       conditional densities, workvx(nx).
c        ywork       :  real vector used to save the variables of,
c                       one observation, ywork(nvar).
c
c=======================================================================

      implicit none 

c+++++Data
      integer nrec,nx,nvar,nmissi,nmiss
      integer missp(nmiss,2)  
      double precision z(nrec,nvar)

c+++++Prediction
      integer cband,npred,ngrid,tband
      double precision xpred(npred,nx)   
      double precision grid(ngrid) 

c+++++Prior 
      integer nuvec(3),nu1,nu2,m1rand
      double precision aa0,ab0,a0b0(2)
      double precision psiinv2(nvar,nvar)
      double precision tau(2),tau1,tau2
      double precision s2inv(nvar,nvar),s2invm2(nvar)

c+++++MCMC parameters
      integer mcmc(5),nburn,nskip,nsave,ndisplay

c+++++Output
      double precision cpo(nrec,2)
      double precision thetasave(nsave,nvar+nvar*(nvar+1)/2+3)
      double precision denspm(npred,ngrid)
      double precision denspl(npred,ngrid)
      double precision densph(npred,ngrid)
      double precision meanfpm(npred)
      double precision meanfpl(npred)
      double precision meanfph(npred)

c+++++Current values of the parameters
      integer ncluster,ss(nrec)
      double precision alpha,k0,m1(nvar),muclus(nrec+100,nvar)
      double precision psi1(nvar,nvar),psiinv1(nvar,nvar)
      double precision sigmaclus(nrec+100,nvar*(nvar+1)/2)

c+++++Seeds
      integer seed(2),seed1,seed2

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++External working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer ccluster(nrec)
      integer cstrt(nrec,nrec)
      integer iflag(nvar)
      integer iflagx(nx)
      double precision muwork(nvar),prob(nrec+100)
      double precision s1(nvar,nvar)
      double precision sigmawork(nvar,nvar),sigworkinv(nvar,nvar)
      double precision theta(nvar)
      double precision workm1(nvar,nvar),workm2(nvar,nvar)
      double precision  workm3(nvar,nvar)
      double precision workmh1(nvar*(nvar+1)/2),workmh2(nvar*(nvar+1)/2)
      double precision workv1(nvar),workv2(nvar),workv3(nvar)
      double precision ywork(nvar)

      double precision workvx(nx) 
      double precision workmx(nx,nx) 

      double precision num(npred,ngrid)
      double precision denom(npred) 
      double precision fs(ngrid) 
      double precision fm(npred) 

      double precision worksam(nsave) 

      double precision numcpo(nrec),denomcpo(nrec)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Internal working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++General
      integer counter,evali
      integer i,ii,ihmssf,j,jj,k,l,nuniqs,nuwork,ns,ok 
      integer since,sprint
      integer isample

      double precision detlog,detlogx
      double precision muc,sigmac  
      double precision tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7
      double precision tpi
      double precision cpotest
      parameter(tpi=6.283185307179586476925286766559d0)
      
c+++++MCMC
      integer dispcount,isave,iscan,nscan,skipcount 

c+++++RNG and distributions
      double precision dnrm,rnorm,rgamma

c+++++DP (functional parameter)
      double precision eps,rbeta,weight
      parameter(eps=0.01)

c+++++CPU time
      double precision sec00,sec0,sec1,sec

c++++ Define parameters

      aa0=a0b0(1)
      ab0=a0b0(2)

      tau1=tau(1)
      tau2=tau(2)
      
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      cband=mcmc(4)
      tband=mcmc(5)

      nuniqs=nvar*(nvar+1)/2
      nu1=nuvec(1)
      nu2=nuvec(2)
      m1rand=nuvec(3)

      cpotest=0.d0
      
c++++ opening files

      open(unit=1,file='dppackage1.out',status='unknown',
     &     form='unformatted')

      open(unit=2,file='dppackage2.out',status='unknown',
     &     form='unformatted')

c++++ set random number generator

      seed1=seed(1)
      seed2=seed(2)

      call setall(seed1,seed2)

c++++ cluster structure

      do i=1,nrec
         ccluster(ss(i))=ccluster(ss(i))+1
         cstrt(ss(i),ccluster(ss(i)))=i
      end do

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      isample=1
      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*(nsave)

      call cpu_time(sec0)
      sec00=0.d0
      
      do iscan=1,nscan

c+++++++++++++++++++++++++++++++++++++++++
c+++++++ Imputing the missing data
c+++++++++++++++++++++++++++++++++++++++++

         if(nmissi.eq.1)then

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            do ii=1,nmiss

               i=missp(ii,1)
               jj=missp(ii,2)   

               do j=1,nvar
                  muwork(j)=muclus(ss(i),j)
                  do k=1,nvar
                     sigmawork(j,k)=sigmaclus(ss(i),ihmssf(j,k,nvar))
                  end do
               end do                

               call condmvn(jj,sigmawork,nvar,workm1,workm2)
               tmp1=sqrt(workm1(1,1))

               do k=1,nvar
                  tmp2=0.d0
                  do l=1,nvar
                     tmp2=tmp2+workm2(k,l)*muwork(l) 
                  end do
                  workv2(k)=tmp2
 
                  tmp2=0.d0
                  do l=1,nvar
                     tmp2=tmp2+workm2(k,l)*z(i,l) 
                  end do
                  workv3(k)=tmp1
               end do

               tmp2=muwork(jj)
               do k=2,nvar
                  tmp2=tmp2-workm1(1,k)*(workv3(k)-workv2(k))
               end do            

               z(i,jj)=rnorm(tmp2,tmp1) 
            end do
         end if 

c++++++++++++++++++++++++++++++++++         
c+++++++ DP part
c++++++++++++++++++++++++++++++++++

c++++++++++++++++++++++++++++++
c+++++++ a) configurations 
c++++++++++++++++++++++++++++++

c         call intpr("Step 1a",-1,0,1)

         do i=1,nrec

            ns=ccluster(ss(i))
            
c++++++++++ observation in cluster with more than 1 element
             
            if(ns.gt.1)then
 
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

               isample=1
            end if

c++++++++++ observation in cluster with only 1 element
             
            if(ns.eq.1)then
                
               since=ss(i)

               if(since.lt.ncluster)then

                   call relabelddr(i,since,nrec,nvar,ncluster,
     &                             ccluster,cstrt,ss,muclus,sigmaclus,
     &                             muwork,sigmawork)                   
               end if

               ccluster(ncluster)=ccluster(ncluster)-1 
               ncluster=ncluster-1

               isample=2
            end if

c++++++++++ sampling

            do j=1,nvar
               ywork(j)=z(i,j)
            end do
               
            do j=1,ncluster
               do k=1,nvar
                  muwork(k)=muclus(j,k)
                  do l=1,nvar
                     sigmawork(k,l)=sigmaclus(j,ihmssf(k,l,nvar))
                  end do
               end do                
                   
               call dmvn(nvar,ywork,muwork,sigmawork,tmp1,
     &                   workv1,workm1,workm2,workv2,iflag)

               prob(j)=dble(ccluster(j))*exp(tmp1)
            end do
            
            if(isample.eq.1)then
               
               do k=1,nvar
                  do l=1,nvar
                     workm3(k,l)=psiinv1(k,l)
                  end do
               end do

               call riwishart(nvar,nu1,workm3,workm1,workm2,workv1,
     &                        workmh1,workmh2,iflag)

               do k=1,nvar
                  do l=1,nvar
                     s1(k,l)=workm3(k,l)/k0
                     sigmaclus(ncluster+1,ihmssf(k,l,nvar))=
     &                      workm3(k,l)
                  end do
               end do

               call rmvnorm(nvar,m1,s1,workmh1,workv1,theta) 
   
               do k=1,nvar
                  muclus(ncluster+1,k)=theta(k)
               end do  
            end if

            do k=1,nvar
                  muwork(k)=muclus(ncluster+1,k)
                  do l=1,nvar
                  sigmawork(k,l)=
     &                     sigmaclus(ncluster+1,ihmssf(k,l,nvar))
               end do
            end do      

            call dmvn(nvar,ywork,muwork,sigmawork,tmp1,
     &                workv1,workm1,workm2,workv2,iflag)
                   
            prob(ncluster+1)=alpha*exp(tmp1)
               
            call simdisc(prob,nrec+100,ncluster+1,evali)
               
               
            if(evali.le.ncluster)then
               ss(i)=evali
               ccluster(evali)=ccluster(evali)+1
               cstrt(evali,ccluster(evali))=i
            end if   
               
            if(evali.gt.ncluster)then
               ncluster=ncluster+1
               ss(i)=ncluster
               ccluster(ncluster)=1
               cstrt(ncluster,ccluster(ncluster))=i
               do j=1,nvar
                  muclus(ncluster,j)=muclus(evali,j)
                  do k=j,nvar
                     sigmaclus(ncluster,ihmssf(j,k,nvar))=
     &                         sigmaclus(evali,ihmssf(j,k,nvar))
                  end do
               end do
            end if               

         end do


c++++++++++++++++++++++++++++++
c+++++++ b) Resampling step
c++++++++++++++++++++++++++++++

c         call intpr("Step 1b",-1,0,1)

         do ii=1,ncluster

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

c++++++++++ cluster's means and variances

            ns=ccluster(ii)

            do i=1,nvar
               workv1(i)=m1(i)*(dble(k0)/(dble(k0)+dble(ns)))
               workv2(i)=0.d0
            end do

            do i=1,ns
               do j=1,nvar 
                  workv2(j)=workv2(j)+z(cstrt(ii,i),j)                  
               end do
            end do

            do i=1,nvar
               workv2(i)=workv2(i)/dble(ns)
               muwork(i)=workv1(i)+workv2(i)*
     &                  (dble(ns)/(dble(k0)+dble(ns)))
            end do

            do i=1,nvar
               do j=1,nvar
                  workm1(i,j)=0.d0
                  workm2(i,j)=0.d0
               end do
            end do

            do i=1,ns
               do j=1,nvar 
                  do k=1,nvar 
                     workm1(j,k)=workm1(j,k)+
     &                     (z(cstrt(ii,i),j)-workv2(j))*                  
     &                     (z(cstrt(ii,i),k)-workv2(k))
                  end do   
               end do
            end do

            do i=1,nvar
               do j=1,nvar
                  workm1(i,j)=workm1(i,j)+psiinv1(i,j)
               end do
            end do
            
            do i=1,nvar
               do j=1,nvar
                  workm2(i,j)=workm2(i,j)+
     &                       (workv2(i)-m1(i))*                  
     &                       (workv2(j)-m1(j))
               end do
            end do
            
            do i=1,nvar
               do j=1,nvar
                  workm1(i,j)=workm1(i,j)+workm2(i,j)
     &                       *(dble(k0*ns)/dble(k0+ns))
               end do
            end do

            call riwishart(nvar,nu1+ns,workm1,workm2,workm3,workv1,
     &                     workmh1,workmh2,iflag)

            do i=1,nvar
               do j=i,nvar
                  sigmaclus(ii,ihmssf(i,j,nvar))=workm1(i,j)
               end do
            end do            

            do i=1,nvar
               do j=1,nvar
                  sigmawork(i,j)=workm1(i,j)/dble(k0+ns)
               end do
            end do            
            
            call rmvnorm(nvar,muwork,sigmawork,workmh1,workv1,theta) 
            
            do i=1,nvar
               muclus(ii,i)=theta(i)
            end do
         end do   


c++++++++++++++++++++++++++++++++++         
c+++++++ Baseline distribution
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++ scale matrix of the inverted-Wishart component

c         call intpr("Step 2a",-1,0,1)
         
         if(nu2.gt.0)then
         
            do i=1,nvar
               do j=1,nvar
                  workm1(i,j)=0.d0
               end do
            end do

            do i=1,ncluster
               do j=1,nvar 
                  do k=1,nvar 
                     sigmawork(j,k)=
     &                           sigmaclus(i,ihmssf(j,k,nvar))
                  end do   
               end do
            
               call invdet(sigmawork,nvar,sigworkinv,detlog,iflag,
     &                     workv3)

               do j=1,nvar 
                  do k=1,nvar 
                     workm1(j,k)=workm1(j,k)+sigworkinv(j,k)
                  end do   
               end do
            end do

            do i=1,nvar
               do j=1,nvar
                  workm1(i,j)=workm1(i,j)+psiinv2(i,j)
               end do
            end do

            nuwork=nu2+ncluster*nu1

            call riwishart(nvar,nuwork,workm1,workm2,workm3,
     &                     workv1,workmh1,workmh2,iflag)

            do i=1,nvar
               do j=1,nvar
                  psi1(i,j)=workm1(i,j)
                  psiinv1(i,j)=workm2(i,j)
               end do
            end do

         end if

c+++++++ mean of the normal component

         if(m1rand.eq.1)then
         
            do i=1,nvar
               workv1(i)=s2invm2(i)
               workv2(i)=0.d0
               do j=1,nvar
                  workm1(i,j)=0.d0
               end do
            end do

            do i=1,ncluster
         
               do j=1,nvar 
                  do k=1,nvar 
                     sigmawork(j,k)=
     &                           sigmaclus(i,ihmssf(j,k,nvar))
                  end do   
               end do
            
               call invdet(sigmawork,nvar,sigworkinv,detlog,iflag,
     &                     workv3)

               do j=1,nvar
                  tmp1=0.d0  
                  do k=1,nvar 
                     workm1(j,k)=workm1(j,k)+sigworkinv(j,k)
                     tmp1=tmp1+dble(k0)*sigworkinv(j,k)*muclus(i,k)
                  end do
                  workv2(j)=workv2(j)+tmp1
               end do         
            end do

            do i=1,nvar
               do j=1,nvar
                  workm1(i,j)=workm1(i,j)*k0
               end do
            end do

            do i=1,nvar
               workv1(i)=workv1(i)+workv2(i)
               do j=1,nvar
                  sigmawork(i,j)=s2inv(i,j)+workm1(i,j)
               end do
            end do
            
            call invdet(sigmawork,nvar,sigworkinv,detlog,iflag,workv3)

            do i=1,nvar
               tmp1=0.d0
               do j=1,nvar
                  tmp1=tmp1+sigworkinv(i,j)*workv1(j)    
                  sigmawork(i,j)=sigworkinv(i,j)
               end do
               muwork(i)=tmp1
            end do

            call rmvnorm(nvar,muwork,sigmawork,workmh1,workv1,theta) 
         
            do i=1,nvar
               m1(i)=theta(i)
            end do

         end if

c+++++++ k0

         if(tau1.gt.0)then

            tmp1=0.d0
            do i=1,ncluster 
               do j=1,nvar
                  ywork(j)=muclus(i,j)-m1(j)
                  do k=1,nvar
                     sigmawork(j,k)=sigmaclus(i,ihmssf(j,k,nvar))
                  end do
               end do
            
               call invdet(sigmawork,nvar,sigworkinv,detlog,iflag,
     &                     workv3)

               do j=1,nvar
                  do k=1,nvar
                     tmp1=tmp1+ywork(j)*sigworkinv(j,k)*ywork(k)
                  end do
               end do
            end do   
         
            k0=rgamma(0.5d0*(dble(ncluster)+tau1),0.5d0*(tmp1+tau2))

         end if

c++++++++++++++++++++++++++++++++++         
c+++++++ Precision parameter
c++++++++++++++++++++++++++++++++++

c         call intpr("Step 3",-1,0,1)
         
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

               counter=0

c+++++++++++++ normal baseline mean

               do i=1,nvar
                  counter=counter+1
                  thetasave(isave,counter)=m1(i)
               end do   

c+++++++++++++ k0 parameter

               counter=counter+1
               thetasave(isave,counter)=k0

c+++++++++++++ IW baseline scale

               do i=1,nvar
                  do j=i,nvar
                     counter=counter+1
                     thetasave(isave,counter)=psi1(i,j)
                  end do
               end do
 
c+++++++++++++ cluster information
               
               counter=counter+1
               thetasave(isave,counter)=ncluster
               counter=counter+1
               thetasave(isave,counter)=alpha   
               
c+++++++++++++ Partially sampling the DP and computing the cond. dist.

               do i=1,ncluster
                  prob(i)=dble(ccluster(i))/(alpha+dble(nrec))
               end do
               prob(ncluster+1)=alpha/(alpha+dble(nrec))
               call simdisc(prob,nrec+100,ncluster+1,evali)

               if(evali.le.ncluster)then
                  do i=1,nvar  
                     muwork(i)=muclus(evali,i)
                     do j=1,nvar
                        sigmawork(i,j)=sigmaclus(evali,
     &                                 ihmssf(i,j,nvar))
                     end do
                  end do  
               end if
               if(evali.eq.ncluster+1)then 
                  do k=1,nvar
                     do l=1,nvar
                        workm3(k,l)=psiinv1(k,l)
                     end do
                  end do

                  call riwishart(nvar,nu1,workm3,workm1,workm2,
     &                           workv1,workmh1,workmh2,iflag)
                  do k=1,nvar
                     do l=1,nvar
                        s1(k,l)=workm3(k,l)/dble(k0)
                        sigmawork(k,l)=workm3(k,l)
                     end do
                  end do
                  call rmvnorm(nvar,m1,s1,workmh1,workv1,theta) 
                  do k=1,nvar
                     muwork(k)=theta(k)
                  end do      
               end if

               do i=2,nvar
                  workvx(i-1)=muwork(i)
                  do j=2,nvar
                     workmx(i-1,j-1)=sigmawork(i,j)         
                  end do  
               end do   
               call inversedet(workmx,nx,iflagx,detlogx)

               call condmvn(1,sigmawork,nvar,workm1,workm2)
               sigmac=sqrt(workm1(1,1))

               tmp1=rbeta(1.d0,alpha+dble(nrec))
               tmp2=tmp1
               weight=(1.d0-tmp1)

               do i=1,npred  
                  call rchkusr()
                  muc=muwork(1)
                  do k=2,nvar
                     muc=muc-workm1(1,k)*(xpred(i,k-1)-muwork(k))
                  end do            

                  tmp4=-(dble(nx)*log(tpi))
                  tmp5=0.d0  
                  do j=1,nx
                     do k=1,nx
                        tmp5=tmp5+(xpred(i,j)-workvx(j))*
     &                       workmx(j,k)*(xpred(i,k)-workvx(k))
                     end do
                  end do
                  tmp6=(tmp4-detlogx-tmp5)/2.d0
                  tmp7=exp(tmp6)
                  
                  do j=1,ngrid  
                     num(i,j)=tmp1*dnrm(grid(j),muc,sigmac,0,0)*tmp7
                  end do 
                  denom(i)=tmp1*tmp7
                  fm(i)=tmp1*tmp7*muc
               end do

               do while((1.d0-tmp2).gt.eps)
                  call rchkusr()

                  tmp3=rbeta(1.d0,alpha+dble(nrec))
                  tmp1=weight*tmp3
                  weight=weight*(1.d0-tmp3)

                  call simdisc(prob,nrec+2,ncluster+1,evali)

                  if(evali.le.ncluster)then
                     do i=1,nvar  
                        muwork(i)=muclus(evali,i)
                        do j=1,nvar
                           sigmawork(i,j)=sigmaclus(evali,
     &                                              ihmssf(i,j,nvar))
                        end do
                     end do  
                  end if

                  if(evali.eq.ncluster+1)then 
                     do k=1,nvar
                        do l=1,nvar
                           workm3(k,l)=psiinv1(k,l)
                        end do
                     end do

                     call riwishart(nvar,nu1,workm3,workm1,workm2,
     &                              workv1,workmh1,workmh2,iflag)

                     do k=1,nvar
                        do l=1,nvar
                           s1(k,l)=workm3(k,l)/dble(k0)
                           sigmawork(k,l)=workm3(k,l)
                        end do
                     end do
                     call rmvnorm(nvar,m1,s1,workmh1,workv1,theta) 
                     do k=1,nvar
                        muwork(k)=theta(k)
                     end do      
                  end if

                  do i=2,nvar
                     workvx(i-1)=muwork(i)
                     do j=2,nvar
                        workmx(i-1,j-1)=sigmawork(i,j)         
                     end do  
                  end do   
                  call inversedet(workmx,nx,iflagx,detlogx)

                  call condmvn(1,sigmawork,nvar,workm1,workm2)
                  sigmac=sqrt(workm1(1,1))

                  do i=1,npred  
                     call rchkusr()
                     muc=muwork(1)
                     do k=2,nvar
                        muc=muc-workm1(1,k)*(xpred(i,k-1)-muwork(k))
                     end do            

                     tmp4=-(dble(nx)*log(tpi))
                     tmp5=0.d0  
                     do j=1,nx
                        do k=1,nx
                           tmp5=tmp5+(xpred(i,j)-workvx(j))*
     &                       workmx(j,k)*(xpred(i,k)-workvx(k))
                        end do
                     end do
                     tmp6=(tmp4-detlogx-tmp5)/2.d0
                     tmp7=exp(tmp6)

                     do j=1,ngrid  
                        num(i,j)=num(i,j)+
     &                           tmp1*dnrm(grid(j),muc,sigmac,0,0)*tmp7
                     end do 
                     denom(i)=denom(i)+tmp1*tmp7
                     fm(i)=fm(i)+tmp1*tmp7*muc
                 end do

                  tmp2=tmp2+tmp1
               end do

               do i=1,npred
                  call rchkusr()

                  fm(i)=fm(i)/denom(i)
                  meanfpm(i)=meanfpm(i)+fm(i) 
                  do j=1,ngrid
                     fs(j)=num(i,j)/denom(i)
                     denspm(i,j)=denspm(i,j)+fs(j)
                  end do
                  write(1) (fs(j),j=1,ngrid)
               end do 
               write(2) (fm(j),j=1,npred)
 
c+++++++++++++ cpo

               do i=1,ncluster
                  prob(i)=dble(ccluster(i))/(alpha+dble(nrec))
               end do
               do i=ncluster+1,ncluster+100
                  prob(i)=alpha/(100.d0*(alpha+dble(nrec)))
               end do   
               call simdisc(prob,nrec+100,ncluster+100,evali)

               do i=ncluster+1,ncluster+100
                  do k=1,nvar
                     do l=1,nvar
                        workm3(k,l)=psiinv1(k,l)
                     end do
                  end do

                  call riwishart(nvar,nu1,workm3,workm1,workm2,
     &                           workv1,workmh1,workmh2,iflag)
                  do k=1,nvar
                     do l=1,nvar
                        s1(k,l)=workm3(k,l)/dble(k0)
                        sigmaclus(i,ihmssf(k,l,nvar))=workm3(k,l)
                     end do
                  end do
                  call rmvnorm(nvar,m1,s1,workmh1,workv1,theta) 
                  do k=1,nvar
                     muclus(i,k)=theta(k)
                  end do      
               end do

               do i=1,nrec
                  numcpo(i)=0.d0
                  denomcpo(i)=0.d0
               end do   

               do ii=1,ncluster+100
     
                  do j=1,nvar
                     muwork(j)=muclus(ii,j)
                     do k=1,nvar
                        sigmawork(j,k)=
     &                         sigmaclus(ii,ihmssf(j,k,nvar))
                     end do
                  end do                

                  do j=2,nvar
                     workvx(j-1)=muwork(j)
                     do k=2,nvar
                        workmx(j-1,k-1)=sigmawork(j,k)         
                     end do  
                  end do   
                  call inversedet(workmx,nx,iflagx,detlogx)

                  call condmvn(1,sigmawork,nvar,workm1,workm2)
                  sigmac=sqrt(workm1(1,1))

                  do i=1,nrec    
                     call rchkusr()
                     muc=muwork(1)
                     do k=2,nvar
                        muc=muc-workm1(1,k)*(z(i,k)-muwork(k))
                     end do            

                     tmp4=-(dble(nx)*log(tpi))
                     tmp5=0.d0  
                     do j=1,nx
                        do k=1,nx
                           tmp5=tmp5+(z(i,j+1)-workvx(j))*
     &                       workmx(j,k)*(z(i,k+1)-workvx(k))
                        end do
                     end do
                     tmp6=(tmp4-detlogx-tmp5)/2.d0
                     tmp7=exp(tmp6)

                     numcpo(i)=numcpo(i)+
     &                   prob(ii)*dnrm(z(i,1),muc,sigmac,0,0)*tmp7
                     denomcpo(i)=denomcpo(i)+prob(ii)*tmp7
                     
                  end do
               end do
                
               tmp2=0.d0
               do i=1,nrec
                  tmp3=numcpo(i)/denomcpo(i)
                  cpo(i,1)=cpo(i,1)+1.0d0/tmp3  
                  cpo(i,2)=cpo(i,2)+tmp3                   

                  tmp2=tmp2+log(dble(isave)/cpo(i,1))
               end do

c               call dblepr("LPML",-1,tmp2,1)

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
         call rchkusr()
         cpo(i,1)=dble(nsave)/cpo(i,1)
         cpo(i,2)=cpo(i,2)/dble(nsave)
      end do
 
      do i=1,npred
         call rchkusr()

         meanfpm(i)=meanfpm(i)/dble(nsave)
         do j=1,ngrid
            denspm(i,j)= denspm(i,j)/dble(nsave)
         end do
      end do

      close(unit=1)
      close(unit=2)

      if(cband.eq.1)then

         call hpddensreg(nsave,npred,ngrid,0.05d0,tband,
     &                   worksam,fs,denspl,densph)

         call hpddensregmf(nsave,npred,0.05d0,tband,
     &                     worksam,meanfpl,meanfph)
      end if

      return
      end


c=======================================================================      
      subroutine hpddensregmf(nsave,npred,alpha,tint,
     &                        worksam,llower,lupper)
c=======================================================================
c     Compute CI for the conditional means.
c
c     Alejandro Jara, 2007-2008
c=======================================================================
      implicit none 

c+++++External parameters
      integer tint
      integer nsave,npred
      double precision alpha

c+++++External working
      double precision worksam(nsave)

c+++++Output      
      double precision llower(npred)
      double precision lupper(npred)

c+++++Internal working
      integer maxnsave,maxnpred
      parameter(maxnsave=30000,maxnpred=500)
      double precision aupp(2),alow(2)
      double precision workm(maxnsave,maxnpred)

c+++++Internal working
      integer i,ii,j   

c+++++algorithm

      if(maxnsave.lt.nsave)then
         call rexit("Increase 'maxnsave' in 'hpddensreg'")
      end if   

      if(maxnpred.lt.npred)then
         call rexit("Increase 'maxnpred' in 'hpddensreg'")
      end if   

      open(unit=2,file='dppackage2.out',status='old',
     &     form='unformatted')

      do i=1,nsave
         read(2) (workm(i,j),j=1,npred)
      end do

      do ii=1,npred
         do i=1,nsave 
            call rchkusr()
            worksam(i)=workm(i,ii)
         end do  
          
         call hpd(nsave,alpha,worksam,alow,aupp)
          
         if(tint.eq.1)then
c               (alow(1),aupp(1)): 100(1-alpha)% HPD interval
c 
            llower(ii)=alow(1)
            lupper(ii)=aupp(1)

           else
c              (alow(2),aupp(2)): 100(1-alpha)% Bayesian credible 
c              interval

            llower(ii)=alow(2)
            lupper(ii)=aupp(2)
         end if

      end do

      close(unit=2)      
      return
      end

c=======================================================================      
      subroutine hpddensreg(nsave,npred,ngrid,alpha,tint,
     &                      workv1,fs,llower,lupper)
c=======================================================================
c     Compute CI for the conditional densities.
c
c     Alejandro Jara, 2007-2008
c=======================================================================
      implicit none 
c+++++External parameters
      integer tint
      integer nsave,npred,ngrid 
      double precision alpha

c+++++External working
      double precision fs(ngrid)
      double precision workv1(nsave)

c+++++Output      
      double precision llower(npred,ngrid)
      double precision lupper(npred,ngrid)

c+++++Internal parameters
      integer maxnsave,maxngrid
      parameter(maxnsave=30000,maxngrid=500)
      double precision aupp(2),alow(2)
      double precision workm(maxnsave,maxngrid)

c+++++Internal working
      integer i,ii,j,l   

c+++++algorithm

      if(maxnsave.lt.nsave)then
         call rexit("Increase 'maxnsave' in 'hpddensreg'")
      end if   

      if(maxngrid.lt.ngrid)then
         call rexit("Increase 'maxngrid' in 'hpddensreg'")
      end if   

      open(unit=1,file='dppackage1.out',status='old',
     &     form='unformatted')

      do ii=1,npred
         do i=1,nsave 
            call rchkusr()
            do j=1,npred
               read(1) (fs(l),l=1,ngrid)
               if(ii.eq.j)then
                  do l=1,ngrid
                     workm(i,l)=fs(l) 
                  end do
               end if
            end do
          end do  
          rewind(unit=1)
          
          do i=1,ngrid
             do j=1,nsave
                workv1(j)=workm(j,i) 
             end do
          
             call hpd(nsave,alpha,workv1,alow,aupp)
          
             if(tint.eq.1)then
c               (alow(1),aupp(1)): 100(1-alpha)% HPD interval
c 
                llower(ii,i)=alow(1)
                lupper(ii,i)=aupp(1)

              else
c              (alow(2),aupp(2)): 100(1-alpha)% Bayesian credible 
c              interval

                llower(ii,i)=alow(2)
                lupper(ii,i)=aupp(2)
             end if

          end do

      end do      

      close(unit=1)      
      return
      end


c=======================================================================      
      subroutine relabelddr(ind,since,nrec,nvar,ncluster,ccluster,cstrt,
     &                      ss,muclus,sigmaclus,muwork,sigmawork)
c=======================================================================
c     relabel the clusters after the elimination of one of them
c     A.J.V., 2006
      implicit none
      integer dimen,i,ihmssf,j
      integer ind,since,nrec,nvar,ncluster,ccluster(nrec)
      integer cstrt(nrec,nrec)
      integer ss(nrec)
      double precision muclus(nrec+100,nvar),
     1  sigmaclus(nrec+100,nvar*(nvar+1)/2)
      double precision muwork(nvar),sigmawork(nvar,nvar)

      integer ns,ii  
      dimen=nvar*(nvar+1)/2

      do i=1,nvar
         muwork(i)=muclus(since,i)
         do j=1,nvar
            sigmawork(i,j)=sigmaclus(since,ihmssf(i,j,nvar))
         end do
      end do
      
      do i=since+1,ncluster

         ns=ccluster(i)    
         
         do j=1,ns
c++++++++++ check if the user has requested an interrupt
            call rchkusr()
            ii=cstrt(i,j) 
            ss(ii)=i-1
         end do

         do j=1,ns
            cstrt(i-1,j)=cstrt(i,j) 
         end do

         do j=1,nvar
            muclus(i-1,j)=muclus(i,j)
         end do
         do j=1,dimen
            sigmaclus(i-1,j)=sigmaclus(i,j)
         end do
         
         ccluster(i-1)=ccluster(i)
      end do
      
      ss(ind)=ncluster
      
      do i=1,nvar
         muclus(ncluster,i)=muwork(i)
         do j=i,nvar
             sigmaclus(ncluster,ihmssf(i,j,nvar))=sigmawork(i,j)
         end do
      end do
      
      ccluster(ncluster)=1
      
      return
      end  
            
