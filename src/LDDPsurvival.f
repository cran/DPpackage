
c=======================================================================
      subroutine lddpsurvival(nrec,p,censor,typed,ymat,z,
     &                        ngrid,npred,grid,zpred,
     &                        a0b0,tau1,taus1,taus2,m0,sbeta0i,nu,
     &                        psiinv,
     &                        ncluster,ss,alpha,betaclus,sigmaclus,
     &                        mub,sb,tau2,y,
     &                        cpo,thetasave,randsave,survsave,
     &                        denspm,denspl,densph,
     &                        hazpm,hazpl,hazph,
     &                        meanfpm,meanfpl,meanfph,
     &                        survpm,survpl,survph,
     &                        mcmc,nsave,seed,
     &                        cstrt,ccluster,iflagp,betam,betawork,
     &                        prob,workmh1,workmh2,workv1,
     &                        xtx,xtx2,xty,xty2,fs,fm,worksam,workcpo,
     &                        workcpo2)
c=======================================================================
c
c     Subroutine `lddpsurvival' to run a Markov chain for a 
c     Linear Dependent Dirichlet Process prior for 
c     survival analysis.
c
c     Copyright: Alejandro Jara, 2009-2012.
c
c     Version 3.0: 
c
c     Last modification: 15-07-2012.
c
c     Changes and Bug fixes: 
c
c     Version 2.0 to Version 3.0:
c          - The subject-specific variances are returned in
c            the randsave object.
c
c     Version 1.0 to Version 2.0:
c          - Computation the cdf functions added.
c          - Now CPO are computed using the epsilon-DP approximation.
c          - Computation useful for ROC curve estimation added.
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
c        nrec        :  integer giving the number of data points. 
c        p           :  integer giving the number of predictors.
c        censor      :  integer indicating whether censored observations
c                       are present (1) or not (0).
c        typed       :  integer vector giving the type of observation,
c                       typed(nrec).
c        ymat        :  real matrix giving the limits of the intervals,
c                       ymat(nrec,2). The data must be in the positive
c                       real line.
c        z           :  real matrix giving the design matrix, z(nrec,p).
c
c-----------------------------------------------------------------------
c
c---- Prediction -------------------------------------------------------
c 
c        ngrid       :  integer giving the number of grid points where 
c                       the density estimates are evaluated.. 
c        npred       :  integer giving the number of predictions..
c        grid        :  real vector giving the grid, grid(ngrid). 
c        zpred       :  real matrix giving the design matrix of the 
c                       predictions, zpred(npred,p).
c
c-----------------------------------------------------------------------
c
c---- Prior information ------------------------------------------------
c 
c        aa0, ab0    :  real giving the hyperparameters of the prior
c                       distribution for the precision parameter,
c                       alpha ~ Gamma(aa0,ab0). If aa0<0 the precision 
c                       parameter is considered as a constant.
c        tau1        :  real giving the value of the tau1 hyperparameter
c                       in the Gamma component of the centering 
c                       distribution. 
c        taus1,taus2 :  real giving the values of hyperparameters of the
c                       Gamma prior for tau2.
c        m0          :  real vector giving the mean of the normal prior
c                       for mub, m0(p).
c        sbeta0i     :  real vector giving the inverse of the covariance
c                       matrix of the normal prior for mub, 
c                       sbeta0i(p,p).
c        nu          :  integer giving the degres of freedom parameter
c                       of the inverted-Wishart prior for sb.
c        psiinv      :  real matrix giving the scale matrix of the
c                       inverted-Wishart prior for sd, psiinv(p,p).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        ncluster    :  integer giving the number of clusters.
c        ss          :  integer vector giving the configurations.
c        alpha       :  real giving the current value of the precision
c                       parameter of the DP.
c        betaclus    :  real matrix giving the value of the regression  
c                       coefficients, betaclus(nrec+100,p).
c        sigmaclus   :  real vector giving the value of the kernel 
c                       variances, sigmaclus(nrec+100).
c        mub         :  real vector giving the mean of the normal 
c                       centering distribution, mub(p).
c        sb          :  real matrix giving the variance of the normal
c                       centering distribution, sb(p,p).
c        tau2        :  real giving the value of the hyperparameter
c                       in the gamma centering distribution.
c        y           :  real vector giving the log latent survival
c                       times, y(nrec). 
c
c-----------------------------------------------------------------------
c
c---- Output -----------------------------------------------------------
c
c        cpo         :  real giving the cpos and fsos, cpo(nrec,2). 
c        randsave    :  real matrix containing the mcmc samples for
c                       the regression coeff and variance, 
c                       randsave(nsave,nrec*(p+1)).
c        survsave    :  real matrix containing the mcmc samples for
c                       the survival curves, randsave(nsave,npred*ngrid).
c        thetasave   :  real matrix containing the mcmc samples for
c                       the parameters, thetasave(nsave,p+p*(p+1)/2+3)
c        denspm      :  real matrix giving the posterior mean of the 
c                       density, denspm(npred,ngrid).
c        denspl      :  real matrix giving the lower limit of the 
c                       HPD of the density, denspl(npred,ngrid).
c        densph      :  real matrix giving the upper limit of the  
c                       HPD of the density, densph(npred,ngrid).
c        hazpm       :  real matrix giving the posterior mean of the 
c                       hazard, hazpm(npred,ngrid).
c        hazpl       :  real matrix giving the lower limit of the 
c                       HPD of the hazard, hazpl(npred,ngrid).
c        hazph       :  real matrix giving the upper limit of the  
c                       HPD of the hazard, hazph(npred,ngrid).
c        meanfpm     :  real vector giving the posterior mean of the 
c                       mean function, meanfpm(npred).
c        meanfpl     :  real vector giving the lower limit of the 
c                       HPD of the mean function meanfpl(npred).
c        meanfph     :  real vector giving the upper limit of the  
c                       HPD of the mean function, meanfph(npred).
c        survpm      :  real matrix giving the posterior mean of the 
c                       survival, survpm(npred,ngrid).
c        survpl      :  real matrix giving the lower limit of the 
c                       HPD of the survival, survpl(npred,ngrid).
c        survph      :  real matrix giving the upper limit of the  
c                       HPD of the survival, survph(npred,ngrid).
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
c        
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        betam       :  real vector used to update the reg. coeff.,
c                       betam(p).
c        betawork    :  real vector used to update the reg. coeff.,
c                       betawork(p).
c        ccluster    :  integer vector indicating the number of
c                       subjects in each cluster, ccluster(nrec).
c        cstrt       :  integer matrix used to save the cluster
c                       structure, cstrt(nrec,nrec).
c        fs          :  real vector used to evaluate the conditional
c                       densities, fs(ngrid).
c        fm          :  real vector used to evaluate the conditional
c                       means, fs(npred).
c        iflagp      :  integer vector used to evaluate reg. coeff.,
c                       iflagp(p).
c        prob        :  real vector used to update the cluster 
c                       structure, prob(nrec+100).
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        workcpo     :  real vector used to compute cpos, 
c                       workcpo(nrec).
c        workcpo2    :  real vector used to compute cpos, 
c                       workcpo2(nrec).
c        workmh1     :  real vector used to update the reg. coeff.,
c                       workmh1(p*(p+1)/2).
c        workmh2     :  real vector used to update the reg. coeff.,
c                       workmh2(p*(p+1)/2).
c        worksam     :  real vector used to comput HPD bands,
c                       worksam(nsave).
c        workv1      :  real vector used to update the reg. coeff.,
c                       workv1(p).
c        xtx         :  real matrix used to update the reg. coeff.,
c                       xtx(p,p).
c        xtx2        :  real matrix used to update the reg. coeff.,
c                       xtx2(p,p).
c        xty         :  real vector used to update the reg. coeff.,
c                       xty(p).
c        xty2        :  real vector used to update the reg. coeff.,
c                       xty2(p).
c
c=======================================================================                  

      implicit none

c++++ data
      integer nrec,p
      integer censor
      integer typed(nrec)
      real*8 ymat(nrec,2)
      real*8 z(nrec,p)

c++++ prediction
      integer ngrid,npred
      real*8 grid(ngrid)
      real*8 zpred(npred,p)

c++++ prior
      integer nu  
      real*8 a0b0(2)
      real*8 tau1
      real*8 taus1,taus2
      real*8 m0(p)
      real*8 sbeta0i(p,p)
      real*8 psiinv(p,p)

c++++ current value
      integer ncluster
      integer ss(nrec)
      real*8 alpha
      real*8 betaclus(nrec+100,p)
      real*8 sigmaclus(nrec+100)
      real*8 mub(p)
      real*8 sb(p,p)
      real*8 tau2
      real*8 y(nrec)

c++++ output
      real*8 cpo(nrec,2)
      real*8 thetasave(nsave,p+p*(p+1)/2+3)
      real*8 randsave(nsave,nrec*(p+1))
      real*8 survsave(nsave,npred*ngrid)
      real*8 denspm(npred,ngrid)
      real*8 denspl(npred,ngrid)
      real*8 densph(npred,ngrid)
      real*8 hazpm(npred,ngrid)
      real*8 hazpl(npred,ngrid)
      real*8 hazph(npred,ngrid)
      real*8 meanfpm(npred)
      real*8 meanfpl(npred)
      real*8 meanfph(npred)
      real*8 survpm(npred,ngrid)
      real*8 survpl(npred,ngrid)
      real*8 survph(npred,ngrid)

c++++ mcmc
      integer mcmc(3),nburn,nskip,nsave,ndisplay
  
c++++ seeds
      integer seed1,seed2,seed(2)

c++++ external working space
      integer cstrt(nrec,nrec)
      integer ccluster(nrec)
      integer iflagp(p)
      real*8 betam(p)
      real*8 betawork(p)
      real*8 prob(nrec+100)
      real*8 workmh1(p*(p+1)/2)
      real*8 workmh2(p*(p+1)/2)
      real*8 workv1(p) 
      real*8 xtx(p,p)
      real*8 xtx2(p,p)
      real*8 xty(p)
      real*8 xty2(p)

      real*8 fs(ngrid) 
      real*8 fm(npred) 
 
      real*8 worksam(nsave) 

      real*8 workcpo(nrec)
      real*8 workcpo2(nrec)
      
c++++ internal working space
      integer count
      integer isample
      integer i,ii,j,k,l,ok
      integer dispcount
      integer evali
      integer iscan
      integer isave
      integer ns
      integer nscan
      integer since
      integer skipcount
      integer sprint
      real*8 cdflnorm
      real*8 cdfnorm
      real*8 dnrm
      real*8 dlnrm
      real*8 rgamma
      real*8 liminf
      real*8 limsup
      real*8 muwork
      real*8 rtnorm
      real*8 sigmawork
      real*8 tmp1,tmp2,tmp3,tmp4
      real*8 yc

      real*8 cpotest
      real*8 cpotest2
      
      logical ainf,asup
   
c++++ DP (functional parameter)
      real*8 eps,rbeta,weight
      parameter(eps=0.01)

c++++ CPU time
      real*8 sec00,sec0,sec1,sec

c++++ opening files

      open(unit=1,file='dppackage1.out',status='unknown',
     &     form='unformatted')
      open(unit=2,file='dppackage2.out',status='unknown',
     &     form='unformatted')
      open(unit=3,file='dppackage3.out',status='unknown',
     &     form='unformatted')
      open(unit=4,file='dppackage4.out',status='unknown',
     &     form='unformatted')

c++++++++++++++++++++++++++
c     initialize variables
c++++++++++++++++++++++++++

      cpotest=0.d0
      cpotest2=0.d0

      sigmawork=0.d0 
c++++ mcmc, priors and "zipped"

      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
   
      seed1=seed(1)
      seed2=seed(2)
      
c++++ set random number generator

      call setall(seed1,seed2)
      
c++++ start the MCMC algorithm

      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*nsave

      call cpu_time(sec0)
      sec00=0.0

c++++ cluster structure

      do i=1,nrec
         ccluster(ss(i))=ccluster(ss(i))+1
         cstrt(ss(i),ccluster(ss(i)))=i
      end do
 
      do iscan=1,nscan


c++++++++++++++++++++++++++++++++++
c+++++++ latent data            +++
c++++++++++++++++++++++++++++++++++

         if(censor.eq.1)then

            do i=1,nrec
               tmp1=0.d0
               do j=1,p
                  tmp1=tmp1+z(i,j)*betaclus(ss(i),j)
               end do
               tmp2=sigmaclus(ss(i))

               if(typed(i).eq.1)then
                  ainf=.true.
                  asup=.false.
                  liminf=0.d0
                  limsup=log(ymat(i,2))
               end if

               if(typed(i).eq.2)then
                  ainf=.false.
                  asup=.false.
                  liminf=log(ymat(i,1))
                  limsup=log(ymat(i,2))
               end if

               if(typed(i).eq.3)then
                  ainf=.false.
                  asup=.true.
                  liminf=log(ymat(i,1))
                  limsup=0.d0
               end if

               if(typed(i).ne.4)then
                  yc=rtnorm(tmp1,dsqrt(tmp2),liminf,limsup,ainf,asup) 
                  y(i)=yc
               end if   

            end do

         end if

c++++++++++++++++++++++++++++++++++
c+++++++ clustering structure   +++
c++++++++++++++++++++++++++++++++++

         do i=1,nrec

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

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
                  call relabellddp(i,since,nrec,p,ncluster,ccluster,
     &                        ss,cstrt,betaclus,sigmaclus,betawork)
               end if
               ccluster(ncluster)=ccluster(ncluster)-1 
               ncluster=ncluster-1
            end if 
         
            do j=1,ncluster

               tmp1=0.d0
               do k=1,p
                  tmp1=tmp1+z(i,k)*betaclus(j,k)
               end do                
               sigmawork=sigmaclus(j)

               tmp2=dnrm(y(i),tmp1,sqrt(sigmawork),0)
               prob(j)=dble(ccluster(j))*tmp2
            end do

            if(isample.eq.1)then
               sigmawork=1.d0/rgamma(0.5d0*tau1,0.5d0*tau2)
               sigmaclus(ncluster+1)=sigmawork
        
               call rmvnorm(p,mub,sb,workmh1,workv1,betawork) 
               do k=1,p
                  betaclus(ncluster+1,k)=betawork(k)
               end do               
            end if         

            sigmawork=sigmaclus(ncluster+1)

            tmp1=0.d0
            do k=1,p
               tmp1=tmp1+z(i,k)*betaclus(ncluster+1,k)
            end do               
         
            tmp2=dnrm(y(i),tmp1,sqrt(sigmawork),0)
            prob(ncluster+1)=alpha*tmp2

            call simdisc(prob,nrec+100,ncluster+1,evali)

            ss(i)=evali
            ccluster(evali)=ccluster(evali)+1
            cstrt(evali,ccluster(evali))=i
                
            if(evali.gt.ncluster)then
               ncluster=ncluster+1
            end if
      
         end do
     
c         call intpr("ss",-1,ss,nrec)
c         call intpr("ncluster",-1,ncluster,1)

c+++++++++++++++++++++++++++++++++++
c+++++++ regression coefficients +++
c+++++++++++++++++++++++++++++++++++

         do k=1,p
            do l=1,p
               xtx2(k,l)=sb(k,l)
            end do
         end do   
         call inverse(xtx2,p,iflagp)

c          xty2=matmul(xtx2,mub)
         do k=1,p
            tmp1=0.d0 
            do l=1,p
               tmp1=tmp1+xtx2(k,l)*mub(l)
            end do
           xty2(k)=tmp1
         end do   

         do i=1,ncluster

            do j=1,p
               do k=1,p
                  xtx(j,k)=xtx2(j,k)
               end do
               xty(j)=xty2(j)
            end do

            ns=ccluster(i)
            do j=1,ns

c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

               ii=cstrt(i,j)

               do k=1,p
                  do l=1,p
                     xtx(k,l)=xtx(k,l)+z(ii,k)*z(ii,l)/sigmaclus(i)
                  end do
                  xty(k)=xty(k)+z(ii,k)*y(ii)/sigmaclus(i)
               end do
            end do

            call inverse(xtx,p,iflagp)      

c            betam=matmul(xtx,xty)
            do k=1,p
               tmp1=0.d0 
               do l=1,p
                  tmp1=tmp1+xtx(k,l)*xty(l)
               end do
               betam(k)=tmp1
            end do   
      
            call rmvnorm(p,betam,xtx,workmh1,workv1,betawork)
  
c            call dblepr("betawork",-1,betawork,p)

            do j=1,p
               betaclus(i,j)=betawork(j)
            end do
         end do


c+++++++++++++++++++++++++++++++++++
c+++++++ variances               +++
c+++++++++++++++++++++++++++++++++++

         do i=1,ncluster

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            ns=ccluster(i)
            tmp2=0.d0
            do j=1,ns

c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

               ii=cstrt(i,j)
               tmp1=0.d0
               do k=1,p
                  tmp1=tmp1+z(ii,k)*betaclus(i,k)
               end do
               tmp1=y(ii)-tmp1
               tmp2=tmp2+tmp1*tmp1
            end do

            sigmaclus(i)=1.d0/rgamma(0.5d0*(tau1+dble(ns)),
     &                   0.5d0*(tau2+tmp2))
         end do

c         call dblepr("sigma",-1,sigmaclus,ncluster)

c+++++++++++++++++++++++++++++++++++
c+++++++ precision parameter     +++
c+++++++++++++++++++++++++++++++++++

         if(a0b0(1).gt.0.0)then
            call samalph(alpha,a0b0(1),a0b0(2),ncluster,nrec)

c            call dblepr("alpha",-1,alpha,1)
         end if 

c+++++++++++++++++++++++++++++++++++
c+++++++ baseline mean           +++
c+++++++++++++++++++++++++++++++++++

c         xty=matmul(sbeta0i,m0)
         do k=1,p
            tmp1=0.d0 
            do l=1,p
               tmp1=tmp1+sbeta0i(k,l)*m0(l)
            end do
           xty(k)=tmp1
         end do   

         do i=1,p
            do j=1,p
               xtx2(i,j)=sb(i,j)
            end do
         end do
         call inverse(xtx2,p,iflagp)

         do i=1,p
            do j=1,p
               xtx(i,j)=sbeta0i(i,j)+dble(ncluster)*xtx2(i,j)
            end do
         end do
         call inverse(xtx,p,iflagp)

         do ii=1,ncluster      
            do i=1,p
               betawork(i)=betaclus(ii,i)
            end do

c            xty=xty+matmul(xtx2,betawork)
            do k=1,p
               tmp1=0.d0 
               do l=1,p
                  tmp1=tmp1+xtx2(k,l)*betawork(l)
               end do
              xty(k)=xty(k)+tmp1
            end do   
         end do

c         betawork=matmul(xtx,xty)
         do k=1,p
            tmp1=0.d0 
            do l=1,p
               tmp1=tmp1+xtx(k,l)*xty(l)
            end do
            betawork(k)=tmp1
         end do   


         call rmvnorm(p,betawork,xtx,workmh1,workv1,mub)

c         call dblepr("mub",-1,mub,p)

c++++++++++++++++++++++++++++++++++++++
c+++++++ baseline covariance matrix +++
c++++++++++++++++++++++++++++++++++++++

         do i=1,p
            do j=1,p
               xtx(i,j)=0.d0
            end do
         end do
  
         do ii=1,ncluster 
            do i=1,p
               betawork(i)=betaclus(ii,i)-mub(i)
            end do
 
            do i=1,p
               do j=1,p
                  xtx(i,j)=xtx(i,j)+betawork(i)*betawork(j)
               end do
            end do
         end do

         do i=1,p
            do j=1,p
               sb(i,j)=xtx(i,j)+psiinv(i,j)
               xtx2(i,j)=0.d0
            end do
            xty(i)=0.d0
            iflagp(i)=0 
         end do

         call riwishart(p,nu+ncluster,sb,xtx,xtx2,xty,workmh1,
     &                  workmh2,iflagp)

c         call dblepr("sb",-1,sb,p*p)

c++++++++++++++++++++++++++++++++++++++
c+++++++ baseline gamma parameter   +++
c++++++++++++++++++++++++++++++++++++++

         tmp1=0.d0
         do i=1,ncluster
            tmp1=tmp1+1.d0/sigmaclus(i)
         end do 

         tau2=rgamma(0.5d0*(dble(ncluster)*tau1+taus1),
     &               0.5d0*(tmp1+taus2))   

c         call dblepr("tau2",-1,tau2,1)

c+++++++++++++++++++++++++++++++++++
c+++++++ saving samples          +++
c+++++++++++++++++++++++++++++++++++

         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1

               count=0

c+++++++++++++ normal baseline mean

               do i=1,p
                  count=count+1
                  thetasave(isave,count)=mub(i)
               end do   

c+++++++++++++ IW baseline scale

               do i=1,p
                  do j=i,p
                     count=count+1
                     thetasave(isave,count)=sb(i,j)
                  end do
               end do
 
c+++++++++++++ tau2
               count=count+1
               thetasave(isave,count)=tau2

c+++++++++++++ cluster information
               
               count=count+1
               thetasave(isave,count)=ncluster
               count=count+1
               thetasave(isave,count)=alpha               

c+++++++++++++ random effects
               count=0
               do i=1,nrec
                  do j=1,p
                     count=count+1
                     randsave(isave,count)=betaclus(ss(i),j) 
                  end do   
                  count=count+1
                  randsave(isave,count)=sigmaclus(ss(i)) 
               end do

c+++++++++++++ Partially sampling the DP (and CPO).

               do i=1,nrec
                  workcpo(i)=0.d0
               end do   

               do i=1,ncluster
                   prob(i)=real(ccluster(i))/(alpha+real(nrec))
               end do
               prob(ncluster+1)=alpha/(alpha+dble(nrec))
               call simdisc(prob,nrec+100,ncluster+1,evali)

               if(evali.le.ncluster)then
                  sigmawork=sigmaclus(evali)
                  do k=1,p
                     betawork(k)=betaclus(evali,k)
                  end do
               end if
               if(evali.eq.ncluster+1)then 
                  sigmawork=1.d0/rgamma(0.5d0*tau1,0.5d0*tau2)
                  sigmaclus(ncluster+1)=sigmawork
        
                  call rmvnorm(p,mub,sb,workmh1,workv1,betawork) 
                  do k=1,p
                     betaclus(ncluster+1,k)=betawork(k)
                  end do               
               end if

               tmp1=rbeta(1.d0,alpha+dble(nrec))
               tmp2=tmp1
               weight=(1.d0-tmp1)
    
               do i=1,npred  
                  call rchkusr()
     
                  muwork=0.d0
                  do j=1,p
                     muwork=muwork+zpred(i,j)*betawork(j)
                  end do
           
                  fm(i)=muwork*tmp1
                  do j=1,ngrid  
                     denspl(i,j)=tmp1*dlnrm(grid(j),muwork,
     &                           sqrt(sigmawork),0)
 
                     survpl(i,j)=tmp1*cdflnorm(grid(j),muwork,
     &                           sqrt(sigmawork),0,0)                
                  end do 
               end do

               do i=1,nrec

                  muwork=0.d0
                  do j=1,p
                     muwork=muwork+z(i,j)*betawork(j)
                  end do                

                  if(typed(i).eq.1)then
                     tmp3=cdfnorm(log(ymat(i,2)),muwork,
     &                             sqrt(sigmawork),1,0)
                     workcpo(i)=workcpo(i)+tmp1*tmp3
                  end if

                  if(typed(i).eq.2)then
                     tmp3=cdfnorm(log(ymat(i,2)),muwork,
     &                       sqrt(sigmawork),1,0)

                     tmp4=cdfnorm(log(ymat(i,1)),muwork,
     &                       sqrt(sigmawork),1,0)

                     workcpo(i)=workcpo(i)+tmp1*(tmp3-tmp4)
                  end if

                  if(typed(i).eq.3)then
                     tmp3=cdfnorm(log(ymat(i,1)),muwork,
     &                       sqrt(sigmawork),0,0)
                     workcpo(i)=workcpo(i)+tmp1*tmp3
                  end if

                  if(typed(i).eq.4)then
                     tmp3=dnrm(log(ymat(i,1)),muwork,
     &                       sqrt(sigmawork),0)*(1.d0/ymat(i,1))
                     workcpo(i)=workcpo(i)+tmp1*tmp3
                  end if
               end do

               do while((1.d0-tmp2).gt.eps)
                  call rchkusr()

                  tmp3=rbeta(1.d0,alpha+real(nrec))
                  tmp1=weight*tmp3
                  weight=weight*(1.d0-tmp3)

                  call simdisc(prob,nrec+100,ncluster+1,evali)

                  if(evali.le.ncluster)then
                     sigmawork=sigmaclus(evali)
                     do k=1,p
                        betawork(k)=betaclus(evali,k)
                     end do
                  end if

                  if(evali.eq.ncluster+1)then 
                     sigmawork=1.d0/rgamma(0.5d0*tau1,0.5d0*tau2)
                     sigmaclus(ncluster+1)=sigmawork
        
                     call rmvnorm(p,mub,sb,workmh1,workv1,betawork) 
                     do k=1,p
                        betaclus(ncluster+1,k)=betawork(k)
                     end do               
                  end if

                  do i=1,npred  
                     call rchkusr()

                     muwork=0.d0
                     do j=1,p
                        muwork=muwork+zpred(i,j)*betawork(j)
                     end do

                     fm(i)=fm(i)+muwork*tmp1
              
                     do j=1,ngrid  
                        denspl(i,j)=denspl(i,j)+tmp1*dlnrm(
     &                              grid(j),
     &                              muwork,sqrt(sigmawork),0)

                        survpl(i,j)=survpl(i,j)+tmp1*
     &                              cdflnorm(grid(j),muwork,
     &                              sqrt(sigmawork),0,0)                
                     end do 
                  end do

                  do i=1,nrec

                     muwork=0.d0
                     do j=1,p
                        muwork=muwork+z(i,j)*betawork(j)
                     end do                

                     if(typed(i).eq.1)then
                       tmp3=cdfnorm(log(ymat(i,2)),muwork,
     &                             sqrt(sigmawork),1,0)
                       workcpo(i)=workcpo(i)+tmp1*tmp3
                     end if

                     if(typed(i).eq.2)then
                       tmp3=cdfnorm(log(ymat(i,2)),muwork,
     &                       sqrt(sigmawork),1,0)

                       tmp4=cdfnorm(log(ymat(i,1)),muwork,
     &                       sqrt(sigmawork),1,0)

                       workcpo(i)=workcpo(i)+tmp1*(tmp3-tmp4)
                     end if

                     if(typed(i).eq.3)then
                       tmp3=cdfnorm(log(ymat(i,1)),muwork,
     &                       sqrt(sigmawork),0,0)
                       workcpo(i)=workcpo(i)+tmp1*tmp3
                     end if

                     if(typed(i).eq.4)then
                       tmp3=dnrm(log(ymat(i,1)),muwork,
     &                       sqrt(sigmawork),0)*(1.d0/ymat(i,1))
                       workcpo(i)=workcpo(i)+tmp1*tmp3
                     end if
                  end do

                  tmp2=tmp2+tmp1
               end do

               call simdisc(prob,nrec+100,ncluster+1,evali)

               if(evali.le.ncluster)then
                  sigmawork=sigmaclus(evali)
                  do k=1,p
                     betawork(k)=betaclus(evali,k)
                  end do
               end if
               if(evali.eq.ncluster+1)then 
                  sigmawork=1.0d0/rgamma(0.5d0*tau1,0.5d0*tau2)
                  sigmaclus(ncluster+1)=sigmawork
        
                  call rmvnorm(p,mub,sb,workmh1,workv1,betawork) 
                  do k=1,p
                     betaclus(ncluster+1,k)=betawork(k)
                  end do               
               end if

               tmp1=1.d0-tmp2
               do i=1,npred  
                  call rchkusr()
     
                  muwork=0.d0
                  do j=1,p
                     muwork=muwork+zpred(i,j)*betawork(j)
                  end do
           
                  fm(i)=fm(i)+muwork*tmp1
                  do j=1,ngrid  
                     denspl(i,j)=denspl(i,j)+tmp1*dlnrm(
     &                           grid(j),
     &                           muwork,sqrt(sigmawork),0)
 
                     survpl(i,j)=survpl(i,j)+tmp1*cdflnorm(grid(j),
     &                           muwork,sqrt(sigmawork),0,0)                
                  end do 
               end do

               do i=1,nrec

                  muwork=0.d0
                  do j=1,p
                     muwork=muwork+z(i,j)*betawork(j)
                  end do                

                  if(typed(i).eq.1)then
                     tmp3=cdfnorm(log(ymat(i,2)),muwork,
     &                             sqrt(sigmawork),1,0)
                     workcpo(i)=workcpo(i)+tmp1*tmp3
                  end if

                  if(typed(i).eq.2)then
                     tmp3=cdfnorm(log(ymat(i,2)),muwork,
     &                       sqrt(sigmawork),1,0)

                     tmp4=cdfnorm(log(ymat(i,1)),muwork,
     &                       sqrt(sigmawork),1,0)

                     workcpo(i)=workcpo(i)+tmp1*(tmp3-tmp4)
                  end if

                  if(typed(i).eq.3)then
                     tmp3=cdfnorm(log(ymat(i,1)),muwork,
     &                       sqrt(sigmawork),0,0)
                     workcpo(i)=workcpo(i)+tmp1*tmp3
                  end if

                  if(typed(i).eq.4)then
                     tmp3=dnrm(log(ymat(i,1)),muwork,
     &                       sqrt(sigmawork),0)*(1.d0/ymat(i,1))
                     workcpo(i)=workcpo(i)+tmp1*tmp3
                  end if
               end do

               ii=0
               do i=1,npred
                  call rchkusr()
                  meanfpm(i)=meanfpm(i)+fm(i)
                  do j=1,ngrid
                     denspm(i,j)=denspm(i,j)+denspl(i,j)
                     survpm(i,j)=survpm(i,j)+survpl(i,j)
                     hazpl(i,j)=denspl(i,j)/survpl(i,j)
                     hazpm(i,j)=hazpm(i,j)+hazpl(i,j)
                     ii=ii+1
                     survsave(isave,ii)=survpl(i,j)
                  end do
                  write(1) (denspl(i,j),j=1,ngrid)
                  write(3) (survpl(i,j),j=1,ngrid)
                  write(4) (hazpl(i,j),j=1,ngrid)
               end do 
               write(2) (fm(i),i=1,npred)

c+++++++++++++ lpml

               tmp2=0.d0
               do i=1,nrec
                  tmp3=workcpo(i)
                  cpo(i,1)=cpo(i,1)+1.d0/tmp3  
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

c      tmp1=log(dble(nsave)/cpotest2)
c      call dblepr("cpotest",-1,tmp1,1)

      do i=1,npred
         call rchkusr()
         meanfpm(i)=meanfpm(i)/dble(nsave)
         do j=1,ngrid
            denspm(i,j)=denspm(i,j)/dble(nsave)
            survpm(i,j)=survpm(i,j)/dble(nsave)
            hazpm(i,j)=hazpm(i,j)/dble(nsave)
         end do
      end do

      close(unit=1)
      close(unit=2)
      close(unit=3)
      close(unit=4)

      call hpddensreg(nsave,npred,ngrid,0.05d0,1,worksam,fs,
     &                denspl,densph)
      call hpddensregmf(nsave,npred,0.05d0,1,worksam,meanfpl,meanfph)

      call hpdslddpsurv(nsave,npred,ngrid,0.05d0,1,worksam,fs,
     &                  survpl,survph)

      call hpdhlddpsurv(nsave,npred,ngrid,0.05d0,1,worksam,fs,
     &                  hazpl,hazph)
      
      return
      end


c=======================================================================      
      subroutine hpdhlddpsurv(nsave,npred,ngrid,alpha,tint,
     &                        workv1,fs,llower,lupper)
c=======================================================================
c     Compute CI for hazard functions.
c
c     Alejandro Jara, 2009.
c=======================================================================
      implicit none 
c+++++External parameters
      integer tint
      integer nsave,npred,ngrid 
      real*8 alpha

c+++++External working
      real*8 fs(ngrid)
      real*8 workv1(nsave)

c+++++Output      
      real*8 llower(npred,ngrid)
      real*8 lupper(npred,ngrid)

c+++++Internal parameters
      integer maxnsave,maxngrid
      parameter(maxnsave=30000,maxngrid=500)
      real*8 aupp(2),alow(2)
      real*8 workm(maxnsave,maxngrid)

c+++++Internal working
      integer i,ii,j,l   

c+++++algorithm

      if(maxnsave.lt.nsave)then
         call rexit("Increase 'maxnsave' in 'hpdhlddpsurv'")
      end if   

      if(maxngrid.lt.ngrid)then
         call rexit("Increase 'maxngrid' in 'hpdhlddpsurv'")
      end if   

      open(unit=4,file='dppackage4.out',status='old',
     &     form='unformatted')

      do ii=1,npred
         do i=1,nsave 
            call rchkusr()
            do j=1,npred
               read(4) (fs(l),l=1,ngrid)
               if(ii.eq.j)then
                  do l=1,ngrid
                     workm(i,l)=fs(l) 
                  end do
               end if
            end do
          end do  
          rewind(unit=4)
          
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

      close(unit=4)      
      return
      end



c=======================================================================      
      subroutine hpdslddpsurv(nsave,npred,ngrid,alpha,tint,
     &                        workv1,fs,llower,lupper)
c=======================================================================
c     Compute CI for survival functions.
c
c     Alejandro Jara, 2009.
c=======================================================================
      implicit none 
c+++++External parameters
      integer tint
      integer nsave,npred,ngrid 
      real*8 alpha

c+++++External working
      real*8 fs(ngrid)
      real*8 workv1(nsave)

c+++++Output      
      real*8 llower(npred,ngrid)
      real*8 lupper(npred,ngrid)

c+++++Internal parameters
      integer maxnsave,maxngrid
      parameter(maxnsave=30000,maxngrid=300)
      real*8 aupp(2),alow(2)
      real*8 workm(maxnsave,maxngrid)

c+++++Internal working
      integer i,ii,j,l   

c+++++algorithm

      if(maxnsave.lt.nsave)then
         call rexit("Increase 'maxnsave' in 'hpdslddpsurv'")
      end if   

      if(maxngrid.lt.ngrid)then
         call rexit("Increase 'maxngrid' in 'hpdslddpsurv'")
      end if   

      open(unit=3,file='dppackage3.out',status='old',
     &     form='unformatted')

      do ii=1,npred
         do i=1,nsave 
            call rchkusr()
            do j=1,npred
               read(3) (fs(l),l=1,ngrid)
               if(ii.eq.j)then
                  do l=1,ngrid
                     workm(i,l)=fs(l) 
                  end do
               end if
            end do
          end do  
          rewind(unit=3)
          
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

      close(unit=3)      
      return
      end

