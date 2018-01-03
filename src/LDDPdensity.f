
c=======================================================================
      subroutine lddpcdensity(nrec,p,y,z,
     &                        ngrid,npred,nroc,grid,rocgrid,zpred,
     &                        a0b0,tau1,taus1,taus2,m0,sbeta0i,nu,
     &                        psiinv,
     &                        ncluster,ss,alpha,betaclus,sigmaclus,
     &                        mub,sb,tau2,
     &                        cpo,thetasave,randsave,aucsave,
     &                        cdfpm,cdfpl,cdfph,
     &                        denspm,denspl,densph,
     &                        meanfpm,meanfpl,meanfph,
     &                        rocpm,rocpl,rocph,
     &                        mcmc,nsave,seed,
     &                        cstrt,ccluster,iflagp,betam,betawork,
     &                        prob,workmh1,workmh2,workv1,
     &                        xtx,xtx2,xty,xty2,fs,fm,worksam,workcpo,
     &                        rocquan,rocqgrid)
c=======================================================================
c
c     Subroutine `lddpcdensity' to run a Markov chain for a 
c     Linear Dependent Dirichlet Process prior for 
c     conditional density estimation.
c
c     Copyright: Alejandro Jara, 2008 - 2012.
c
c     Version 4.0: 
c
c     Last modification: 15-07-2012.
c     
c     Changes and Bug fixes: 
c
c     Version 3.0 to Version 4.0:
c          - The subject-specific variances are returned in
c            the randsave object.
c
c     Version 2.0 to Version 3.0:
c          - Computation the cdf functions added.
c          - Now CPO are computed using the epsilon-DP approximation.
c          - Computation useful for ROC curve estimation added.
c
c     Version 1.0 to Version 2.0:
c          - Save number of clusters, cluster indicators, and cluster
c            parameters in files.
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
c        y           :  real vector giving the responses, y(nrec). 
c        z           :  real matrix giving the design matrix, z(nrec,p).
c
c-----------------------------------------------------------------------
c
c---- Prediction -------------------------------------------------------
c 
c        cband       :  integer value indicating whether the 
c                       credible bands need to be computed or not.
c        ngrid       :  integer giving the number of grid points where 
c                       the density estimates are evaluated.. 
c        npred       :  integer giving the number of predictions.
c        tband       :  integer indicating the type of credible 
c                       band that need to be computed.
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
c
c-----------------------------------------------------------------------
c
c---- Output -----------------------------------------------------------
c
c        cpo         :  real giving the cpos and fsos, cpo(nrec,2). 
c        randsave    :  real matrix containing the mcmc samples for
c                       the regression coeff and variance, 
c                       randsave(nsave,nrec*(p+1)).
c        thetasave   :  real matrix containing the mcmc samples for
c                       the parameters, thetasave(nsave,p+p*(p+1)/2+3)
c        cdfpm       :  real matrix giving the posterior mean of the 
c                       cdf, cdfpm(npred,ngrid).
c        cdfpl       :  real matrix giving the lower limit of the 
c                       HPD of the cdf, cdfpl(npred,ngrid).
c        cdfph       :  real matrix giving the upper limit of the  
c                       HPD of the cdf, cdfph(npred,ngrid).
c        denspm      :  real matrix giving the posterior mean of the 
c                       density, denspm(npred,ngrid).
c        denspl      :  real matrix giving the lower limit of the 
c                       HPD of the density, denspl(npred,ngrid).
c        densph      :  real matrix giving the upper limit of the  
c                       HPD of the density, densph(npred,ngrid).
c        meanfpm     :  real vector giving the posterior mean of the 
c                       mean function, meanfpm(npred).
c        meanfpl     :  real vector giving the lower limit of the 
c                       HPD of the mean function meanfpl(npred).
c        meanfph     :  real vector giving the upper limit of the  
c                       HPD of the mean function, meanfph(npred).
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
      double precision y(nrec)
      double precision z(nrec,p)

c++++ prediction
      integer cband,ngrid,npred,tband,rocc,nroc
      double precision grid(ngrid)
      double precision rocgrid(ngrid)
      double precision zpred(npred,p)

c++++ prior
      integer nu  
      double precision a0b0(2)
      double precision tau1
      double precision taus1,taus2
      double precision m0(p)
      double precision sbeta0i(p,p)
      double precision psiinv(p,p)

c++++ current value
      integer ncluster
      integer ss(nrec)
      double precision alpha
      double precision betaclus(nrec+100,p)
      double precision sigmaclus(nrec+100)
      double precision mub(p)
      double precision sb(p,p)
      double precision tau2

c++++ output
      integer nsave
      double precision cpo(nrec,2)
      double precision thetasave(nsave,p+p*(p+1)/2+3)
      double precision randsave(nsave,nrec*(p+1))
      double precision aucsave(nsave,npred)
      double precision cdfpm(npred,ngrid)
      double precision cdfpl(npred,ngrid)
      double precision cdfph(npred,ngrid)
      double precision denspm(npred,ngrid)
      double precision denspl(npred,ngrid)
      double precision densph(npred,ngrid)
      double precision meanfpm(npred)
      double precision meanfpl(npred)
      double precision meanfph(npred)

      double precision rocpm(npred,nroc)
      double precision rocpl(npred,nroc)
      double precision rocph(npred,nroc)

c++++ mcmc
      integer mcmc(6),nburn,nskip,ndisplay

c++++ roc computation
      double precision rocquan(nroc)
      double precision rocqgrid(npred,nroc)
  
c++++ seeds
      integer seed1,seed2,seed(2)

c++++ external working space
      integer cstrt(nrec,nrec)
      integer ccluster(nrec)
      integer iflagp(p)
      double precision betam(p)
      double precision betawork(p)
      double precision prob(nrec+100)
      double precision workmh1(p*(p+1)/2)
      double precision workmh2(p*(p+1)/2)
      double precision workv1(p) 
      double precision xtx(p,p)
      double precision xtx2(p,p)
      double precision xty(p)
      double precision xty2(p)

      double precision fs(ngrid) 
      double precision fm(npred) 
 
      double precision worksam(nsave) 

      double precision workcpo(nrec)
      
c++++ internal working space
      integer count
      integer isample
      integer i,ii,j,k,l,ok
      integer j1,j2,j3
      integer dispcount
      integer evali
      integer iscan
      integer isave
      integer ns
      integer nscan
      integer since
      integer skipcount
      integer sprint
      double precision cdfnorm
      double precision dnrm
      double precision rgamma
      double precision muwork
      double precision sigmawork
      double precision tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7
      double precision rocacum
   
c++++ DP (functional parameter)
      double precision eps,rbeta,weight
      parameter(eps=0.01)

c++++ CPU time
      double precision sec00,sec0,sec1,sec

c++++++++++++++++++++++++++
c     initialize variables
c++++++++++++++++++++++++++

c++++ mcmc, priors and "zipped"

      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      cband=mcmc(4)
      tband=mcmc(5)
      rocc=mcmc(6)
   
      seed1=seed(1)
      seed2=seed(2)
      sigmawork=0.d0


c++++ set random number generator

      call setall(seed1,seed2)

c++++ opening files

      open(unit=1,file='dppackage1.out',status='unknown',
     &     form='unformatted')

      open(unit=2,file='dppackage2.out',status='unknown',
     &     form='unformatted')

      open(unit=3,file='dppackage3.out',status='unknown',
     &     form='unformatted')

      open(unit=4,file='dppackage4.out',status='unknown',
     &     form='unformatted')

      open(unit=5,file='dppackage5.out',status='unknown',
     &     form='unformatted')

      open(unit=6,file='dppackage6.out',status='unknown',
     &     form='unformatted')

      open(unit=7,file='dppackage7.out',status='unknown',
     &     form='unformatted')


      if(rocc.eq.1)then
         open(unit=8,file='dppackage8.out',status='unknown',
     &        form='unformatted')
       else if(rocc.eq.2)then
         open(unit=8,file='dppackage8.out',status='old',
     &        form='unformatted')
       end if

      open(unit=9,file='dppackage9.out',status='unknown',
     &     form='unformatted')

      
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
 
      write(3) nrec,p

      do iscan=1,nscan

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

               tmp1=0.0
               do k=1,p
                  tmp1=tmp1+z(i,k)*betaclus(j,k)
               end do                
               sigmawork=sigmaclus(j)

               tmp2=dnrm(y(i),tmp1,sqrt(sigmawork),1)
                        
               prob(j)=dble(ccluster(j))*exp(tmp2)
            end do

            if(isample.eq.1)then
               sigmawork=1.d0/rgamma(0.5*tau1,0.5*tau2)
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
         
            tmp2=dnrm(y(i),tmp1,sqrt(sigmawork),1)

            prob(ncluster+1)=alpha*exp(tmp2)

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

         do i=1,p
            do j=1,p
               xtx2(i,j)=sb(i,j)
            end do
         end do   
         call inverse(xtx2,p,iflagp)

         do i=1,p
            tmp1=0.d0
            do j=1,p
               tmp1=tmp1+xtx2(i,j)*mub(j)
            end do
            xty2(i)=tmp1
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

         do i=1,p
            tmp1=0.d0
            do j=1,p
               tmp1=tmp1+sbeta0i(i,j)*m0(j)
            end do
            xty(i)=tmp1
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

            do i=1,p
               tmp1=0.d0
               do j=1,p 
                  tmp1=tmp1+xtx2(i,j)*betawork(j)
               end do
               xty(i)=xty(i)+tmp1
            end do
         end do

         do i=1,p
            tmp1=0.d0
            do j=1,p 
               tmp1=tmp1+xtx(i,j)*xty(j)
            end do
            betawork(i)=tmp1
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

c+++++++++++++ Partially sampling the DP and CPO computation.

               if(rocc.eq.2)then
                 do i=1,npred
                     read(8) (rocquan(j),j=1,nroc)
                    do j=1,nroc
                        rocqgrid(i,j)=rocquan(j)
                      end do
                 end do
               end if

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
                  sigmawork=1.0d0/rgamma(0.5d0*tau1,0.5d0*tau2)
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
                     denspl(i,j)=tmp1*dnrm(grid(j),muwork,
     &                           sqrt(sigmawork),0)

                     cdfpl(i,j)=tmp1*cdfnorm(grid(j),muwork,
     &                          sqrt(sigmawork),1,0)
                  end do 

                  if(rocc.eq.2)then
                     do j=1,nroc 
                        rocpl(i,j)=tmp1*cdfnorm(rocqgrid(i,j),
     &                             muwork,sqrt(sigmawork),1,0)
                     end do 
                  end if

               end do

               do i=1,nrec
                  muwork=0.d0
                  do j=1,p
                     muwork=muwork+z(i,j)*betawork(j)
                  end do           

                  tmp3=dnrm(y(i),muwork,sqrt(sigmawork),0)
                  workcpo(i)=tmp1*tmp3
               end do

               do while((1.0-tmp2).gt.eps)
                  call rchkusr()

                  tmp3=rbeta(1.d0,alpha+real(nrec))
                  tmp1=weight*tmp3
                  weight=weight*(1.0-tmp3)

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

                  do i=1,npred  
                     call rchkusr()

                     muwork=0.0
                     do j=1,p
                        muwork=muwork+zpred(i,j)*betawork(j)
                     end do

                     fm(i)=fm(i)+muwork*tmp1
              
                     do j=1,ngrid  
                        denspl(i,j)=denspl(i,j)+tmp1*dnrm(grid(j),
     &                              muwork,sqrt(sigmawork),0)

                        cdfpl(i,j)=cdfpl(i,j)+tmp1*cdfnorm(grid(j),
     &                             muwork,sqrt(sigmawork),1,0)
                     end do 

                     if(rocc.eq.2)then
                        do j=1,nroc 
                           rocpl(i,j)=rocpl(i,j)+
     &                             tmp1*cdfnorm(rocqgrid(i,j),
     &                             muwork,sqrt(sigmawork),1,0)
                        end do 
                     end if
                  end do


                  do i=1,nrec
                     muwork=0.d0
                     do j=1,p
                        muwork=muwork+z(i,j)*betawork(j)
                     end do           

                     tmp3=dnrm(y(i),muwork,sqrt(sigmawork),0)
                     workcpo(i)=workcpo(i)+tmp1*tmp3
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

               tmp1=weight
           
               do i=1,npred  
                  call rchkusr()
     
                  muwork=0.0
                  do j=1,p
                     muwork=muwork+zpred(i,j)*betawork(j)
                  end do
           
                  fm(i)=fm(i)+muwork*tmp1
                  do j=1,ngrid  
                     denspl(i,j)=denspl(i,j)+tmp1*dnrm(grid(j),
     &                           muwork,sqrt(sigmawork),0)

                     cdfpl(i,j)=cdfpl(i,j)+tmp1*cdfnorm(grid(j),
     &                          muwork,sqrt(sigmawork),1,0)

                  end do 

                  if(rocc.eq.2)then
                     do j=1,nroc 
                        rocpl(i,j)=rocpl(i,j)+
     &                             tmp1*cdfnorm(rocqgrid(i,j),
     &                             muwork,sqrt(sigmawork),1,0)
                     end do 
                  end if
               end do

               do i=1,nrec
                  muwork=0.d0
                  do j=1,p
                     muwork=muwork+z(i,j)*betawork(j)
                  end do           

                  tmp3=dnrm(y(i),muwork,sqrt(sigmawork),0)
                  workcpo(i)=workcpo(i)+tmp1*tmp3
               end do

               do i=1,npred
                  call rchkusr()
                  meanfpm(i)=meanfpm(i)+fm(i)
                  do j=1,ngrid
                     denspm(i,j)=denspm(i,j)+denspl(i,j)
                     cdfpm(i,j)=cdfpm(i,j)+cdfpl(i,j)
                  end do
                  write(1) (denspl(i,j),j=1,ngrid)
                  write(3) (cdfpl(i,j),j=1,ngrid)

                  if(rocc.eq.2)then

                     rocacum=0.0
                     do j=1,nroc 
                        rocpl(i,j)=1.0-rocpl(i,j)
                        rocacum=rocacum+rocpl(i,j) 
                        rocpm(i,j)=rocpm(i,j)+rocpl(i,j)
                     end do 
                     aucsave(isave,i)=rocacum/dble(nroc)
                  end if
                  write(9) (rocpl(i,j),j=1,nroc)
               end do 
               write(2) (fm(i),i=1,npred)


               if(rocc.eq.1)then
                  do i=1,npred
                     call rchkusr()
                     do j=1,nroc
                        tmp1=1.0-rocgrid(j)
                        tmp2=cdfpl(i,1)
                        tmp3=cdfpl(i,ngrid)

                        if(tmp1.le.tmp2)then

                           tmp5=cdfpl(i,1)  
                           tmp6=cdfpl(i,2) 
                           ii=2
                           
                           do while((tmp5-tmp6).eq.0.0)
                              ii=ii+1
                              tmp6=cdfpl(i,ii) 
                           end do

                           tmp7=grid(ii)+(grid(1)-grid(ii))*
     &                          (tmp1-tmp6)/(tmp5-tmp6)

                           if(tmp7.gt.grid(1))then
                              call dblepr("tmp1",-1,tmp1,1)
                              call dblepr("tmp2",-1,tmp2,1)
                              call dblepr("tmp3",-1,tmp3,1)
                              call dblepr("tmp5",-1,tmp5,1)
                              call dblepr("tmp6",-1,tmp6,1)
                              call dblepr("tmp7",-1,tmp7,1)
                              call dblepr("grid1",-1,grid(1),1)

                              do ii=1,ngrid
                                 call dblepr("cdf",-1,cdfpl(i,ii),1)
                              end do

                              call rexit("Quantile outside range 1")
                           end if

                         else if(tmp1.gt.tmp3)then

                           tmp5=cdfpl(i,ngrid-1)  
                           tmp6=cdfpl(i,ngrid) 
 
                           tmp7=grid(ngrid-1)+
     &                              (grid(ngrid)-grid(ngrid-1))*
     &                              (tmp1-tmp5)/(tmp6-tmp5)

                           if(tmp7.lt.grid(ngrid))then
                              call rexit("Quantile outside range 2")
                           end if

                         else
                           j1=1
                           j2=ngrid
                           ok=0
                           do while(ok.eq.0)
                              j3=(j1+j2)/2
                              tmp4=cdfpl(i,j3)

                              if(tmp1.le.tmp4)then
                                 j2=j3
                                else
                                 j1=j3
                              end if
                              if((j1+1).ge.j2)ok=1
                           end do
                           
                          tmp5=cdfpl(i,j1)  
                          tmp6=cdfpl(i,j1+1)  
                          tmp7=grid(j1)+(grid(j1+1)-grid(j1))*
     &                                  (tmp1-tmp5)/(tmp6-tmp5)

                          if(tmp7.lt.grid(j1))then
                             call rexit("Quantile outside range 3")
                          end if
                          if(tmp7.gt.grid(j1+1))then
                             call dblepr("tmp1",-1,tmp1,1)  
                             call dblepr("tmp5",-1,tmp5,1)  
                             call dblepr("tmp6",-1,tmp6,1)  
                             call dblepr("tmp7",-1,tmp7,1)  
                             call dblepr("grid1",-1,grid(j1),1)  
                             call dblepr("grid2",-1,grid(j1+1),1)  
                             call rexit("Quantile outside range 4")
                          end if

                        end if
                        rocquan(j)=tmp7

c                        call dblepr("quan",-1,tmp7,1)
                     end do

                     write(8) (rocquan(j),j=1,nroc)
                  end do
               end if

c+++++++++++++ save elements in files

c               write(4) ncluster
c               write(5) (ss(i),i=1,nrec)
c               do i=1,ncluster 
c                  write(6) (betaclus(i,j),j=1,p)
c               end do
c               write(7) (sigmaclus(i),i=1,ncluster)

c+++++++++++++ lpml

               tmp2=0.d0
               do i=1,nrec
                  tmp3=workcpo(i)
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
            denspm(i,j)=denspm(i,j)/dble(nsave)
            cdfpm(i,j)=cdfpm(i,j)/dble(nsave)
         end do

         if(rocc.eq.2)then
            do j=1,nroc 
               rocpm(i,j)=rocpm(i,j)/dble(nsave)
            end do 
         end if

      end do

      close(unit=1)
      close(unit=2)
      close(unit=3)
      close(unit=4)
      close(unit=5)
      close(unit=6)
      close(unit=7)
      close(unit=8)
      close(unit=9)

      if(cband.eq.1)then

         call hpddensreg(nsave,npred,ngrid,0.05d0,tband,worksam,fs,
     &                   denspl,densph)

         call hpddensregmf(nsave,npred,0.05d0,tband,worksam,
     &                     meanfpl,meanfph)

         call hpdlddpcdf(nsave,npred,ngrid,0.05d0,tband,worksam,fs,
     &                   cdfpl,cdfph)
      
      end if

      if(rocc.eq.2)then
         call hpdlddproc(nsave,npred,nroc,0.05d0,tband,worksam,
     &                   rocquan,rocpl,rocph)
      end if


      return
      end

c=======================================================================      
      subroutine hpdlddproc(nsave,npred,nroc,alpha,tint,
     &                      workv1,fs,llower,lupper)
c=======================================================================
c     Compute CI for survival functions.
c
c     Alejandro Jara, 2009.
c=======================================================================
      implicit none 
c+++++External parameters
      integer tint
      integer nsave,npred,nroc
      double precision alpha

c+++++External working
      double precision fs(nroc)
      double precision workv1(nsave)

c+++++Output      
      double precision llower(npred,nroc)
      double precision lupper(npred,nroc)

c+++++Internal parameters
      integer maxnsave,maxngrid
      parameter(maxnsave=30000,maxngrid=300)
      double precision aupp(2),alow(2)
      double precision workm(maxnsave,maxngrid)

c+++++Internal working
      integer i,ii,j,l   

c+++++algorithm

      if(maxnsave.lt.nsave)then
         call rexit("Increase 'maxnsave' in 'hpdslddproc'")
      end if   

      if(maxngrid.lt.nroc)then
         call rexit("Increase 'maxngrid' in 'hpdslddproc'")
      end if   

      open(unit=9,file='dppackage9.out',status='old',
     &     form='unformatted')

      do ii=1,npred
         do i=1,nsave 
            call rchkusr()
            do j=1,npred
               read(9) (fs(l),l=1,nroc)
               if(ii.eq.j)then
                  do l=1,nroc
                     workm(i,l)=fs(l) 
                  end do
               end if
            end do
          end do  
          rewind(unit=9)
          
          do i=1,nroc
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

      close(unit=9)      
      return
      end



c=======================================================================      
      subroutine hpdlddpcdf(nsave,npred,ngrid,alpha,tint,
     &                       workv1,fs,llower,lupper)
c=======================================================================
c     Compute CI for survival functions.
c
c     Alejandro Jara, 2009.
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
      parameter(maxnsave=30000,maxngrid=300)
      double precision aupp(2),alow(2)
      double precision workm(maxnsave,maxngrid)

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


c=======================================================================
      subroutine readlddpdens(nsave,nrec,p,
     &                        ss,betaw,sigmaw,
     &                        nclus,clusind,clusreg,
     &                        clussig)
c=======================================================================
      implicit none

c++++ input
      integer nsave,nrec,p

c++++ working
      integer ss(nrec)
      double precision betaw(p)
      double precision sigmaw(nrec)

c++++ output
      integer nclus(nsave)
      integer clusind(nsave,nrec)
      double precision clusreg(nsave,nrec*p)
      double precision clussig(nsave,nrec)

c++++ internal working
      integer count,i,j,k,ncluster

      open(unit=4,file='dppackage4.out',status='unknown',
     &     form='unformatted')

      open(unit=5,file='dppackage5.out',status='unknown',
     &     form='unformatted')

      open(unit=6,file='dppackage6.out',status='unknown',
     &     form='unformatted')

      open(unit=7,file='dppackage7.out',status='unknown',
     &     form='unformatted')


c++++ algorithm

      do i=1,nsave
c+++++++ check if the user has requested an interrupt
         call rchkusr()

         read(4) ncluster
         nclus(i)=ncluster

         read(5) (ss(j),j=1,nrec)
         do j=1,nrec
            clusind(i,j)=ss(j)
         end do  

         count=0
         do j=1,ncluster
            read(6) (betaw(k),k=1,p)
            do k=1,p
               count=count+1
               clusreg(i,count)=betaw(k) 
            end do
         end do

         read(7) (sigmaw(j),j=1,ncluster)
         do j=1,ncluster
            clussig(i,j)=sigmaw(j)
         end do
      end do 
     
      close(unit=4)
      close(unit=5)
      close(unit=6)
      close(unit=7)

      return
      end

c=======================================================================
      subroutine relabellddp(ind,since,nsubject,q,ncluster,ccluster,
     &                       ss,cstrt,
     &                       betaclus,sigmaclus,betawork)
c=======================================================================
      implicit none
      integer i,j,ind,since,nsubject,q,ncluster,ccluster(nsubject)
      integer ss(nsubject),cstrt(nsubject,nsubject)
      double precision betaclus(nsubject+100,q)
      double precision betawork(q)
      double precision sigmaclus(nsubject+100)
      double precision sigmawork

      integer ns,ii
      
      do i=1,q
         betawork(i)=betaclus(since,i)
      end do
      sigmawork=sigmaclus(since)

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

         do j=1,q
            betaclus(i-1,j)=betaclus(i,j)
         end do
         sigmaclus(i-1)=sigmaclus(i)
         ccluster(i-1)=ccluster(i)
      end do

      ss(ind)=ncluster
      
      do i=1,q
         betaclus(ncluster,i)=betawork(i)
      end do
      sigmaclus(ncluster)=sigmawork
  
      ccluster(ncluster)=1
      
      return
      end
      
