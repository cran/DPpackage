
c=======================================================================
      subroutine lddpcdensity(nrec,p,y,z,
     &                        ngrid,npred,grid,zpred,
     &                        a0b0,tau1,taus1,taus2,m0,sbeta0i,nu,
     &                        psiinv,
     &                        ncluster,ss,alpha,betaclus,sigmaclus,
     &                        mub,sb,tau2,
     &                        cpo,thetasave,randsave,
     &                        denspm,denspl,densph,
     &                        meanfpm,meanfpl,meanfph,
     &                        mcmc,nsave,seed,
     &                        cstrt,ccluster,iflagp,betam,betawork,
     &                        prob,workmh1,workmh2,workv1,
     &                        xtx,xtx2,xty,xty2,fs,fm,worksam,workcpo)
c=======================================================================
c
c     Subroutine `lddpcdensity' to run a Markov chain for a 
c     Linear Dependent Dirichlet Process prior for 
c     conditional density estimation.
c
c     Copyright: Alejandro Jara, 2008 - 2009
c
c     Version 1.0: 
c
c     Last modification: 14-03-2009.
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
c     Department of Statistics
c     Facultad de Ciencias Físicas y Matemáticas
c     Universidad de Concepción
c     Avenida Esteban Iturra S/N
c     Barrio Universitario
c     Concepción
c     Chile
c     Voice: +56-41-2203163  URL  : http://www2.udec.cl/~ajarav
c     Fax  : +56-41-2251529  Email: ajarav@udec.cl
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
c
c-----------------------------------------------------------------------
c
c---- Output -----------------------------------------------------------
c
c        cpo         :  real giving the cpos and fsos, cpo(nrec,2). 
c        randsave    :  real matrix containing the mcmc samples for
c                       the regression coeff, randsave(nsave,nrec*p).
c        thetasave   :  real matrix containing the mcmc samples for
c                       the parameters, thetasave(nsave,p+p*(p+1)/2+3)
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
      real*8 y(nrec)
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

c++++ output
      real*8 cpo(nrec,2)
      real*8 thetasave(nsave,p+p*(p+1)/2+3)
      real*8 randsave(nsave,p*nrec)
      real*8 denspm(npred,ngrid)
      real*8 denspl(npred,ngrid)
      real*8 densph(npred,ngrid)
      real*8 meanfpm(npred)
      real*8 meanfpl(npred)
      real*8 meanfph(npred)

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
      real*8 dnrm
      real*8 rgamma
      real*8 muwork
      real*8 sigmawork
      real*8 tmp1,tmp2,tmp3
   
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

c++++++++++++++++++++++++++
c     initialize variables
c++++++++++++++++++++++++++

c++++ mcmc, priors and "zipped"

      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
   
      seed1=seed(1)
      seed2=seed(2)
      sigmawork=0.d0

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
                        
               prob(j)=exp(log(dble(ccluster(j)))+tmp2)
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

            tmp1=0.0
            do k=1,p
               tmp1=tmp1+z(i,k)*betaclus(ncluster+1,k)
            end do               
         
            tmp2=dnrm(y(i),tmp1,sqrt(sigmawork),1)

            prob(ncluster+1)=exp(log(alpha)+tmp2)

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
         xty2=matmul(xtx2,mub)

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
            betam=matmul(xtx,xty)
      
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
            tmp2=0.0
            do j=1,ns

c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

               ii=cstrt(i,j)
               tmp1=0.0
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

         xty=matmul(sbeta0i,m0)

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
            xty=xty+matmul(xtx2,betawork)
         end do

         betawork=matmul(xtx,xty)
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
               end do

c+++++++++++++ Partially sampling the DP.

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
    
c               call dblepr("V1",-1,tmp1,1)
           
               do i=1,npred  
                  call rchkusr()
     
                  muwork=0.0
                  do j=1,p
                     muwork=muwork+zpred(i,j)*betawork(j)
                  end do
           
                  fm(i)=muwork*tmp1
                  do j=1,ngrid  
                     denspl(i,j)=tmp1*dnrm(grid(j),muwork,
     &                           sqrt(sigmawork),0)
                  end do 
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
                     end do 
                  end do

                  tmp2=tmp2+tmp1
               end do

               do i=1,npred
                  call rchkusr()
                  meanfpm(i)=meanfpm(i)+fm(i)
                  do j=1,ngrid
                     denspm(i,j)=denspm(i,j)+denspl(i,j)
                  end do
                  write(1) (denspl(i,j),j=1,ngrid)
               end do 
               write(2) (fm(i),i=1,npred)

c+++++++++++++ cpo


               do i=1,ncluster
                  prob(i)=dble(ccluster(i))/(alpha+dble(nrec))
               end do
               do i=ncluster+1,ncluster+100
                  prob(i)=alpha/(100.d0*(alpha+dble(nrec)))
               end do   
               call simdisc(prob,nrec+100,ncluster+100,evali)

               do i=ncluster+1,ncluster+100

                  sigmawork=1.0d0/rgamma(0.5d0*tau1,0.5d0*tau2)
                  sigmaclus(i)=sigmawork
        
                  call rmvnorm(p,mub,sb,workmh1,workv1,betawork) 
                  do k=1,p
                     betaclus(i,k)=betawork(k)
                  end do 
               end do
               
               do i=1,nrec
                  workcpo(i)=0.d0
               end do   
         
               do ii=1,ncluster+100
   
                  do i=1,nrec
                     muwork=0.0
                     do k=1,p
                        muwork=muwork+z(i,k)*betaclus(ii,k)
                     end do                
                     sigmawork=sigmaclus(ii)
                     tmp3=dnrm(y(i),muwork,sqrt(sigmawork),0)
                     workcpo(i)=workcpo(i)+prob(ii)*tmp3
                  end do   
               end do
                
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
         end do
      end do

      close(unit=1)
      close(unit=2)

      call hpddensreg(nsave,npred,ngrid,0.05d0,1,worksam,fs,
     &                denspl,densph)
      call hpddensregmf(nsave,npred,0.05d0,1,worksam,meanfpl,meanfph)
      
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
      real*8 betaclus(nsubject+100,q)
      real*8 betawork(q)
      real*8 sigmaclus(nsubject+100)
      real*8 sigmawork

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
      