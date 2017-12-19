c=======================================================================                      
      subroutine dpmdenscens(nrec,nvar,tint,llower,lupper,
     &                       ngrid,grid,ngridb,
     &                       a0b0,m1rand,nuvec,psiinv2,tau,s2inv,
     &                       s2invm2,                                            
     &                       mcmc,nsave,                                       
     &                       alpha,k0,m1,muclus,ncluster,psi1,psiinv1,
     &                       ss,sigmaclus,y,
     &                       randsave,thetasave,fbiv,funi,                             
     &                       seed,
     &                       muwork,sigmawork,workm1,workm2,workm3,
     &                       workv1,workv2,
     &                       ccluster,cstrt,iflag,
     &                       prob,workmh1,workmh2,ywork)              

c=======================================================================                      
c
c     # of arguments = 46.
c
c     Subroutine `dpmdensces' to generate a Markov chain in a 
c     DPM of normals model for multivariate interval-censored data.
c     In this routine, inference is based on the Polya urn representation 
c     of the Dirichlet process on the augmented posterior. 
c     The algorithm 8 of Neal (2000) is used, with m=1, for the DP part
c     and the latent variables are updated using Geweke (1991)'s 
c     algorithm. 
c
c     Copyright: Alejandro Jara, 2010.
c
c     Version 1.0: 
c
c     Last modification: 19-07-2010.
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
c        nvar        :  integer giving the number of variables.
c        tint        :  integer matrix giving the type of response:
c                       (1) left, (2) interval, (3) right censored,
c                       (4) missing and (5) known.
c        llower      :  real matrix giving the lower limits of 
c                       intervals, llower(nrec,nvar).
c        llupper     :  real matrix giving the upper limits of 
c                       intervals, lupper(nrec,nvar).
c
c-----------------------------------------------------------------------
c
c---- Grid -------------------------------------------------------------
c 
c        ngrid       :  integer giving the number of grid 
c                       points where the marginal density is evaluated.
c        grid        :  real matrix giving the grid points for each
c                       coordinate, grid(ngrid,nvar).
c        ngridb      :  integer giving the lentgh of the vectorization
c                       of the bivariate grids. ngridb=nvar*ngrid*ngrid
c                       if nvar > 1 and 1 otherwise.
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
c---- Current value of the parameters ----------------------------------
c
c        alpha       :  real giving the current value of the precision
c                       parameter of the Dirichlet process.
c        k0          :  real giving the precision parameter for the 
c                       normal baseline
c        m1          :  real vector giving the mean of the normal 
c                       component of the baseline distribution, m1(nvar)
c        muclus      :  real matrix giving the current value of the 
c                       means, muclus(nrec+100,nvar).
c        ncluster    :  integer giving the number of clusters in the
c                       data.
c        psi1        :  real matrix giving the scale matrix for the
c                       inverted-Wishart component of the baseline
c                       distribution, psi1(nvar,nvar).
c        psiinv1     :  real matrix giving the inverse of the scale 
c                       matrix for the inverted-Wishart component of 
c                       the baseline distribution, psiinv1(nvar,nvar).
c        sigmaclus   :  real matrix giving the current value of the
c                       variances, sigmaclus(nrec+100,nvar*(nvar+1)/2) .
c        ss          :  integer vector giving the cluster label for 
c                       each record, ss(nrec).
c        y           :  real matrix giving the current value of the
c                       imputed latent variables, y(nvar,nrec).
c
c-----------------------------------------------------------------------
c
c---- Output -----------------------------------------------------------
c
c        randsave    :  real matrix containing the mcmc samples for
c                       the parameters of the density, 
c                       randsave(nsave,(nrec+2)*nvar+
c                       (nrec+1)*nvar*(nvar+1)/2+nvar).
c        thetasave   :  real matrix containing the mcmc samples for
c                       the parameters, 
c                       thetasave(nsave,nvar+nvar*(nvar+1)/2+3).
c        fbiv        :  real vector giving the posterior mean of the
c                       bivariate densities, fbiv(ngridb).
c        funi        :  real matrix giving the posterior mean of the
c                       marginal densities, funi(ngrid,nvar). 
c
c-----------------------------------------------------------------------
c
c---- External working space -------------------------------------------
c
c        ccluster    :  integer vector indicating the number of
c                       subjects in each cluster, ccluster(nrec).
c        cstrt       :  integer matrix used to save the cluster
c                       structure, cstrt(nrec,nrec).
c        iflag       :  integer vector used to evaluate the mvn density,
c                       iflag(nvar).
c        muwork      :  real vector used to save the mean,
c                       one observation, muwork(nvar).
c        prob        :  real vector used to update the cluster 
c                       structure, prob(nrec+100).
c        seed        :  seeds for random number generation.
c        sigmawork   :  real matrix used to save the variance of
c                       one observation, sigmawork(nvar,nvar).
c        workm1      :  real matrix used to update the cluster 
c                       structure, workm1(nvar,nvar).
c        workm2      :  real matrix used to update the cluster 
c                       structure, workm2(nvar,nvar).
c        workm3      :  real matrix used to update the cluster 
c                       structure, workm3(nvar,nvar).
c        workmh1     :  real vector used to update the cluster
c                       structure, workmh1(nvar*(nvar+1)/2).
c        workmh2     :  real vector used to update the cluster
c                       structure, workmh2(nvar*(nvar+1)/2).
c        workv1      :  real vector used to update the cluster
c                       structure, workv1(nvar).
c        workv2      :  real vector used to update the cluster
c                       structure, workv2(nvar).
c        ywork       :  real vector used to save the variables of,
c                       one observation, ywork(nvar).
c
c=======================================================================                      

      implicit none

c+++++Data
      integer nrec,nvar
      integer tint(nrec,nvar)
      double precision llower(nrec,nvar)
      double precision lupper(nrec,nvar)

c+++++Functionals
      integer ngrid,ngridb
      double precision grid(ngrid,nvar)
 
c+++++Prior 
      integer nuvec(2),nu1,nu2,m1rand
      double precision aa0,ab0,a0b0(2)
      double precision psiinv2(nvar,nvar)
      double precision tau(2),tau1,tau2
      double precision s2inv(nvar,nvar),s2invm2(nvar)

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Current values of the parameters
      integer ncluster,ss(nrec)
      double precision alpha,k0,m1(nvar)
      double precision muclus(nrec+100,nvar)
      double precision psi1(nvar,nvar),psiinv1(nvar,nvar)
      double precision sigmaclus(nrec+100,nvar*(nvar+1)/2)
      double precision y(nrec,nvar)
 
c+++++Output
      double precision randsave(nsave,
     1 (nrec+2)*nvar+(nrec+1)*nvar*(nvar+1)/2)
      double precision thetasave(nsave,nvar+nvar*(nvar+1)/2+3)
      double precision fbiv(ngridb)
      double precision funi(ngrid,nvar)

c+++++Seeds
      integer seed(2),seed1,seed2

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++External working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++Latent 
      double precision muwork(nvar)
      double precision sigmawork(nvar,nvar)
      double precision workm1(nvar,nvar)
      double precision workm2(nvar,nvar)
      double precision workm3(nvar,nvar)
      double precision workv1(nvar)
      double precision workv2(nvar)

c+++++DP
      integer ccluster(nrec)
      integer cstrt(nrec,nrec)
      integer iflag(nvar)
      double precision prob(nrec+100)
      double precision ywork(nvar)
      double precision workmh1(nvar*(nvar+1)/2)
      double precision workmh2(nvar*(nvar+1)/2)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Internal working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++General
      integer count,count2
      integer i1,j1
      integer evali,i,ii,j,jj,k,l,nuniqs,nuwork,ns,ok,since,sprint
      integer ihmssf
      double precision tmp1,tmp2,tmp3

c+++++Bivariate density
      integer iflagd(2)
      double precision gridd(2)
      double precision mud(2)
      double precision sigmad1(2,2)
      double precision tpi 

c+++++Latent variables
      double precision slow,supp  
      logical ainf,binf

c+++++MCMC
      integer dispcount,isave,iscan,nscan,skipcount 

c+++++RNG and distributions
      double precision dnrm
      double precision rgamma,rtnorm

c+++++CPU time
      double precision sec00,sec0,sec1,sec

c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)

      aa0=a0b0(1)
      ab0=a0b0(2)

      tau1=tau(1)
      tau2=tau(2)

      nuniqs=nvar*(nvar+1)/2
      nu1=nuvec(1)
      nu2=nuvec(2)

c++++ set random number generator
      seed1=seed(1)
      seed2=seed(2)
      call setall(seed1,seed2)

c++++ opening files

c++++ set configurations

      do i=1,nrec
         ccluster(ss(i))=ccluster(ss(i))+1
         cstrt(ss(i),ccluster(ss(i)))=i
      end do

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      tpi=6.283185307179586476925286766559d0

      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*nsave

      call cpu_time(sec0)
      sec00=0.d0
      
      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Latent data
c+++++++++++++++++++++++++++++++++++++++++++++++++

         do i=1,nrec
            do j=1,nvar
               muwork(j)=muclus(ss(i),j)
               do k=1,nvar
                  sigmawork(j,k)=sigmaclus(ss(i),ihmssf(j,k,nvar))  
               end do
            end do

            if(nvar.gt.1)then

               do j=1,nvar

                  if(tint(i,j).ne.5)then

                     call condmvn(j,sigmawork,nvar,workm1,workm2)
                     tmp1=sqrt(workm1(1,1))

                     do k=1,nvar
                        tmp2=0.d0
                        do l=1,nvar
                           tmp2=tmp2+workm2(k,l)*muwork(l) 
                        end do
                        workv1(k)=tmp2

                        tmp2=0.d0
                        do l=1,nvar
                           tmp2=tmp2+workm2(k,l)*y(i,l) 
                        end do
                        workv2(k)=tmp2
                     end do

                     tmp2=muwork(j)
                     do k=2,nvar
                        tmp2=tmp2-workm1(1,k)*(workv2(k)-workv1(k))
                     end do            

                     if(tint(i,j).eq.1)then
                        ainf=.true.
                        binf=.false. 
                        slow=0.d0
                        supp=lupper(i,j)

                      else if(tint(i,j).eq.2)then
                        ainf=.false.
                        binf=.false. 
                        slow=llower(i,j)
                        supp=lupper(i,j)

                      else if(tint(i,j).eq.3)then
                        ainf=.false.
                        binf=.true.
                        slow=llower(i,j)
                        supp=0.d0
                      else
                        ainf=.true.
                        binf=.true.
                        slow=0.d0
                        supp=0.d0
                     end if
                    y(i,j)=rtnorm(tmp2,tmp1,slow,supp,ainf,binf)
                  end if
               end do

             else

               if(tint(i,1).ne.5)then

                  if(tint(i,1).eq.1)then
                     ainf=.true.
                     binf=.false. 
                     slow=0.d0
                     supp=lupper(i,1)

                   else if(tint(i,1).eq.2)then
                     ainf=.false.
                     binf=.false. 
                     slow=llower(i,1)
                     supp=lupper(i,1)

                   else if(tint(i,1).eq.3)then
                     ainf=.false.
                     binf=.true.
                     slow=llower(i,1)
                     supp=0.d0
                   else
                     ainf=.true.
                     binf=.true.
                     slow=0.d0
                     supp=0.d0
                  end if
                  tmp2=muwork(1)
                  tmp1=sqrt(sigmawork(1,1))
                   
                  y(i,1)=rtnorm(tmp2,tmp1,slow,supp,ainf,binf)
               end if

            end if

         end do


c++++++++++++++++++++++++++++++++++         
c+++++++ DP part
c++++++++++++++++++++++++++++++++++

c++++++++++++++++++++++++++++++
c+++++++ a) configurations 
c++++++++++++++++++++++++++++++

         do i=1,nrec

            ns=ccluster(ss(i))
            
c++++++++++ observation in cluster with more than 1 element
             
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
               
               do j=1,nvar
                  ywork(j)=y(i,j)
               end do
               
               do j=1,ncluster
                  do k=1,nvar
                     muwork(k)=muclus(j,k)
                     do l=1,nvar
                        sigmawork(k,l)=sigmaclus(j,ihmssf(k,l,nvar))
                     end do
                  end do                
                  call dmvnd(nvar,ywork,muwork,sigmawork,tmp1,iflag)        
                  prob(j)=dble(ccluster(j))*exp(tmp1)
               end do
               
               do k=1,nvar
                  do l=1,nvar
                     workm1(k,l)=psiinv1(k,l)
                  end do
               end do

               call riwishart(nvar,nu1,workm1,workm2,workm3,workv1,
     &                        workmh1,workmh2,iflag)

               do k=1,nvar
                  do l=1,nvar
                     workm2(k,l)=workm1(k,l)/dble(k0)
                     sigmaclus(ncluster+1,ihmssf(k,l,nvar))=
     &                         workm1(k,l)
                  end do
               end do

               call rmvnorm(nvar,m1,workm2,workmh1,workv1,workv2) 
                  
               do k=1,nvar
                  muclus(ncluster+1,k)=workv2(k)
                  muwork(k)=muclus(ncluster+1,k)
                  do l=1,nvar
                     sigmawork(k,l)=
     &                         sigmaclus(ncluster+1,ihmssf(k,l,nvar))
                  end do
               end do      
               call dmvnd(nvar,ywork,muwork,sigmawork,tmp1,iflag)        
               prob(ncluster+1)=alpha*exp(tmp1)
               ccluster(ncluster+1)=0
               
               call simdisc(prob,nrec+100,ncluster+1,evali)

               ss(i)=evali
               ccluster(evali)=ccluster(evali)+1
               cstrt(evali,ccluster(evali))=i

               if(evali.gt.ncluster)then
                  ncluster=ncluster+1
               end if               
            end if

c++++++++++ observation in cluster with only 1 element
             
            if(ns.eq.1)then
                
               since=ss(i)

               if(since.lt.ncluster)then
                  call relabeldpms(i,since,nrec,nvar,ncluster,
     &                             cstrt,ccluster,ss,muclus,
     &                             sigmaclus,workv1,workmh1)
               end if
               ccluster(ncluster)=ccluster(ncluster)-1 
               ncluster=ncluster-1

               do j=1,nvar
                  ywork(j)=y(i,j)
               end do
 
               do j=1,ncluster
                  do k=1,nvar
                     muwork(k)=muclus(j,k)
                     do l=1,nvar
                        sigmawork(k,l)=sigmaclus(j,ihmssf(k,l,nvar))
                     end do
                  end do                
                  call dmvnd(nvar,ywork,muwork,sigmawork,tmp1,iflag)        
                  prob(j)=dble(ccluster(j))*exp(tmp1)
               end do

               do k=1,nvar
                  muwork(k)=muclus(ncluster+1,k)
                  do l=1,nvar
                     sigmawork(k,l)=
     &                         sigmaclus(ncluster+1,ihmssf(k,l,nvar))
                  end do
               end do      
               call dmvnd(nvar,ywork,muwork,sigmawork,tmp1,iflag)        
               prob(ncluster+1)=alpha*exp(tmp1)
               ccluster(ncluster+1)=0

               call simdisc(prob,nrec+100,ncluster+1,evali)

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

c++++++++++ cluster's means and variances

            ns=ccluster(ii)

            do i=1,nvar
               workv1(i)=m1(i)*(dble(k0)/(dble(k0)+dble(ns)))
               workv2(i)=0.d0
            end do

            do i=1,ns
               do j=1,nvar 
                  workv2(j)=workv2(j)+y(cstrt(ii,i),j)                  
               end do
            end do

            do i=1,nvar
               workv2(i)=workv2(i)/dble(ns)
               muwork(i)=workv1(i)+workv2(i)*
     &                   dble(ns)/(dble(k0)+dble(ns))
            end do

            do i=1,nvar
               do j=1,nvar
                  workm1(i,j)=0.d0
                  workm2(i,j)=0.d0
               end do
            end do

            do i=1,ns
               jj=cstrt(ii,i)
               do j=1,nvar 
                  do k=1,nvar 
                     workm1(j,k)=workm1(j,k)+
     &                 (y(jj,j)-workv2(j))*                  
     &                 (y(jj,k)-workv2(k))
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
            
            call rmvnorm(nvar,muwork,sigmawork,workmh1,workv1,workv2) 
            
            do i=1,nvar
               muclus(ii,i)=workv2(i)
            end do

         end do   

c++++++++++++++++++++++++++++++++++         
c+++++++ Base line distribution
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++ scale matrix of the inverted-Wishart component

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
               call inverse(sigmawork,nvar,iflag)
           
               do j=1,nvar 
                  do k=1,nvar 
                     workm1(j,k)=workm1(j,k)+sigmawork(j,k)
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
     &                        sigmaclus(i,ihmssf(j,k,nvar))
                  end do   
               end do
               call inverse(sigmawork,nvar,iflag)
            
               do j=1,nvar
                  tmp1=0.d0  
                  do k=1,nvar 
                     workm1(j,k)=workm1(j,k)+sigmawork(j,k)
                     tmp1=tmp1+dble(k0)*sigmawork(j,k)*muclus(i,k)
                  end do
                  workv2(j)=workv2(j)+tmp1
               end do         
            end do

            do i=1,nvar
               do j=1,nvar
                  workm1(i,j)=workm1(i,j)*dble(k0)
               end do
            end do

            do i=1,nvar
               workv1(i)=workv1(i)+workv2(i)
               do j=1,nvar
                  sigmawork(i,j)=s2inv(i,j)+workm1(i,j)
               end do
            end do
            call inverse(sigmawork,nvar,iflag)
            
            do i=1,nvar
               tmp1=0.d0
               do j=1,nvar
                  tmp1=tmp1+sigmawork(i,j)*workv1(j)    
               end do
               muwork(i)=tmp1
            end do
            call rmvnorm(nvar,muwork,sigmawork,workmh1,workv1,workv2) 
         
            do i=1,nvar
               m1(i)=workv2(i)
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
               call inverse(sigmawork,nvar,iflag)
            
               do j=1,nvar
                  do k=1,nvar
                     tmp1=tmp1+ywork(j)*sigmawork(j,k)*ywork(k)
                  end do
               end do
            end do   
         
            k0=rgamma(0.5d0*(dble(ncluster)+tau1),0.5d0*(tmp1+tau2))

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

               count=0

c+++++++++++++ normal baseline mean

               do i=1,nvar
                  count=count+1
                  thetasave(isave,count)=m1(i)
               end do   

c+++++++++++++ k0 parameter

               count=count+1
               thetasave(isave,count)=k0

c+++++++++++++ IW baseline scale

               do i=1,nvar
                  do j=i,nvar
                     count=count+1
                     thetasave(isave,count)=psi1(i,j)
                  end do
               end do
 
c+++++++++++++ cluster information
               
               count=count+1
               thetasave(isave,count)=ncluster
               count=count+1
               thetasave(isave,count)=alpha               

c+++++++++++++ mixture model parameters
               count=0
               do i=1,nrec
                  do j=1,nvar
                     muwork(j)=muclus(ss(i),j) 
                     do k=1,nvar
                        sigmawork(j,k)=sigmaclus(ss(i),ihmssf(j,k,nvar))
                     end do
                  end do

                  do j=1,nvar
                     count=count+1
                     randsave(isave,count)=muwork(j)
                  end do   
                  
                  do j=1,nvar
                     do k=j,nvar
                        count=count+1
                        randsave(isave,count)=sigmawork(j,k)
                     end do
                  end do               
               end do

c+++++++++++++ predictive information
       
               do i=1,ncluster
                  prob(i)=dble(ccluster(i))/(alpha+dble(nrec))
               end do
               prob(ncluster+1)=alpha/(alpha+dble(nrec))
         
               call simdisc(prob,nrec+100,ncluster+1,evali)

               do k=1,nvar
                  do l=1,nvar
                     workm3(k,l)=psiinv1(k,l)
                  end do
               end do

               call riwishart(nvar,nu1,workm3,workm1,workm2,workv1,
     &                        workmh1,workmh2,iflag)

               do k=1,nvar
                  do l=1,nvar
                     workm2(k,l)=workm3(k,l)/dble(k0)
                     sigmawork(k,l)=workm3(k,l)
                  end do
               end do

               do k=1,nvar
                   do l=k,nvar
                      sigmaclus(ncluster+1,ihmssf(k,l,nvar))=
     &                    sigmawork(k,l)
                   end do
               end do            

               call rmvnorm(nvar,m1,workm2,workmh1,workv1,workv2) 

               do k=1,nvar
                  muclus(ncluster+1,k)=workv2(k)
               end do      

               do i=1,nvar
                  do j=1,ngrid 
                     do k=1,ncluster+1
                        tmp1=muclus(k,i)
                        tmp2=sigmaclus(k,ihmssf(i,i,nvar))
                        tmp3=dnrm(grid(j,i),tmp1,sqrt(tmp2),0)
                        funi(j,i)=funi(j,i)+prob(k)*tmp3   
                     end do
                  end do
               end do


c++++++++++++++only if nvar bigger than 1+++++++++++++++++++++++++++++++++
     
               if(nvar.gt.1)then

               count2=0
               do ii=1,nvar-1
                  do jj=ii+1,nvar 
                     count2=count2+1
             
                     do k=1,ncluster+1

                        mud(1)=muclus(k,ii)
                        mud(2)=muclus(k,jj)
                        sigmad1(1,1)=sigmaclus(k,ihmssf(ii,ii,nvar))
                        sigmad1(1,2)=sigmaclus(k,ihmssf(ii,jj,nvar))
                        sigmad1(2,1)=sigmaclus(k,ihmssf(ii,jj,nvar))
                        sigmad1(2,2)=sigmaclus(k,ihmssf(jj,jj,nvar))

                        call inversedet(sigmad1,2,iflagd,tmp1)
 
                        do i=1,ngrid
                            do j=1,ngrid
                               gridd(1)=grid(i,ii)-mud(1)
                               gridd(2)=grid(j,jj)-mud(2)

                               tmp2=0.d0
                               do i1=1,2
                                  do j1=1,2
                                     tmp2=tmp2+gridd(i1)*
     &                                    sigmad1(i1,j1)*
     &                                    gridd(j1)          
                                 end do
                               end do

                              tmp3=(-(2.d0*log(tpi))-tmp1-tmp2)/2.d0
                                   
                              i1=(count2-1)*ngrid*ngrid+(i-1)*ngrid+j
                              fbiv(i1)=fbiv(i1)+prob(k)*exp(tmp3)   
                           end do
                         end do
                     end do
                  end do
               end do
               end if
 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

               do i=1,nvar  
                  muwork(i)=muclus(evali,i)
                  do j=1,nvar
                     sigmawork(i,j)=
     &                  sigmaclus(evali,ihmssf(i,j,nvar))
                  end do
               end do  
               
               do i=1,nvar
                  count=count+1 
                  randsave(isave,count)=muwork(i)
               end do
               
               do i=1,nvar
                  do j=i,nvar
                     count=count+1
                     randsave(isave,count)=sigmawork(i,j)
                  end do
               end do
               
               call rmvnorm(nvar,muwork,sigmawork,workmh1,workv1,workv2) 

               do i=1,nvar
                  count=count+1 
                  randsave(isave,count)=workv2(i)
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

      do i=1,ngrid
         do j=1,nvar
            funi(i,j)=funi(i,j)/dble(nsave) 
         end do
      end do

      if(nvar.gt.1)then
         do i=1,ngridb
            fbiv(i)=fbiv(i)/dble(nsave) 
         end do
      end if

      return
      end


c=======================================================================      
      subroutine relabeldpms(ind,since,nrec,nvar,ncluster,cstrt,
     &                       ccluster,ss,muclus,sigmaclus,
     &                       workv1,workmh1)
c=======================================================================
c     relabel the clusters after the elimination of one of them
c     A.J.V., 2005
c=======================================================================
      implicit none
      integer i,j,ind,since,nrec,nvar,ncluster,ccluster(nrec)
      integer ss(nrec),cstrt(nrec,nrec)
      double precision muclus(nrec+100,nvar),
     1  sigmaclus(nrec+100,nvar*(nvar+1)/2)
      double precision workv1(nvar),workmh1(nvar*(nvar+1)/2)

      integer ns,ii,nt

      nt=nvar*(nvar+1)/2 


      do i=1,nvar
         workv1(i)=muclus(since,i)
      end do
     
      do i=1,nt
         workmh1(i)=sigmaclus(since,i)
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

         do j=1,nt
            sigmaclus(i-1,j)=sigmaclus(i,j)
         end do

         ccluster(i-1)=ccluster(i)
      end do
      
      ss(ind)=ncluster
      
      do i=1,nvar
         muclus(ncluster,i)=workv1(i)
      end do
      do i=1,nt
         sigmaclus(ncluster,i)=workmh1(i)
      end do
      ccluster(ncluster)=1
      
      return
      end  


