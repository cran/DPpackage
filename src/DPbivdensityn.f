c=======================================================================                      
      subroutine bivspdenn(ngrid,nrec,nvar,y,a0b0,k0,nuvec,m1rand,s2inv,
     &                 s2invm2,psiinv2,tau,mcmc,nsave,cpo,f,fun1,fun2,
     &                 randsave,thetasave,alpha,m1,muclus,ncluster,psi1,
     &                 psiinv1,s1,sigmaclus,ss,ccluster,grid1,grid2,
     &                 iflag,muwork,muwork2,
     &                 prob,seed,sigmawork,sigmawork2,sigworkinv,theta,
     &                 workm1,workm2,workm3,workmh1,workmh2,workv1,
     &                 workv2,workv3,ywork,workcpo)
c=======================================================================                      
c
c     Subroutine `bivspdenn' to run a Markov chain in the DP mixture of  
c     bivariate normals model. In this routine, inference is based on the 
c     Polya urn representation of the Dirichlet process. The algorithm
c     8 of Neal (2000) is used with m=1. The difference with respect to
c     `spdenn` is that this includes a grid of points where the density
c     estimate is evaluated. This fucntion can be used only for dimension
c     <=2.
c
c     Copyright: Alejandro Jara, 2006-2010.
c
c     Version 1.0: 
c
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
c        ngrid       :  integer giving the size of the grid where
c                       the density estimate is evaluated.
c        nrec        :  integer giving the number of observations.
c        nvar        :  integer giving the number of variables.
c        y           :  real matrix giving the response variables,
c                       y(nrec,nvar).
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
c        f           :  real matrix giving the density estimate at the
c                       grid, f(ngrid,ngrid).
c        fun1        :  real giving the marginal density estimate at the
c                       grid, fun1(ngrid).
c        fun2        :  real giving the marginal density estimate at the
c                       grid, fun2(ngrid).
c        cpo         :  real giving the cpo. 
c        randsave    :  real matrix containing the mcmc samples for
c                       the parameters of the density, 
c                       randsave(nsave,(nrec+2)*nvar+
c                       (nrec+1)*nvar*(nvar+1)/2+nvar).
c        thetasave   :  real matrix containing the mcmc samples for
c                       the parameters, 
c                       thetasave(nsave,nvar+nvar*(nvar+1)/2+3).
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
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        ccluster    :  integer vector indicating the number of
c                       subjects in each cluster, ccluster(nrec).
c        count       :  index.           
c        detlog      :  real used to save the log-determinant in an
c                       matrix inversion process.
c        dispcount   :  index. 
c        dnrm        :  density of a normal distribution.
c        evali       :  integer indicator used in updating the state.
c        grid1       :  real vector giving the grid where the density
c                       estimate is evaluated for first coordinate,
c                       grid1(ngrid).
c        grid2       :  real vector giving the grid where the density
c                       estimate is evaluated for second coordinate,
c                       grid2(ngrid).
c        i           :  index. 
c        ii          :  index. 
c        ihmssf      :  integer function to determine the position of a
c                       half-stored matrix.
c        iflag       :  integer vector used to evaluate the mvn density,
c                       iflag(nvar).
c        isave       :  index. 
c        iscan       :  index.
c        j           :  index. 
c        k           :  index. 
c        l           :  index. 
c        l1          :  index. 
c        l2          :  index. 
c        ns          :  integer indicator used in updating the state. 
c        nscan       :  integer indicating the total number of MCMC
c                       scans.
c        nuniqs      :  integer giving the dimension of the half-stored
c                       covariance matrix.
c        muwork      :  real vector used to save the mean,
c                       one observation, muwork(nvar).
c        muwork2     :  real vector used to save the mean,
c                       one observation, muwork2(nvar).
c        nuwork      :  index.
c        prob        :  real vector used to update the cluster 
c                       structure, prob(nrec+100).
c        rgamma      :  gamma random number generator
c        s1          :  real matrix giving the covariance matrix of 
c                       the normal component of the baseline 
c                       distribution, s1(nvar,nvar).
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        seed3       :  seed for random number generation.
c        sigmawork   :  real matrix used to save the variance of
c                       one observation, sigmawork(nvar,nvar).
c        sigmawork2  :  real matrix used to save the variance of
c                       one observation, sigmawork2(nvar,nvar).
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
c        workmh1     :  real vector used to update the cluster
c                       structure, workmh1(nvar*(nvar+1)/2).
c        workmh2     :  real vector used to update the cluster
c                       structure, workmh2(nvar*(nvar+1)/2).
c        workv1      :  real vector used to update the cluster
c                       structure, workv1(nvar).
c        workv2      :  real vector used to update the cluster
c                       structure, workv2(nvar).
c        workv3      :  real vector used to update the cluster
c                       structure, workv3(nvar).
c        ywork       :  real vector used to save the variables of,
c                       one observation, ywork(nvar).
c
c=======================================================================

      implicit none 

c+++++Data
      integer ngrid,nrec,nvar
      double precision y(nrec,nvar)

c+++++Prior 
      integer nuvec(2),nu1,nu2,m1rand
      double precision aa0,ab0,a0b0(2)
      double precision psiinv2(nvar,nvar)
      double precision tau(2),tau1,tau2
      double precision s2inv(nvar,nvar),s2invm2(nvar)

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Output
      double precision cpo(nrec)
      double precision f(ngrid,ngrid),fun1(ngrid),fun2(ngrid)
      double precision randsave(nsave,
     1  (nrec+2)*nvar+(nrec+1)*nvar*(nvar+1)/2)
      double precision thetasave(nsave,nvar+nvar*(nvar+1)/2+3)

c+++++Current values of the parameters
      integer ncluster,ss(nrec)
      double precision alpha,k0,m1(nvar)
      double precision muclus(nrec+100,nvar)
      double precision psi1(nvar,nvar),psiinv1(nvar,nvar)
      double precision sigmaclus(nrec+100,nvar*(nvar+1)/2)
      
c+++++Working space
      integer ccluster(nrec),count,dispcount,evali
      integer i,ii,iflag(nvar),ihmssf,isave,iscan
      integer j,k,l,l1,l2
      integer ns,nscan,nuniqs,nuwork,sprint
      integer seed(3),seed1,seed2,seed3,since,skipcount
      double precision detlog,dnrm
      double precision muwork(nvar),muwork2(nvar),prob(nrec+100),rgamma
      double precision s1(nvar,nvar)
      double precision sigmawork(nvar,nvar),sigmawork2(nvar,nvar)
      double precision sigworkinv(nvar,nvar)
      double precision theta(nvar),tmp1,tmp2
      double precision workm1(nvar,nvar),workm2(nvar,nvar),
     1  workm3(nvar,nvar)
      double precision workmh1(nvar*(nvar+1)/2),workmh2(nvar*(nvar+1)/2)
      double precision workv1(nvar),workv2(nvar),workv3(nvar)
      double precision ywork(nvar)
      double precision workcpo(nrec)

c+++++Working space - Density
      double precision grid1(ngrid),grid2(ngrid)

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
      
      nuniqs=nvar*(nvar+1)/2
      nu1=nuvec(1)
      nu2=nuvec(2)

c++++ set random number generator

      seed1=seed(1)
      seed2=seed(2)
      seed3=seed(3)

      call setall(seed1,seed2)
     
c++++ cluster structure

      do i=1,nrec
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
c+++++++ DP part
c++++++++++++++++++++++++++++++++++

c++++++++++++++++++++++++++++++
c+++++++ a) configurations 
c++++++++++++++++++++++++++++++

         do i=1,nrec

            ns=ccluster(ss(i))
            
c++++++++++ observation in cluster with more than 1 element
             
            if(ns.gt.1)then
 
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
                   
                  call dmvn(nvar,ywork,muwork,sigmawork,tmp1,
     &                      workv1,workm1,workm2,workv2,iflag)

                  prob(j)=dble(ccluster(j))*exp(tmp1)
               end do
               
               do k=1,nvar
                  do l=1,nvar
                     workm3(k,l)=psiinv1(k,l)
                  end do
               end do

               call riwishart(nvar,nu1,workm3,workm1,workm2,workv1,
     &                           workmh1,workmh2,iflag)

               do k=1,nvar
                  do l=1,nvar
                     s1(k,l)=workm3(k,l)/dble(k0)
                     sigmaclus(ncluster+1,ihmssf(k,l,nvar))=
     &                         workm3(k,l)
                  end do
               end do

               call rmvnorm(nvar,m1,s1,workmh1,workv1,theta) 
                  
               do k=1,nvar
                  muclus(ncluster+1,k)=theta(k)
                  muwork(k)=muclus(ncluster+1,k)
                  do l=1,nvar
                     sigmawork(k,l)=
     &                         sigmaclus(ncluster+1,ihmssf(k,l,nvar))
                  end do
               end do      

               call dmvn(nvar,ywork,muwork,sigmawork,tmp1,
     &                      workv1,workm1,workm2,workv2,iflag)
                   
               prob(ncluster+1)=alpha*exp(tmp1)
               
               call simdisc(prob,nrec+100,ncluster+1,evali)
               
               
               if(evali.le.ncluster)then
                  ss(i)=evali
                  ccluster(evali)=ccluster(evali)+1
               end if   
               
               if(evali.gt.ncluster)then
                  ncluster=ncluster+1
                  ss(i)=ncluster
                  ccluster(ncluster)=1

                  do j=1,nvar
                     muclus(ncluster,j)=muclus(evali,j)
                     do k=j,nvar
                        sigmaclus(ncluster,ihmssf(j,k,nvar))=
     &                            sigmaclus(evali,ihmssf(j,k,nvar))
                     end do
                  end do
               end if               
            end if

c++++++++++ observation in cluster with only 1 element
             
            if(ns.eq.1)then
                
               since=ss(i)

               if(since.lt.ncluster)then
                   call relabeldd(i,since,nrec,nvar,ncluster,
     &                           ccluster,ss,muclus,sigmaclus,
     &                           muwork,sigmawork)                   
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
                   
                  call dmvn(nvar,ywork,muwork,sigmawork,tmp1,
     &                      workv1,workm1,workm2,workv2,iflag)


                  prob(j)=dble(ccluster(j))*exp(tmp1)
               end do

               do k=1,nvar
                  muwork(k)=muclus(ncluster+1,k)
                  do l=1,nvar
                     sigmawork(k,l)=
     &                         sigmaclus(ncluster+1,ihmssf(k,l,nvar))
                  end do
               end do      

               call dmvn(nvar,ywork,muwork,sigmawork,tmp1,
     &                      workv1,workm1,workm2,workv2,iflag)
                   
               prob(ncluster+1)=alpha*exp(tmp1)


               call simdisc(prob,nrec+100,ncluster+1,evali)


               if(evali.le.ncluster)then
                  ss(i)=evali
                  ccluster(evali)=ccluster(evali)+1
               end if   
               
               if(evali.gt.ncluster)then
                  ncluster=ncluster+1
                  ss(i)=ncluster
                  ccluster(ncluster)=1

                  do j=1,nvar
                     muclus(ncluster,j)=muclus(evali,j)
                     do k=j,nvar
                        sigmaclus(ncluster,ihmssf(j,k,nvar))=
     &                            sigmaclus(evali,ihmssf(j,k,nvar))
                     end do
                  end do
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

            do i=1,nrec
               if(ss(i).eq.ii)then
                 do j=1,nvar 
                    workv2(j)=workv2(j)+y(i,j)                  
                 end do
               end if          
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

            do i=1,nrec
               if(ss(i).eq.ii)then
                 do j=1,nvar 
                    do k=1,nvar 
                       workm1(j,k)=workm1(j,k)+
     &                       (y(i,j)-workv2(j))*                  
     &                       (y(i,k)-workv2(k))
                    end do   
                 end do
               end if          
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
     &                        sigmaclus(i,ihmssf(j,k,nvar))
               end do   
            end do
            
            call invdet(sigmawork,nvar,sigworkinv,detlog,iflag,workv3)

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
     &                  workv1,workmh1,workmh2,iflag)


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
            
            call invdet(sigmawork,nvar,sigworkinv,detlog,iflag,workv3)

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
               workm1(i,j)=workm1(i,j)*dble(k0)
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
            
            call invdet(sigmawork,nvar,sigworkinv,detlog,iflag,workv3)

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

c+++++++++++++ random effects
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
                     s1(k,l)=workm3(k,l)/dble(k0)
                     sigmawork(k,l)=workm3(k,l)
                  end do
               end do

               call rmvnorm(nvar,m1,s1,workmh1,workv1,muwork) 

               do i=1,ngrid
                  tmp1=0.d0
                  do j=1,ncluster
                     tmp2=dnrm(grid1(i),muclus(j,1),
     &                    sqrt(sigmaclus(j,1)),0)
                     tmp1=tmp1+tmp2*prob(j)
                  end do

                  tmp2=dnrm(grid1(i),muwork(1),
     &                 sqrt(sigmawork(1,1)),0)
                  tmp1=tmp1+tmp2*prob(ncluster+1)
                  
                  fun1(i)=fun1(i)+tmp1                    
               end do

               if(nvar.eq.2)then
                 do i=1,ngrid
                    tmp1=0.d0
                    do j=1,ncluster
                       tmp2=dnrm(grid2(i),muclus(j,2),
     &                      sqrt(sigmaclus(j,3)),0)
                       tmp1=tmp1+tmp2*prob(j)
                    end do

                    tmp2=dnrm(grid2(i),muwork(2),
     &                   sqrt(sigmawork(2,2)),0)
                    tmp1=tmp1+tmp2*prob(ncluster+1)
                  
                    fun2(i)=fun2(i)+tmp1                    
                 end do
                     
                 do i=1,ngrid
                    ywork(1)=grid1(i)
                    do j=1,ngrid
                       ywork(2)=grid2(j)
                       tmp1=0.d0
                       do k=1,ncluster
                          do l1=1,nvar
                             muwork2(l1)=muclus(k,l1)
                             do l2=1,nvar
                                sigmawork2(l1,l2)=
     &                          sigmaclus(k,ihmssf(l1,l2,nvar))
                             end do
                          end do                
                   
                          call dmvn(nvar,ywork,muwork2,sigmawork2,
     &                      tmp2,workv1,workm1,workm2,workv2,iflag)

                          tmp1=tmp1+exp(tmp2)*prob(k)
                       end do
                      
                       call dmvn(nvar,ywork,muwork,sigmawork,
     &                      tmp2,workv1,workm1,workm2,workv2,iflag)
                       
                       tmp1=tmp1+exp(tmp2)*prob(ncluster+1)
                       
                       f(i,j)=f(i,j)+tmp1                    
                    end do
                 end do
               end if

               if(evali.le.ncluster)then
                  do i=1,nvar  
                     muwork(i)=muclus(evali,i)
                     do j=1,nvar
                        sigmawork(i,j)=sigmaclus(evali,ihmssf(i,j,nvar))
                     end do
                  end do  
               end if
               
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
               
               call rmvnorm(nvar,muwork,sigmawork,workmh1,workv1,theta) 
               
               do i=1,nvar
                  count=count+1 
                  randsave(isave,count)=theta(i)
               end do

               
c+++++++++++++ cpo and save samples

               do i=1,nrec
                  workcpo(i)=0.d0
               end do   

               do i=ncluster+1,ncluster+100
                  do k=1,nvar
                     do l=1,nvar
                        workm3(k,l)=psiinv1(k,l)
                     end do
                  end do

                  call riwishart(nvar,nu1,workm3,workm1,workm2,workv1,
     &                           workmh1,workmh2,iflag)

                 do k=1,nvar
                    do l=1,nvar
                       s1(k,l)=workm3(k,l)/dble(k0)
                       sigmaclus(i,ihmssf(k,l,nvar))=
     &                         workm3(k,l)
                    end do
                 end do

                 call rmvnorm(nvar,m1,s1,workmh1,workv1,theta) 
                 do k=1,nvar
                     muclus(i,k)=theta(k)
                 end do      
               end do

               do ii=1,ncluster+100
                  do i=1,nrec
                     if(ii.le.ncluster)then
                        ns=ccluster(ii)
                        if(ss(i).eq.ii)ns=ns-1
                        prob(ii)=dble(ns)/(alpha+dble(nrec-1))
                      else
                        prob(ii)=alpha/(100.d0*(alpha+dble(nrec-1)))
                     end if   

                     do j=1,nvar
                        ywork(j)=y(i,j)
                        muwork(j)=muclus(ii,j) 
                        do k=1,nvar
                           sigmawork(j,k)=sigmaclus(ii,ihmssf(j,k,nvar))
                        end do
                     end do

                     call dmvn(nvar,ywork,muwork,sigmawork,tmp1,
     &                      workv1,workm1,workm2,workv2,iflag)
 
                     workcpo(i)=workcpo(i)+prob(ii)*exp(tmp1)
                  end do
               end do

               tmp2=0.d0
               do i=1,nrec 
                  tmp1=workcpo(i)
                  cpo(i)=cpo(i)+1.0d0/tmp1
                  tmp2=tmp2+log(dble(isave))/cpo(i)  
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
         cpo(i)=dble(nsave)/cpo(i)
      end do

      do i=1,ngrid
         fun1(i)=fun1(i)/dble(nsave)       
         fun2(i)=fun2(i)/dble(nsave)       
         do j=1,ngrid      
            f(i,j)=f(i,j)/dble(nsave)
         end do   
      end do


      return
      end
