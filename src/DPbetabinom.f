c=======================================================================
      subroutine dpbetabinom(nrec,y,
     &                       ngrid,grid, 
     &                       a0,b0,a1,b1,
     &                       ncluster,ss,alpha,p,
     &                       mcmc,nsave,
     &                       cpo,densm,thetasave,randsave,
     &                       ccluster,cstrt,prob,lprob,workcpo,
     &                       seed)
c=======================================================================
c     # of arguments = 24.
c
c     Subroutine `dpbetabinom' to run a Markov chain for a 
c     semiparametric Beta-Binomial model using a DP prior.
c
c     Copyright: Alejandro Jara, 2009-2010.
c
c     Version 1.0: 
c
c     Last modification: 05-07-2007.
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
c        y           :  real vector giving the transformed data, y(nrec). 
c
c-----------------------------------------------------------------------
c
c---- Prediction -------------------------------------------------------
c 
c        ngrid       :  integer giving the number of grid points. 
c        grid        :  real vector giving the grid points, grid(ngrid). 
c
c-----------------------------------------------------------------------
c
c---- Prior information ------------------------------------------------
c 
c        a0, b0      :  real giving the hyperparameters of the prior
c                       distribution for the precision parameter,
c                       alpha ~ Gamma(a0,b0). If a0<0 the precision 
c                       parameter is considered as a constant.
c        a1,b1       :  real giving the parameters of the Beta centering
c                       distribution.
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        ncluster    :  integer giving the number of clusters.
c        ss          :  integer vector giving the configurations.
c        alpha       :  real giving the current value of the precision
c                       parameter of the DP.
c        p           :  real vector giving the value of the binomial 
c                       probabilities, p(nrec+1).
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
c---- Output -----------------------------------------------------------
c
c        cpo         :  real matrix giving the cpo's, cpo(nrec,2). 
c        densm       :  real vector giving the posterior mean of the
c                       density estimate, densm(ngrid).
c        thetasave   :  real matrix containing the mcmc samples for
c                       k, the number of clusters and alpha, 
c                       thetsave(nsave,3).
c        randsave    :  real matrix containing the mcmc samples for
c                       the latent variables and prediction,
c                       randsave(nsave,nrec+1).
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        ccluster    :  integer vector indicating the number of
c                       subjects in each cluster, ccluster(nrec).
c        cstrt       :  integer matrix used to save the cluster
c                       structure, cstrt(nrec,nrec).
c        prob        :  real vector used to update the cluster 
c                       structure, prob(nrec+1).
c        lprob       :  real vector used to update the cluster 
c                       structure, prob(nrec+1).
c        workcpo     :  real vector used to compute the cpo, 
c                       workcpo(nrec).
c        seed        :  integer vector giving the seeds for the random
c                       number generator, seed(2).
c
c=======================================================================                  
      implicit none
c++++ data
      integer nrec
      double precision y(nrec,2)

c++++ prediction
      integer ngrid
      double precision grid(ngrid)

c++++ prior
      double precision a0,b0
      double precision a1,b1
      
c++++ current value
      integer ncluster
      integer ss(nrec)
      double precision alpha
      double precision p(nrec+1)

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++output
      double precision cpo(nrec,2)
      double precision densm(ngrid)
      double precision thetasave(nsave,3)
      double precision randsave(nsave,nrec+1)

c+++++External Working space
      integer ccluster(nrec)
      integer cstrt(nrec,nrec)
      double precision prob(nrec+1)
      double precision lprob(nrec+1)
      double precision workcpo(nrec)

c+++++External Working space - RNG
      integer seed(2),seed1,seed2

c+++++Internal Working space
      integer dispcount
      integer evali
      integer i,ii,j
      integer iscan,isave
      integer ns,nscan
      integer ok
      integer since
      integer skipcount
      integer sprint
      double precision dbin,dbet
      double precision mbetabin
      double precision pwork
      double precision tmp1,tmp2,tmp3
      
c+++++CPU time
      double precision sec00,sec0,sec1,sec


c++++ DP (functional parameter)
      double precision eps,rbeta
      parameter(eps=0.01)

c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)

c++++ set random number generator
      seed1=seed(1)
      seed2=seed(2)
      
      call setall(seed1,seed2)

c++++ cluster structure
      do i=1,nrec
         ccluster(ss(i))=ccluster(ss(i))+1
         cstrt(ss(i),ccluster(ss(i)))=i
         prob(i)=0.d0
         lprob(i)=0.d0
      end do
      pwork=0.d0
     
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

c++++++++++++++++++++++++++++++
c+++++++ a) Polya Urn 
c++++++++++++++++++++++++++++++

         do i=1,nrec

c++++++++++ check if the user has requested an interrupt
            call rchkusr()
         
            ns=ccluster(ss(i))

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

             else

               since=ss(i)

               if(since.lt.ncluster)then

                  call relabeldpbetabin(i,since,nrec,ncluster,
     &                                  ccluster,ss,cstrt,p)
               end if
               ccluster(ncluster)=ccluster(ncluster)-1 
               ncluster=ncluster-1
            end if 

            do j=1,ncluster
               tmp2=dbin(dble(y(i,1)),dble(y(i,1)+y(i,2)),p(j),0)
               prob(j)=dble(ccluster(j))*tmp2
            end do

            tmp1=a1+dble(y(i,1))
            tmp2=b1+dble(y(i,2))
            p(ncluster+1)=rbeta(tmp1,tmp2)

            tmp1=dble(y(i,1))
            tmp2=dble(y(i,2))
            prob(ncluster+1)=alpha*exp(mbetabin(a1,b1,tmp1,tmp2)) 

            call simdisc(prob,nrec+1,ncluster+1,evali)

            ss(i)=evali
            ccluster(evali)=ccluster(evali)+1
            cstrt(evali,ccluster(evali))=i
                
            if(evali.gt.ncluster)then
               ncluster=ncluster+1
            end if

         end do
         
c++++++++++++++++++++++++++++++
c+++++++ b) Resampling step 
c++++++++++++++++++++++++++++++

         do i=1,ncluster

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            ns=ccluster(i)
            tmp1=a1
            tmp2=b1
            
            do j=1,ns

c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

               ii=cstrt(i,j)
               tmp1=tmp1+dble(y(ii,1))
               tmp2=tmp2+dble(y(ii,2))
            end do

            p(i)=rbeta(tmp1,tmp2)
         end do

c++++++++++++++++++++++++++++++++++         
c+++++++ Precision parameter
c++++++++++++++++++++++++++++++++++
         if(a0.gt.0.d0)then
            call samalph(alpha,a0,b0,ncluster,nrec)
         end if 
          

c++++++++++++++++++++++++++++++++++         
c+++++++ save samples
c++++++++++++++++++++++++++++++++++         
         
         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1

c+++++++++++++ cluster information

               thetasave(isave,1)=ncluster
               thetasave(isave,2)=alpha

c+++++++++++++ subject specific information

               do i=1,nrec
                  randsave(isave,i)=p(ss(i))
               end do

c+++++++++++++ predictive information

               do i=1,ncluster
                   prob(i)=dble(ccluster(i))/(alpha+dble(nrec))
               end do
               prob(ncluster+1)=alpha/(alpha+dble(nrec))
               call simdisc(prob,nrec+1,ncluster+1,evali)

               pwork=0.d0
               if(evali.le.ncluster)then
                  pwork=p(evali)
               end if
               if(evali.eq.ncluster+1)then 
                  pwork=rbeta(a1,b1)
               end if
               randsave(isave,nrec+1)=pwork

               do i=1,ncluster
                  ns=ccluster(i)
                  tmp1=a1
                  tmp2=b1
            
                  do j=1,ns
                     ii=cstrt(i,j)
                     tmp1=tmp1+dble(y(ii,1))
                     tmp2=tmp2+dble(y(ii,2))
                  end do

                  do j=1,ngrid
                     densm(j)=densm(j)+
     &                        prob(i)*dbet(grid(j),tmp1,tmp2,0) 
                  end do
               end do
               do j=1,ngrid
                  densm(j)=densm(j)+
     &                     prob(ncluster+1)*dbet(grid(j),a1,b1,0) 
               end do

c+++++++++++++ cpo

               do i=1,nrec
                  workcpo(i)=0.d0
               end do   

               do i=1,nrec
                  do j=1,ncluster
                     prob(j)=dble(ccluster(i))/(alpha+dble(nrec))
                     tmp2=dbin(dble(y(i,1)),dble(y(i,1)+y(i,2)),p(j),0)
                     workcpo(i)=workcpo(i)+prob(j)*tmp2
                  end do
                  prob(ncluster+1)=alpha/(alpha+dble(nrec))

                  tmp1=dble(y(i,1))
                  tmp2=dble(y(i,2))
                  workcpo(i)=workcpo(i)+prob(ncluster+1)*
     &                       exp(mbetabin(a1,b1,tmp1,tmp2)) 
               end do

               tmp2=0.d0
               do i=1,nrec
                  tmp3=workcpo(i)
                  cpo(i,1)=cpo(i,1)+1.0d0/tmp3  
                  cpo(i,2)=cpo(i,2)+tmp3                   
                  tmp2=tmp2+log(dble(isave)/cpo(i,1))
               end do
               thetasave(isave,3)=tmp2

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

      do j=1,ngrid
         densm(j)=densm(j)/dble(nsave)
      end do
      

      return
      end
      
c=======================================================================
      subroutine relabeldpbetabin(ind,since,nrec,ncluster,ccluster,
     &                            ss,cstrt,p)
c=======================================================================
      implicit none
      integer ind,since,nrec,ncluster,ccluster(nrec)
      integer ss(nrec),cstrt(nrec,nrec)
      double precision p(nrec+1)
 
      integer i,ii,j,ns
      double precision pwork
      
      pwork=p(since) 

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
         p(i-1)=p(i)
         ccluster(i-1)=ccluster(i)
      end do

      ss(ind)=ncluster

      p(ncluster)=pwork         
      ccluster(ncluster)=1
      
      return
      end

c=======================================================================
      double precision function mbetabin(a1,b1,x1,x2)
c=======================================================================
      implicit none
      double precision a1,b1
      double precision x1,x2
      double precision dgamlog
      mbetabin=dgamlog(a1+x1)+dgamlog(b1+x2)+dgamlog(a1+b1)
      mbetabin=mbetabin-dgamlog(a1+b1+x1+x2)
      mbetabin=mbetabin-dgamlog(a1)-dgamlog(b1)
      mbetabin=mbetabin+dgamlog(x1+x2+1.d0)
      mbetabin=mbetabin-dgamlog(x1+1.d0)-dgamlog(x2+1.d0)
      return
      end

c=======================================================================
      subroutine sindidp(maxn,n,lprob,aprob,evali)
c=======================================================================
      implicit none
      integer maxn,n
      double precision lprob(maxn),aprob(maxn)

c+++++internal working space 
      integer evali
      
c+++++internal working space 
      integer cond,i
      double precision mmax,tmp1
      real runif
      
      mmax=lprob(1)
      do i=2,n
         if(lprob(i).gt.mmax) mmax=lprob(i)
      end do
      do i=1,n
         lprob(i)=dexp(lprob(i)-mmax)
      end do

      aprob(1)=lprob(1)
      do i=2,n
         aprob(i)=aprob(i-1)+lprob(i)
      enddo
      do i=1,n
         aprob(i)=aprob(i)/aprob(n)
      enddo

      tmp1=dble(runif())  
      i=0
      cond=0
      do while((cond.eq.0).and.(i.le.n))
         i=i+1
         if (aprob(i).ge.tmp1) cond=1
      end do
      evali=i
      return
      end
      
            
