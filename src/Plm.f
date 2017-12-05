c=======================================================================
      subroutine plm(nrec,p,x,y,
     &               tau1,tau2,ms,sbeta0i,
     &               beta,sigma2,
     &               mcmc,nsave,
     &               cpo,thetasave,
     &               iflagp,betam,workmh,xtx,xty,
     &               seed)
c=======================================================================
c     Subroutine `plm' to run a Markov chain for a 
c     parametric linear model.
c
c     Copyright: Alejandro Jara, 2009-2010.
c
c     Version 1.0: 
c
c     Last modification: 15-03-2009.
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
c        x           :  real matrix giving the design matrix, x(nrec,p).
c        y           :  real vector giving the responses, y(nrec). 
c
c-----------------------------------------------------------------------
c
c---- Prior information ------------------------------------------------
c 
c        tau1,tau2   :  real giving the value of the hyperparameters
c                       in the inverted-Gamma prior for sigma2. 
c        ms          :  real vector giving the product of the prior
c                       precision and prior mean for beta, ms(p).
c        sbeta0i     :  real vector giving the inverse of the covariance
c                       matrix of the normal prior for beta, 
c                       sbeta0i(p,p).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        beta        :  real vector giving the value of the regression  
c                       coefficients, beta(p).
c        sigma2      :  real giving the value of the error variance.
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
c        cpo         :  real giving the cpos and fsos, cpo(nrec,2). 
c        thetasave   :  real matrix containing the mcmc samples for
c                       the parameters, thetasave(nsave,p+2)
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        betam       :  real vector used to update the reg. coeff.,
c                       betam(p).
c        iflagp      :  integer vector used to evaluate reg. coeff.,
c                       iflagp(p).
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        workmh      :  real vector used to update the reg. coeff.,
c                       workmh(p*(p+1)/2).
c        xtx         :  real matrix used to update the reg. coeff.,
c                       xtx(p,p).
c        xty         :  real vector used to update the reg. coeff.,
c                       xty(p).
c
c=======================================================================                  

      implicit none
c++++ data
      integer nrec,p
      double precision x(nrec,p)
      double precision y(nrec)

c++++ prior
      double precision tau1,tau2
      double precision ms(p)
      double precision sbeta0i(p,p)
     
c++++ current value
      double precision beta(p)
      double precision sigma2

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++output
      double precision cpo(nrec,2)
      double precision thetasave(nsave,p+2)

c+++++External Working space
      integer iflagp(p)
      double precision betam(p)
      double precision workmh(p*(p+1)/2)
      double precision xtx(p,p)
      double precision xty(p)

c+++++External Working space - RNG
      integer seed(2),seed1,seed2

c+++++Internal Working space
      integer dispcount
      integer i,j,k
      integer iscan,isave
      integer nscan
      integer skipcount
      integer sprint
      double precision dnrm
      double precision rgamma
      double precision tmp1,tmp2,tmp3
      
c+++++CPU time
      double precision sec00,sec0,sec1,sec

c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)

c++++ set random number generator
      seed1=seed(1)
      seed2=seed(2)
      
      call setall(seed1,seed2)

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

c++++++++++++++++++++++++++++++++++
c+++++++ regression coefficients 
c+++++++++++++++++++++++++++++++++

         do i=1,p
            do j=1,p
               xtx(i,j)=sbeta0i(i,j)
            end do
            xty(i)=ms(i)
         end do   
         
         do i=1,nrec

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            do j=1,p
               do k=1,p
                  xtx(j,k)=xtx(j,k)+x(i,j)*x(i,k)/sigma2
               end do
               xty(j)=xty(j)+x(i,j)*y(i)/sigma2
            end do
 
         end do
 
         call inverse(xtx,p,iflagp)
         betam=matmul(xtx,xty)
         call rmvnorm(p,betam,xtx,workmh,xty,beta)
 
c++++++++++++++++++++++++++++++
c+++++++ error variance
c++++++++++++++++++++++++++++++
         tmp1=0.d0
         do i=1,nrec
            tmp2=0.d0
            do j=1,p
               tmp2=tmp2+x(i,j)*beta(j)
            end do
            tmp1=tmp1+(y(i)-tmp2)*(y(i)-tmp2)
         end do

         sigma2=1.d0/rgamma(0.5d0*(tau1+dble(nrec)),
     &                      0.5d0*(tau2+tmp1))

c++++++++++++++++++++++++++++++++++         
c+++++++ save samples
c++++++++++++++++++++++++++++++++++         
         
         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1

c+++++++++++++ regression coefficients

               do i=1,p 
                  thetasave(isave,i)=beta(i)
               end do

c+++++++++++++ error variance

               thetasave(isave,p+1)=sigma2

c+++++++++++++ cpo

               tmp3=0.d0
               do i=1,nrec
                  tmp1=0.d0
                  do j=1,p
                     tmp1=tmp1+x(i,j)*beta(j)
                  end do
                  tmp2=dnrm(y(i),tmp1,sqrt(sigma2),0)

                  cpo(i,1)=cpo(i,1)+1.0d0/tmp2  
                  cpo(i,2)=cpo(i,2)+tmp2                   
                  tmp3=tmp3+log(dble(isave)/cpo(i,1))
               end do
               thetasave(isave,p+2)=tmp3
               
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

     

      return
      end
      
