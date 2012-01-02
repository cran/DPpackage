c=======================================================================
      subroutine rsbamdensity(nrec,y,
     &                        ngrid,grid,
     &                        maxn,maxr,maxc,mu0,sigma20,mu1,sigma21,
     &                        tau1,taus1,taus2,m1,s1,tauh1,tauh2,
     &                        mcmc,nsave,
     &                        ss,sba,sigma2,tau2,
     &                        cpo,densave,sbasave,thetasave,
     &                        seed,
     &                        ccluster,cstrt,possv,fw,prob,theta,w)
c=======================================================================
c
c     Subroutine `rsbamdensity' to run a Markov chain for a 
c     Random Sequential Barycenter Array (RSBA) mixture of
c     normals model for density estimation.
c
c     Copyright: Alejandro Jara, 2011.
c
c     Version 1.0: 
c
c     Last modification: 22-07-2011.
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
c=======================================================================
      implicit none

c+++++Data
      integer nrec
      real*8 y(nrec)

c+++++Density
      integer ngrid
      real*8 grid(ngrid)
 
c+++++Prior
      integer maxn,maxr,maxc
      real*8 mu0,sigma20
      real*8 mu1,sigma21
      real*8 tau1
      real*8 taus1,taus2
      real*8 m1,s1
      real*8 tauh1,tauh2

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Current value of the parameters
      integer ss(nrec)
      real*8 sba(maxn,maxr)
      real*8 sigma2(maxc)
      real*8 tau2

c+++++Output
      real*8 cpo(nrec,2)
      real*8 densave(nsave,ngrid) 
      real*8 sbasave(nsave,maxr+maxc)
      real*8 thetasave(nsave,3)

c+++++Seeds
      integer seed(2),seed1,seed2

c+++++External working space
      integer ccluster(maxc)
      integer cstrt(maxc,nrec)
      integer possv(maxn)
      real*8 fw(maxn-1,maxc-1)
      real*8 prob(maxc)
      real*8 theta(maxc)
      real*8 w(maxc)

c+++++Internal working space
      integer i,j,count
      integer dispcount
      integer iscan
      integer isave
      integer nscan
      integer since
      integer skipcount
      integer sprint
      real*8 dnrm
      real*8 tmp1

      real runif

c++++ Working space slice sampling
      integer evali
      real*8 rexpo,re,uwork
      real*8 logy,xx0,xx1,llim,rlim
      real*8 grlim,gllim,gxx0,gxx1

c+++++CPU time
      real*8 sec00,sec0,sec1,sec

c++++++++++++++++++++++++++
c     initialize variables
c++++++++++++++++++++++++++

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

c++++ compute weights and atoms from SBA

      call discrsba(maxn,maxr,maxc,sba,fw,theta,w)

      do iscan=1,nscan

c++++++++++++++++++++++++++++++++++
c+++++++ clustering structure   +++
c++++++++++++++++++++++++++++++++++

         call rsbauclus(nrec,y,maxn,maxc,w,theta,sigma2,prob,
     &                  ccluster,cstrt,ss)

c         call intpr("ss",-1,ss,nrec)

c++++++++++++++++++++++++++++++++++
c+++++++ barycenter parameters  +++
c++++++++++++++++++++++++++++++++++

         call rsbauba(nrec,y,
     &                maxn,maxr,maxc,
     &                mu0,mu1,sigma20,sigma21,
     &                ccluster,cstrt,sigma2,
     &                fw,possv,
     &                sba,theta,w)

c         call dblepr("theta",-1,theta,maxc)
c         call dblepr("w",-1,w,maxc)

c++++++++++++++++++++++++++++++++++
c+++++++ kernel variances       +++
c++++++++++++++++++++++++++++++++++

         call rsbaukvar(nrec,y,maxn,maxc,ccluster,cstrt,
     &                  ss,theta,
     &                  tau1,tau2,
     &                  sigma2)

c         call dblepr("sigma2",-1,sigma2,maxc)

c++++++++++++++++++++++++++++++++++
c+++++++ gamma parameter        +++
c++++++++++++++++++++++++++++++++++

         call rsbaugamma(maxc,sigma2,
     &                   tau1,taus1,taus2,
     &                   tau2)

c         call dblepr("tau2",-1,tau2,1)


c++++++++++++++++++++++++++++++++++
c+++++++ SS Update mu1          +++
c++++++++++++++++++++++++++++++++++
         if(s1.gt.0.d0)then

            evali=1
            xx0=mu1
            re=rexpo(1.d0)
            logy=-re

            call lposth1msba(maxn,maxr,xx0,
     &                       sigma21,sba,m1,
     &                       s1,tmp1)
            logy=logy+tmp1
            uwork=dble(runif())*0.25d0  
            llim=xx0-uwork
            rlim=llim+0.25d0

            evali=evali+1
            call lposth1msba(maxn,maxr,llim,
     &                       sigma21,sba,m1,
     &                       s1,gllim)

            evali=evali+1
            call lposth1msba(maxn,maxr,rlim,
     &                       sigma21,sba,m1,
     &                       s1,grlim)

            do while(gllim.gt.logy)
               llim=llim-0.25d0
               call lposth1msba(maxn,maxr,llim,
     &                          sigma21,sba,m1,
     &                          s1,gllim)
            end do 

            do while(grlim.gt.logy)
               rlim=rlim+0.25d0
               call lposth1msba(maxn,maxr,rlim,
     &                          sigma21,sba,m1,
     &                          s1,grlim)

            end do 

            xx1=llim+(rlim-llim)*dble(runif())
            evali=evali+1
            call lposth1msba(maxn,maxr,xx1,
     &                       sigma21,sba,m1,
     &                       s1,gxx1)


            do while(gxx1.lt.logy)
               if(xx1.gt.xx0)rlim=xx1
               if(xx1.le.xx0)llim=xx1

               evali=evali+1
               xx1=llim+(rlim-llim)*dble(runif())
               call lposth1msba(maxn,maxr,xx1,
     &                          sigma21,sba,m1,
     &                          s1,gxx1)
            end do

            mu1=xx1

c            call dblepr("mu1",-1,mu1,1)

         end if


c++++++++++++++++++++++++++++++++++
c+++++++ SS Update sd1          +++
c++++++++++++++++++++++++++++++++++
         if(tauh1.gt.0.d0)then

            evali=1
            xx0=sigma21
            re=rexpo(1.d0)
            logy=-re
            call lposth1ssba(maxn,maxr,mu1,
     &                       sigma21,sba,tauh1,
     &                       tauh2,tmp1)
            logy=logy+tmp1
            uwork=dble(runif())*0.25d0  
            llim=xx0-uwork
            rlim=llim+0.25d0

            if(llim.lt.0.01d0)llim=0.01d0 

            evali=evali+1
            tmp1=llim
            call lposth1ssba(maxn,maxr,mu1,
     &                       tmp1,sba,tauh1,
     &                       tauh2,gllim)

            evali=evali+1
            tmp1=rlim
            call lposth1ssba(maxn,maxr,mu1,
     &                       tmp1,sba,tauh1,
     &                       tauh2,grlim)

            do while(gllim.gt.logy)
               llim=llim-0.25d0
               if(llim.lt.0.01d0)then
                  llim=0.01d0
                  gllim=logy-1.d0
                 else
                  evali=evali+1
                  tmp1=llim
                  call lposth1ssba(maxn,maxr,
     &                           mu1,tmp1,sba,
     &                           tauh1,tauh2,
     &                           gllim)
               end if
            end do 

            do while(grlim.gt.logy)
               rlim=rlim+0.25d0
               tmp1=rlim
               call lposth1ssba(maxn,maxr,mu1,
     &                          tmp1,sba,tauh1,
     &                          tauh2,grlim)
            end do 

            xx1=llim+(rlim-llim)*dble(runif())
            evali=evali+1
            tmp1=xx1
            call lposth1ssba(maxn,maxr,mu1,
     &                       tmp1,sba,tauh1,
     &                       tauh2,gxx1)


            do while(gxx1.lt.logy)
               if(xx1.gt.xx0)rlim=xx1
               if(xx1.lt.xx0)llim=xx1
               if(llim.lt.0.01d0)llim=0.01d0 

               evali=evali+1
               xx1=llim+(rlim-llim)*dble(runif())
               tmp1=xx1
               call lposth1ssba(maxn,maxr,mu1,
     &                          tmp1,sba,tauh1,
     &                          tauh2,gxx1)
            end do

            sigma21=xx1

c            call dblepr("sigma21",-1,sigma21,1)

         end if

c+++++++++++++++++++++++++++++++++++
c+++++++ saving samples          +++
c+++++++++++++++++++++++++++++++++++

         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1

               count=0
c+++++++++++++ sba parameters
               do i=1,maxr
                  count=count+1
                  sbasave(isave,count)=sba(maxn,i)
               end do          

c+++++++++++++ kernel variances
               do i=1,maxc
                  count=count+1
                  sbasave(isave,count)=sigma2(i)
               end do          

c+++++++++++++ gamma parameter
               thetasave(isave,1)=tau2

c+++++++++++++ H1 parameters
               thetasave(isave,2)=mu1
               thetasave(isave,3)=sigma21

c+++++++++++++ density
               do i=1,ngrid
                  tmp1=0.d0
                  do j=1,maxc
                     tmp1=tmp1+w(j)*dnrm(grid(i),theta(j),
     &                                   sqrt(sigma2(j)),0)
                  end do
                  densave(isave,i)=tmp1
               end do
     
c+++++++++++++ cpo's
               do i=1,nrec
                  tmp1=0.d0
                  do j=1,maxc
                     tmp1=tmp1+w(j)*dnrm(y(i),theta(j),
     &                                   sqrt(sigma2(j)),0)
                  end do
                  cpo(i,1)=cpo(i,1)+1.0d0/tmp1  
                  cpo(i,2)=cpo(i,2)+tmp1                   
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

      do i=1,nrec
         call rchkusr()
         cpo(i,1)=dble(nsave)/cpo(i,1)
         cpo(i,2)=cpo(i,2)/dble(nsave)
      end do

      return
      end
       

