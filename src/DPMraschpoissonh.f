
c=======================================================================                      
      subroutine dpmraschph(datastrm,imiss,nmissing,nsubject,p,y,
     &                      roffset,ngrid,grid,
     &                      maxn,a0b0,m0,s0,prec,sb,tau,
     &                      mcmc,nsave,
     &                      acrate,cpo,cpov,
     &                      randsave,thetasave,densave,
     &                      cdfsave,
     &                      alpha,b,beta,ss,ncluster,mub,sigmab,tauk2,
     &                      wdp,vdp,muclus,sigmaclus,
     &                      betac,workvp,workmhp,xtx,xty,iflagp,
     &                      cstrt,ccluster,prob,workv,
     &                      seed)
c=======================================================================                      
c
c     Copyright: Alejandro Jara, 2007-2010.
c
c     Version 1.0:
c
c     Last modification: 25-09-2009.
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
      integer imiss,nmissing,nsubject,p
      integer datastrm(nmissing,2),y(nsubject,p)
      double precision roffset(nsubject,p)

c+++++Density estimation
      integer ngrid
      double precision grid(ngrid)

c+++++Prior 
      integer maxn
      double precision aa0,ab0,a0b0(2),sb(p-1,2),prec(p-1,p-1)
      double precision m0,s0
      double precision tau(5),tauk1,taub1,taub2,taus1,taus2

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Output
      double precision acrate
      double precision cpo(nsubject,p)
      double precision cpov(nsubject)
      double precision randsave(nsave,nsubject+1)
      double precision thetasave(nsave,p+6)
      double precision densave(nsave,ngrid)
      double precision cdfsave(nsave,ngrid)

c+++++Current values of the parameters
      integer ncluster,ss(nsubject)     
      double precision alpha,beta(p-1),b(nsubject)
      double precision muclus(maxn),sigmaclus(maxn)
      double precision mub,sigmab
      double precision tauk2
      double precision wdp(maxn),vdp(maxn)

c+++++Working space - General
      integer i,ii,j
      integer sprint
      integer ns,ok
      double precision ztz,zty
      double precision dpoiss,dnrm,cdfnorm
      double precision tmp1,tmp2
      double precision muwork,sigmawork
      double precision mureal,sigma2real

c+++++Working space - RNG
      integer evali,seed(2),seed1,seed2
      integer rpois
      double precision rnorm,rgamma
      real runif

c+++++Working space - MCMC
      integer iscan,isave,nscan
      integer skipcount,dispcount

c+++++Working space - Difficulty parameters
      integer iflagp(p-1)
      double precision betac(p-1)
      double precision xtx(p-1,p-1)
      double precision xty(p-1)
      double precision workmhp((p-1)*p/2)
      double precision workvp(p-1)

c+++++DPM
      integer cstrt(maxn,nsubject)
      integer ccluster(maxn)
      double precision prob(maxn)
      double precision workv(maxn+1)

c+++++Working space - GLM part
      integer yij
      double precision eta,offset,gprime,ytilde,mean
       
c+++++Working space - MH 
      double precision logcgkn,logcgko,logliko,loglikn,ratio
      double precision logpriorn,logprioro 

c+++++Working space slice sampling
      double precision rexpo,re,uwork
      double precision logy,xx0,xx1,llim,rlim
      double precision grlim,gllim,gxx1

c+++++CPU time
      double precision sec00,sec0,sec1,sec

c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      aa0=a0b0(1)
      ab0=a0b0(2)
      
      tauk1=tau(1)
      taub1=tau(2)
      taub2=tau(3)
      taus1=tau(4)
      taus2=tau(5)

c++++ set random number generator
      seed1=seed(1)
      seed2=seed(2)
      
      call setall(seed1,seed2)
     
c++++ set configurations
      do i=1,nsubject
         ccluster(ss(i))=ccluster(ss(i))+1
         cstrt(ss(i),ccluster(ss(i)))=i
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

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            do ii=1,nmissing
               i=datastrm(ii,1)
               j=datastrm(ii,2)
               if(j.eq.1)then
                 eta=b(i)+roffset(i,j)
                else
                 eta=b(i)-beta(j-1)+roffset(i,j)
               end if
               
               mean=exp(eta)
               evali=rpois(mean)
               y(i,j)=evali
            end do
         end if

c++++++++++++++++++++++++++++++++++
c+++++++ difficulty parameters
c++++++++++++++++++++++++++++++++++

c+++++++ generating the candidate

         do i=1,p-1
            do j=1,p-1
               xtx(i,j)=prec(i,j)
            end do
            xty(i)=sb(i,1)
         end do

         logliko=0.d0         
         
         do i=1,nsubject

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            do j=1,p-1
               yij=y(i,j+1)            
               eta=b(i)-beta(j)+roffset(i,j+1) 
               offset=b(i)+roffset(i,j+1) 
               mean=exp(eta)
               logliko=logliko+dpoiss(dble(yij),mean,1)
               gprime=exp(-eta)

               ytilde=eta+(dble(yij)-mean)*gprime-offset
               xtx(j,j)=xtx(j,j)+1.d0/gprime
               xty(j)=xty(j)-ytilde/gprime

            end do
         end do

         call inverse(xtx,p-1,iflagp)      
         do i=1,p-1
            tmp1=0.d0
            do j=1,p-1
               tmp1=tmp1+xtx(i,j)*xty(j) 
            end do
            workvp(i)=tmp1
         end do

         call rmvnorm(p-1,workvp,xtx,workmhp,xty,betac)

c+++++++ evaluating the candidate generating kernel

         call dmvnd(p-1,betac,workvp,xtx,logcgko,iflagp)
  
c+++++++ prior ratio

         logprioro=0.d0
         logpriorn=0.d0
         
         do i=1,p-1
            do j=1,p-1
               logpriorn=logpriorn+(betac(i)-sb(i,2))* 
     &                    prec(i,j)      *
     &                   (betac(j)-sb(j,2))

               logprioro=logprioro+(beta(i) -sb(i,2))* 
     &                    prec(i,j)      *
     &                   (beta(j) -sb(j,2))
            end do
         end do
         
         logpriorn=-0.5d0*logpriorn
         logprioro=-0.5d0*logprioro
            
c+++++++ candidate generating kernel contribution

         do i=1,p-1
            do j=1,p-1
               xtx(i,j)=prec(i,j)
            end do
            xty(i)=sb(i,1)
         end do

         loglikn=0.d0         
         
         do i=1,nsubject

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            do j=1,p-1
               yij=y(i,j+1)            
               eta=b(i)-betac(j)+roffset(i,j+1) 
               offset=b(i)+roffset(i,j+1) 
               mean=exp(eta)
               loglikn=loglikn+dpoiss(dble(yij),mean,1)

               gprime=exp(-eta)

               ytilde=eta+(dble(yij)-mean)*gprime-offset
               xtx(j,j)=xtx(j,j)+1.d0/gprime
               xty(j)=xty(j)-ytilde/gprime
            end do
         end do

         call inverse(xtx,p-1,iflagp)      
         do i=1,p-1
            tmp1=0.d0
            do j=1,p-1
               tmp1=tmp1+xtx(i,j)*xty(j) 
            end do
            workvp(i)=tmp1
         end do

c+++++++ evaluating the candidate generating kernel
            
         call dmvnd(p-1,beta,workvp,xtx,logcgkn,iflagp)

c+++++++ mh step
           
         ratio=loglikn-logliko+logcgkn-logcgko+
     &         logpriorn-logprioro

         if(log(runif()).lt.ratio)then
            acrate=acrate+1.d0
            do i=1,p-1
               beta(i)=betac(i) 
            end do
         end if

c         call dblepr("loglikn",-1,loglikn,1)
c         call dblepr("logliko",-1,logliko,1)
c         call dblepr("logcgkn",-1,logcgkn,1)
c         call dblepr("logcgko",-1,logcgko,1)
c         call dblepr("logpriorn",-1,logpriorn,1)
c         call dblepr("logprioro",-1,logprioro,1)
c         call dblepr("lratio",-1,ratio,1)
c         call dblepr("beta",-1,beta,p-1)

c++++++++++++++++++++++++++++++++++         
c+++++++ random effects 
c++++++++++++++++++++++++++++++++++

         do i=1,nsubject

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            muwork=muclus(ss(i))
            sigmawork=sigmaclus(1)
 
            evali=1
            xx0=b(i)

            call targetreslicedpm(i,nsubject,p,y,roffset,
     &                            muwork,sigmawork,xx0,beta,tmp1)
            re=rexpo(1.d0)
            logy=tmp1-re

            uwork=dble(runif())  
            llim=xx0-uwork
            rlim=llim+1.d0

            evali=evali+1
            call targetreslicedpm(i,nsubject,p,y,roffset,
     &                            muwork,sigmawork,llim,beta,gllim)

            evali=evali+1
            call targetreslicedpm(i,nsubject,p,y,roffset,
     &                            muwork,sigmawork,rlim,beta,grlim)

            do while(gllim.gt.logy)
               llim=llim-1.d0

               evali=evali+1
               call targetreslicedpm(i,nsubject,p,y,roffset,
     &                               muwork,sigmawork,llim,beta,gllim)
            end do 

            do while(grlim.gt.logy)
               rlim=rlim+1.d0

               evali=evali+1
               call targetreslicedpm(i,nsubject,p,y,roffset,
     &                               muwork,sigmawork,rlim,beta,grlim)
            end do 

            xx1=llim+(rlim-llim)*dble(runif())

            evali=evali+1
            call targetreslicedpm(i,nsubject,p,y,roffset,
     &                            muwork,sigmawork,xx1,beta,gxx1)

            do while(gxx1.lt.logy)
               if(xx1.gt.xx0)rlim=xx1
               if(xx1.lt.xx0)llim=xx1

               xx1=llim+(rlim-llim)*dble(runif())

               evali=evali+1
               call targetreslicedpm(i,nsubject,p,y,roffset,
     &                               muwork,sigmawork,xx1,beta,gxx1)
            end do

            b(i)=xx1

         end do

c         call dblepr("b",-1,b,nsubject)

c+++++++++++++++++++++++++++++++++++++         
c+++++++ clustering structure 
c+++++++++++++++++++++++++++++++++++++

         do i=1,nsubject

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
               ccluster(ss(i))=ccluster(ss(i))-1 
               ncluster=ncluster-1
            end if 
         
            do j=1,maxn
               muwork=muclus(j)                
               sigmawork=sigmaclus(1)
               tmp2=dnrm(b(i),muwork,sqrt(sigmawork),0)
               prob(j)=wdp(j)*tmp2
            end do
            call simdisc(prob,maxn,maxn,evali)

            ss(i)=evali
            ccluster(evali)=ccluster(evali)+1
            cstrt(evali,ccluster(evali))=i
            if(ccluster(evali).eq.1)ncluster=ncluster+1
         end do

c         call intpr("ss",-1,ss,nsubject)

c+++++++++++++++++++++++++++++++++++++         
c+++++++ DP weights 
c+++++++++++++++++++++++++++++++++++++
         call dpweightsimbl(maxn,ccluster,alpha,workv,wdp,vdp)

c         call dblepr("wdp",-1,wdp,maxn)
     

c+++++++++++++++++++++++++++++++++++
c+++++++ kernel means            +++
c+++++++++++++++++++++++++++++++++++

         do i=1,maxn

            ztz=1.d0/sigmab
            zty=mub/sigmab

            ns=ccluster(i)
            if(ns.gt.0)then 
               ztz=ztz+dble(ns)/sigmaclus(1)
               do j=1,ns
c++++++++++++++++ check if the user has requested an interrupt
                  call rchkusr()

                  ii=cstrt(i,j)
                  zty=zty+b(ii)/sigmaclus(1)
               end do
            end if  
 
            tmp1=zty/ztz
            tmp2=1.d0/ztz

            muclus(i)=rnorm(tmp1,sqrt(tmp2))

         end do

c         call dblepr("mu",-1,muclus,maxn)

c+++++++++++++++++++++++++++++++++++
c+++++++ kernel variance         +++
c+++++++++++++++++++++++++++++++++++

         do i=1,nsubject

c++++++++++ check if the user has requested an interrupt
            call rchkusr()
            tmp1=b(i)-muclus(ss(i))
                 tmp2=tmp2+tmp1*tmp1
         end do
         sigmaclus(1)=1.d0/rgamma(0.5d0*(tauk1+dble(nsubject)),
     &                0.5d0*(tauk2+tmp2))

c         call dblepr("sigma",-1,sigmaclus,1)

c+++++++++++++++++++++++++++++++++++
c+++++++ precision parameter     +++
c+++++++++++++++++++++++++++++++++++

         if(a0b0(1).gt.0.d0)then

            call samalph(alpha,a0b0(1),a0b0(2),ncluster,nsubject)

c            tmp1=log(wdp(maxn))
c            alpha=rgamma(dble(maxn)+a0b0(1)-1.d0,a0b0(2)-tmp1)
c            call dblepr("alpha",-1,alpha,1)
         end if 


c++++++++++++++++++++++++++++++++++++++         
c+++++++ baseline mean
c++++++++++++++++++++++++++++++++++++++

         if(s0.gt.0.d0)then
            ztz=1.d0/s0
            zty=m0/s0

            do i=1,maxn
c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()
               ztz=ztz+1.d0
               zty=zty+muclus(i)/sigmab
            end do

            tmp1=zty/ztz
            tmp2=1.d0/ztz

            mub=rnorm(tmp1,sqrt(tmp2))
         end if


c++++++++++++++++++++++++++++++++++++++         
c+++++++ baseline variance
c++++++++++++++++++++++++++++++++++++++

         if(taub1.gt.0.d0)then

            tmp1=0.d0
            do i=1,maxn
               tmp2=muclus(i)-mub
               tmp1=tmp1+tmp2*tmp2
            end do
 
            sigmab=1.d0/
     &           rgamma(0.5d0*(dble(maxn)+taub1),0.5d0*(tmp1+taub2))
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

c+++++++++++++ computing mean and variance

               mureal=0.d0
               sigma2real=0.d0

               do ii=1,maxn
                  sigmawork=sigmaclus(1)
                  muwork=muclus(ii)
                  mureal=mureal+wdp(ii)*muwork
                  sigma2real=sigma2real+wdp(ii)*
     &                (muwork*muwork+sigmawork)
               end do
               sigma2real=sigma2real-mureal*mureal

c+++++++++++++ real mean and variance
               thetasave(isave,p)=mureal
               thetasave(isave,p+1)=sigma2real

c+++++++++++++ precision parameter
               thetasave(isave,p+2)=ncluster
               thetasave(isave,p+3)=alpha

c+++++++++++++ baseline parameters
               thetasave(isave,p+4)=mub
               thetasave(isave,p+5)=sigmab
               thetasave(isave,p+6)=tauk2

c+++++++++++++ random effects
               do i=1,nsubject
                  randsave(isave,i)=b(i)
               end do

c+++++++++++++ predictive

               call simdisc(wdp,maxn,maxn,evali)
               muwork=muclus(evali)
               sigmawork=sigmaclus(1)
  
               randsave(isave,nsubject+1)=rnorm(muwork,
     &                                    sqrt(sigmawork))

c+++++++++++++ density and cdf

               do i=1,ngrid
                  ytilde=grid(i)
                  tmp1=0.d0
                  tmp2=0.d0
                  do j=1,maxn
                     muwork=muclus(j)
                     sigmawork=sigmaclus(1)
                     tmp1=tmp1+wdp(j)*dnrm(ytilde,muwork,
     &                                     sqrt(sigmawork),0) 
                     tmp2=tmp2+wdp(j)*cdfnorm(ytilde,muwork,
     &                                        sqrt(sigmawork),1,0) 
                  end do
                  densave(isave,i)=tmp1
                  cdfsave(isave,i)=tmp2
               end do

c+++++++++++++ cpo

               do i=1,nsubject
                  tmp2=0.d0
                  do j=1,p
                     yij=y(i,j)
                     if(j.eq.1)then
                       eta=b(i)+roffset(i,j)
                      else
                       eta=b(i)-beta(j-1)+roffset(i,j)
                     end if  
                     mean=exp(eta)
                     tmp1=dpoiss(dble(yij),mean,1)
                     cpo(i,j)=cpo(i,j)+1.0d0/exp(tmp1)   
                     tmp2=tmp2+tmp1
                  end do
                  cpov(i)=cpov(i)+1.d0/exp(tmp2)
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
      
      acrate=acrate/dble(nscan)    
      
      do i=1,nsubject
         do j=1,p
            cpo(i,j)=dble(nsave)/cpo(i,j)
         end do 
         cpov(i)=dble(nsave)/cpov(i)  
      end do

      return
      end
