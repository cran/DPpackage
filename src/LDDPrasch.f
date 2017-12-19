c=======================================================================                      
      subroutine lddprasch(ngrid,npred,nsubject,p,q,y,roffset,
     &                     grid,z,zpred,     
     &                     murand,a0b0,b0,prec1,sb,mu0,prec2,smu,tinv,   
     &                     acrate,cpo,
     &                     denspm,randsave,thetasave,densave,          
     &                     ncluster,ss,beta,b,
     &                     alphaclus,sigmaclus,
     &                     alpha,mu,sigma,tau2,        
     &                     mcmc,nsave,
     &                     iflagp,betac,xtx,xty,workmhp1,workvp1,       
     &                     cstrt,ccluster,iflagq,
     &                     alphawork,densw,
     &                     prob,quadf,sigmainv,            
     &                     workmhq1,workmhq2,                            
     &                     workvq1,workvq2,      
     &                     ztz,zty,                                      
     &                     seed)                                          
c=======================================================================                  
c   # 58 arguments
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
      integer ngrid,npred,nsubject,p,q
      integer y(nsubject,p)
      double precision roffset(nsubject,p)
      double precision grid(ngrid),z(nsubject,q),zpred(npred,q)

c+++++Prior
      integer murand
      double precision a0b0(6)
      double precision aa0,ab0
      double precision b0(p-1)
      double precision prec1(p-1,p-1)
      double precision sb(p-1)
      double precision mu0(q)
      double precision prec2(q,q)
      double precision smu(q)
      double precision tau1
      double precision taus1,taus2
      double precision nu
      double precision tinv(q,q) 

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Output
      double precision acrate(2)
      double precision cpo(nsubject,p)
      double precision denspm(npred,ngrid)
      double precision randsave(nsave,nsubject+npred)
      double precision thetasave(nsave,p+q+(q*(q+1)/2)+2)
      double precision densave(nsave,npred*ngrid)

c+++++Current values of the parameters
      integer ncluster,ss(nsubject) 
      double precision beta(p-1)
      double precision b(nsubject)
      double precision alphaclus(nsubject+100,q)
      double precision sigmaclus(nsubject+100)
      double precision alpha
      double precision mu(q)
      double precision sigma(q,q)
      double precision tau2

c+++++External Working space - Difficulty parameters
      integer iflagp(p-1)
      double precision betac(p-1)
      double precision xtx(p-1,p-1),xty(p-1)
      double precision workmhp1((p-1)*p/2)
      double precision workvp1(p-1)

c+++++External Working space - DDP part
      integer cstrt(nsubject,nsubject)
      integer ccluster(nsubject)
      integer iflagq(q)
      double precision alphawork(q)
      double precision densw(npred,ngrid)
      double precision prob(nsubject+100)
      double precision quadf(q,q)
      double precision sigmainv(q,q)
      double precision workmhq1(q*(q+1)/2)
      double precision workmhq2(q*(q+1)/2)
      double precision workvq1(q),workvq2(q)
      double precision ztz(q,q),zty(q)

c+++++External Working space - RNG
      integer seed(2),seed1,seed2

c+++++Internal Working space
      integer counter
      integer dispcount
      integer evali
      integer i,ii
      integer iscan,isave,isample
      integer j,k,l
      integer ns,nscan
      integer ok
      integer since
      integer skipcount
      integer sprint
      integer yij
      double precision acrate2
      double precision dbin,dnrm
      double precision eta,gprime
      double precision logcgkn,logcgko
      double precision loglikn,logliko
      double precision logpriorn,logprioro
      double precision mean
      double precision muwork
      double precision offset
      double precision ratio,rgamma,rnorm,runif
      double precision thetac
      double precision tmp1,tmp2,tmp3
      double precision ytilde
      double precision wtw,wtwinv,wty
      double precision sigmawork

c+++++CPU time
      double precision sec00,sec0,sec1,sec

c++++ DP (functional parameter)
      double precision eps,rbeta,weight
      parameter(eps=0.01)

c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)

      aa0=a0b0(1)
      ab0=a0b0(2)
      tau1=a0b0(3)
      taus1=a0b0(4)
      taus2=a0b0(5)
      nu=a0b0(6)

c++++ set random number generator
      seed1=seed(1)
      seed2=seed(2)
      
      call setall(seed1,seed2)

c++++ cluster structure
      do i=1,nsubject
         ccluster(ss(i))=ccluster(ss(i))+1
         cstrt(ss(i),ccluster(ss(i)))=i
      end do
     
c++++ inverse of sigma

      do i=1,q
         do j=1,q
            sigmainv(i,j)=sigma(i,j)
         end do
      end do
      call inverse(sigmainv,q,iflagq) 

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

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++++
c+++++++ difficulty parameters
c++++++++++++++++++++++++++++++++++

         do i=1,p-1
            do j=1,p-1
               xtx(i,j)=prec1(i,j)
            end do
            xty(i)=sb(i)
         end do
         
         logliko=0.d0
         
         do i=1,nsubject
            do j=1,p-1
               eta=b(i)-beta(j)+roffset(i,j+1) 
               offset=b(i)+roffset(i,j+1)
               mean=exp(eta)/(1.d0+dexp(eta))

               tmp1=mean*(1.d0-mean)
               gprime=1.0d0/(tmp1)           

               ytilde=eta+(dble(y(i,j+1))-mean)*gprime-offset

               xtx(j,j)=xtx(j,j)+tmp1
               xty(j)=xty(j)-ytilde*tmp1

               logliko=logliko+dbin(dble(y(i,j+1)),1.d0,mean,1)               
            end do
         end do

         call inverse(xtx,p-1,iflagp) 

         do i=1,p-1
            tmp1=0.d0
            do j=1,p-1
               tmp1=tmp1+xtx(i,j)*xty(j) 
            end do
            workvp1(i)=tmp1
         end do

         call rmvnorm(p-1,workvp1,xtx,workmhp1,xty,betac)

         call dmvnd(p-1,betac,workvp1,xtx,logcgko,iflagp)

c+++++++ prior ratio

         tmp1=0.d0
         tmp2=0.d0
         
         do i=1,p-1
            do j=1,p-1
               tmp1=tmp1+(betac(i)-b0(i))* 
     &                    prec1(i,j)      *
     &                   (betac(j)-b0(j))

               tmp2=tmp2+(beta(i) -b0(i))* 
     &                    prec1(i,j)      *
     &                   (beta(j) -b0(j))
            end do
         end do

         logpriorn=-0.5d0*tmp1
         logprioro=-0.5d0*tmp2

            
c+++++++ candidate generating kernel contribution

         do i=1,p-1
            do j=1,p-1
               xtx(i,j)=prec1(i,j)
            end do
            xty(i)=sb(i)
         end do

         loglikn=0.d0

         do i=1,nsubject
            do j=1,p-1
               eta=b(i)-betac(j)+roffset(i,j+1) 
               offset=b(i)+roffset(i,j+1)
               mean=exp(eta)/(1.d0+dexp(eta))

               tmp1=mean*(1.d0-mean)
               gprime=1.0d0/(tmp1)       
    
               ytilde=eta+(dble(y(i,j+1))-mean)*gprime-offset

               loglikn=loglikn+dbin(dble(y(i,j+1)),1.d0,mean,1)

               xtx(j,j)=xtx(j,j)+1.d0*tmp1
               xty(j)=xty(j)-ytilde*tmp1
            end do
         end do

         call inverse(xtx,p-1,iflagp)      
         do i=1,p-1
            tmp1=0.d0
            do j=1,p-1
               tmp1=tmp1+xtx(i,j)*xty(j) 
            end do
            workvp1(i)=tmp1
         end do
            
         call dmvnd(p-1,beta,workvp1,xtx,logcgkn,iflagp)
 
c+++++++ mh step

         ratio=loglikn-logliko+logcgkn-logcgko+
     &         logpriorn-logprioro

         if(log(runif()).lt.ratio)then
            acrate(1)=acrate(1)+1.d0
            do i=1,p-1
               beta(i)=betac(i) 
            end do
         end if

c         call dblepr("beta",-1,beta,p-1)

c++++++++++++++++++++++++++++++++++         
c+++++++ random effects 
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         acrate2=0.d0

         do i=1,nsubject

            tmp1=0.d0
            do j=1,q
               tmp1=tmp1+z(i,j)*alphaclus(ss(i),j)
            end do
            sigmawork=sigmaclus(ss(i))

            logprioro=dnrm(b(i),tmp1,sqrt(sigmawork),1)

            wtw=1.d0/sigmawork
            wty=tmp1/sigmawork

            logliko=0.d0
            
            do j=1,p
               if(j.eq.1)then
                  eta=b(i)+roffset(i,j)
                  offset=roffset(i,j)
                 else
                  eta=b(i)-beta(j-1)+roffset(i,j)
                  offset=-beta(j-1)+roffset(i,j)
               end if
               
               yij=y(i,j)

               mean=exp(eta)/(1.d0+dexp(eta))
               tmp1=mean*(1.d0-mean)
               gprime=1.0d0/(tmp1)           

               ytilde=eta+(dble(yij)-mean)*gprime-offset    

               logliko=logliko+dbin(dble(yij),1.d0,mean,1)
               
               wtw=wtw+tmp1
               wty=wty+ytilde*tmp1
            end do

            wtwinv=1.d0/wtw
            tmp1=wtwinv*wty
 
            thetac=rnorm(tmp1,sqrt(wtwinv))

            logcgko=dnrm(thetac,tmp1,sqrt(wtwinv),1)

c++++++++++ candidate generating kernel contribution

            tmp1=0.d0
            do j=1,q
               tmp1=tmp1+z(i,j)*alphaclus(ss(i),j)
            end do
            sigmawork=sigmaclus(ss(i))

            logpriorn=dnrm(thetac,tmp1,sqrt(sigmawork),1)

            wtw=1.d0/sigmawork
            wty=tmp1/sigmawork

            loglikn=0.d0
            
            do j=1,p
               if(j.eq.1)then
                  eta=thetac+roffset(i,j)
                  offset=roffset(i,j)
                 else
                  eta=thetac-beta(j-1)+roffset(i,j)
                  offset=-beta(j-1)+roffset(i,j)
               end if
               
               yij=y(i,j)
               
               mean=exp(eta)/(1.d0+dexp(eta))
               tmp1=mean*(1.d0-mean)
               gprime=1.0d0/(tmp1)           
               
               ytilde=eta+(dble(yij)-mean)*gprime-offset    

               loglikn=loglikn+dbin(dble(yij),1.d0,mean,1)
               
               wtw=wtw+tmp1
               wty=wty+ytilde*tmp1
            end do

            wtwinv=1.d0/wtw
            tmp1=wtwinv*wty
 
            logcgkn=dnrm(b(i),tmp1,sqrt(wtwinv),1)

c++++++++++ mh step
            ratio=loglikn-logliko+logcgkn-logcgko+
     &            logpriorn-logprioro

            if(log(runif()).lt.ratio)then
               acrate2=acrate2+1.d0
               b(i)=thetac
            end if
         end do

         acrate(2)=acrate(2)+acrate2/dble(nsubject)

c         call dblepr("b",-1,b,nsubject)

c++++++++++++++++++++++++++++++++++
c+++++++ clustering structure   +++
c++++++++++++++++++++++++++++++++++

         do i=1,nsubject

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
                  call relabellddp(i,since,nsubject,q,ncluster,ccluster,
     &                             ss,cstrt,alphaclus,sigmaclus,
     &                             alphawork)
               end if
               ccluster(ncluster)=ccluster(ncluster)-1 
               ncluster=ncluster-1
            end if 
         
            do j=1,ncluster

               tmp1=0.0
               do k=1,q
                  tmp1=tmp1+z(i,k)*alphaclus(j,k)
               end do                
               sigmawork=sigmaclus(j)

               tmp2=dnrm(b(i),tmp1,sqrt(sigmawork),1)
                        
               prob(j)=dble(ccluster(j))*exp(tmp2)
            end do

            if(isample.eq.1)then
               sigmawork=1.d0/rgamma(0.5*tau1,0.5*tau2)
               sigmaclus(ncluster+1)=sigmawork
        
               call rmvnorm(q,mu,sigma,workmhq1,workvq1,alphawork) 
               do k=1,q
                  alphaclus(ncluster+1,k)=alphawork(k)
               end do               
            end if         
            sigmawork=sigmaclus(ncluster+1)

            tmp1=0.0
            do k=1,q
               tmp1=tmp1+z(i,k)*alphaclus(ncluster+1,k)
            end do               
            tmp2=dnrm(b(i),tmp1,sqrt(sigmawork),1)

            prob(ncluster+1)=alpha*exp(tmp2)

            call simdisc(prob,nsubject+100,ncluster+1,evali)

            ss(i)=evali
            ccluster(evali)=ccluster(evali)+1
            cstrt(evali,ccluster(evali))=i
                
            if(evali.gt.ncluster)then
               ncluster=ncluster+1
            end if
      
         end do
     
c         call intpr("ncluster",-1,ncluster,1)
c         call intpr("ss",-1,ss,nsubject)

c+++++++++++++++++++++++++++++++++++
c+++++++ regression coefficients +++
c+++++++++++++++++++++++++++++++++++

         do i=1,ncluster

            do k=1,q
               tmp1=0.d0
               do l=1,q
                  ztz(k,l)=sigmainv(k,l)
                  tmp1=tmp1+sigmainv(k,l)*mu(l)
               end do
               zty(k)=tmp1
            end do   

            ns=ccluster(i)
            do j=1,ns

c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

               ii=cstrt(i,j)

               do k=1,q
                  do l=1,q
                     ztz(k,l)=ztz(k,l)+z(ii,k)*z(ii,l)/sigmaclus(i)
                  end do
                  zty(k)=zty(k)+z(ii,k)*b(ii)/sigmaclus(i)
               end do
            end do

            call inverse(ztz,q,iflagq)

            do k=1,q 
               tmp1=0.d0
               do l=1,q
                  tmp1=tmp1+ztz(k,l)*zty(l)
               end do
               workvq1(k)=tmp1
            end do  
      
            call rmvnorm(q,workvq1,ztz,workmhq1,workvq2,alphawork)
  
            do j=1,q
               alphaclus(i,j)=alphawork(j)
            end do

c            call dblepr("alphawork",-1,alphawork,q)

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
               tmp1=0.d0
               do k=1,q
                  tmp1=tmp1+z(ii,k)*alphaclus(i,k)
               end do
               tmp1=b(ii)-tmp1
               tmp2=tmp2+tmp1*tmp1
            end do

            sigmaclus(i)=1.d0/rgamma(0.5d0*(tau1+dble(ns)),
     &                   0.5d0*(tau2+tmp2))
         end do

c         call dblepr("sigma",-1,sigmaclus,ncluster)

c++++++++++++++++++++++++++++++++++         
c+++++++ Precision parameter
c++++++++++++++++++++++++++++++++++
         if(aa0.gt.0.d0)then
            call samalph(alpha,aa0,ab0,ncluster,nsubject)
         end if 

c+++++++++++++++++++++++++++++++++++
c+++++++ baseline mean           +++
c+++++++++++++++++++++++++++++++++++

         if(murand.eq.1)then
            do i=1,q
               do j=1,q
                  ztz(i,j)=sigmainv(i,j)+dble(ncluster)*prec2(i,j)
               end do
               zty(i)=smu(i)
            end do
            call inverse(ztz,q,iflagq)

            do ii=1,ncluster      
               do i=1,q
                  alphawork(i)=alphaclus(ii,i)
               end do
               do i=1,q
                  tmp1=0.d0
                  do j=1,q
                     tmp1=tmp1+sigmainv(i,j)*alphawork(j)
                  end do 
                  zty(i)=zty(i)+tmp1
               end do
            end do

            do i=1,q
               tmp1=0.d0
               do j=1,q
                  tmp1=tmp1+ztz(i,j)*zty(j)
               end do
               workvq1(i)=tmp1
            end do

            call rmvnorm(q,workvq1,ztz,workmhq1,workvq2,mu)

c            call dblepr("mu",-1,mu,q)
         end if

c++++++++++++++++++++++++++++++++++++++
c+++++++ baseline covariance matrix +++
c++++++++++++++++++++++++++++++++++++++

c         call dblepr("nu",-1,nu,1)

         if(nu.gt.0)then
            do i=1,q
               do j=1,q
                  ztz(i,j)=0.d0
               end do
            end do
  
            do ii=1,ncluster 
               do i=1,q
                  alphawork(i)=alphaclus(ii,i)-mu(i)
               end do
 
               do i=1,q
                  do j=1,q
                     ztz(i,j)=ztz(i,j)+alphawork(i)*alphawork(j)
                  end do
               end do
            end do

            do i=1,q	
               do j=1,q
                  sigma(i,j)=ztz(i,j)+tinv(i,j)
               end do
            end do

            call riwishart(q,int(nu)+ncluster,sigma,sigmainv,
     &                     quadf,zty,workmhq1,workmhq2,iflagq)

c            call dblepr("sigma",-1,sigma,q*q)

         end if

c++++++++++++++++++++++++++++++++++++++
c+++++++ baseline gamma parameter   +++
c++++++++++++++++++++++++++++++++++++++

         if(taus1.gt.0.d0)then
            tmp1=0.d0
            do i=1,ncluster
               tmp1=tmp1+1.d0/sigmaclus(i)
            end do 

            tau2=rgamma(0.5d0*(dble(ncluster)*tau1+taus1),
     &                  0.5d0*(tmp1+taus2))   

c           call dblepr("tau2",-1,tau2,1)
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

c+++++++++++++ random effects variance
               thetasave(isave,p)=tau2

c+++++++++++++ baseline information
               
               do i=1,q
                  thetasave(isave,p+i)=mu(i)
               end do   
              
               counter=0
               do i=1,q
                  do j=i,q	
                     counter=counter+1
                     thetasave(isave,p+q+counter)=sigma(i,j)
                  end do   
               end do

c+++++++++++++ cluster information

               thetasave(isave,p+q+counter+1)=ncluster
               thetasave(isave,p+q+counter+2)=alpha

c+++++++++++++ random effects

               do i=1,nsubject
                  randsave(isave,i)=b(i)
               end do

c+++++++++++++ predictive information

               do i=1,npred

                  call simdisc(prob,nsubject+100,ncluster+1,evali)

                  if(evali.le.ncluster)then
                     sigmawork=sigmaclus(evali)
                     do k=1,q
                        alphawork(k)=alphaclus(evali,k)
                     end do
                  end if
                  if(evali.eq.ncluster+1)then 
                     sigmawork=1.0d0/rgamma(0.5d0*tau1,0.5d0*tau2)
                     sigmaclus(ncluster+1)=sigmawork
        
                     call rmvnorm(q,mu,sigma,workmhq1,workvq1,alphawork) 
                     do k=1,q
                        alphaclus(ncluster+1,k)=alphawork(k)
                     end do               
                  end if

                  tmp1=0.d0 
                  do k=1,q
                     tmp1=tmp1+zpred(i,k)*alphaclus(evali,k) 
                  end do
                  sigmawork=sigmaclus(evali)
                  thetac=rnorm(tmp1,sqrt(sigmawork))
                  randsave(isave,nsubject+i)=thetac
               end do
               
c+++++++++++++ Partially sampling the DP

               do i=1,ncluster
                   prob(i)=real(ccluster(i))/(alpha+real(nsubject))
               end do
               prob(ncluster+1)=alpha/(alpha+dble(nsubject))
               call simdisc(prob,nsubject+100,ncluster+1,evali)

               if(evali.le.ncluster)then
                  sigmawork=sigmaclus(evali)
                  do k=1,q
                     alphawork(k)=alphaclus(evali,k)
                  end do
               end if
               if(evali.eq.ncluster+1)then 
                  sigmawork=1.0d0/rgamma(0.5d0*tau1,0.5d0*tau2)
                  sigmaclus(ncluster+1)=sigmawork
        
                  call rmvnorm(q,mu,sigma,workmhq1,workvq1,alphawork) 
                  do k=1,q
                     alphaclus(ncluster+1,k)=alphawork(k)
                  end do               
               end if

               tmp1=rbeta(1.d0,alpha+dble(nsubject))
               tmp2=tmp1
               weight=(1.d0-tmp1)
    
               do i=1,npred  
                  call rchkusr()
     
                  muwork=0.0
                  do j=1,q
                     muwork=muwork+zpred(i,j)*alphawork(j)
                  end do
           
                  do j=1,ngrid  
                     densw(i,j)=tmp1*dnrm(grid(j),muwork,
     &                          sqrt(sigmawork),0)
                 
                  end do 
               end do

               do while((1.0-tmp2).gt.eps)
                  call rchkusr()

                  tmp3=rbeta(1.d0,alpha+real(nsubject))
                  tmp1=weight*tmp3
                  weight=weight*(1.0-tmp3)

                  call simdisc(prob,nsubject+100,ncluster+1,evali)

                  if(evali.le.ncluster)then
                     sigmawork=sigmaclus(evali)
                     do k=1,q
                        alphawork(k)=alphaclus(evali,k)
                     end do
                  end if

                  if(evali.eq.ncluster+1)then 
                     sigmawork=1.0d0/rgamma(0.5d0*tau1,0.5d0*tau2)
                     sigmaclus(ncluster+1)=sigmawork
        
                     call rmvnorm(q,mu,sigma,workmhq1,workvq1,alphawork) 
                     do k=1,q
                        alphaclus(ncluster+1,k)=alphawork(k)
                     end do               
                  end if

                  do i=1,npred  
                     call rchkusr()

                     muwork=0.0
                     do j=1,q
                        muwork=muwork+zpred(i,j)*alphawork(j)
                     end do

                     do j=1,ngrid  
                        densw(i,j)=densw(i,j)+tmp1*dnrm(grid(j),
     &                             muwork,sqrt(sigmawork),0)

                     end do 
                  end do

                  tmp2=tmp2+tmp1
               end do

               call simdisc(prob,nsubject+100,ncluster+1,evali)

               if(evali.le.ncluster)then
                  sigmawork=sigmaclus(evali)
                  do k=1,q	
                     alphawork(k)=alphaclus(evali,k)
                  end do
               end if
               if(evali.eq.ncluster+1)then 
                  sigmawork=1.0d0/rgamma(0.5d0*tau1,0.5d0*tau2)
                  sigmaclus(ncluster+1)=sigmawork
        
                  call rmvnorm(q,mu,sigma,workmhq1,workvq1,alphawork) 
                  do k=1,q
                     alphaclus(ncluster+1,k)=alphawork(k)
                  end do               
               end if

               tmp1=weight
               do i=1,npred  
                  call rchkusr()
     
                  muwork=0.0
                  do j=1,q
                     muwork=muwork+zpred(i,j)*alphawork(j)
                  end do
           
                  do j=1,ngrid  
                     densw(i,j)=densw(i,j)+tmp1*dnrm(grid(j),
     &                           muwork,sqrt(sigmawork),0)
                   end do 
               end do

               ii=0
               do i=1,npred
                  call rchkusr()
                  do j=1,ngrid
                     denspm(i,j)=denspm(i,j)+densw(i,j)
                     ii=ii+1
                     densave(isave,ii)=densw(i,j)
                  end do
               end do 

c+++++++++++++ cpo

               do i=1,nsubject
                  do j=1,p
                     yij=y(i,j)
                     if(j.eq.1)then
                       eta=b(i)+roffset(i,j)
                      else
                       eta=b(i)-beta(j-1)+roffset(i,j)
                     end if
                     mean=exp(eta)/(1.d0+dexp(eta))
                     tmp1=dbin(dble(yij),1.d0,mean,0)
                     cpo(i,j)=cpo(i,j)+1.d0/tmp1
                  end do
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
      
      do i=1,2
         acrate(i)=acrate(i)/dble(nscan)    
      end do   
      
      do i=1,nsubject
         do j=1,p
            cpo(i,j)=dble(nsave)/cpo(i,j)
         end do   
      end do

      do i=1,npred
         do j=1,ngrid
            denspm(i,j)=denspm(i,j)/dble(nsave)
         end do   
      end do

      return
      end

