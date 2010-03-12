c=======================================================================                      
      subroutine ptmdensityf(nrec,nvar,y,
     &                       ngrid,grid,
     &                       frstl,maxm,a0b0,
     &                       mcmc,nsave,tune1,tune2,
     &                       acrate,cpo,f,randsave,thetasave,
     &                       alpha,mu,sigma,ortho,
     &                       seed,iflag,kmat,kphi,kvec,
     &                       omat,propvr,sigmac,
     &                       sigmachol,sigmacholc,uinv,uinvc,
     &                       workm,workm2,workmh,workmh2,workv)
c======================================================================= 
c
c     Copyright: Alejandro Jara, 2008-2010.
c
c     Version 1.0:
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
      integer nrec,nvar
      real*8 y(nrec,nvar)

c+++++Prediction 
      integer ngrid
      real*8 grid(ngrid,2)

c+++++Prior information
      integer frstl
      integer maxm
      real*8 a0b0(2)

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay
      real*8 tune1,tune2

c+++++Stored output
      real*8 acrate(2)
      real*8 cpo(nrec)
      real*8 f(ngrid,ngrid)
      real*8 randsave(nsave,nvar)
      real*8 thetasave(nsave,nvar+nvar*(nvar+1)/2+1)

c+++++Current values of the parameters
      real*8 alpha
      real*8 mu(nvar)
      real*8 sigma(nvar,nvar)
      real*8 ortho(nvar,nvar)

c+++++Random number generator
      integer seed(2),seed1,seed2

c+++++External working space
      integer iflag(nvar)
      integer kmat(nrec,maxm)
      integer kphi(nvar)
      integer kvec(maxm)
      real*8 omat(nvar,nvar)
      real*8 propvr(nvar,nvar)
      real*8 sigmac(nvar,nvar) 
      real*8 sigmachol(nvar,nvar) 
      real*8 sigmacholc(nvar,nvar) 
      real*8 uinv(nvar,nvar)
      real*8 uinvc(nvar,nvar)
      real*8 workm(nvar,nvar) 
      real*8 workm2(nvar,nvar) 
      real*8 workmh(nvar*(nvar+1)/2)
      real*8 workmh2(nvar*(nvar+1)/2)
      real*8 workv(nvar)

c+++++Working space - CPU time
      real*8 sec00,sec0,sec1,sec

c+++++Working space - Distributions
      real*8 dnrm,dlnrm,rtlnorm 

c+++++Working space - General
      integer i,i1,j,j1,k,k1
      integer ihmssf
      integer pprn,sprint
      real*8 alphac
      real*8 ldet,ldetc
      real*8 tmp1      

c+++++Working space - MCMC scans
      integer dispcount,isave,iscan,nscan,skipcount

c+++++Working space - MH steps
      real*8 logcgkn,logcgko
      real*8 loglikn,logliko
      real*8 logpriorn,logprioro
      real*8 ratio

c+++++Working space - Random numbers
      real runif

c++++ Working space - slice sampling
      integer evali
      real*8 rexpo,re,uwork
      real*8 logy,xx0,xx1,llim,rlim
      real*8 grlim,gllim,gxx0,gxx1

c+++++Adaptive MH for sigma
      integer nadaptive,nu
      parameter(nadaptive=2000)
      integer adaptives,sigmaskip
      real*8 aratesigma

c+++++Adaptive MH for c
      integer adaptivec,cskip
      real*8 aratec

c++++ initialize variables
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      
c++++ set random number generator
      seed1=seed(1)
      seed2=seed(2)
      call setall(seed1,seed2)

c++++ parameters for adaptive MH for sigma

      adaptives=0
      aratesigma=0.d0
      sigmaskip=0
      if(tune1.lt.0.d0)then
         adaptives=1
         tune1=10.0
         nburn=nburn+nadaptive
      end if  

      adaptivec=0
      aratec=0.d0
      cskip=0
      if(a0b0(1).gt.0.d0.and.tune2.lt.0.d0)then
         adaptivec=1
         tune2=1.0
         nburn=nburn+nadaptive
      end if  

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*(nsave)
      
      call cpu_time(sec0)
      sec00=0.d0
  
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ First computations
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      do i=1,nvar
         do j=1,nvar   
            workm2(i,j)=ortho(i,j)
         end do
      end do
      call rhaar2(workm,workm2,nvar,omat)

      do i=1,nvar
         do j=1,nvar
            sigmac(i,j)=sigma(i,j)
         end do
      end do
      call inversedet(sigmac,nvar,iflag,ldet)
      call cholesky(nvar,sigma,workmh)
      do i=1,nvar
         do j=1,i
            sigmachol(i,j)=workmh(ihmssf(i,j,nvar))
         end do
      end do

      do i=1,nvar
         do j=1,nvar
            tmp1=0.d0
            do k=1,nvar 
               tmp1=tmp1+sigmachol(i,k)*omat(k,j)
            end do
            uinv(i,j)=tmp1
         end do
      end do
      call inverse(uinv,nvar,iflag)      


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ Scanning the posterior distribution
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      do iscan=1,nscan
  
c+++++++ check if the user has requested an interrupt
         call rchkusr()
 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Updating mu using Slice sampling                            +++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         
         do i=1,nvar
            evali=1
            xx0=mu(i)
            re=rexpo(1.d0)
            logy=-re

c++++++++++ evaluates the log-posterior

            call loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,uinv,
     &                       ldet,kphi,kmat,tmp1)
            logy=logy+tmp1
            uwork=dble(runif())*0.5  
            llim=xx0-uwork
            rlim=llim+0.5

            evali=evali+1
            mu(i)=llim
            call loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,uinv,
     &                       ldet,kphi,kmat,gllim)


            evali=evali+1
            mu(i)=rlim
            call loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,uinv,
     &                       ldet,kphi,kmat,grlim)

            do while(gllim.gt.logy)
               llim=llim-0.5
               evali=evali+1
               mu(i)=llim
               call loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,uinv,
     &                          ldet,kphi,kmat,gllim)
            end do 

            do while(grlim.gt.logy)
               rlim=rlim+0.5
               evali=evali+1
               mu(i)=rlim
               call loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,uinv,
     &                          ldet,kphi,kmat,grlim)
            end do 
 
            xx1=llim+(rlim-llim)*dble(runif())
            evali=evali+1
            mu(i)=xx1
            call loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,uinv,
     &                       ldet,kphi,kmat,gxx1)

            do while(gxx1.lt.logy)
               if(xx1.gt.xx0)rlim=xx1
               if(xx1.lt.xx0)llim=xx1
               xx1=llim+(rlim-llim)*dble(runif())
               evali=evali+1
               mu(i)=xx1
               call loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,uinv,
     &                          ldet,kphi,kmat,gxx1)
            end do
            mu(i)=xx1
         end do
         
         call dblepr("mu",-1,mu,nvar)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Updating sigma using a MH step                              +++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(adaptives.eq.1)then  
            sigmaskip = sigmaskip + 1
            if(sigmaskip.eq.100)then
               aratesigma=aratesigma/dble(100)

               if(iscan.le.nadaptive)then  
                  if(nvar.eq.1)then
                     if(aratesigma.lt.0.44)then
                        tune1=exp(log(tune1)+(0.44-aratesigma))
                      else
                        tune1=exp(log(tune1)-(aratesigma-0.44))
                     end if  
                    else
                     if(aratesigma.lt.0.234)then
                        tune1=exp(log(tune1)+(0.234-aratesigma))
                      else
                        tune1=exp(log(tune1)-(aratesigma-0.234))
                     end if  
                  end if  

                else 
                  if(nvar.eq.1)then
                     if(aratesigma.lt.0.44)then
                        tune1=exp(log(tune1)+
     &                    min(0.01,1.d0/sqrt(dble(iscan-nadaptive))))
                       else
                        tune1=exp(log(tune1)-
     &                    min(0.01,1.d0/sqrt(dble(iscan-nadaptive))))
                     end if  
                    else
                     if(aratesigma.lt.0.234)then
                        tune1=exp(log(tune1)+
     &                   min(0.01,1.d0/sqrt(dble(iscan-nadaptive))))
                       else
                        tune1=exp(log(tune1)-
     &                   min(0.01,1.d0/sqrt(dble(iscan-nadaptive))))
                     end if  
                  end if  
               end if
                  
               nu=(dble(nrec))*tune1
               if(nu.le.(nvar+1))tune1=dble(nvar+2)/dble(nrec)
                  
               sigmaskip=0
               aratesigma=0.d0
            end if
         end if


c+++++++ generating the candidate value

         nu=(dble(nrec))*tune1
         
         do i=1,nvar
            do j=1,nvar
               sigmac(i,j)=dble(nu-nvar-1)*sigma(i,j)
            end do
         end do

         call riwishart(nvar,nu,sigmac,workm2,workm,workv,
     &                  workmh,workmh2,iflag)


c+++++++ evaluating the candidate generating kernel

         do i=1,nvar
            do j=1,nvar
               propvr(i,j)=dble(nu-nvar-1)*sigma(i,j)
            end do
         end do

         call diwishart(nvar,nu,sigmac,propvr,workm,workm2,workv,
     &                  iflag,logcgko)        

         do i=1,nvar
            do j=1,nvar
               propvr(i,j)=dble(nu-nvar-1)*sigmac(i,j)
            end do
         end do

         call diwishart(nvar,nu,sigma,propvr,workm,workm2,workv,
     &                  iflag,logcgkn)        


c+++++++ evaluating likelihood contribution

         call loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,uinv,
     &                    ldet,kphi,kmat,logliko)

         do i=1,nvar
            do j=1,nvar
               workm(i,j)=sigmac(i,j)
            end do
         end do
         call inversedet(workm,nvar,iflag,ldetc)
         call cholesky(nvar,sigmac,workmh)
         do i=1,nvar
            do j=1,i
               sigmacholc(i,j)=workmh(ihmssf(i,j,nvar))
            end do
         end do

         do i=1,nvar
            do j=1,nvar
               tmp1=0.d0
               do k=1,nvar 
                  tmp1=tmp1+sigmacholc(i,k)*omat(k,j)
               end do
               uinvc(i,j)=tmp1
            end do
         end do
         call inverse(uinvc,nvar,iflag)      

         call loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,uinvc,
     &                    ldetc,kphi,kmat,loglikn)

c+++++++ evaluating the prior

         logpriorn=-0.5d0*dble(nvar+1)*ldetc  
         logprioro=-0.5d0*dble(nvar+1)*ldet  

c+++++++ acceptance step
         
         ratio=loglikn-logliko+logcgkn-logcgko+
     &         logpriorn-logprioro

         if(log(dble(runif())).lt.ratio)then
            do i=1,nvar
               do j=1,nvar
                  sigma(i,j)=sigmac(i,j)
                  uinv(i,j)=uinvc(i,j)
                  sigmachol(i,j)=sigmacholc(i,j)
               end do
            end do
            ldet=ldetc
            acrate(1)=acrate(1)+1.d0
               
            if(adaptives.eq.1)aratesigma=aratesigma+1.d0

          else
            call loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,uinv,
     &                       ldet,kphi,kmat,loglikn)
              
         end if


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Updating partition - ortho matrix                           +++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         do i=1,nvar
            do j=1,nvar

                evali=1
                xx0=ortho(i,j)
                re=rexpo(1.d0)
                logy=-re

c++++++++++++++ evaluates the log-posterior
                do i1=1,nvar
                   do j1=1,nvar   
                       workm2(i1,j1)=ortho(i1,j1)
                   end do
                end do
                call rhaar2(workm,workm2,nvar,omat)

                do i1=1,nvar
                   do j1=1,nvar
                      tmp1=0.d0
                      do k1=1,nvar 
                         tmp1=tmp1+sigmachol(i1,k1)*omat(k1,j1)
                      end do
                      uinv(i1,j1)=tmp1
                   end do
                end do
                call inverse(uinv,nvar,iflag)      

                call loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,uinv,
     &                           ldet,kphi,kmat,tmp1)

                logy=logy+tmp1+dnrm(xx0,0.d0,1.d0,1)

                uwork=dble(runif())*0.5  
                llim=xx0-uwork
                rlim=llim+0.5

                evali=evali+1
                ortho(i,j)=llim
                do i1=1,nvar
                   do j1=1,nvar   
                       workm2(i1,j1)=ortho(i1,j1)
                   end do
                end do
                call rhaar2(workm,workm2,nvar,omat)

                do i1=1,nvar
                   do j1=1,nvar
                      tmp1=0.d0
                      do k1=1,nvar 
                         tmp1=tmp1+sigmachol(i1,k1)*omat(k1,j1)
                      end do
                      uinv(i1,j1)=tmp1
                   end do
                end do
                call inverse(uinv,nvar,iflag)      

                call loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,uinv,
     &                           ldet,kphi,kmat,gllim)
                gllim=gllim+dnrm(llim,0.d0,1.d0,1)


                evali=evali+1
                ortho(i,j)=rlim
                do i1=1,nvar
                   do j1=1,nvar   
                       workm2(i1,j1)=ortho(i1,j1)
                   end do
                end do
                call rhaar2(workm,workm2,nvar,omat)

                do i1=1,nvar
                   do j1=1,nvar
                      tmp1=0.d0
                      do k1=1,nvar 
                         tmp1=tmp1+sigmachol(i1,k1)*omat(k1,j1)
                      end do
                      uinv(i1,j1)=tmp1
                   end do
                end do
                call inverse(uinv,nvar,iflag)      

                call loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,uinv,
     &                           ldet,kphi,kmat,grlim)
                grlim=grlim+dnrm(rlim,0.d0,1.d0,1)


                do while(gllim.gt.logy)
                   llim=llim-0.5
                   evali=evali+1
                   ortho(i,j)=llim
                   do i1=1,nvar
                      do j1=1,nvar   
                          workm2(i1,j1)=ortho(i1,j1)
                      end do
                   end do
                   call rhaar2(workm,workm2,nvar,omat)

                   do i1=1,nvar
                      do j1=1,nvar
                         tmp1=0.d0
                         do k1=1,nvar 
                            tmp1=tmp1+sigmachol(i1,k1)*omat(k1,j1)
                         end do
                         uinv(i1,j1)=tmp1
                      end do
                   end do
                   call inverse(uinv,nvar,iflag)      

                   call loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,
     &                              uinv,ldet,kphi,kmat,gllim)
                   gllim=gllim+dnrm(llim,0.d0,1.d0,1)

                end do 

                do while(grlim.gt.logy)
                   rlim=rlim+0.5
                   evali=evali+1
                   ortho(i,j)=rlim
                   do i1=1,nvar
                      do j1=1,nvar   
                          workm2(i1,j1)=ortho(i1,j1)
                      end do
                   end do
                   call rhaar2(workm,workm2,nvar,omat)

                   do i1=1,nvar
                      do j1=1,nvar
                         tmp1=0.d0
                         do k1=1,nvar 
                            tmp1=tmp1+sigmachol(i1,k1)*omat(k1,j1)
                         end do
                         uinv(i1,j1)=tmp1
                      end do
                   end do
                   call inverse(uinv,nvar,iflag)      

                   call loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,
     &                              uinv,ldet,kphi,kmat,grlim)
                   grlim=grlim+dnrm(rlim,0.d0,1.d0,1)
                end do 
 
                xx1=llim+(rlim-llim)*dble(runif())
                evali=evali+1
                ortho(i,j)=xx1
                do i1=1,nvar
                   do j1=1,nvar   
                       workm2(i1,j1)=ortho(i1,j1)
                   end do
                end do
                call rhaar2(workm,workm2,nvar,omat)

                do i1=1,nvar
                   do j1=1,nvar
                      tmp1=0.d0
                      do k1=1,nvar 
                         tmp1=tmp1+sigmachol(i1,k1)*omat(k1,j1)
                      end do
                      uinv(i1,j1)=tmp1
                   end do
                end do
                call inverse(uinv,nvar,iflag)      

                call loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,uinv,
     &                           ldet,kphi,kmat,gxx1)
                gxx1=gxx1+dnrm(xx1,0.d0,1.d0,1)


                do while(gxx1.lt.logy)
                   if(xx1.gt.xx0)rlim=xx1
                   if(xx1.lt.xx0)llim=xx1
                   xx1=llim+(rlim-llim)*dble(runif())
                   evali=evali+1
                   ortho(i,j)=xx1
                   do i1=1,nvar
                      do j1=1,nvar   
                          workm2(i1,j1)=ortho(i1,j1)
                      end do
                   end do
                   call rhaar2(workm,workm2,nvar,omat)

                   do i1=1,nvar
                      do j1=1,nvar
                         tmp1=0.d0
                         do k1=1,nvar 
                            tmp1=tmp1+sigmachol(i1,k1)*omat(k1,j1)
                         end do
                         uinv(i1,j1)=tmp1
                      end do
                   end do
                   call inverse(uinv,nvar,iflag)      

                   call loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,
     &                              uinv,ldet,kphi,kmat,gxx1)
                   gxx1=gxx1+dnrm(xx1,0.d0,1.d0,1)

                end do
                ortho(i,j)=xx1
 
            end do
         end do

         do i1=1,nvar
            do j1=1,nvar   
               workm2(i1,j1)=ortho(i1,j1)
               workm(i1,j1)=0.d0
            end do
         end do
         call rhaar2(workm,workm2,nvar,omat)

         call dblepr("ortho",-1,ortho,nvar*nvar)


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Updating the precision parameter using a MH step            +++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(a0b0(1).gt.0.d0)then

c++++++++++ Addaptive MH

            if(adaptivec.eq.1)then  
               cskip = cskip + 1
               if(cskip.eq.100)then
                  aratec=aratec/dble(100)
                  if(iscan.le.nadaptive)then  
                     if(aratec.lt.0.44)then
                        tune2=exp(log(tune2)+(0.44-aratec))
                      else
                        tune2=exp(log(tune2)-(aratec-0.44))
                     end if  
                   else 
                     if(aratec.gt.0.44)then
                        tune2=exp(log(tune2)+
     &                        min(0.01,1.d0/sqrt(dble(iscan))))
                       else
                        tune2=exp(log(tune2)-
     &                        min(0.01,1.d0/sqrt(dble(iscan))))
                     end if 
                  end if    
                  cskip=0
                  aratec=0.d0
               end if
            end if

c++++++++++ sample candidates

            alphac=rtlnorm(log(alpha),tune2,0,0,.true.,.true.)
            logcgkn=dlnrm(alpha ,log(alphac),tune2,1) 
            logcgko=dlnrm(alphac,log(alpha) ,tune2,1) 

c++++++++++ logpost

            call loglikefmptcp(nrec,nvar,y,maxm,frstl,alpha,kmat,a0b0,
     &                         logliko)

            call loglikefmptcp(nrec,nvar,y,maxm,frstl,alphac,kmat,a0b0,
     &                         loglikn)

c++++++++++ acceptance step
            ratio=loglikn+logliko+
     &            logcgkn-logcgko

            if(log(dble(runif())).lt.ratio)then
               alpha=alphac
               acrate(2)=acrate(2)+1.d0

               if(adaptivec.eq.1)aratec=aratec+1.d0
            end if            
 
            call dblepr("alpha",-1,alpha,1)
         end if



c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Save samples                                 +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1
              
c+++++++++++++ mu
               do i=1,nvar
                  thetasave(isave,i)=mu(i)
               end do   

c+++++++++++++ sigma
               k1=0
               do i=1,nvar
                  do j=i,nvar
                     k1=k1+1
                     thetasave(isave,nvar+k1)=sigma(i,j)
                  end do
               end do

c+++++++++++++ c parameter
               
               thetasave(isave,nvar+k1+1)=alpha

c+++++++++++++ cpo


c+++++++++++++ density 

               if(nvar.eq.2)then
                  call pdensfmpt(nrec,nvar,maxm,frstl,alpha,
     &                     mu,uinv,ldet,kphi,kvec,kmat,ngrid,grid,f)
               end if 

c+++++++++++++ print
               skipcount = 0
               if(dispcount.ge.ndisplay)then
                  call cpu_time(sec1)
                  sec00=sec00+(sec1-sec0)
                  sec=sec00
                  sec0=sec1
                  pprn=sprint(isave,nsave,sec)
                  dispcount=0
               end if   
            end if        
         end if

      end do

      do i=1,2
         acrate(i)=acrate(i)/dble(nscan)      
      end do     

      do i=1,nrec
         cpo(i)=dble(nsave)/cpo(i)
      end do

      do i=1,ngrid
         do j=1,ngrid
            f(i,j)=f(i,j)/dble(nsave)
         end do
      end do

      return
      end
      

c=======================================================================                      
      subroutine loglikefmpt(nrec,nvar,y,maxm,frstl,alpha,mu,uinv,ldet,
     &                       kphi,kmat,out)
c=======================================================================
c     Alejandro Jara, 2009                  
c=======================================================================
      implicit none
c+++++Input
      integer nrec,nvar
      real*8 y(nrec,nvar)

      integer maxm
      integer frstl
      real*8 alpha
      real*8 mu(nvar)
      real*8 uinv(nvar,nvar)
      real*8 ldet       

c+++++External working space
      integer kphi(nvar)

c+++++Output
      integer kmat(nrec,maxm)
      real*8 out

c+++++Internal working space
      integer i,j,j1,j2,k,k1,ind
      integer n1,n2
      real*8 cdfnorm,dnrm
      real*8 tmp1,tmp2
   
c+++++Algorithm

      out=-0.5*real(nrec)*ldet

      do i=1,nrec

c+++++++ check if the user has requested an interrupt
         call rchkusr()
          
         do j=1,nvar
            tmp1=0.d0
            do k=1,nvar
               tmp1=tmp1+uinv(j,k)*(y(i,k)-mu(k))
            end do

            if(tmp1.gt.4.0d0)then
               tmp2=0.999968d0
             else if(tmp1.lt.-4.0d0)then
               tmp2=0.000032d0
             else 
               tmp2=cdfnorm(tmp1,0.d0,1.d0,1,0)
            end if
            kphi(j)=int(real(2**maxm)*tmp2)
            out=out+dnrm(tmp1,0.d0,1.d0,1)     
         end do

         do j1=1,maxm
            j2=maxm-j1+1

            ind=0
            do k1=1,nvar
               ind=ind+kphi(k1)*(2**((k1-1)*j2))
            end do 
            kmat(i,j2)=ind
 
            do k1=1,nvar
               kphi(k1)=int(real(kphi(k1))/2.0)
            end do
         end do

c+++++++ Add q contribution

         if(i.gt.1)then
            do j1=2,maxm
               n1=count(kmat(1:i-1,j1)==kmat(i,j1))   
               n2=count(kmat(1:i-1,j1-1)==kmat(i,j1-1))
               out=out+real(nvar)*log(2.0)+
     &                 log(alpha*(j1**2)+n1)-
     &                 log((2.0**nvar)*alpha*(j1**2)+n2)
            end do
            if(frstl.ne.0)then
               n1=count(kmat(1:i-1,1)==kmat(i,1))
               out=out+real(nvar)*log(2.0)+
     &                 log(alpha+n1)-
     &                 log((2.0**nvar)*alpha+i-1)
              else
               out=out-real(nvar)*log(2.d0)
            end if   
         end if 

      end do


      return
      end



c=======================================================================                      
      subroutine loglikefmptcp(nrec,nvar,y,maxm,frstl,alpha,kmat,a0b0,
     &                         out)
c=======================================================================
c     Alejandro Jara, 2009                  
c=======================================================================
      implicit none
c+++++Input
      integer nrec,nvar
      real*8 y(nrec,nvar)

      integer maxm
      integer frstl
      real*8 alpha

      integer kmat(nrec,maxm)
      
      real*8 a0b0(2)
 
c+++++Output
      real*8 out

c+++++Internal working space
      integer i,j,j1,j2,k,k1,ind
      integer n1,n2
      real*8 tmp1,tmp2
   
c+++++Algorithm

      out=(a0b0(1)-1.d0)*log(alpha)-
     &     a0b0(2)*alpha

      do i=2,nrec

c+++++++ check if the user has requested an interrupt
         call rchkusr()
          
         do j1=2,maxm
            n1=count(kmat(1:i-1,j1)==kmat(i,j1))   
            n2=count(kmat(1:1-1,j1-1)==kmat(i,j1-1))
            out=out+real(nvar)*log(2.0)+
     &              log(alpha*(j1**2)+n1)-
     &              log((2.0**nvar)*alpha*(j1**2)+n2)
         end do
         if(frstl.ne.0)then
            n1=count(kmat(1:i-1,1)==kmat(i,1))
            out=out+real(nvar)*log(2.0)+
     &              log(alpha+n1)-
     &              log((2.0**nvar)*alpha+i-1)
           else
            out=out-real(nvar)*log(2.d0)
         end if   

      end do

      return
      end

      

c=======================================================================                      
      subroutine pdensfmpt(nrec,nvar,maxm,frstl,alpha,mu,uinv,ldet,
     &                     kphi,kvec,kmat,ngrid,grid,f)
c=======================================================================
c     Alejandro Jara, 2009                  
c=======================================================================
      implicit none
c+++++Input
      integer nrec,nvar
      integer maxm
      integer frstl
      real*8 alpha
      real*8 mu(nvar)
      real*8 uinv(nvar,nvar)
      real*8 ldet       

      integer kmat(nrec,maxm)

      integer ngrid
      real*8 grid(ngrid,nvar)

c+++++External working space
      integer kphi(nvar)
      integer kvec(maxm)

c+++++Output
      real*8 f(ngrid,ngrid)

c+++++Internal working space
      integer i,j,j1,j2,k,k1,ind,l1,l2
      integer n1,n2
      real*8 cdfnorm,dnrm
      real*8 tmp1,tmp2 
      real*8 out
      real*8 ywork(2)
   
c+++++Algorithm

      do l1=1,ngrid
         do l2=1,ngrid

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            out=-0.5d0*ldet
            ywork(1)=grid(l1,1)         
            ywork(2)=grid(l2,2)         

            do j=1,2
               tmp1=0.d0
               do k=1,2
                  tmp1=tmp1+uinv(j,k)*(ywork(k)-mu(k))
               end do

               if(tmp1.gt.4.0d0)then
                  tmp2=0.999968d0
                else if(tmp1.lt.-4.0d0)then
                  tmp2=0.000032d0
                else 
                  tmp2=cdfnorm(tmp1,0.d0,1.d0,1,0)
               end if
               kphi(j)=int(real(2**maxm)*tmp2)
               out=out+dnrm(tmp1,0.d0,1.d0,1)     
            end do

            do j1=1,maxm
               j2=maxm-j1+1

               ind=0
               do k1=1,nvar
                  ind=ind+kphi(k1)*(2**((k1-1)*j2))
               end do 
               kvec(j2)=ind
 
               do k1=1,nvar
                  kphi(k1)=int(real(kphi(k1))/2.0)
               end do
            end do

c++++++++++ Add q contribution

            do j1=2,maxm
               n1=count(kmat(1:nrec,j1)==kvec(j1))   
               n2=count(kmat(1:nrec,j1-1)==kvec(j1-1))
               out=out+real(nvar)*log(2.0)+
     &                 log(alpha*(j1**2)+n1)-
     &                 log((2.0**nvar)*alpha*(j1**2)+n2)
            end do
            if(frstl.ne.0)then
               n1=count(kmat(1:nrec,1)==kvec(1))
               out=out+real(nvar)*log(2.0)+
     &                 log(alpha+n1)-
     &                 log((2.0**nvar)*alpha+i-1)
              else
               out=out-real(nvar)*log(2.d0)
            end if   

            f(l1,l2)=f(l1,l2)+exp(out)

         end do
      end do

      return
      end




