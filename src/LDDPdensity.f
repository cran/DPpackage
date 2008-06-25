c=======================================================================                      
      subroutine lddpmnormals(nrec,p,y,x,                               #4
     &                        npred,xpred,ngrid,grid,                   #4
     &                        a0b0,tau,murand,sm,sinv,nu,psiinv,        #7
     &                        mcmc,nsave,                               #2
     &                        cpo,densr,randsave,thetasave,             #4
     &                        alpha,b,mub,sigmab,sigma2k,               #5
     &                        ncluster,ss,ccluster,cstrt,               #4
     &                        iflagc,prob,seed,                         #3
     &                        theta,ztz,zty,sigmabinv,                  #4
     &                        workmc,workmc2,workvc,workmhc,workmhc2)   #5
c=======================================================================                      
c     # of arguments = 42.
c 
c
c=======================================================================

      implicit none 

c+++++Data
      integer nrec,p
      real*8 y(nrec)
      real*8 x(nrec,p)

c+++++Prediction
      integer npred,ngrid
      real*8 xpred(npred,p)
      real*8 grid(ngrid)

c+++++Prior 
      integer murand,nu 
      real*8 aa0,ab0,a0b0(2)
      real*8 tau(2),tau1,tau2
      real*8 sm(p),sinv(p,p)
      real*8 psiinv(p,p) 

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Output
      real*8 cpo(nrec,2)
      real*8 densr(npred,ngrid)
      real*8 randsave(nsave,nrec*p)
      real*8 thetasave(nsave,p+p*(p+1)/2+3)

c+++++Current values of the parameters
      integer ncluster,ss(nrec)
      real*8 alpha,b(nrec,p),mub(p),sigmab(p,p)
      real*8 sigma2k
      
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++External working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer iflagc(p)
      integer seed(2)
      real*8 sigmabinv(p,p)   
      real*8 theta(p)
      real*8 workmc(p,p),workmc2(p,p)  
      real*8 workvc(p),workmhc(p*(p+1)/2),workmhc2(p*(p+1)/2)
      real*8 ztz(p,p),zty(p)
     
c+++++DPM
      integer cstrt(nrec,nrec)
      integer ccluster(nrec)
      real*8 prob(nrec+1)


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Internal working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++General
      integer count,dispcount,evali
      integer i,ii,ihmssf,isave,iscan
      integer j,k,l,j1,j2
      integer ns,nscan,nuwork,ok,sprint
      integer seed1,seed2,since,skipcount
      real*8 detlog
      real*8 tmp1,tmp2
      real*8 ssb

c+++++RNG and functions
      real*8 dnrm 
      real*8 rgamma

c+++++CPU time
      real*8 sec00,sec0,sec1,sec
      
c++++ Define parameters

      aa0=a0b0(1)
      ab0=a0b0(2)

      tau1=tau(1)
      tau2=tau(2)

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
      end do

      do i=1,p
         do j=1,p
            sigmabinv(i,j)=sigmab(i,j)
         end do
      end do 
      call inverse(sigmabinv,p,iflagc) 

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

c+++++++ a) Polya Urn based on a collapsed state

         call lddpmregconfuni(nrec,p,
     &                        ncluster,ccluster,cstrt,ss,
     &                        y,x,sigma2k,
     &                        mub,sigmab,sigmabinv,alpha,
     &                        iflagc,prob,theta,
     &                        ztz,zty)


c+++++++ b) Resampling step

         call lddpmreglocuni(nrec,p,
     &                       ncluster,ccluster,cstrt,
     &                       y,x,sigma2k,
     &                       mub,sigmabinv,b,
     &                       iflagc,theta,
     &                       workvc,workmhc,
     &                       ztz,zty)


c++++++++++++++++++++++++++++++++++         
c+++++++ Kernel variance
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(tau1.gt.0.d0)then
            ssb=0.d0 
            do i=1,nrec
               tmp1=0.d0
               do j=1,p
                  tmp1=tmp1+x(i,j)*b(ss(i),j)
               end do
               tmp2=y(i)-tmp1
               ssb=ssb+tmp2*tmp2
            end do

            sigma2k=1.d0/
     &           rgamma(0.5d0*(dble(nrec)+tau1),0.5d0*(ssb+tau2))
     
         end if

c++++++++++++++++++++++++++++++++++         
c+++++++ Baseline mean
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()
         
         if(murand.eq.1)then
            do i=1,p
               zty(i)=sm(i)
               do j=1,p
                  ztz(i,j)=(sigmabinv(i,j)*dble(ncluster))+sinv(i,j)
               end do
            end do
 
            call inverse(ztz,p,iflagc) 

            do i=1,ncluster
               do j=1,p
                  tmp1=0.d0
                  do k=1,p
                     tmp1=tmp1+sigmabinv(j,k)*b(i,k)
                  end do
                  zty(j)=zty(j)+tmp1
               end do
            end do

            do i=1,p
               tmp1=0.d0
               do j=1,p
                  tmp1=tmp1+ztz(i,j)*zty(j)
               end do
               workvc(i)=tmp1
            end do

            call rmvnorm(p,workvc,ztz,workmhc,zty,mub)
         end if
        
c++++++++++++++++++++++++++++++++++         
c+++++++ Baseline variance
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(nu.gt.0)then

            do i=1,p
               do j=1,p
                  ztz(i,j)=psiinv(i,j)
               end do
            end do

            do i=1,ncluster
               do j=1,p
                  do k=1,p
                     ztz(j,k)=ztz(j,k)+               
     &                    (b(i,j)-mub(j))*(b(i,k)-mub(k))
                  end do
               end do
            end do

            call riwishart(p,nu+ncluster,ztz,workmc,workmc2,
     &                     workvc,workmhc,workmhc2,iflagc)

            do i=1,p
               do j=1,p
                  sigmab(i,j)=ztz(i,j)
                  sigmabinv(i,j)=workmc(i,j)
               end do
            end do
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

               do i=1,p
                  count=count+1
                  thetasave(isave,count)=mub(i)
               end do   

c+++++++++++++ IW baseline scale

               do i=1,p
                  do j=i,p
                     count=count+1
                     thetasave(isave,count)=sigmab(i,j)
                  end do
               end do
 
c+++++++++++++ kernel variance
               count=count+1
               thetasave(isave,count)=sigma2k

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
                     randsave(isave,count)=b(ss(i),j) 
                  end do   
               end do

c+++++++++++++ predictive information

               do i=1,ncluster
                  prob(i)=dble(ccluster(i))/(alpha+dble(nrec))
               end do
               prob(ncluster+1)=alpha/(alpha+dble(nrec))
               call rmvnorm(p,mub,sigmab,workmhc,workvc,theta)

               do i=1,npred

                  call simdisc(prob,nrec+1,ncluster+1,evali)

                  do j=1,ncluster
                     tmp1=0.d0 
                     do k=1,p
                        tmp1=tmp1+xpred(i,k)*b(j,k) 
                     end do
                     
                     do k=1,ngrid
                        tmp2=prob(j)*
     &                       dnrm(grid(k),tmp1,sqrt(sigma2k),0)
                        densr(i,k)=densr(i,k)+tmp2
                     end do   
                  end do

                  tmp1=0.d0 
                  do k=1,p
                     tmp1=tmp1+xpred(i,k)*theta(k) 
                  end do                
                  
                  do k=1,ngrid
                     tmp2=prob(ncluster+1)*
     &                    dnrm(grid(k),tmp1,sqrt(sigma2k),0)
                     densr(i,k)=densr(i,k)+tmp2
                  end do   
               end do
               
c+++++++++++++ cpo and save samples

               do i=1,nrec
                  tmp1=0.d0
                  do l=1,p
                     tmp1=tmp1+x(i,l)*b(ss(i),l)
                  end do
                  tmp2=dnrm(y(i),tmp1,sqrt(sigma2k),0)
                  cpo(i,1)=cpo(i,1)+1.0d0/tmp2 
                  cpo(i,2)=cpo(i,2)+tmp2                     
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
         cpo(i,1)=dble(nsave)/cpo(i,1)
         cpo(i,2)=cpo(i,2)/dble(nsave)                                    
      end do

      do i=1,npred
         do j=1,ngrid  
            densr(i,j)=densr(i,j)/dble(nsave)
         end do
      end do
 
      return
      end
