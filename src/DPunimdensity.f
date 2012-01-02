c=======================================================================                      
      subroutine dpunimdensity(nrec,y,
     &                         ngrid,grid,
     &                         a0b0,epsilon,
     &                         m0,s0,tau1,taus1,taus2,mb,sb,
     &                         taub1,taub2,
     &                         mcmc,nsave,seed,
     &                         densm,randsave,
     &                         ncluster,ss,alpha,mu,muclus,sigmaclus,
     &                         uvec,
     &                         m1,s1,tau2,
     &                         cstrt,ccluster,prob)
c=======================================================================                      
c
c=======================================================================

      implicit none 

c+++++Data
      integer nrec
      real*8 y(nrec)

c+++++Grid
      integer ngrid
      real*8 grid(ngrid)

c+++++Prior 
      real*8 aa0,ab0,a0b0(2)
      real*8 epsilon
      real*8 m0,s0
      real*8 tau1
      real*8 taus1,taus2
      real*8 mb,sb
      real*8 taub1,taub2

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay
      integer seed(2)

c+++++Output
      real*8 densm(ngrid)
      real*8 randsave(nsave)

c+++++Current values of the parameters
      integer ncluster
      integer ss(nrec)
      real*8 alpha
      real*8 mu
      real*8 muclus(nrec+100)
      real*8 sigmaclus(nrec+100)
      real*8 uvec(nrec)
      real*8 m1,s1
      real*8 tau2

c+++++External working space
      integer cstrt(nrec,nrec)
      integer ccluster(nrec)
      real*8 prob(nrec+100)

c+++++Working space - Distributions and RNG
      integer seed1,seed2
      real*8 dnrm
      real*8 rnorm,rgamma
      real runif

c+++++Working space slice sampling
      integer evali
      real*8 rexpo,re,uwork
      real*8 logy,xx0,xx1,llim,rlim
      real*8 grlim,gllim,gxx0,gxx1

c+++++General working space
      integer i,ii,j,k
      integer isample,ok
      integer ns,nscan,sprint
      integer since,skipcount
      real*8 muwork,sigmawork,uniwork
      real*8 tmp1,tmp2,tmp3,tmp4

c+++++MCMC working space
      integer dispcount
      integer isave,iscan

c+++++CPU time
      real*8 sec00,sec0,sec1,sec
      
c++++ Define parameters

      aa0=a0b0(1)
      ab0=a0b0(2)

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
c+++++++ Mode
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         tmp1=m0/s0
         tmp2=1/s0

         do i=1,nrec
            muwork=uvec(i)*muclus(ss(i))
            sigmawork=uvec(i)*uvec(i)*sigmaclus(ss(i))

            tmp1=tmp1+(y(i)-muwork)/sigmawork
            tmp2=tmp2+1.d0/sigmawork
         end do
         tmp3=tmp1/tmp2
         tmp4=1.d0/tmp2

         mu=rnorm(tmp3,dsqrt(tmp4))

c         call dblepr("mu",-1,mu,1)

c++++++++++++++++++++++++++++++++++         
c+++++++ Uniforms (epsilon,1)
c++++++++++++++++++++++++++++++++++

c         call dblepr("uvec",-1,uvec,nrec)

         do i=1,nrec

            muwork=muclus(ss(i))
            sigmawork=sigmaclus(ss(i))

            evali=1
            xx0=uvec(i)

            tmp3=mu+xx0*muwork
            tmp4=dsqrt(sigmawork)*xx0       
            tmp1=dnrm(y(i),tmp3,tmp4,1)

            re=rexpo(1.d0)
            logy=tmp1-re

            uwork=dble(runif())*0.05d0  
            llim=xx0-uwork
            rlim=xx0+(0.05d0-uwork)

            if(llim.lt.epsilon)llim=epsilon 
            if(rlim.gt.1.0)rlim=1.0 

            evali=evali+1
            tmp3=mu+llim*muwork
            tmp4=dsqrt(sigmawork)*llim       
            gllim=dnrm(y(i),tmp3,tmp4,1)

            evali=evali+1
            tmp3=mu+rlim*muwork
            tmp4=dsqrt(sigmawork)*rlim       
            grlim=dnrm(y(i),tmp3,tmp4,1)

            do while(gllim.gt.logy)
               llim=llim-0.05d0

               if(llim.lt.epsilon)then
                  llim=epsilon
                  gllim=logy-1.d0
                 else   
                  evali=evali+1

                  tmp3=mu+llim*muwork
                  tmp4=dsqrt(sigmawork)*llim       
                  gllim=dnrm(y(i),tmp3,tmp4,1)
               end if
            end do 

            do while(grlim.gt.logy)
               rlim=rlim+0.05d0

               if(rlim.gt.1.0)then
                  rlim=1.0
                  grlim=logy-1.d0
                 else   
                  evali=evali+1

                  tmp3=mu+rlim*muwork
                  tmp4=dsqrt(sigmawork)*rlim       
                  grlim=dnrm(y(i),tmp3,tmp4,1)
               end if
            end do 

            xx1=llim+(rlim-llim)*dble(runif())

            evali=evali+1
            tmp3=mu+xx1*muwork
            tmp4=dsqrt(sigmawork)*xx1       
            gxx1=dnrm(y(i),tmp3,tmp4,1)

            do while(gxx1.lt.logy)
               if(xx1.gt.xx0)rlim=xx1
               if(xx1.lt.xx0)llim=xx1

               if(llim.lt.epsilon)llim=epsilon 
               if(rlim.gt.1.0)rlim=1.0 

               xx1=llim+(rlim-llim)*dble(runif())

               evali=evali+1
               tmp3=mu+xx1*muwork
               tmp4=dsqrt(sigmawork)*xx1       
               gxx1=dnrm(y(i),tmp3,tmp4,1)

c               call intpr("eavli",-1,evali,1)
c               call dblepr("llim",-1,llim,1)
c               call dblepr("rlim",-1,rlim,1)
c               call dblepr("xx1",-1,xx1,1)

            end do


            uvec(i)=xx1

c            call dblepr("uvec",-1,uvec(i),1)
         end do

c++++++++++++++++++++++++++++++++++         
c+++++++ DP part
c++++++++++++++++++++++++++++++++++

c++++++++++++++++++++++++++++++
c+++++++ a) configurations 
c++++++++++++++++++++++++++++++

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
                  call relabelunim(i,since,nrec,ncluster,ccluster,
     &                       ss,cstrt,muclus,sigmaclus)

               end if
               ccluster(ncluster)=ccluster(ncluster)-1 
               ncluster=ncluster-1
            end if 
         
            do j=1,ncluster
               muwork=mu+uvec(i)*muclus(j)
               sigmawork=dsqrt(sigmaclus(j))*uvec(i)
               tmp1=dnrm(y(i),muwork,sigmawork,1)
               prob(j)=exp(log(dble(ccluster(j)))+tmp1)
            end do

            if(isample.eq.1)then
               sigmawork=1.d0/rgamma(0.5*tau1,0.5*tau2)
               sigmaclus(ncluster+1)=sigmawork

               muwork=rnorm(m1,dsqrt(s1))
               muclus(ncluster+1)=muwork
            end if         


            muwork=mu+uvec(i)*muclus(ncluster+1)
            sigmawork=dsqrt(sigmaclus(ncluster+1))*uvec(i)

            tmp1=dnrm(y(i),muwork,sigmawork,1)

            prob(ncluster+1)=exp(log(alpha)+tmp1)

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


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ b) Resampling step
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++++++++++++++++++++++++++++++
c+++++++ means                   +++
c+++++++++++++++++++++++++++++++++++

         do i=1,ncluster
            
            tmp1=m1/s1
            tmp2=1/s1
            sigmawork=sigmaclus(i)

            ns=ccluster(i)
            do j=1,ns

c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

               ii=cstrt(i,j)

               tmp1=tmp1+(y(ii)-mu)/sigmawork
               tmp2=tmp2+uvec(ii)/sigmawork

            end do 

            tmp3=tmp1/tmp2
            tmp4=1.d0/tmp2

            muclus(i)=rnorm(tmp3,dsqrt(tmp4))

         end do

c         call dblepr("mu",-1,muclus,ncluster)

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
               tmp1=(y(ii)-mu-uvec(ii)*muclus(i))/uvec(ii)
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
            call samalph(alpha,aa0,ab0,ncluster,nrec)
         end if 

c++++++++++++++++++++++++++++++++++         
c+++++++ baseline normal parameters
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         tmp1=mb/sb
         tmp2=1/sb

         do i=1,ncluster
            tmp1=tmp1+muclus(i)/s1
            tmp2=tmp2+1.d0/s1
         end do
         tmp3=tmp1/tmp2
         tmp4=1.d0/tmp2

         m1=rnorm(tmp3,dsqrt(tmp4))

         tmp1=0.d0
         do i=1,ncluster
            tmp1=tmp1+(muclus(i)-m1)*(muclus(i)-m1)
         end do

         s1=1.d0/
     &         rgamma(0.5d0*(dble(ncluster)+taub1),0.5d0*(tmp1+taub2))


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


c++++++++++++++++++++++++++++++++++         
c+++++++ save samples
c++++++++++++++++++++++++++++++++++         

         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1

c+++++++++++++ sampling from the predictive

               sigmawork=1.d0/rgamma(0.5*tau1,0.5*tau2)
               muwork=rnorm(m1,dsqrt(s1))
               muclus(ncluster+1)=muwork
               sigmaclus(ncluster+1)=sigmawork

               do i=1,ncluster
                  ns=ccluster(i)
                  prob(i)=dble(ns)/(dble(nrec)+alpha)
               end do
               prob(ncluster+1)=alpha/(dble(nrec)+alpha)

               call simdisc(prob,nrec+100,ncluster+1,evali)
               muwork=muclus(evali)
               sigmawork=sigmaclus(evali)
               uniwork=epsilon+(1.d0-epsilon)*dble(runif())

               tmp1=mu+uniwork*muwork
               tmp2=dsqrt(sigmawork)*uniwork
               tmp3=rnorm(tmp1,tmp2)
               randsave(isave)=tmp3


c+++++++++++++ evaluating the density in a grid

               uniwork=epsilon+(1.d0-epsilon)*dble(runif())

               do i=1,ngrid 
                  tmp1=grid(i)
                  tmp2=0.d0
                  do j=1,ncluster+1
                     muwork=mu+uniwork*muclus(j)
                     sigmawork=dsqrt(sigmaclus(j))*uniwork 
                     tmp2=tmp2+prob(j)*dnrm(tmp1,muwork,sigmawork,0)
                  end do
                  densm(i)=densm(i)+tmp2
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
         densm(i)=densm(i)/dble(nsave)
      end do

      return
      end

c=======================================================================
      subroutine relabelunim(ind,since,nrec,ncluster,ccluster,
     &                       ss,cstrt,
     &                       muclus,sigmaclus)
c=======================================================================
      implicit none
      integer i,j,ind,since,nrec,ncluster,ccluster(nrec)
      integer ss(nrec),cstrt(nrec,nrec)
      real*8 muclus(nrec+100)
      real*8 sigmaclus(nrec+100)
      real*8 muwork
      real*8 sigmawork

      integer ns,ii
      
      muwork=muclus(since)
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

         muclus(i-1)=muclus(i)
         sigmaclus(i-1)=sigmaclus(i)
         ccluster(i-1)=ccluster(i)
      end do

      ss(ind)=ncluster
      
      muclus(ncluster)=muwork
      sigmaclus(ncluster)=sigmawork
  
      ccluster(ncluster)=1
      
      return
      end
      
