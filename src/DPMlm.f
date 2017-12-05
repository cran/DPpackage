c=======================================================================                      
      subroutine dpmlm(nfixed,nrec,p,x,y,                               
     &                 a0b0,sb,prec,smu,psiinv,                         
     &                 tau,murand,sigmarand,                            
     &                 ncluster,ss,alpha,beta,mu,mub,sigmab,sigma,      
     &                 mcmc,nsave,                                      
     &                 cpo,musave,                                      
     &                 randsave,thetasave,clustsave,                     
     &                 ngrid,grid,fun,                                   
     &                 seed,                                            
     &                 iflagp,xtx,xty,workmhp1,workvp1,                 
     &                 res,                                             
     &                 cstrt,ccluster,prob,                             
     &                 betasave,bsave,mc)
c=======================================================================                      
c     # of arguments = 44.
c
c     Copyright: Alejandro Jara, 2007-2010.
c
c     Version 1.0:
c
c     Last modification: 20-04-2007.
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
      integer nfixed,nrec,p
      double precision x(nrec,p),y(nrec)

c+++++Prior 
      integer murand,sigmarand
      double precision aa0,ab0,a0b0(2)
      double precision sb(p),prec(p,p)
      double precision smu,psiinv
      double precision tau(4),tau01,tau02,tau11,tau12

c+++++Current values of the parameters
      integer ncluster,ss(nrec)
      double precision alpha
      double precision beta(p)
      double precision sigma
      double precision mub,sigmab
      double precision mu(nrec)

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Output
      integer ngrid
      integer clustsave(nsave,nrec) 
      double precision cpo(nrec,2)
      double precision grid(ngrid),fun(ngrid)
      double precision musave(nsave,nrec)
      double precision randsave(nsave,nrec+1)
      double precision thetasave(nsave,nfixed+6)

c+++++Seeds
      integer seed(2),seed1,seed2

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++External working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++DP
      integer cstrt(nrec,nrec)
      integer ccluster(nrec)
      double precision prob(nrec+1)
      
c+++++Residuals
      double precision res(nrec)

c+++++Regression
      integer iflagp(p) 
      double precision xtx(p,p)
      double precision xty(p)
      double precision workmhp1(p*(p+1)/2)
      double precision workvp1(p)

c++++ model's performance
      double precision mc(5)
      double precision betasave(p+1),bsave(nrec)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Internal working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++General
      integer evali,ii,i,j,k,l,ns 
      integer ok
      integer since,sprint 
      double precision betar
      double precision theta,tmp1,tmp2,tmp3
      double precision ztz,zty
      double precision sigmainv,sigmabinv

c+++++MCMC
      integer dispcount,isave,iscan,nscan,skipcount 

c+++++CPU time
      double precision sec00,sec0,sec1,sec

c+++++RNG and distributions
      double precision dnrm,rnorm,rgamma

c++++ model's performance
      double precision dbarc,dbar,dhat,pd,lpml

c+++++DP (functional parameter)
      double precision eps,rbeta,weight
      parameter(eps=0.01)

c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      
      aa0=a0b0(1)
      ab0=a0b0(2)
      
      tau01=tau(1)
      tau02=tau(2)
      tau11=tau(3)
      tau12=tau(4)

c++++ set random number generator
      seed1=seed(1)
      seed2=seed(2)
      call setall(seed1,seed2)

c++++ set configurations
      do i=1,nrec
         ccluster(ss(i))=ccluster(ss(i))+1
         cstrt(ss(i),ccluster(ss(i)))=i
      end do


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      dbar=0.d0
      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*(nsave)


      sigmainv=1.d0/sigma
      sigmabinv=1.d0/sigmab

      call cpu_time(sec0)
      sec00=0.d0
      
      do iscan=1,nscan

c+++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Regression coefficients
c+++++++++++++++++++++++++++++++++++++++++++++++++

         if(nfixed.gt.0)then

            do i=1,p
               xty(i)=sb(i)            
               do j=1,p
                  xtx(i,j)=prec(i,j)
               end do
            end do

            do i=1,nrec
               tmp1=y(i)-mu(ss(i))

               do j=1,p
                  do l=1,p
                     xtx(j,l)=xtx(j,l)+x(i,j)*x(i,l)/sqrt(sigma)
                  end do
                  xty(j)=xty(j)+x(i,j)*tmp1/sqrt(sigma)
               end do
            end do

            call inverse(xtx,p,iflagp)      

            do i=1,p
               tmp1=0.d0
               do j=1,p
                  tmp1=tmp1+xtx(i,j)*xty(j) 
               end do
               workvp1(i)=tmp1
            end do

            call rmvnorm(p,workvp1,xtx,workmhp1,xty,beta)

            do i=1,nrec
               tmp1=0.d0
               do j=1,p
                  tmp1=tmp1+x(i,j)*beta(j)
               end do
               res(i)=y(i)-tmp1
            end do
            
         end if          

c+++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ a) Polya Urn based on a collapsed state
c+++++++++++++++++++++++++++++++++++++++++++++++++

         do i=1,nrec
         
            ns=ccluster(ss(i))

c++++++++++ subject in cluster with more than 1 observations
             
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

               do j=1,ncluster
                  ztz=sigmabinv+dble(ccluster(j))*sigmainv
                  zty=sigmabinv*mub
                  ztz=1.d0/ztz
                  do k=1,ccluster(j)
                     zty=zty+sigmainv*res(cstrt(j,k))   
                  end do 
                  tmp1=ztz*zty
                  ztz=ztz+sigma
                  tmp2=dnrm(res(i),tmp1,sqrt(ztz),1)
                  prob(j)=exp(log(dble(ccluster(j)))+tmp2)
               end do
               
               tmp1=mub
               ztz=sigma+sigmab
               tmp2=dnrm(res(i),tmp1,sqrt(ztz),1)
               prob(ncluster+1)=exp(log(alpha)+tmp2)

               call simdisc(prob,nrec+1,ncluster+1,evali)

               ss(i)=evali
               
               ccluster(evali)=ccluster(evali)+1
               
               cstrt(evali,ccluster(evali))=i
               
               if(evali.gt.ncluster)then
                  ncluster=ncluster+1
               end if
            end if


c++++++++++ subject in cluster with only 1 observation
             
            if(ns.eq.1)then
                
               since=ss(i)
                
               if(since.lt.ncluster)then
                   call relabeldpm(i,since,nrec,1,ncluster,
     &                             ccluster,ss,cstrt)                   
               end if

               ccluster(ncluster)=ccluster(ncluster)-1 
               ncluster=ncluster-1

               do j=1,ncluster
                  ztz=sigmabinv+dble(ccluster(j))*sigmainv
                  zty=sigmabinv*mub
                  ztz=1.d0/ztz
                  do k=1,ccluster(j)
                     zty=zty+sigmainv*res(cstrt(j,k))   
                  end do 
                  tmp1=ztz*zty
                  ztz=ztz+sigma
                  tmp2=dnrm(res(i),tmp1,sqrt(ztz),1)
                  prob(j)=exp(log(dble(ccluster(j)))+tmp2)
               end do
               
               tmp1=mub
               ztz=sigma+sigmab
               tmp2=dnrm(res(i),tmp1,sqrt(ztz),1)
               prob(ncluster+1)=exp(log(alpha)+tmp2)

               call simdisc(prob,nrec+1,ncluster+1,evali)

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

            ns=ccluster(ii)
            
            ztz=sigmabinv+dble(ns)*sigmainv
            zty=sigmabinv*mub
            ztz=1.d0/ztz
            
            do k=1,ns
               zty=zty+sigmainv*res(cstrt(ii,k))    
            end do 

            tmp1=ztz*zty

            theta=rnorm(tmp1,sqrt(ztz))            
            
            mu(ii)=theta
         end do

c++++++++++++++++++++++++++++++
c+++++++ Kernel variance
c++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         tmp1=0.d0
         do i=1,nrec
            tmp1=tmp1+(res(i)-mu(ss(i)))*(res(i)-mu(ss(i)))
         end do

         sigma=1.d0/
     &         rgamma(0.5d0*(dble(nrec)+tau01),0.5d0*(tmp1+tau02))

         sigmainv=1.d0/sigma 


c++++++++++++++++++++++++++++++++++         
c+++++++ Base line distribution
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(murand.eq.1)then
            zty=smu
            ztz=(sigmabinv*dble(ncluster))+psiinv
            ztz=1.d0/ztz

            do i=1,ncluster
               zty=zty+sigmabinv*mu(i) 
            end do
            tmp1=ztz*zty

            theta=rnorm(tmp1,sqrt(ztz))            
            mub=theta
         end if   

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(sigmarand.eq.1)then
            tmp1=0.d0
            do i=1,ncluster
               tmp1=tmp1+(mu(i)-mub)*(mu(i)-mub)
            end do

            sigmab=1.d0/
     &          rgamma(0.5d0*(dble(ncluster)+tau11),0.5d0*(tmp1+tau12))

            sigmabinv=1.d0/sigmab 
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

c+++++++++++++ regression coefficients

c+++++++++++++ random effects
               do i=1,ncluster
                  musave(isave,i)=mu(i)
                  clustsave(isave,i)=ccluster(i)
               end do

               do i=1,nrec
                  bsave(i)=bsave(i)+mu(ss(i))
                  randsave(isave,i)=res(i)
               end do

c+++++++++++++ predictive information

               do i=1,ncluster
                  prob(i)=dble(ccluster(i))/(alpha+dble(nrec))
               end do
               prob(ncluster+1)=alpha/(alpha+dble(nrec))

               call simdisc(prob,nrec+1,ncluster+1,evali)
               
               if(evali.le.ncluster)then
                  theta=mu(evali)
                else
                  theta=rnorm(mub,sqrt(sigmab))  
               end if
               tmp1=rnorm(theta,sqrt(sigma))  
               randsave(isave,nrec+1)=tmp1

               do i=1,ngrid
                  tmp2=0.d0
                  do j=1,ncluster
                     tmp2=tmp2+prob(j)*dnrm(grid(i),mu(j),
     &                                      sqrt(sigma),0)
                  end do
                  tmp2=tmp2+prob(ncluster+1)*dnrm(grid(i),tmp1,
     &                                      sqrt(sigma),0)
                  fun(i)=fun(i)+tmp2                    
               end do

c+++++++++++++ functional parameters
               
               tmp1=rbeta(1.d0,alpha+dble(nrec))
               betar=tmp1*theta
               tmp2=tmp1
               weight=(1.d0-tmp1)
               
               do while((1.d0-tmp2).gt.eps)
                  tmp3=rbeta(1.d0,alpha+dble(nrec))
                  tmp1=weight*tmp3
                  weight=weight*(1.d0-tmp3)

                  do i=1,ncluster
                     prob(i)=dble(ccluster(i))/(alpha+dble(nrec))
                  end do
                  prob(ncluster+1)=alpha/(alpha+dble(nrec))

                  call simdisc(prob,nrec+1,ncluster+1,evali)
               
                  if(evali.le.ncluster)then
                     theta=mu(evali)
                   else
                     theta=rnorm(mub,sqrt(sigmab))  
                  end if
                  betar=betar+tmp1*theta
                  tmp2=tmp2+tmp1
               end do

               do i=1,ncluster
                  prob(i)=dble(ccluster(i))/(alpha+dble(nrec))
               end do
               prob(ncluster+1)=alpha/(alpha+dble(nrec))

               call simdisc(prob,nrec+1,ncluster+1,evali)
               
               if(evali.le.ncluster)then
                  theta=mu(evali)
                else
                  theta=rnorm(mub,sqrt(sigmab))
               end if
               
               tmp1=weight
               betar=betar+tmp1*theta

c+++++++++++++ regression coefficients

               thetasave(isave,1)=betar

               if(nfixed.gt.0)then
                  do i=1,p
                     thetasave(isave,1+i)=beta(i)
                     betasave(i)=betasave(i)+beta(i)
                  end do
               end if
               betasave(p+1)=betasave(p+1)+sigma

c+++++++++++++ kernel variance
               thetasave(isave,1+nfixed+1)=sigma

c+++++++++++++ baseline mean
               thetasave(isave,1+nfixed+2)=mub

c+++++++++++++ baseline covariance
               thetasave(isave,1+nfixed+3)=sigmab

c+++++++++++++ cluster information
               thetasave(isave,1+nfixed+4)=ncluster
               thetasave(isave,1+nfixed+5)=alpha

c+++++++++++++ cpo
               dbarc=0.d0
               do i=1,nrec
                  tmp1=0.d0
                  if(nfixed.gt.0)then
                     do j=1,p
                        tmp1=tmp1+x(i,j)*beta(j)
                     end do   
                  end if
                  tmp1=tmp1+mu(ss(i)) 
                  tmp2=dnrm(y(i),tmp1,sqrt(sigma),0)
                  cpo(i,1)=cpo(i,1)+1.0d0/tmp2  
                  cpo(i,2)=cpo(i,2)+tmp2                    
                  tmp2=dnrm(y(i),tmp1,sqrt(sigma),1)
                  dbarc=dbarc+tmp2
               end do

c+++++++++++++ dic
               dbar=dbar-2.d0*dbarc

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
         fun(i)=fun(i)/dble(nsave)       
      end do     

      do i=1,nrec
         cpo(i,1)=dble(nsave)/cpo(i,1)
         cpo(i,2)=cpo(i,2)/dble(nsave)
      end do

      do i=1,p+1
         betasave(i)=betasave(i)/dble(nsave)
      end do

      do i=1,nrec
         bsave(i)=bsave(i)/dble(nsave)
      end do   

      dhat=0.d0
      lpml=0.d0
      do i=1,nrec
         tmp1=0.d0
         if(nfixed.gt.0)then
            do j=1,p
               tmp1=tmp1+x(i,j)*betasave(j)
            end do   
         end if
         tmp1=tmp1+bsave(i) 
         dhat=dhat+dnrm(y(i),tmp1,sqrt(betasave(p+1)),1)
         lpml=lpml+log(cpo(i,1))
      end do
      dhat=-2.d0*dhat

      dbar=dbar/dble(nsave)
      pd=dbar-dhat
      
      mc(1)=dbar
      mc(2)=dhat
      mc(3)=pd
      mc(4)=dbar+pd
      mc(5)=lpml


      return
      end
