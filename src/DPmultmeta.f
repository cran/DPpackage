c=======================================================================                      
      subroutine dpmultmeta(nrec,nvar,y,sigma2e,                        
     &                      a0b0,m1rand,s2inv,s2invm2,                  
     &                      nu,s1ra nd,psiinv,                          
     &                      mcmc,nsave,                                 
     &                      cpo,randsave,thetasave,                     
     &                      alpha,m1,s1,                                
     &                      ncluster,muclus,ss,ccluster,cstrt,          
     &                      iflag,prob,sigma2ei,seed,                   
     &                      s1inv,s1invm1,                              
     &                      theta,workm1,workm2,workm3,                 
     &                      workmh1,workmh2,                            
     &                      workv1,workv2,ywork)                        
c=======================================================================                      
c 
c     Copyright: Alejandro Jara, 2008-2010.
c
c     Version 1.0:
c
c     Last modification: 20-04-2008.
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
      double precision y(nrec,nvar)
      double precision sigma2e(nrec,nvar*(nvar+1)/2)

c+++++Prior 
      integer nu,m1rand,s1rand
      double precision aa0,ab0,a0b0(2)
      double precision psiinv(nvar,nvar)
      double precision s2inv(nvar,nvar),s2invm2(nvar)

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Output
      double precision cpo(nrec,2)
      double precision randsave(nsave,(nrec+1)*nvar)
      double precision thetasave(nsave,nvar+nvar*(nvar+1)/2+2)

c+++++Current values of the parameters
      integer ncluster,ss(nrec)
      double precision alpha,m1(nvar),s1(nvar,nvar)
      double precision muclus(nrec+1,nvar)
      
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++External working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer iflag(nvar)
      integer seed(2)
      double precision s1inv(nvar,nvar)
      double precision s1invm1(nvar) 
      double precision sigma2ei(nrec,nvar*(nvar+1)/2)
      double precision theta(nvar)
      double precision workm1(nvar,nvar),workm2(nvar,nvar)
      double precision workm3(nvar,nvar)
      double precision workmh1(nvar*(nvar+1)/2) 
      double precision workmh2(nvar*(nvar+1)/2) 
      double precision workv1(nvar),workv2(nvar)
      double precision ywork(nvar)

c+++++DPM
      integer cstrt(nrec,nrec)
      integer ccluster(nrec)
      double precision prob(nrec+1)


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Internal working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++General
      integer count,dispcount,evali
      integer i,ii,ihmssf,isave,iscan
      integer j,k,l,j1,j2
      integer ns,nscan,nuniqs,ok,sprint
      integer seed1,seed2,since,skipcount
      double precision tmp1,tmp2

c+++++CPU time
      double precision sec00,sec0,sec1,sec
      
c++++ Define parameters

      aa0=a0b0(1)
      ab0=a0b0(2)

      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      
      nuniqs=nvar*(nvar+1)/2

c++++ set random number generator

      seed1=seed(1)
      seed2=seed(2)
      call setall(seed1,seed2)
     
c++++ cluster structure

      do i=1,nrec
         ccluster(ss(i))=ccluster(ss(i))+1
         cstrt(ss(i),ccluster(ss(i)))=i
         do j1=1,nvar
            do j2=1,nvar
               workm1(j1,j2)=sigma2e(i,ihmssf(j1,j2,nvar))
            end do
         end do
         call inverse(workm1,nvar,iflag)
         do j1=1,nvar
            do j2=j1,nvar
               sigma2ei(i,ihmssf(j1,j2,nvar))=workm1(j1,j2)
            end do
         end do
      end do
       
      do j1=1,nvar
         do j2=1,nvar
            s1inv(j1,j2)=s1(j1,j2)
         end do
      end do
      call inverse(s1inv,nvar,iflag)
          
      do j1=1,nvar
         tmp1=0.d0
         do j2=1,nvar
            tmp1=tmp1+s1inv(j1,j2)*m1(j2)
         end do
         s1invm1(j1)=tmp1
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

c+++++++ a) Polya Urn based on a collapsed state

         do i=1,nrec
         
            ns=ccluster(ss(i))

            do j=1,nvar
               ywork(j)=y(i,j)
            end do

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
                  do k=1,nvar
                     do l=1,nvar 
                        workm1(k,l)=s1inv(k,l)
                     end do
                     workv1(k)=s1invm1(k) 
                  end do  

                  do k=1,ccluster(j)
                     do j1=1,nvar
                        tmp1=0.d0
                        do j2=1,nvar
                           tmp1=tmp1+sigma2ei(cstrt(j,k),
     &                          ihmssf(j1,j2,nvar))*y(cstrt(j,k),j2)
                           workm1(j1,j2)=workm1(j1,j2)+
     &                          sigma2ei(cstrt(j,k),ihmssf(j1,j2,nvar))
                        end do
                        workv1(j1)=workv1(j1)+tmp1
                     end do
                  end do 

                  call inverse(workm1,nvar,iflag)

                  do j1=1,nvar
                     tmp1=0.d0
                     do j2=1,nvar
                        tmp1=tmp1+workm1(j1,j2)*workv1(j2)
                        workm2(j1,j2)=workm1(j1,j2)+
     &                                sigma2e(i,ihmssf(j1,j2,nvar))
                     end do
                     workv2(j1)=tmp1
                  end do

                  call dmvnd(nvar,ywork,workv2,workm2,tmp2,iflag)
                  prob(j)=exp(log(dble(ccluster(j)))+tmp2)
               end do
               
               do j1=1,nvar
                  workv2(j1)=m1(j1)
                  do j2=1,nvar
                     workm2(j1,j2)=s1(j1,j2)+
     &                             sigma2e(i,ihmssf(j1,j2,nvar))  
                  end do
               end do 
               call dmvnd(nvar,ywork,workv2,workm2,tmp2,iflag)        
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
                  do k=1,nvar
                     do l=1,nvar 
                        workm1(k,l)=s1inv(k,l)
                     end do
                     workv1(k)=s1invm1(k) 
                  end do  

                  do k=1,ccluster(j)
                     do j1=1,nvar
                        tmp1=0.d0
                        do j2=1,nvar
                           tmp1=tmp1+sigma2ei(cstrt(j,k),
     &                          ihmssf(j1,j2,nvar))*y(cstrt(j,k),j2)
                           workm1(j1,j2)=workm1(j1,j2)+
     &                          sigma2ei(cstrt(j,k),ihmssf(j1,j2,nvar))
                        end do
                        workv1(j1)=workv1(j1)+tmp1
                     end do
                  end do 

                  call inverse(workm1,nvar,iflag)

                  do j1=1,nvar
                     tmp1=0.d0
                     do j2=1,nvar
                        tmp1=tmp1+workm1(j1,j2)*workv1(j2)
                        workm2(j1,j2)=workm1(j1,j2)+
     &                                sigma2e(i,ihmssf(j1,j2,nvar))
                     end do
                     workv2(j1)=tmp1
                  end do

                  call dmvnd(nvar,ywork,workv2,workm2,tmp2,iflag)
                  prob(j)=exp(log(dble(ccluster(j)))+tmp2)
               end do
               
               do j1=1,nvar
                  workv2(j1)=m1(j1)
                  do j2=1,nvar
                     workm2(j1,j2)=s1(j1,j2)+
     &                  sigma2e(i,ihmssf(j1,j2,nvar))
                  end do
               end do 
               call dmvnd(nvar,ywork,workv2,workm2,tmp2,iflag)        
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


c+++++++ b) Resampling step

         do ii=1,ncluster

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            ns=ccluster(ii)

            do k=1,nvar
               do l=1,nvar 
                  workm1(k,l)=s1inv(k,l)
               end do
               workv1(k)=s1invm1(k) 
            end do  

            do k=1,ns
               do j1=1,nvar
                  tmp1=0.d0
                  do j2=1,nvar
                     tmp1=tmp1+sigma2ei(cstrt(ii,k),ihmssf(j1,j2,nvar))*
     &                    y(cstrt(ii,k),j2)
                     workm1(j1,j2)=workm1(j1,j2)+sigma2ei(cstrt(ii,k),
     &                             ihmssf(j1,j2,nvar))
                  end do
                  workv1(j1)=workv1(j1)+tmp1
               end do
            end do 

            call inverse(workm1,nvar,iflag)

            do j1=1,nvar
               tmp1=0.d0
               do j2=1,nvar
                  tmp1=tmp1+workm1(j1,j2)*workv1(j2)
               end do
               workv2(j1)=tmp1
            end do

            call rmvnorm(nvar,workv2,workm1,workmh1,workv1,theta)
            
            do i=1,nvar
               muclus(ii,i)=theta(i)
            end do
         end do


c++++++++++++++++++++++++++++++++++         
c+++++++ Base line distribution
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++ mean of the normal component

         if(m1rand.eq.1)then
         
            do i=1,nvar
               workv1(i)=s2invm2(i)
               do j=1,nvar
                  workm1(i,j)=s2inv(i,j)+dble(ncluster)*s1inv(i,j)
               end do
            end do
            call inverse(workm1,nvar,iflag)

            do i=1,ncluster
               do j=1,nvar
                  tmp1=0.d0  
                  do k=1,nvar 
                     tmp1=tmp1+s1inv(j,k)*muclus(i,k)
                  end do
                  workv1(j)=workv1(j)+tmp1
               end do         
            end do

            do j1=1,nvar
               tmp1=0.d0
               do j2=1,nvar
                  tmp1=tmp1+workm1(j1,j2)*workv1(j2)
               end do
               workv2(j1)=tmp1
            end do

            call rmvnorm(nvar,workv2,workm1,workmh1,workv1,theta)
                  
            do i=1,nvar
               m1(i)=theta(i)
            end do

         end if

c+++++++ scale matrix of the inverted-Wishart component

         if(s1rand.eq.1)then
         
            do i=1,nvar
               do j=1,nvar
                  workm1(i,j)=psiinv(i,j)
               end do
            end do

            do i=1,ncluster
               do j=1,nvar 
                  do k=1,nvar 
                     workm1(j,k)=workm1(j,k)+(muclus(i,j)-m1(j))*
     &                                       (muclus(i,k)-m1(k))
                  end do   
               end do
            end do

            call riwishart(nvar,nu+ncluster,workm1,workm2,workm3,
     &                     workv1,workmh1,workmh2,iflag)

            do i=1,nvar
               do j=1,nvar
                  s1(i,j)=workm1(i,j)
                  s1inv(i,j)=workm2(i,j)
               end do
            end do
         end if

         do j1=1,nvar
            tmp1=0.d0
            do j2=1,nvar
               tmp1=tmp1+s1inv(j1,j2)*m1(j2)
            end do
            s1invm1(j1)=tmp1
         end do

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

c+++++++++++++ IW baseline scale

               do i=1,nvar
                  do j=i,nvar
                     count=count+1
                     thetasave(isave,count)=s1(i,j)
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
                     count=count+1
                     randsave(isave,count)=muclus(ss(i),j) 
                  end do   
               end do

c+++++++++++++ predictive information
       
               do i=1,ncluster
                  prob(i)=dble(ccluster(i))/(alpha+dble(nrec))
               end do
               prob(ncluster+1)=alpha/(alpha+dble(nrec))
         
               call simdisc(prob,nrec+1,ncluster+1,evali)

               if(evali.le.ncluster)then
                  do i=1,nvar  
                     theta(i)=muclus(evali,i)
                  end do  
                else
                  call rmvnorm(nvar,m1,s1,workmh1,workv1,theta) 
               end if
               
               do i=1,nvar
                  count=count+1 
                  randsave(isave,count)=theta(i)
               end do
               
c+++++++++++++ cpo and save samples

               do i=1,nrec
                  do j=1,nvar
                     ywork(j)=y(i,j)
                     workv1(j)=muclus(ss(i),j) 
                     do k=1,nvar
                        workm1(j,k)=sigma2e(i,ihmssf(j,k,nvar))
                     end do
                  end do

                  call dmvnd(nvar,ywork,workv1,workm1,tmp1,iflag)

                  tmp1=exp(tmp1)
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
         cpo(i,1)=dble(nsave)/cpo(i,1)
         cpo(i,2)=cpo(i,2)/dble(nsave)
      end do

      return
      end
