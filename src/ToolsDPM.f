c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR DPM
c=======================================================================                  
c=======================================================================                  

c=======================================================================                  
c=======================================================================                  
c=======================================================================                  
c     AUXILIARY FUNCTIONS/COMMON FUNCTIONS
c=======================================================================                  
c=======================================================================                  
c=======================================================================                  

c=======================================================================      
      subroutine relabeldpmc(ind,since,nsubject,ncluster,ccluster,ss,
     &                       cstrt)
c=======================================================================
c     relabel the clusters after the elimination of one of them
c     A.J.V., 2007
      implicit none
      integer i,j,ind,since,nsubject,ncluster,ccluster(nsubject)
      integer ss(nsubject),cstrt(nsubject,nsubject)

      integer ns,ii
      
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

         ccluster(i-1)=ccluster(i)
      end do
      
      ss(ind)=ncluster
      ccluster(ncluster)=1
      
      return
      end  


c=======================================================================      
      subroutine relabeldpm(ind,since,nsubject,q,ncluster,ccluster,ss,
     &                      cstrt)
c=======================================================================
c     relabel the clusters after the elimination of one of them
c     A.J.V., 2007
      implicit none
      integer i,j,ind,since,nsubject,q,ncluster,ccluster(nsubject)
      integer ss(nsubject),cstrt(nsubject,nsubject)

      integer ns,ii
      ii=q        

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

         ccluster(i-1)=ccluster(i)
      end do
      
      ss(ind)=ncluster
      ccluster(ncluster)=1
      
      return
      end  

c=======================================================================                  
c=======================================================================                  
c=======================================================================                  
c     EVALUATION OF PREDICTIVE DISTRIBUTION IN A GRID
c=======================================================================                  
c=======================================================================                  
c=======================================================================                  

c=======================================================================            
      subroutine predictivedpmu(nsave,nsubject,q,musave,sigmakmat,
     &                          clustsave,  
     &                          nclusvec,alphavec,mubmat,sigmabmat,
     &                          ngrid,grid,fs,seed1,seed2)
c=======================================================================                 
      implicit none
c+++++Input      
      integer ngrid,nsave,nsubject,q
      integer clustsave(nsave,nsubject),nclusvec(nsave)
      integer seed1,seed2
      double precision grid(ngrid)
      double precision sigmakmat(nsave,q*(q+1)/2)
      double precision alphavec(nsave),mubmat(nsave,q),
     1  sigmabmat(nsave,q*(q+1)/2)
      double precision musave(nsave,q*nsubject)
      
c+++++Output
      double precision fs(ngrid)
      
c+++++Working
      integer i,j,k,ncluster,ns 
      double precision alpha,denom,dnrm,meanc,meanb,meanw,rnorm,sdb,sdk
      double precision tmp1

c++++ Set the random number generator
      call setall(seed1,seed2)

c+++++Algorithm
      do i=1,nsave
c+++++++ check if the user has requested an interrupt
         call rchkusr()

         alpha=alphavec(i)
         sdk=sqrt(sigmakmat(i,1))
         meanb=mubmat(i,1)
         sdb=sqrt(sigmabmat(i,1))
         ncluster=nclusvec(i)
         meanc=rnorm(meanb,sdb)
         denom=dble(nsubject)+alpha

         do j=1,ngrid
c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            tmp1=0.d0
            do k=1,ncluster
               ns=clustsave(i,k)
               meanw=musave(i,k)

               tmp1=tmp1+(dble(ns)/denom)*dnrm(grid(j),meanw,sdk,0)
            end do
            tmp1=tmp1+(alpha/denom)*dnrm(grid(j),meanc,sdk,0)
            fs(j)=fs(j)+tmp1
         end do
      end do
      
      do i=1,ngrid
         fs(i)=fs(i)/dble(nsave) 
      end do
      return
      end
      

c=======================================================================            
      subroutine predictivedpmb(nsave,nsubject,q,musave,sigmakmat,
     &                          clustsave,  
     &                          nclusvec,alphavec,mubmat,sigmabmat,
     &                          ngrid1,ngrid2,grid1,grid2,
     &                          fs,f1,f2,seed1,seed2,iflagr,
     &                          meanc,meanb,meanw,sigmab,sigmak,theta,
     &                          workmhr,workvr)
c=======================================================================                 
      implicit none
c+++++Input      
      integer ngrid1,ngrid2,nsave,nsubject,q
      integer clustsave(nsave,nsubject)
      integer nclusvec(nsave)
      integer seed1,seed2
      double precision grid1(ngrid1),grid2(ngrid2)
      double precision sigmakmat(nsave,q*(q+1)/2)
      double precision alphavec(nsave),mubmat(nsave,q),
     1  sigmabmat(nsave,q*(q+1)/2)
      double precision musave(nsave,q*nsubject)
      
      integer iflagr(q)
      double precision meanc(q),meanb(q),meanw(q),sigmab(q,q),
     1  sigmak(q,q),theta(q)
      double precision workmhr(q*(q+1)/2),workvr(q)
      
c+++++Output
      double precision fs(ngrid1,ngrid2),f1(ngrid1),f2(ngrid2)
      
c+++++Working
      integer counter,i,ii,j,jj,k,l,m,ncluster,ns 
      double precision alpha,denom,det,dnrm
      double precision meanuc1,meanuc2
      double precision meanub1,meanub2
      double precision meanuw1,meanuw2
      double precision sdb1,sdb2
      double precision sdk1,sdk2
      double precision tmp1,tpi,work1,work2,work3

c++++ Set the random number generator
      call setall(seed1,seed2)
      
      tpi=6.283185307179586476925286766559d0
      work1=-(dble(q)*log(tpi))       
      

c+++++Algorithm
      do i=1,nsave

c+++++++ check if the user has requested an interrupt
         call rchkusr()
      
         alpha=alphavec(i)

         sigmak(1,1)=sigmakmat(i,1)
         sigmak(1,2)=sigmakmat(i,2)
         sigmak(2,1)=sigmakmat(i,2)
         sigmak(2,2)=sigmakmat(i,3)
         sdk1=sqrt(sigmak(1,1))
         sdk2=sqrt(sigmak(2,2))
         
         meanb(1)=mubmat(i,1)
         meanb(2)=mubmat(i,2)
         meanub1=meanb(1)           
         meanub2=meanb(2)           
         
         sigmab(1,1)=sigmabmat(i,1)
         sigmab(1,2)=sigmabmat(i,2)
         sigmab(2,1)=sigmabmat(i,2)
         sigmab(2,2)=sigmabmat(i,3)
         sdb1=sqrt(sigmab(1,1))
         sdb2=sqrt(sigmab(2,2))

         ncluster=nclusvec(i)
         
         call rmvnorm(q,meanb,sigmab,workmhr,workvr,meanc)
         meanuc1=meanc(1)           
         meanuc2=meanc(2)           

         call inversedet(sigmak,q,iflagr,det) 
         work2=det

         denom=dble(nsubject)+alpha

c         call intpr("ncluster",-1,ncluster,1) 
c         call dblepr("meanb",-1,meanb,q) 
c         call dblepr("meanc",-1,meanc,q) 
c         call dblepr("sigmak",-1,sigmak,q*q) 
c         call dblepr("sigmab",-1,sigmab,q*q) 
         
c++++++++Bivariate

         do ii=1,ngrid1 
            theta(1)=grid1(ii)
            do jj=1,ngrid2
               theta(2)=grid2(jj)

c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

               tmp1=0.d0
               counter=0
               do k=1,ncluster
                  do l=1,q
                     counter=counter+1
                     meanw(l)=musave(i,counter)
                  end do
                  
                  work3=0.d0
                  do l=1,q
                     do m=1,q
                        work3=work3+(theta(l)-meanw(l))*
     &                               sigmak(l,m)*
     &                              (theta(m)-meanw(m))
                     end do
                  end do
                 
                  ns=clustsave(i,k)
               
                  tmp1=tmp1+(dble(ns)/denom)*
     &                      exp((work1-work2-work3)/2.d0)  
               end do

               work3=0.d0
               do l=1,q
                  do m=1,q
                     work3=work3+(theta(l)-meanc(l))*
     &                            sigmak(l,m)*
     &                           (theta(m)-meanc(m))
                  end do
               end do
               tmp1=tmp1+(alpha/denom)*exp((work1-work2-work3)/2.d0)  
               
               fs(ii,jj)=fs(ii,jj)+tmp1
            end do
         end do


c++++++++Univariates

         do j=1,ngrid1
c++++++++++ check if the user has requested an interrupt
            call rchkusr()
         
            tmp1=0.d0
            counter=0 
            do k=1,ncluster
               ns=clustsave(i,k)
               do l=1,q
                  counter=counter+1
                  meanw(l)=musave(i,counter)
               end do
               meanuw1=meanw(1) 

               tmp1=tmp1+(dble(ns)/denom)*dnrm(grid1(j),meanuw1,sdk1,0)
            end do
            tmp1=tmp1+(alpha/denom)*dnrm(grid1(j),meanuc1,sdk1,0)
            f1(j)=f1(j)+tmp1
         end do

         do j=1,ngrid2
c++++++++++ check if the user has requested an interrupt
            call rchkusr()
        
            tmp1=0.d0
            counter=0 
            do k=1,ncluster
               ns=clustsave(i,k)
               do l=1,q
                  counter=counter+1
                  meanw(l)=musave(i,counter)
               end do
               meanuw2=meanw(2) 

               tmp1=tmp1+(dble(ns)/denom)*dnrm(grid2(j),meanuw2,sdk2,0)
            end do
            tmp1=tmp1+(alpha/denom)*dnrm(grid2(j),meanuc2,sdk2,0)
            f2(j)=f2(j)+tmp1
         end do
       
      end do
      
      do i=1,ngrid1
         f1(i)=f1(i)/dble(nsave) 
         do j=1,ngrid2
            fs(i,j)=fs(i,j)/dble(nsave) 
         end do   
      end do
      
      do i=1,ngrid2
         f2(i)=f2(i)/dble(nsave)       
      end do
      return
      end
      
            
