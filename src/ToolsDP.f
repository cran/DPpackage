c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR DP Models
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
      subroutine relabel(ind,since,nsubject,q,ncluster,cstrt,ccluster,
     &                   ss,bclus,theta)
c=======================================================================
c     relabel the clusters after the elimination of one of them
c     A.J.V., 2005
      implicit none
      integer i,j,ind,since,nsubject,q,ncluster,ccluster(nsubject)
      integer ss(nsubject),cstrt(nsubject,nsubject)
      real*8 bclus(nsubject,q),theta(q)

      do i=1,q
         theta(i)=bclus(since,i)
      end do
      
      do i=since+1,ncluster
         do j=1,nsubject
            if(ss(j).eq.i)then
               ss(j)=i-1 
            end if
         end do
         do j=1,ccluster(i)
            cstrt(i-1,j)=cstrt(i,j) 
         end do
         do j=1,q
            bclus(i-1,j)=bclus(i,j)
         end do
         ccluster(i-1)=ccluster(i)
      end do
      
      ss(ind)=ncluster
      
      do i=1,q
         bclus(ncluster,i)=theta(i)
      end do
      ccluster(ncluster)=1
      
      return
      end  


c=======================================================================      
      subroutine relabelmeta(ind,since,nrec,ncluster,cstrt,ccluster,
     &                       ss,bclus)
c=======================================================================
c     relabel the clusters after the elimination of one of them
c     A.J.V., 2005
      implicit none
      integer i,j,ind,since,nrec,ncluster,ccluster(nrec)
      integer ss(nrec),cstrt(nrec,nrec)
      real*8 bclus(nrec),theta

      theta=bclus(since)
      
      do i=since+1,ncluster
         do j=1,nrec
            if(ss(j).eq.i)then
               ss(j)=i-1 
            end if
         end do
         do j=1,ccluster(i)
            cstrt(i-1,j)=cstrt(i,j) 
         end do
         bclus(i-1)=bclus(i)
         ccluster(i-1)=ccluster(i)
      end do
      
      ss(ind)=ncluster
      
      bclus(ncluster)=theta
      
      ccluster(ncluster)=1
      
      return
      end  
            