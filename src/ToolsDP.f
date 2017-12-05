c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR DP Models
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
      double precision bclus(nsubject,q),theta(q)

      integer ns,ii

      do i=1,q
         theta(i)=bclus(since,i)
      end do
      
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
      double precision bclus(nrec),theta

      integer ns,ii

      theta=bclus(since)
      
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

         bclus(i-1)=bclus(i)
         ccluster(i-1)=ccluster(i)
      end do
      
      ss(ind)=ncluster
      
      bclus(ncluster)=theta
      
      ccluster(ncluster)=1
      
      return
      end  
            
