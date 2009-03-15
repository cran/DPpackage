c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR STICK-BREAKING
c=======================================================================                  
c=======================================================================                  

c=======================================================================
      subroutine sickbreaksim(nr,n,alpha,beta,workv,w,v)
c=======================================================================
c     Subroutine to compute the stick breaking weights up to
c     finite level n.
c
c     A.J.V. 2007 
c=======================================================================
      implicit none

c+++++Input
      integer nr,n
      real*8 alpha(nr),beta(nr)

c+++++External working space
      real*8 workv(nr+1)

c+++++Output
      real*8 w(nr),v(nr)

c+++++Internal working space
      integer i 
      real*8 rbeta
      real*8 tmp1,tmp2

c+++++Algorihtm

      workv(1)=0.d0
      tmp1=0.d0
      do i=1,n
         v(i)=rbeta(alpha(i),beta(i)) 
         tmp2=1.d0-v(i)  
         if(i.ne.n)then
            tmp1=tmp1-dlog(tmp2) 
            workv(i+1)=tmp1 
         end if
      end do 
      workv(n+1)=1.0D300

      tmp1=0.d0
      do i=1,n
         w(i)=dexp(-workv(i))-dexp(-workv(i+1))
         tmp1=tmp1+w(i)
      end do  

      do i=1,n
         w(i)=w(i)/tmp1
      end do 

      return
      end

