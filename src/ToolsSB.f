c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR STICK-BREAKING
c=======================================================================                  
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
      double precision alpha(nr),beta(nr)

c+++++External working space
      double precision workv(nr+1)

c+++++Output
      double precision w(nr),v(nr)

c+++++Internal working space
      integer i 
      double precision rbeta
      double precision tmp1,tmp2

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

c=======================================================================
      subroutine dpweightsimbl(n,ccluster,alpha,workv,w,v)
c=======================================================================
c     Subroutine to compute the DP weights up to
c     finite level n.
c
c     A.J.V. 2007 
c=======================================================================
      implicit none

c+++++Input
      integer n
      integer ccluster(n)
      double precision alpha

c+++++External working space
      double precision workv(n+1)

c+++++Output
      double precision w(n),v(n)

c+++++Internal working space
      integer i,j,ns 
      double precision rbeta
      double precision tmp1,tmp2
      double precision tmp3
      double precision tmp5,tmp6

c+++++Algorihtm

      workv(1)=0.d0
      tmp1=0.d0
      do i=1,n-1
         tmp3=dble(ccluster(i))
         ns=0
         do j=i+1,n
            ns=ns+ccluster(j)
         end do
         tmp5=1.d0+tmp3
         tmp6=alpha+dble(ns) 

         v(i)=rbeta(tmp5,tmp6) 
         tmp2=1.d0-v(i)  
         tmp1=tmp1-dlog(tmp2) 
         workv(i+1)=tmp1 
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


c=======================================================================
      subroutine sickbreak(n,v,workv,w)
c=======================================================================
c     Subroutine to compute the stick breaking weights.
c=======================================================================
      implicit none

c+++++Input
      integer n
      double precision v(n)

c+++++External working space
      double precision workv(n+1)

c+++++Output
      double precision w(n)

c+++++Internal working space
      integer i 
      double precision tmp1,tmp2

c+++++Algorihtm
    
c      w(1)=v(1)
c      tmp1=1.d0-v(1)
c      do i=2,n
c	 w(i)=tmp1*v(i)
c         tmp1=tmp1*(1.d0-v(i))	
c      end do


      workv(1)=0.d0
      tmp1=0.d0
      do i=1,n
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


 
