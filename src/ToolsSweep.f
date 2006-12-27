c=======================================================================
      subroutine sweep(a,dim,m,k1,k2,ier)
c=======================================================================
c     Subroutine to sweep the mxm matrix a on its k1 st thru k2 th
c     diagonals. ier is 1 if a is singular. dim is the real dimension
c     of a.
c     
c     Typically in regression: sweep(a,dim,m+1,1,m,ier) 
c     The output is:
c                   (xtx)^(-1)     : (xtx)^(-1) xty
c                  -ytx (xtx)^(-1) : yty-ytx (xtx)^(-1) xty
c
c     A.J.V. 2007 
c=======================================================================
      implicit none
      integer dim,i,ier,j,k,k1,k2,m
      real*8 a(dim,dim),d

      ier=1
 
      do k=k1,k2
         if(abs(a(k,k)).lt.1.e-20) return
         d=1./a(k,k)
         a(k,k)=1.
         do i=1,m
            a(k,i)=d*a(k,i)
            if(i.ne.k) a(i,k)=-a(i,k)*d
         end do

         do i=1,m
            do j=1,m
               if((i.ne.k).and.(j.ne.k)) a(i,j)=a(i,j)+a(i,k)*a(k,j)/d
            end do
         end do      
       end do

       ier=0
       return
       end
