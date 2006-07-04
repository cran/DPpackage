c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR MATRIX OPERATIONS
c=======================================================================                  
c=======================================================================                  

c=======================================================================
      integer function efind(x,maxend,endp,npoints)
c=======================================================================
c     this function finds the position of x in the sorted vector 
c     endp (of dimension maxend) using a binary search. The search is
c     performed in the first npoints elements.
c     A.J.V., 2006
      integer maxend,i,j,k
      real*8 x,endp(maxend)
      
      i=1
      j=npoints
      efind=0

1     continue
      
      k=(i+j)/2

      if(x.eq.endp(k))then
          efind=k
          return
c       value is higher than list position
        else if(x.gt.endp(k))then
          i=k+1
c       value is lower than list position
        else
          j=k-1
      end if
      go to 1          
      return
      end
      
c=======================================================================
      subroutine sortvec(n,vec,i1,i2)
c=======================================================================
c     this routine sort the vector vec of dimension n in the positions
c     i1,...,i2.
c
c     A.J.V., 2005
      implicit none
      integer n,i1,i2
      real*8 vec(n)
      
      if(i1.ge.i2)then
	call rexit("i1 >= i2 in subroutine sortvec")
      end if
      
      call qsort3(vec,i1,i2)
      
      return
      end

c=======================================================================
      subroutine invdet2(a,np,n,ainv,detlog,indx,vv)
c=======================================================================
c     this routine returns the inverse of A in AINV and the log of the
c     determinate of A in DETLOG. N.B: log of determinate only exists
c     (as a real number) if A is a postive definate matrix.  This routine
c     assume that the matrix is a covariance matrix, and therefore it is
c     the dabs of a(j,j) is used.
c     This is from numerical recipies, pages 38 and 39. (The fortran
c     version).
c     This program calls the lu decomposition routines LUDCMP and LUBKSB,
c     see numerical recipies pages 35-37.
c
c     This is a modified version where np es the dimension of the objects
c     and n es the used dimension
c
c     A.J.V., 2005
      implicit none
      integer i,j,n,np,indx(np)
      double precision a(np,np),ainv(np,np),vv(np),detlog
c
c  set up identity matrix
c
      do 6 i=1,n
	  do 5 j=1,n
	    ainv(i,j) = 0.0d0
 5	    continue
	  ainv(i,i) = 1.0d0
 6	  continue

c
c  decompose the matrix
c

      call dludcmp2(a,n,np,indx,detlog,vv)
c
c  calculate the determinate
c
      detlog = 0.d0
      do 11 j=1,n
	    detlog = detlog + dlog(dabs(a(j,j)))
 11	    continue
c
c  find inverse by columns
c
      do 13 j=1,n
	    call dlubksb2(a,n,np,indx,ainv(1,j))
 13	    continue
      return
      end

c=======================================================================
      subroutine invdet(a,n,ainv,detlog,indx,vv)
c=======================================================================
c     this routine returns the inverse of A in AINV and the log of the
c     determinate of A in DETLOG. N.B: log of determinate only exists
c     (as a real number) if A is a postive definate matrix.  This routine
c     assume that the matrix is a covariance matrix, and therefore it is
c     the dabs of a(j,j) is used.
c     This is from numerical recipies, pages 38 and 39. (The fortran
c     version).
c     This program calls the lu decomposition routines LUDCMP and LUBKSB,
c     see numerical recipies pages 35-37.
c     A.J.V., 2005
      implicit none   
      integer i,j,n,indx(n)
      real*8 a(n,n),ainv(n,n),vv(n),detlog
c
c  set up identity matrix
c
      do 6 i=1,n
	  do 5 j=1,n
	    ainv(i,j) = 0.0d0
 5	    continue
	  ainv(i,i) = 1.0d0
 6	  continue
c
c  decompose the matrix
c
      call dludcmp(a,n,indx,detlog,vv)
c
c  calculate the determinate
c
      detlog = 0.d0
      do 11 j=1,n
	    detlog = detlog + dlog(dabs(a(j,j)))
 11	    continue
c
c  find inverse by columns
c
      do 13 j=1,n
	    call dlubksb(a,n,indx,ainv(1,j))
 13	    continue
      return
      end


c=======================================================================
      subroutine dludcmp2(a,n,np,indx,d,vv)
c=======================================================================
c     This is copied from Press, et. al. {\em Numerical Recipes:
c     Fortran version} pg 35-36
c
c     Given an N x N matrix A, with physical dimension NP, this routine
c     replaces it by the LU decomposition of a rowwise permutation of
c     itself. A and N are input. A is output, arranged as in equation
c     (2.3.14, in Press et. al., pg 34.); INDX is an output vector which
c     records the row permutation effected by the partial pivoting; D is
c     output as +/-1.d0 depending on wheter the number of row interchanges
c     was even or odd, respectively.  This routine is used in combination
c     with DLUBDSB to solve linear equatoin or invert a matrix.
c
c     This is a modified version where np es the dimension of the objects
c     and n es the used dimension
c
c     A.J.V., 2005
      implicit none 
      integer i,j,k,n,np
      real*8 tiny
      parameter (tiny=1.0d-22)
      real*8 a(np,np),vv(np),d,aamax,sum,dum
      integer indx(np),iimax

      iimax=1
      d=1.d0
      do 12 i=1,n
	    aamax=0.d0
	    do 11 j=1,n
	       if(dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
 11	       continue
	       if (aamax.eq.0.d0) then
	          call rexit("matrix is not pd in dludcmp2 subroutine")
               end if
	    vv(i)=1.d0/aamax
 12	    continue
	  do 19 j=1,n
	     do 14 i=1,j-1
	        sum=a(i,j)
	        do 13 k=1,i-1
		       sum=sum-a(i,k)*a(k,j)
 13		       continue
	         a(i,j)=sum
 14	         continue
	     aamax = 0.d0
	     do 16 i=j,n
	        sum=a(i,j)
	        do 15 k=1,j-1
		       sum = sum - a(i,k)*a(k,j)
 15		       continue
	        a(i,j)=sum
	        dum=vv(i)*dabs(sum)
	        if(dum.ge.aamax)then
		       iimax=i
		       aamax=dum
		       endif
 16	        continue
	    if(j.ne.iimax)then
	       do 17 k=1,n
		      dum=a(iimax,k)
		      a(iimax,k)=a(j,k)
		      a(j,k)=dum
 17		      continue
	       d=(-1.d0)*d
	       vv(iimax)=vv(j)
	       endif
	    indx(j)=iimax
	    if(a(j,j).eq.0.d0)a(j,j)=tiny
	    if(j.ne.n)then
	       dum=1.d0/a(j,j)
	       do 18 i=j+1,n
		      a(i,j)=a(i,j)*dum
 18		      continue
	       endif
 19	    continue
	  return
	  end

c=======================================================================
      subroutine dludcmp(a,n,indx,d,vv)
c=======================================================================
c     This is copied from Press, et. al. {\em Numerical Recipes:
c     Fortran version} pg 35-36
c
c     Given an N x N matrix A, with physical dimension NP, this routine
c     replaces it by the LU decomposition of a rowwise permutation of
c     itself. A and N are input. A is output, arranged as in equation
c     (2.3.14, in Press et. al., pg 34.);	INDX is an output vector which
c     records the row permutation effected by the partial pivoting; D is
c     output as +/-1.d0 depending on wheter the number of row interchanges
c     was even or odd, respectively.  This routine is used in combination
c     with DLUBDSB to solve linear equatoin or invert a matrix.
c     A.J.V., 2005
      implicit none 
      integer i,j,k,n
      real*8 tiny
      parameter (tiny=1.0d-22)
      real*8 a(n,n), vv(n),d,aamax,sum,dum
      integer indx(n),iimax
      
      iimax=1
      d=1.d0
      do 12 i=1,n
	    aamax=0.d0
	    do 11 j=1,n
	       if(dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
 11	       continue
	       if (aamax.eq.0.d0) then
	          call rexit("matrix is not pd in dludcmp subroutine")
               end if
	    vv(i)=1.d0/aamax
 12	    continue
	  do 19 j=1,n
	     do 14 i=1,j-1
	        sum=a(i,j)
	        do 13 k=1,i-1
		       sum=sum-a(i,k)*a(k,j)
 13		       continue
	         a(i,j)=sum
 14	         continue
	     aamax = 0.d0
	     do 16 i=j,n
	        sum=a(i,j)
	        do 15 k=1,j-1
		       sum = sum - a(i,k)*a(k,j)
 15		       continue
	        a(i,j)=sum
	        dum=vv(i)*dabs(sum)
	        if(dum.ge.aamax)then
		       iimax=i
		       aamax=dum
		       endif
 16	        continue
	    if(j.ne.iimax)then
	       do 17 k=1,n
		      dum=a(iimax,k)
		      a(iimax,k)=a(j,k)
		      a(j,k)=dum
 17		      continue
	       d=(-1.d0)*d
	       vv(iimax)=vv(j)
	       endif
	    indx(j)=iimax
	    if(a(j,j).eq.0.d0)a(j,j)=tiny
	    if(j.ne.n)then
	       dum=1.d0/a(j,j)
	       do 18 i=j+1,n
		      a(i,j)=a(i,j)*dum
 18		      continue
	       endif
 19	    continue
	  return
	  end


c=======================================================================
       subroutine dlubksb2(a,n,np,indx,b)
c=======================================================================
c      This is copied from Press, et. al. {\em Numerical Recipes:
c      Fortran version} pg 35-36
c
c      Solves the set of N linear equations A*X=B.	
c      Here A is input, not as
c      the matrix A but rather as its LU decomposition, determined by the
c      routine LUDCMP. INDX is input as the permutation vector returned by
c      LUDCMP. B is input as the right-hand side vector B, and returns with
c      the solution vector X. A,N,NP and INDX are not modified by this
c      routine and can be left in place for successive calls with different
c      right hand sides B. This routine tades into account the possibility
c      that B will begin with many zero elements, so it is efficient for use
c      in matrix inversion.
c
c     This is a modified version where np es the dimension of the objects
c     and n es the used dimension
c
c      A.J.V., 2005
       implicit none
       integer n,np,i,j,ii,ll
       real*8 a(np,np), b(np), sum
       integer indx(np)

       ii=0
        do 12 i=1,n
	    ll = indx(i)
	    sum=b(ll)
	    b(ll)=b(i)
	    if(ii.ne.0) then
	       do 11 j=ii,i-1
		  sum=sum-a(i,j)*b(j)
 11		  continue
	    else if (sum.ne.0.d0)then
	       ii=i
	    endif
	    b(i)=sum
 12	    continue
	  do 14 i=n,1,-1
	    sum=b(i)
	    do 13 j=i+1,n
	       sum=sum-a(i,j)*b(j)
 13	       continue
	    b(i)=sum/a(i,i)
 14	    continue
       return
       end


c=======================================================================
       subroutine dlubksb(a,n,indx,b)
c=======================================================================
c      This is copied from Press, et. al. {\em Numerical Recipes:
c      Fortran version} pg 35-36
c
c      Solves the set of N linear equations A*X=B.	
c      Here A is input, not as
c      the matrix A but rather as its LU decomposition, determined by the
c      routine LUDCMP. INDX is input as the permutation vector returned by
c      LUDCMP. B is input as the right-hand side vector B, and returns with
c      the solution vector X. A,N,NP and INDX are not modified by this
c      routine and can be left in place for successive calls with different
c      right hand sides B. This routine tades into account the possibility
c      that B will begin with many zero elements, so it is efficient for use
c      in matrix inversion.
c      A.J.V., 2005
       implicit none
       integer n,i,j,ii,ll
       real*8 a(n,n), b(n), sum
       integer indx(n)

       ii=0
        do 12 i=1,n
	    ll = indx(i)
	    sum=b(ll)
	    b(ll)=b(i)
	    if(ii.ne.0) then
	       do 11 j=ii,i-1
		  sum=sum-a(i,j)*b(j)
 11		  continue
	    else if (sum.ne.0.d0)then
	       ii=i
	    endif
	    b(i)=sum
 12	    continue
	  do 14 i=n,1,-1
	    sum=b(i)
	    do 13 j=i+1,n
	       sum=sum-a(i,j)*b(j)
 13	       continue
	    b(i)=sum/a(i,i)
 14	    continue
       return
       end


c=======================================================================      
      subroutine cholesky(n,a,l)
c=======================================================================      
c     Subroutine to do a Double precision Half stored CHOLesky
c     decomposition.  Calculate Cholesky decomposition of symmetric,
c     positive definite matrix A which is LOWER HALF-STORED.  The matrix
c     l is the output.
c     A.J.V., 2005
      implicit none
      integer n,i,ii,j,jj,k,kk
      real*8 a(n,n),l(n*(n+1)/2)
      real*8 aii,scal,rtemp

      jj=0
      do i=1,n
         do j=i,n
            jj=jj+1
            l(jj)=a(i,j)
         end do
      end do   
      
      ii = 1
      do i = 1,n-1
         aii = l(ii)
         if (aii.le.0.d0) then
            call rexit("matrix is not pd in chol subroutine")
            return
         end if
         aii = sqrt(aii)
         l(ii) = aii
         scal = 1.d0/aii
         do j = ii+1,ii+n-i
            l(j) = scal*l(j)
         end do

         jj = 1
         do j = ii+1,ii+n-i
            if (l(j).ne.0.d0) then
               rtemp = -1.d0 * l(j)
               kk = ii + jj + n - i
               do k = j,ii+n-i
                  l(kk) = l(kk) + l(k)*rtemp
                  kk = kk + 1
               end do
            end if
            jj = jj + n - i - j + ii + 1
         end do
         ii = ii + n - i + 1
      end do
      aii = l(ii)
      if (aii.le.0.d0) then
            call rexit("matrix is not pd in chol subroutine")
          return 
      end if
      aii = sqrt(aii)
      l(ii) = aii
      return
      end

