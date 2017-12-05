c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR MATRIX OPERATIONS
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

      subroutine sortqr(a,n,m,ix,iy,col,icol)
c sorts matrix of reals by icol columns given by col, in ascending
c order. matrix a is n x m, only first ix rows and iy columns are sorted
      integer n,m,ix,iy,icol
      integer i,j,j1,j2,k,l,r,s,col(icol)
      double precision a(n,m),x(20),a1,a2
      integer stack(50,2)
      s=1
      stack(1,1)=1
      stack(1,2)=ix
10    l=stack(s,1)
      r=stack(s,2)
      s=s-1
20    i=l
      j=r
      j1=(l+r)/2
      do 25 k=1,icol
25    x(k)=a(j1,col(k))
30    do 35 k=1,icol
      a1=a(i,col(k))
      a2=x(k)
      if (a1.lt.a2) then
      i=i+1
      goto 30
      else if(a1.gt.a2) then
      goto 40
      endif
35    continue
40    do 45 k=1,icol
      a1=a(j,col(k))
      a2=x(k)
      if (a1.gt.a2) then
      j=j-1
      goto 40
      else if(a1.lt.a2) then
      goto 47
      endif
45    continue
47    if (i.le.j) then
      do 50 j2=1,iy
      a1=a(i,j2)
      a(i,j2)=a(j,j2)
50    a(j,j2)=a1
      i=i+1
      j=j-1
      endif
      if (i.le.j) goto 30
      if (i.lt.r) then
      s=s+1
      stack(s,1)=i
      stack(s,2)=r
      endif
      r=j
      if (l.lt.r) goto 20
      if (s.ne.0) goto 10
      return
      end
 

c=======================================================================
       double precision function inprod(x,y,n)
c=======================================================================
c      computes xty 
c
c      A.J.V. 2007
c=======================================================================
       implicit none
       integer i,n 
       double precision x(n),y(n)
       
       inprod=0.d0
       do i=1,n
          inprod=inprod+x(i)*y(i)
       end do   
       return
       end

c=======================================================================
      integer function efind(x,maxend,endp,npoints)
c=======================================================================
c     this function finds the position of x in the sorted vector 
c     endp (of dimension maxend) using a binary search. The search is
c     performed in the first npoints elements.
c     A.J.V., 2006
      integer maxend,i,j,k
      double precision x,endp(maxend)
      
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
      
c      if(efind.eq.0)then
c        call rexit("Error in efind")
c      end if
c      return
      end

c=======================================================================
      subroutine sortvecind(n,index,vec1,vec2,i1,i2)
c=======================================================================
c     this routine sort the vector vec of dimension n in the positions
c     i1,...,i2.
c
c     A.J.V., 2006
      implicit none
      integer n,i1,i2
      integer index(n)
      integer i
      double precision vec1(n),vec2(n),vwork(10000000)

      if(n.gt.10000000)then
         call rexit("increase dimension subroutine sortvecind")
      end if
      
      if(i1.ge.i2)then
         call rexit("i1 >= i2 in subroutine sortvecind")
      end if
      
      do i=i1,i2
         vwork(i)=vec2(i)
      end do
      
      call qsort4(vec1,index,i1,i2)
      
      do i=i1,i2
         vec2(i)=vwork(index(i))
      end do
      
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
      double precision vec(n)
      
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
 5       continue
         ainv(i,i) = 1.0d0
 6    continue

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
 11   continue
c
c  find inverse by columns
c
      do 13 j=1,n
         call dlubksb2(a,n,np,indx,ainv(1,j))
 13   continue
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
      double precision a(n,n),ainv(n,n),vv(n),detlog
c
c  set up identity matrix
c
      do 6 i=1,n
         do 5 j=1,n
            ainv(i,j) = 0.0d0
 5       continue
         ainv(i,i) = 1.0d0
 6    continue
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
 11   continue
c
c  find inverse by columns
c
      do 13 j=1,n
         call dlubksb(a,n,indx,ainv(1,j))
 13   continue
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
      double precision tiny
      parameter (tiny=1.0d-20)
      double precision a(np,np),vv(np),d,aamax,sum,dum
      integer indx(np),iimax

      iimax=1
      d=1.d0
      do 12 i=1,n
         aamax=0.d0
         do 11 j=1,n
            if(dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
 11      continue
         if(aamax.eq.0.d0) then
            call rexit("matrix is not pd in dludcmp2 subroutine")
         end if
         vv(i)=1.d0/aamax
 12   continue
      do 19 j=1,n
         do 14 i=1,j-1
            sum=a(i,j)
            do 13 k=1,i-1
               sum=sum-a(i,k)*a(k,j)
 13         continue
            a(i,j)=sum
 14      continue
         aamax = 0.d0
         do 16 i=j,n
            sum=a(i,j)
            do 15 k=1,j-1
               sum = sum - a(i,k)*a(k,j)
 15         continue
            a(i,j)=sum
            dum=vv(i)*dabs(sum)
            if(dum.ge.aamax)then
               iimax=i
               aamax=dum
            endif
 16      continue
         if(j.ne.iimax)then
            do 17 k=1,n
               dum=a(iimax,k)
               a(iimax,k)=a(j,k)
               a(j,k)=dum
 17         continue
            d=(-1.d0)*d
            vv(iimax)=vv(j)
         endif
         indx(j)=iimax
         if(a(j,j).eq.0.d0)a(j,j)=tiny
            if(j.ne.n)then
               dum=1.d0/a(j,j)
               do 18 i=j+1,n
                  a(i,j)=a(i,j)*dum
 18            continue
            endif
 19    continue
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
      double precision tiny
      parameter (tiny=1.0d-20)
      double precision a(n,n), vv(n),d,aamax,sum,dum
      integer indx(n),iimax
      
      iimax=1
      d=1.d0
      do 12 i=1,n
         aamax=0.d0
         do 11 j=1,n
            if(dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
 11      continue
         if(aamax.eq.0.d0) then
            call rexit("matrix is not pd in dludcmp subroutine")
         end if
         vv(i)=1.d0/aamax
 12      continue
         do 19 j=1,n
            do 14 i=1,j-1
               sum=a(i,j)
               do 13 k=1,i-1
                  sum=sum-a(i,k)*a(k,j)
 13            continue
               a(i,j)=sum
 14         continue
            aamax = 0.d0
            do 16 i=j,n
               sum=a(i,j)
               do 15 k=1,j-1
                  sum = sum - a(i,k)*a(k,j)
 15            continue
               a(i,j)=sum
               dum=vv(i)*dabs(sum)
               if(dum.ge.aamax)then
                  iimax=i
                  aamax=dum
               endif
 16         continue
            if(j.ne.iimax)then
               do 17 k=1,n
                  dum=a(iimax,k)
                  a(iimax,k)=a(j,k)
                  a(j,k)=dum
 17            continue
               d=(-1.d0)*d
               vv(iimax)=vv(j)
            end if
            indx(j)=iimax
            if(a(j,j).eq.0.d0)a(j,j)=tiny
            if(j.ne.n)then
               dum=1.d0/a(j,j)
               do 18 i=j+1,n
                  a(i,j)=a(i,j)*dum
 18            continue
            end if
 19   continue
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
       double precision a(np,np), b(np), sum
       integer indx(np)

       ii=0
       do 12 i=1,n
          ll = indx(i)
          sum=b(ll)
          b(ll)=b(i)
          if(ii.ne.0) then
             do 11 j=ii,i-1
                 sum=sum-a(i,j)*b(j)
 11          continue
          else if (sum.ne.0.d0)then
             ii=i
          endif
          b(i)=sum
 12    continue
       do 14 i=n,1,-1
          sum=b(i)
          do 13 j=i+1,n
             sum=sum-a(i,j)*b(j)
 13       continue
       b(i)=sum/a(i,i)
 14    continue
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
       double precision a(n,n), b(n), sum
       integer indx(n)

       ii=0
       do 12 i=1,n
          ll = indx(i)
          sum=b(ll)
          b(ll)=b(i)
          if(ii.ne.0) then
             do 11 j=ii,i-1
                sum=sum-a(i,j)*b(j)
 11          continue
           else if (sum.ne.0.d0)then
             ii=i
          endif
          b(i)=sum
 12    continue
       do 14 i=n,1,-1
          sum=b(i)
          do 13 j=i+1,n
             sum=sum-a(i,j)*b(j)
 13       continue
          b(i)=sum/a(i,i)
 14    continue
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
      double precision a(n,n),l(n*(n+1)/2)
      double precision aii,scal,rtemp

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


c=======================================================================      
      subroutine inverse(ainv,n,ipiv)
c=======================================================================      
c     A.J.V., 2007
      integer lwork,n,ipiv(n),info
      parameter(lwork=6400)
      double precision ainv(n,n),work(lwork)

      if(n.gt.100)then
         call rexit("error in inversion, increase n")         
      end if

      call dgetrf(n,n,ainv,n,ipiv,info)
      if(info.ne.0)then
         call rexit("error in inversion")         
      end if
      call dgetri(n,ainv,n,ipiv,work,lwork,info)      
      if(info.ne.0)then
         call rexit("error in inversion")         
      end if
      return
      end

    
c=======================================================================      
      subroutine inversedet(ainv,n,ipiv,detlog)
c=======================================================================      
c     A.J.V., 2007
      implicit none
      integer i,lwork,n,ipiv(n),info
      parameter(lwork=6400)
      double precision ainv(n,n),detlog,work(lwork)

      if(n.gt.100)then
         call rexit("error in inversion, increase n")         
      end if

      call dgetrf(n,n,ainv,n,ipiv,info)
      if(info.ne.0)then
         call rexit("error in inversion")         
      end if
      
      detlog = 0.d0
      do i=1,n
         detlog = detlog + dlog(dabs(ainv(i,i)))
      end do
      
      call dgetri(n,ainv,n,ipiv,work,lwork,info)      
      if(info.ne.0)then
         call rexit("error in inversion")         
      end if
      return
      end
      

C=======================================================================
      SUBROUTINE EIGENV(NM,N,A,D,E,Z)
C=======================================================================
C                   SUBROUTINE TO FIND EIGENVALUES AND VECTORS
C                   OF THE SYMMETRIC MATRIX A (FULL STORED)
C              NM,N - DECLARED AND ACTUAL DIMENSIONS OF A
C              D    - WILL CONTAIN EIGENVALUES OF A
C              E    - WORK VECTOR OF LENGTH N
C              Z    - WILL CONTAIN EIGENVECTORS OF A (AS COLUMNS),
C                     Z AND A MAY BE THE SAME TO SAVE SPACE
C                   THIS SUBROUTINE IS A COMBINATION OF SUBROUTINES
C                   TRED2 AND TQL2 FROM EISPACK, ARGONNE NAT. LAB.
C                   TRED2 REDUCES A TO SYMMETRIC TRIDIAGONAL FORM
C                   TQL2 GETS EIGENVALUES AND VECTORS BY QL METHOD
C
C             Obtained from PVR's EM-REML sire program
C
C    NOTE AJV: This will be used provisorily       
C-----------------------------------------------------------------------
      DOUBLE PRECISION A(NM,N),D(N),E(N),Z(NM,N)
      DOUBLE PRECISION C,C2,C3,DL1,EL1,F,G,H,HH,P,R
      DOUBLE PRECISION  SCALE,S,S2,TST1,TST2,PPP,RRR,SSS,TTT,UUU
      
      c3=0.d0
      s2=0.d0

      DO 10 I = 1, N
        DO 8 J = I, N
   8      Z(J,I) = A(J,I)
  10    D(I) = A(N,I)
      IF (N .EQ. 1) GO TO 51
      DO 30 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 2) GO TO 13
         DO 12 K = 1, L
  12     SCALE = SCALE + DABS(D(K))
         IF (SCALE .NE. 0.0D0) GO TO 14
  13     E(I) = D(L)
         DO 135 J = 1, L
            D(J) = Z(L,J)
            Z(I,J) = 0.0D0
  135       Z(J,I) = 0.0D0
         GO TO 29
  14     DO 15 K = 1, L
            D(K) = D(K) / SCALE
  15        H = H + D(K) * D(K)
         F = D(L)
         G = -DSIGN(DSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         DO 17 J = 1, L
  17       E(J) = 0.0D0
         DO 24 J = 1, L
            F = D(J)
            Z(J,I) = F
            G = E(J) + Z(J,J) * F
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 22
            DO 20 K = JP1, L
               G = G + Z(K,J) * D(K)
  20           E(K) = E(K) + Z(K,J) * F
  22        E(J) = G
  24     CONTINUE
         F = 0.0D0
         DO 245 J = 1, L
            E(J) = E(J) / H
  245       F = F + E(J) * D(J)
         HH = F / (H + H)
         DO 25 J = 1, L
  25       E(J) = E(J) - HH * D(J)
         DO 28 J = 1, L
            F = D(J)
            G = E(J)
            DO 26 K = J, L
  26          Z(K,J) = Z(K,J) - F * E(K) - G * D(K)
            D(J) = Z(L,J)
            Z(I,J) = 0.0D0
  28     CONTINUE
  29     D(I) = H
  30  CONTINUE
      DO 50 I = 2, N
         L = I - 1
         Z(N,L) = Z(L,L)
         Z(L,L) = 1.0D0
         H = D(I)
         IF (H .EQ. 0.0D0) GO TO 38
         DO 33 K = 1, L
  33     D(K) = Z(K,I) / H
         DO 36 J = 1, L
            G = 0.0D0
            DO 34 K = 1, L
  34        G = G + Z(K,I) * Z(K,J)
            DO 36 K = 1, L
               Z(K,J) = Z(K,J) - G * D(K)
  36     CONTINUE
  38     DO 40 K = 1, L
  40     Z(K,I) = 0.0D0
  50  CONTINUE
  51  DO 52 I = 1, N
         D(I) = Z(N,I)
  52     Z(N,I) = 0.0D0
      Z(N,N) = 1.0D0
      E(1) = 0.0D0

      IF (N .EQ. 1) GO TO 1001
      DO 100 I = 2, N
  100 E(I-1) = E(I)
      F = 0.0D0
      TST1 = 0.0D0
      E(N) = 0.0D0
      DO 240 L = 1, N
         J = 0
         H = DABS(D(L)) + DABS(E(L))
         IF (TST1 .LT. H) TST1 = H
         DO 110 M = L, N
            TST2 = TST1 + DABS(E(M))
            IF (TST2 .EQ. TST1) GO TO 120
  110    CONTINUE
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * E(L))
      PPP = DMAX1(DABS(P),DABS(1.0D0))
      IF (PPP .EQ. 0.0D0) GO TO 137
      RRR = (DMIN1(DABS(P),DABS(1.0D0))/PPP)**2
  136    TTT = 4.0D0 + RRR
         IF (TTT .EQ. 4.0D0) GO TO 137
         SSS = RRR/TTT
         UUU = 1.0D0 + 2.0D0*SSS
         PPP = UUU*PPP
         RRR = (SSS/UUU)**2 * RRR
         GO TO 136
  137 R = PPP

         D(L) = E(L) / (P + DSIGN(R,P))
         D(L1) = E(L) * (P + DSIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 150
         DO 140 I = L2, N
  140    D(I) = D(I) - H
  150    F = F + H
         P = D(M)
         C = 1.0D0
         C2 = C
         EL1 = E(L1)
         S = 0.0D0
         MML = M - L
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P

      PPP = DMAX1(DABS(P),DABS(E(I)))
      IF (PPP .EQ. 0.0D0) GO TO 157
      RRR = (DMIN1(DABS(P),DABS(E(I)))/PPP)**2
  156    TTT = 4.0D0 + RRR
         IF (TTT .EQ. 4.0D0) GO TO 157
         SSS = RRR/TTT
         UUU = 1.0D0 + 2.0D0*SSS
         PPP = UUU*PPP
         RRR = (SSS/UUU)**2 * RRR
         GO TO 156
  157 R = PPP

            E(I+1) = S * R
            S = E(I) / R
            C = P / R
            P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
  180          Z(K,I) = C * Z(K,I) - S * H
  200    CONTINUE
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         TST2 = TST1 + DABS(E(L))
         IF (TST2 .GT. TST1) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE

      GO TO 1001
 1000 call rexit("PROBLEMS IN SUBROUTINE EIGVEC")         
 
 1001 RETURN
      END


c=======================================================================      
      subroutine cholesky2(nr,n,a,l)
c=======================================================================      
c     Subroutine to do a Double precision Half stored CHOLesky
c     decomposition.  Calculate Cholesky decomposition of symmetric,
c     positive definite matrix A which is LOWER HALF-STORED.  The matrix
c     l is the output.
c     A.J.V., 2007
      implicit none
      integer nr,n,i,ii,j,jj,k,kk
      double precision a(nr,nr),l(nr*(nr+1)/2)
      double precision aii,scal,rtemp

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


c=======================================================================      
      subroutine cumsum(nr,n,vec,vecsum)
c=======================================================================      
c     Subroutine to compute the cumulative sum. 
c     A.J.V., 2007
      implicit none
      integer i,nr,n  
      double precision vec(nr),vecsum(nr)
      double precision tmp1

      tmp1=vec(1)
      tmp1=vecsum(1)

      do i=2,n
         tmp1=tmp1+vec(i)
         vecsum(i)=tmp1
      end do
      return
      end

c=======================================================================      
      subroutine permut(l,n,r)
c=======================================================================      
c     Subroutine to give the l-th row of a permutation matrix. 
c     A.J.V., 2007
      implicit integer (a-z)
      dimension r(n)

c+++++Init
      h=1
      do 5 k=1,n
           h=h*k
5     r(k)=k

c+++++Body
      do 10 k=1,n
         p=n-k+1
         h=h/p
         q=mod((l-1)/h,p)+1            
         if(q.eq.1) goto 10

c++++++++Rotate data
         s=r(k+q-1)
         do m=k+q-1,k+1,-1
            r(m)=r(m-1)
         end do 
         r(k)=s
10    continue

      return
      end


c=======================================================================      
      subroutine genpermmat(l,n,r,q)
c=======================================================================      
c     Subroutine to genetare the l-th permutation matrix. 
c     A.J.V., 2007
      implicit none

c+++++Input
      integer l,n

c+++++Output
      double precision q(n,n)

c+++++External working space
      integer r(n)   

c+++++Internal working space
      integer i,j

c+++++Body
      do i=1,n
         do j=1,n
            q(i,j)=0.d0
         end do
      end do

      call permut(l,n,r)
      do i=1,n
         q(i,r(i))=1.d0
      end do

      return
      end



