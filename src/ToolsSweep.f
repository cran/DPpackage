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
      integer dim,i,ier,j,k,k1,k2,m
      double precision a(dim,dim),d

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

c=======================================================================
      subroutine condmvn(ind,a,m,aw,perm)
c=======================================================================
c     computes the conditional covariance matrix 
c     A1.2 = A11 - A12*inv(A22)*A21 and the
c     matrix for conditional mean vector A12*inv(A22)
c
c     Note: perm = Q permutation matrix. The ind coordinate is
c           placed in the first entry by using Q*A*Q^T
c
c     Alejandro Jara, 2007 
c=======================================================================
      implicit none

c+++++Input
      integer ind,m
      double precision a(m,m)

c+++++Working
      integer maxm
      parameter(maxm=100)
      integer i,j,k,ier
      integer counter
      double precision aw2(maxm,maxm),tmp1

c+++++Output
      double precision aw(m,m),perm(m,m) 

c+++++Algorithm

      if(maxm.lt.m)then
         call rexit("Increase 'maxm' in 'condmvn'")
      end if   

      if(ind.eq.1)then  
         do i=1,m
            do j=1,m
               perm(i,j)=0.0d0  
            end do
         end do   
         do i=1,m
            perm(i,i)=1.d0
         end do
         do i=1,m
            do j=1,m
               aw(i,j)=a(i,j)    
            end do
         end do   
       else         
         do i=1,m
            do j=1,m
               perm(i,j)=0.0d0  
            end do
         end do   
         perm(1,ind)=1.0d0
         counter=1
         do i=1,m
            if(i.ne.ind)then
               counter=counter+1
               perm(counter,i)=1.0d0
            end if
         end do
         do i=1,m
            do j=1,m
               tmp1=0.d0 
               do k=1,m
                  tmp1=tmp1+perm(i,k)*a(k,j)    
               end do
               aw2(i,j)=tmp1
            end do   
         end do   

         do i=1,m
            do j=1,m
               tmp1=0.d0 
               do k=1,m
                  tmp1=tmp1+aw2(i,k)*perm(j,k)   
               end do
               aw(i,j)=tmp1
            end do   
         end do   
      end if 

c      do i=1,m
c         do j=1,m
c            call dblepr("aw",-1,aw(i,j),1) 
c         end do
c      end do   

      call sweep(aw,m,m,2,m,ier)
      if(ier.eq.1)then
         call rexit("error in computing conditionals")
      end if   
     
      return
      end
      
      

