c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR THE COMPUTAION OF POLYNOMIALS
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
      integer function iceil(x)
c=======================================================================
c     Function to get ceil(x).
c
c     A.J.V. 2007 
c=======================================================================
      implicit none 
      double precision x
      iceil=int(x)
      if(x.le.0.0)return
      if(float(iceil).ne.x)iceil=iceil+1
      return
      end

c=======================================================================
       subroutine legendrepld(n,x,pn,pd)
c=======================================================================
c      compute Legendre polynomials pn(x) and their derivatives pn'(x).
c      In the input x is the argument, n is the degree of pn(x) 
c      ( n = 0,1,...). The output is pn(n) and pd(n)=pn'(x)
c
c      Alejandro Jara, 2007 
c=======================================================================
       implicit none

c+++++ Input
       integer n
       double precision x,pn,pd
       dimension pn(0:n),pd(0:n)

c+++++ Working
       integer k
       double precision p0,p1,pf

c+++++ Algorithm
       pn(0)=1.0d0
       pn(1)=x
       pd(0)=0.0d0
       pd(1)=1.0d0
       p0=1.d0
       p1=x
       
       do k=2,n
           pf=(2.0d0*k-1.0d0)/k*X*p1-(k-1.0d0)/k*p0
           pn(k)=pf
           if(dabs(x).eq.1.0d0)then
              pd(k)=0.5d0*x**(k+1)*k*(k+1.0d0)
            else
              pd(k)=k*(p1-x*pf)/(1.0d0-x*x)
           end if
           p0=p1
           p1=pf
       end do
       return
       end

c=======================================================================
       subroutine legendres(n,x,pn)
c=======================================================================
c      compute 'normalized' Legendre polynomials pn(x) 
c      In the input x is the argument, n is the degree of pn(x) 
c      ( n = 0,1,...).
c
c      Alejandro Jara, 2007 
c=======================================================================
       implicit none

c+++++ Input
       integer n
       double precision x,pn
       dimension pn(0:n)

c+++++ Working
       integer k
       double precision p0,p1,pf

c+++++ Algorithm
       pn(0)=1.0d0
       pn(1)=x
       p0=1.d0
       p1=x
       
       do k=2,n
           pf=(2.0d0*k-1.0d0)/k*X*p1-(k-1.0d0)/k*p0
           pn(k)=pf
           p0=p1
           p1=pf
       end do
       
       do k=0,n
          pn(k)=dsqrt((2.0d0*k+1.d0)/2.0d0)*pn(k)
       end do
       
       return
       end


c=======================================================================
       subroutine legendremat(nrec,x,n,pn,z)
c=======================================================================
c      compute the design matrix of 'normalized' Legendre polynomials of
c      degree n, giving a vector of covariates x(nrec).
c
c      Alejandro Jara, 2007 
c=======================================================================
       implicit none

c+++++ Input
       integer nrec,n
       double precision x(nrec),z(nrec,n+1),pn
       dimension pn(0:n)

c+++++ Working
       integer i,j
       double precision xmin,xmax,zwork          

c+++++ Algorithm
   
       xmin=x(1)
       xmax=x(1)
       do i=2,nrec
          xmin=min(x(i),xmin)          
          xmax=max(x(i),xmax)       
       end do


       do i=1,nrec
          zwork=-1.0d0+2.0d0*((x(i)-xmin)/(xmax-xmin))
          call legendres(n,zwork,pn)
          do j=0,n
             z(i,j+1)=pn(j)
          end do
       end do
      
       return
       end

