c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR THE COMPUTAION OF POLYNOMIALS
c=======================================================================                  
c=======================================================================                  

c=======================================================================
       subroutine jcomponentbd(y,k,j)
c=======================================================================
c      return the component j in {1,...,k} corresponding to the y
c      in [0,1] random variable in a Bernstein-Dirichlet prior.
c
c      Alejandro Jara, 2007 
c=======================================================================
       implicit none
       integer k,j
       real*8 y
       j=int(y*dble(k))+1
       
       if(y.gt.dble(j)/dble(k))then
          call dblepr("y",-1,y,1)
          call intpr("k",-1,k,1)
          call intpr("j",-1,j,1)          
          call rexit("Error in ´jcomponent´")      
       end if

       if(y.lt.dble(j-1)/dble(k))then
          call dblepr("y",-1,y,1)
          call intpr("k",-1,k,1)
          call intpr("j",-1,j,1)                    
          call rexit("Error in ´jcomponent´")      
       end if
      
       return
       end
       
c=======================================================================
       subroutine baseevalbd(x,k,a0,b0,eval)
c=======================================================================
c      evaluates the Bernstein polinomial at the 
c      baseline distribution, Beta(a0,b0), in a Bernstein-Dirichlet 
c      prior.
c
c      Alejandro Jara, 2007 
c=======================================================================
       implicit none
       integer i,k
       real*8 a0,b0,cdfbetas,dbet,x,eval

       eval=0.d0
       do i=1,k
          eval=eval+
     &        (cdfbetas(dble(i)/dble(k),a0,b0,1,0)-
     &         cdfbetas(dble(i-1)/dble(k),a0,b0,1,0))*
     &         dbet(x,dble(i),dble(k-i+1),0) 
       end do 
       return
       end
       
c=======================================================================
       subroutine clustevalbd(x,k,y,eval)
c=======================================================================
c      evaluates the cluster contribution for the cluster
c      defined by y in a Bernstein-Dirichlet prior.
c
c      Alejandro Jara, 2007 
c=======================================================================
       implicit none
       integer j,k
       real*8 dbet,x,y,eval

       call jcomponentbd(y,k,j)       
       eval=dbet(x,dble(j),dble(k-j+1),0) 

       return
       end

c=======================================================================
       subroutine samplekbd(nrec,x,y,prob,kmax,k)
c=======================================================================
c      generate k from its conditional distribution in a 
c      Bernstein-Dirichlet prior. This function assumes that the prior 
c      for k is a uniform U(1,kmax).
c
c      Alejandro Jara, 2007.
c=======================================================================
       implicit none
       integer i,j,kmax,k,nrec
       real*8 eval,prob(kmax),y(nrec),x(nrec),tmp1

       do i=1,kmax
          tmp1=0.d0         
          do j=1,nrec
             call clustevalbd(x(j),i,y(j),eval)
             tmp1=tmp1+log(eval)
          end do
          prob(i)=exp(tmp1)*(1.d0/dble(kmax))
       end do

       call simdisc(prob,kmax,kmax,k)

       return
       end


c=======================================================================
       subroutine sampleybd(x,kmax,a0,b0,k,y)
c=======================================================================
c      generate y from the baseline in a 
c      Bernstein-Dirichlet prior. 
c
c      Alejandro Jara, 2007.
c=======================================================================
       implicit none
       integer i,j,kmax,k,status
       real*8 a0,b0,bound,dbet,prob(kmax),tmp1,tmp2,tmp3,x
       real*8 y,y2
       real*8 tt1,tt2,tt3,tt4
       real runif

       do i=1,k
          prob(i)=dbet(x,dble(i),dble(k-i+1),0)
       end do

       call simdisc(prob,kmax,k,j)

       if(a0.eq.1.d0.and.b0.eq.1)then
          y=(dble(j-1)+dble(runif()))/dble(k)
        else
          tt3=dble(j-1)/dble(k)
          tt4=1.d0-dble(j-1)/dble(k)
          call cdfbet(1,tt1,tt2,tt3,tt4,a0,b0,status,bound)
          if(status.ne.0)then
             call rexit("Error in ´sampleybd´")      
          end if
          tmp1=tt1
       
          tt3=dble(j)/dble(k)
          tt4=1.d0-dble(j)/dble(k)
          call cdfbet(1,tt1,tt2,tt3,tt4,a0,b0,status,bound)
          if(status.ne.0)then
             call rexit("Error in ´sampleybd´")      
          end if
          tmp2=tt1
 
          tmp3=tmp1+dble(runif())*(tmp2-tmp1) 
       
          call cdfbet(2,tmp3,1.d0-tmp3,y,y2,a0,b0,status,bound)
          if(status.ne.0)then
             call rexit("Error in ´sampleybd´")      
          end if

       end if
        
       if(y.gt.dble(j)/dble(k))then
          call rexit("Error in ´sampleybd´")      
       end if

       if(y.le.dble(j-1)/dble(k))then
          call rexit("Error in ´sampleybd´")      
       end if

       if(y.eq.0.d0)then
          call rexit("Error in ´sampleybd´ (0)")      
       end if


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
       real*8 x,pn,pd
       dimension pn(0:n),pd(0:n)

c+++++ Working
       integer k
       real*8 p0,p1,pf

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
       real*8 x,pn
       dimension pn(0:n)

c+++++ Working
       integer k
       real*8 p0,p1,pf

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
       real*8 x(nrec),z(nrec,n+1),pn
       dimension pn(0:n)

c+++++ Working
       integer i,j
       real*8 xmin,xmax,zwork          

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

