
c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES AND FUNCTIONS TO EVALUATE THE DENSITY AND CDF OF RVs
c=======================================================================                  
c=======================================================================                  


c=======================================================================            	  
      double precision function cdfslogistic(x)
c=======================================================================            
c     This function evaluate the cdf of a standard logistic distribution.
c     A.J.V., 2005
      implicit none 
      real*8 x
      
c      if(x.lt.-50.d0)then
c         cdfslogistic=0.d0
c        else if(x.gt.50.d0)then
c         cdfslogistic=1.d0
c        else 
         cdfslogistic=exp(x-log(1+exp(x)))
c      end if   
      return
      end
     
c=======================================================================            	  
      double precision function invcdfslogistic(p)
c=======================================================================            
c     This function evaluate the cdf of a standard logistic distribution.
c     A.J.V., 2005
      implicit none 
      real*8 p
      invcdfslogistic=log(p)-log(1.d0-p)
      return
      end

c=======================================================================            
      real function ppnda(p)
c=======================================================================            
c     This function calculates the inverse normal distribution function
c     usin the algorithm of Beasley and Springer (1977)
c     A.J.V., 2005
      implicit none
      integer ifault
      real zero,split,half,one
      real a0,a1,a2,a3,b1,b2,b3,b4,c0,c1,c2,c3,d1,d2
      real p,q,r
      zero=0.0e0
      half=0.5e0
      one=1.0e0
      split=0.42e0
      a0=2.50662823884e0
      a1=-18.61500062529e0
      a2=41.39119773534e0
      a3=-25.44106049637e0
      b1=-8.47351093090e0
      b2=23.08336743743e0
      b3=-21.06224101826e0
      b4=3.13082909833e0
      c0=-2.78718931138e0
      c1=-2.29796479134e0
      c2=4.85014127135e0
      c3=2.32121276858e0
      d1=3.54388924762e0
      d2=1.63706781897e0
      ifault=0
      q=p-half
      if (abs(q).gt.split) go to 101
      r=q*q
      ppnda=q*(((a3*r+a2)*r+a1)*r+a0)/
     +	  ((((b4*r+b3)*r+b2)*r+b1)*r+one)
      return
101   r=p
	  if (q.gt.zero) r=one-p
	  if (r.le.zero) go to 102
	  r=sqrt(-alog(r))
	  ppnda=(((c3*r+c2)*r+c1)*r+c0)/
     +	  ((d2*r+d1)*r+one)
	  if (q.lt.zero) ppnda=-ppnda
      return
102	  ifault=1
	  ppnda=99.9
      return
      end        

c=======================================================================
      subroutine dmvn(n,x,mu,sigma,eval,vv,a,sigma2,res,iflag)        
c=======================================================================
c     return the log of a multivariate normal density
c     A.J.V., 2005
      implicit none
      integer n,i,j
      real*8 mu(n),sigma(n,n),x(n),vv(n)
      real*8 a(n,n),sigma2(n,n),res(n),det,sse,eval
      integer iflag(n)
      real*8 work1,work2,work3,tpi
     
      work1=0.d0
      work2=0.d0
      work3=0.d0
      det=0.d0
      sse=0.d0
	  
      det=0.d0
      
      tpi=6.283185307179586476925286766559d0
       
      work1=-(dble(n)*log(tpi))
    
      do i=1,n
         do j=1,n
            a(i,j)=sigma(i,j)
         end do
      end do

      call invdet(a,n,sigma2,det,iflag,vv)

      work2=det
   
      do i=1,n
         res(i)=x(i)-mu(i)
      end do
   
      do i=1,n
         do j=1,n
            sse=sse+res(i)*sigma2(i,j)*res(j)          
         end do
      end do

      work3=sse
     
      eval=(work1-work2-work3)/2.d0
      
      return
      end

