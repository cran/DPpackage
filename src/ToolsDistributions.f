
c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES AND FUNCTIONS TO EVALUATE THE DENSITY AND CDF OF RVs
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
      subroutine dgamma(y,mu,v,eval)        
c=======================================================================
c     return the log of a gamma distribution
c     A.J.V., 2006
      double precision y,mu,v,eval
      double precision dgamlog
      eval=(v-1)*log(y)-y*(v/mu)-v*log(mu/v)-dgamlog(v)
      return
      end


c=======================================================================
      subroutine dgamma2(y,alpha,beta,eval)        
c=======================================================================
c     return the log of a gamma distribution
c     A.J.V., 2006
      double precision y,alpha,beta,eval
      double precision dgamlog
      eval=alpha*log(beta)+(alpha-1.d0)*log(y) - beta*y - dgamlog(alpha)
      return
      end


c=======================================================================
      subroutine dgammai(y,alpha,beta,eval)        
c=======================================================================
c     return the log of a inverted gamma distribution
c     A.J.V., 2006
      double precision y,alpha,beta,eval
      double precision dgamlog
      eval=alpha*log(beta)-(alpha+1.d0)*log(y) - beta/y - dgamlog(alpha)
      return
      end


c=======================================================================            	  
      double precision function cdfslogistic(x)
c=======================================================================            
c     This function evaluate the cdf of a standard logistic distribution.
c     A.J.V., 2005
      implicit none 
      double precision x
      double precision cdflogis
      cdfslogistic=cdflogis(x,0.d0,1.d0,1,0)
      return
      end
     
c=======================================================================            	  
      double precision function invcdfslogistic(p)
c=======================================================================            
c     This function evaluate the cdf of a standard logistic distribution.
c     A.J.V., 2005
      implicit none 
      double precision p
      double precision invcdflogis
      invcdfslogistic=invcdflogis(p,0.d0,1.d0,1,0)
      return
      end



c=======================================================================            	  
      double precision function cdfbeta(x,alpha,beta)
c=======================================================================            
c     This function evaluate the cdf of a Beta(alpha,beta) distribution.
c     A.J.V., 2005
      implicit none 
      double precision x,alpha,beta
      double precision cdfbetas
      cdfbeta=cdfbetas(x,alpha,beta,1,0)
      return
      end

     
c=======================================================================            	  
      double precision function invcdfbeta(p,alpha,beta)
c=======================================================================            
c     This function evaluate the cdf of a Beta(alpha,beta) distribution.
c     A.J.V., 2005
      implicit none 
      double precision p,alpha,beta
      double precision invcdfbetas
      invcdfbeta=invcdfbetas(p,alpha,beta,1,0)
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
     +      ((((b4*r+b3)*r+b2)*r+b1)*r+one)
      return
101   r=p
      if (q.gt.zero) r=one-p
      if (r.le.zero) go to 102
      r=sqrt(-alog(r))
      ppnda=(((c3*r+c2)*r+c1)*r+c0)/
     +      ((d2*r+d1)*r+one)
      if (q.lt.zero) ppnda=-ppnda
      return
102   ifault=1
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
      double precision mu(n),sigma(n,n),x(n),vv(n)
      double precision a(n,n),sigma2(n,n),res(n),det,sse,eval
      integer iflag(n)
      double precision work1,work2,work3,tpi
     
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


c=======================================================================
      subroutine dmvn2(n,x,mu,sigma,eval,vv,a,sigma2,iflag)        
c=======================================================================
c     return the log of a multivariate normal density
c     A.J.V., 2005
      implicit none
      integer n,i,j
      double precision mu(n),sigma(n,n),x(n),vv(n)
      double precision a(n,n),sigma2(n,n),det,sse,eval
      integer iflag(n)
      double precision work1,work2,work3,tpi
     
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
         vv(i)=x(i)-mu(i)
      end do
   
      do i=1,n
         do j=1,n
            sse=sse+vv(i)*sigma2(i,j)*vv(j)          
         end do
      end do

      work3=sse
     
      eval=(work1-work2-work3)/2.d0
      
      return
      end


c=======================================================================
      subroutine dmvn3(nr,n,x,mu,sigma,eval,vv,a,sigma2,iflag)        
c=======================================================================
c     return the log of a multivariate normal density
c     A.J.V., 2005
      implicit none
      integer nr,n,i,j
      double precision mu(nr),sigma(nr,nr),x(nr),vv(nr)
      double precision a(nr,nr),sigma2(nr,nr),det,sse,eval
      integer iflag(nr)
      double precision work1,work2,work3,tpi
     
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

      call invdet2(a,nr,n,sigma2,det,iflag,vv)      
      
      work2=det
   
      do i=1,n
         vv(i)=x(i)-mu(i)
      end do
   
      do i=1,n
         do j=1,n
            sse=sse+vv(i)*sigma2(i,j)*vv(j)          
         end do
      end do

      work3=sse
     
      eval=(work1-work2-work3)/2.d0
      
      return
      end

c=======================================================================
      subroutine diwishart(kk,nu,sigma,tinv,workm1,workm2,workv,iflag,
     &                     eval)        
c=======================================================================
c     return the log of a inverted wishart density (proportional only)
c     W~InvWishart(nu,T) 
c     p(W) \propto |Tinv|^{nu/2} |Sigma|^{(nu+kk+1)/2} 
c                  exp(-0.5*tr(Tinv Sigmainv))
c     A.J.V., 2007

      implicit none
c-----input      
      integer kk,nu
      double precision sigma(kk,kk),tinv(kk,kk)

c-----output      
      double precision eval
      
c-----input-working       
      integer iflag(kk)      
      double precision workm1(kk,kk),workm2(kk,kk)
      double precision workv(kk)

c-----working
      integer i,j,k  
      double precision detlogt,detlogs
      double precision trace,tmp1

      do i=1,kk
         do j=1,kk
            workm1(i,j)=sigma(i,j)
         end do
      end do
      call invdet(workm1,kk,workm2,detlogs,iflag,workv)
      
      do i=1,kk
         do j=1,kk
            workm1(i,j)=0.d0
         end do
      end do
      
      do i=1,kk
         do j=1,kk
            tmp1=0.d0
            do k=1,kk
               tmp1=tmp1+tinv(i,k)*workm2(k,j)
            end do
            workm1(i,j)=tmp1
         end do
      end do
      
      trace=0.d0
      do i=1,kk
         trace=trace+workm1(i,i)
      end do

      do i=1,kk
         do j=1,kk
            workm1(i,j)=tinv(i,j)
         end do
      end do
      call invdet(workm1,kk,workm2,detlogt,iflag,workv)

       
      eval=  0.5d0*dble(nu)*detlogt
     &      -0.5d0*(dble(nu+kk+1))*detlogs
     &      -0.5d0*trace
      
      return
      end


c=======================================================================
      subroutine dmvnd(n,x,mu,sigma,eval,iflag)        
c=======================================================================
c     return the log of a multivariate normal density
c     in this routine sigma is destroyed!!!.
c     A.J.V., 2007
      implicit none
      integer n,i,j
      double precision mu(n),sigma(n,n),x(n)
      double precision det,sse,eval
      integer iflag(n)
      double precision work1,work2,work3,tpi
     
      work1=0.d0
      work2=0.d0
      work3=0.d0
      det=0.d0
      sse=0.d0

      det=0.d0
      
      tpi=6.283185307179586476925286766559d0
       
      work1=-(dble(n)*log(tpi))
    
      call inversedet(sigma,n,iflag,det)

      work2=det
   
      do i=1,n
         do j=1,n
            sse=sse+(x(i)-mu(i))*sigma(i,j)*(x(j)-mu(j))
         end do
      end do

      work3=sse
     
      eval=(work1-work2-work3)/2.d0
      
      return
      end
      
c=======================================================================
      subroutine dtmvn(n,mu,sigma,lower,upper,work1,work2,typeint,
     &                 ynew,yold,logcgkn)
c=======================================================================
c     evaluate the canditate generating distribution in the 
c     multivariate truncated normal generation using a 
c     Gibbs sampler. 
c
c     typeint=1 (left)
c     typeint=2 (interval)
c     typeint=3 (right)
c     typeint=4 (unconstrained support)    
c     typeint=5 known, not sampled
c
c     A.J.V. 2007
c=======================================================================
      implicit none

c+++++Input
      integer n
      integer typeint(n)
      double precision mu(n),sigma(n,n)
      double precision lower(n),upper(n)
      double precision work1(n,n),work2(n,n)
      double precision ynew(n),yold(n)

c+++++Output
      double precision logcgkn

c+++++Working
      integer i,j,k,maxn
      parameter(maxn=50)
      double precision dnrm,cdfnorm
      double precision muv(maxn)
      double precision muc,sigmac
      double precision tmp1
      double precision slow,supp
      double precision yv(maxn)
      double precision yt(maxn)
      logical ainf,binf      

c+++++Algorithm

      if(maxn.lt.n)then
         call rexit("Increase 'maxn' in 'dtmvn'")
      end if   

      logcgkn=0.d0 
      
      do i=1,n
         yt(i)=ynew(i)
      end do

      do i=1,n
         yt(i)=yold(i)

         call condmvn(i,sigma,n,work1,work2)
         sigmac=sqrt(work1(1,1))
      
         do j=1,n
            tmp1=0.d0
            do k=1,n
               tmp1=tmp1+work2(j,k)*mu(k) 
            end do
            muv(j)=tmp1

            tmp1=0.d0
            do k=1,n
               tmp1=tmp1+work2(j,k)*yt(k) 
            end do
            yv(j)=tmp1
         end do

         muc=mu(i)
         do j=2,n
             muc=muc-work1(1,j)*(yv(j)-muv(j))
         end do
      
         if(typeint(i).eq.1)then
            ainf=.true.
            binf=.false. 
          else if(typeint(i).eq.2)then
            ainf=.false.
            binf=.false. 
          else if(typeint(i).eq.3)then
            ainf=.false.
            binf=.true. 
          else 
            ainf=.true.
            binf=.true. 
         end if
         
         slow=lower(i)
         supp=upper(i)
         
         if(typeint(i).ne.5)then
            logcgkn=logcgkn+dnrm(yt(i),muc,sigmac,1)

            if(typeint(i).eq.1)then
               logcgkn=logcgkn-cdfnorm(supp,muc,sigmac,1,1)

             else if(typeint(i).eq.2)then
               logcgkn=logcgkn-
     &             log(               
     &                  cdfnorm(supp,muc,sigmac,1,0)-
     &                  cdfnorm(slow,muc,sigmac,1,0)
     &                ) 
             else if(typeint(i).eq.3)then
               logcgkn=logcgkn-cdfnorm(slow,muc,sigmac,0,1)
            end if

         end if  
      end do
      return
      end

c=======================================================================
      subroutine dtriang(x,a,b,c,dens)        
c=======================================================================
c     return the density of a triangular(a,b,c) distribution
c     a = lim inf
c     b = lim sup
c     c = mode
c     A.J.V., 2007
c=======================================================================
      implicit none
      double precision a,b,c,x
      double precision dens
      
      dens=0.d0
      if(a.le.x.and.x.le.c)then
         dens=2.d0*(x-a)/((b-a)*(c-a))
      end if
      if(c.le.x.and.x.le.b)then
         dens=2.d0*(b-x)/((b-a)*(b-c))
      end if

      return
      end
     
      
        
