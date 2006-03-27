
c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES AND FUNCTIONS FOR RN GENERATION
c=======================================================================                  
c=======================================================================                  

c=======================================================================            
      subroutine setrand(seed1,seed2,seed3)
c=======================================================================            
c     This routine initialize the random number generator.
c     A.J.V., 2005
      integer seed1,seed2,seed3
      common/rrunif/iix,iiy,iiz
      iix=seed1
      iiy=seed2
      iiz=seed3
      return
      end

c=======================================================================                  
      real function runif()
c=======================================================================                  
c     This function generates a uniform random variable using the
c     algorithm of Wichman and Hill (1982) 
c     A.J.V., 2005
      integer ix,iy,iz
      common/rrunif/ix,iy,iz
      ix=171*mod(ix,177)-2*(ix/177)
      iy=172*mod(iy,176)-35*(iy/176)
      iz=170*mod(iz,178)-63*(iz/178)
      if (ix.lt.0) ix=ix+30269
      if (iy.lt.0) iy=iy+30307
      if (iz.lt.0) iz=iz+30323
      runif=amod(float(ix)/30269.0+float(iy)/30307.0+
     +		    float(iz)/30323.0,1.0)
      return
      end


c=======================================================================                        
      subroutine rdisc(imin,imax,evali)
c=======================================================================                  
c     This subroutine generates a discrite uniform random variable in 
c     {imin,...,imax}
c     A.J.V., 2005
      implicit none 
      integer imin,imax,evali
      real runif
      evali=imin+((imax-imin)*runif())
      return
      end    


c=======================================================================                        
      double precision function rbeta(a0,b0)
c=======================================================================                  
c     This function generates a beta random variable  
c     A.J.V., 2005
      implicit none
      real*8 a0,b0,rgamma
      rbeta=0.d0
      if(a0.gt.0.d0)then
         if(b0.gt.0.d0)then
            rbeta=rgamma(a0,1.d0) 
            rbeta=rbeta/(rbeta+rgamma(b0,1.d0))
          else
            rbeta=1.d0
         end if   
      end if   
      return
      end         

c=======================================================================                        
      subroutine dirichlet(alpha,kreal,k,x)
c=======================================================================                  
c     This subroutine generates a dirichlet random vector x  
c     A.J.V., 2005
      implicit none
      integer i,k,kreal
      real*8 alpha(kreal),x(kreal),temp,rgamma,a0

      temp=0.d0

      do i=1,k
         x(i)=0.d0 
         a0=alpha(i)
         if(a0.gt.0.d0)then
            x(i)=rgamma(a0,1.d0) 
            temp=temp+x(i)
         end if   
      end do
      
      do i=1,k
         x(i)=x(i)/temp
      end do   
      return
      end         

c=======================================================================                  
      double precision function rgamma(alpha,beta)
c=======================================================================                  
c     This function generates a random gamma value.
c     It call gamdv which is a modified version of Mike West's gamma
c     generator.
c     gamdv(alpha)  ~ gamma(alpha,1)
c     rangam(alpha,beta) = gamdv(alpha)/beta ~ gamma(alpha,beta)
c     A.J.V., 2005 
      implicit none 
      real*8 beta,alpha,gamdv
      rgamma = gamdv(alpha)
      rgamma = rgamma/beta
      return
      end      

c=======================================================================                  
      double precision function gamdv(a)
c=======================================================================            
c     This subrountine is from Mike West's code. I slightly modified it
c     to use runif as the uniform random number generator.
c     Generates a random gamma variable with shape a>0 and scale=1
c     ix : random seed
c     requires uniform random generator, runif
c     A.J.V., 2005
      implicit none
      real*8 a,aa,ea,u0,u1,u2,c1,c2,c3,c4,c5,w
      real runif
      
      aa=a
      if(aa.eq.1.d0)then
         gamdv = -log(dble(runif()))
	 return
      endif
      if(aa.lt.1.d0)then
         ea=2.7182818d0/(aa+2.7182818d0)
11       u0=dble(runif())
         u1=dble(runif())
         if(u0.le.ea)then
            gamdv=(u0/ea)**(1.d0/aa)
            if(u1.gt.exp(-gamdv))then
               go to 11
              else
               return
            end if
           else
            gamdv=-log((1.d0-u0)/(ea*aa))
            if(u1.gt.gamdv**(aa-1.d0))then
               go to 11
              else
               return
              end if
            end if
	   else
            c1=aa-1.d0
            c2=(aa-1.d0/(6.d0*aa))/c1
            c3=2.d0/c1
            c4=c3+2.d0
            c5=1.d0/sqrt(aa)
12          u1=dble(runif())
            u2=dble(runif())
            u1=u2+c5*(1.d0-1.86d0*u1)
            if(abs(u1-0.5d0).ge.0.5d0)go to 12
            w=c2*u2/u1
            if((c3*log(u1)-log(w)+w).ge.1.d0) then
               go to 12
              else
               gamdv=c1*w
	     end if
        end if
        return
        end

c=======================================================================
      subroutine simdisc(prob,n,m,val)
c=======================================================================
c     generates a sample from a discrete distribution
c     n= real dimension
c     m= used dimension      
c     A.J.V., 2005
      implicit none 
      integer n,m,val,i1,ok
      real*8 prob(n),temp1,u,total
      real runif
      
      total=0.d0
      do i1=1,m
         total=total+prob(i1)
      end do

      if(total.eq.0.d0)then
c        call intpr('error in simdisc',-1,0,0)
c        call dblepr('probability',-1,prob,m)
        call rdisc(1,m,val) 
        return
      end if  

c++++ Generating the rn
      temp1=0.d0
      
      u=dble(runif())
      i1=1
      ok=1
      do while(ok.eq.1.and.i1.le.m)
         temp1=temp1+(prob(i1)/total)
         if(u.lt.temp1)then
            val=i1
            ok=0 
         end if
         i1=i1+1
      end do
      
      return
      end           
      
c=======================================================================            	        
      double precision function rtslogistic(ind,eta)
c=======================================================================            	  
c     This function gerenates a truncated logistic(alpha=0,beta=1) 
c     random variable. The truncation region is (-Inf,eta] if ind=1 and
c     (eta,+Inf) if ind=0.
c     A.J.V., 2005
      implicit none 
      integer ind
      real*8 uni,eta,cdfslogistic,invcdfslogistic
      real runif
      uni=runif()
      rtslogistic=0.d0
      if(ind.eq.1)then
         rtslogistic=invcdfslogistic(uni*cdfslogistic(eta))
      end if
      if(ind.eq.0)then
         rtslogistic=invcdfslogistic(uni+(1-uni)*cdfslogistic(eta))
      end if
      return
      end
      
c=======================================================================            	        
      double precision function rtslogistic2(ainf,binf,a,b)
c=======================================================================            	  
c     This function gerenates a truncated logistic(alpha=0,beta=1) 
c     random variable. The truncation region is (a,b) 
c     a,b  = end points of interval (ainf = binf = .false.)   
c     ainf = true, if left endpoint is infinite; otherwise = false
c     binf = true, if right endpoint is infinite; otherwise = false      
c     A.J.V., 2005
      implicit none 
      real*8 uni,cdfslogistic,invcdfslogistic,a,b
      real runif
      logical ainf,binf
      
      uni=dble(runif())
      rtslogistic2=0.d0
      if(ainf.and.binf) go to 110
      if(ainf.or.binf) go to 100
      
      if(a.gt.b)then
        call intpr('error in limits',-1,0,0)
        rtslogistic2=a
        return
      end if  

      rtslogistic2=invcdfslogistic(cdfslogistic(a)+
     &             uni*(cdfslogistic(b)-cdfslogistic(a)))
      go to 120

100   if(ainf)then
         rtslogistic2=invcdfslogistic(uni*cdfslogistic(b))
         go to 120
      end if
      if(binf)then
         rtslogistic2=invcdfslogistic(uni+(1-uni)*cdfslogistic(a))
         go to 120
      end if

110   rtslogistic2=invcdfslogistic(uni)

120   continue
      return
      end      
      
c=======================================================================            
      double precision function rnorm2(mu,sd)
c=======================================================================            
c     This function generates a N(mu,sd^2) random values.
c     A.J.V., 2005
      implicit none
      real*8 mu,sd
      real u,runif,ppnda
      u=runif()
      rnorm2 = mu+sd*ppnda(u)
      return
      end

c=======================================================================            
      double precision function rnorm(mu,sd)
c=======================================================================            
c     This function generates a N(mu,sd^2) random values.
c     A.J.V., 2006
      implicit none
      real*8 mu,sd,stdnrm

      rnorm = mu+stdnrm()*sd
      return
      end


c=======================================================================            
      double precision function stdnrm()
c=======================================================================            
c     Function to return a single standard normal (0,1) deviate.  
c     This function calls the NORMAL subroutine for a value.
c     A.J.V., 2006
      implicit none
      dimension temp(1)
      real*8 temp

      call normal(temp,1)
      stdnrm= temp(1)
      return
      end
      
c=======================================================================            
      subroutine normal(sndvec,n)
c=======================================================================            
c     Draw n random numbers from a normal(0,1) distribution using the
c     algorithm described by Kennedy & Gentle 'Statistical Computing'
c     Marcel Dekker, New York, 1980 pp 205-207.
c     A.J.V., 2006
      implicit none
      integer N
      dimension sndvec(n)
      real*8 sndvec
      real runif
      real*8 a,asq
      parameter(a=2.216035867166471d0,asq=a*a)
      real*8 u1,u2,u3,t,g
      integer i

      do i=1,n
         u1 = dble(runif())
         if (u1.lt.0.884070402298758d0) then
            u2 = dble(runif())
            sndvec(i) = a * (1.13113163544180d0 * u1 + u2 - 1.d0)
            
         else if (u1.ge.0.973310954173898d0) then
            u2 = dble(runif())
            u3 = dble(runif())
            do while (u2*u2.ge.(asq/(asq-2.d0*dlog(U3))))
               u2 = dble(runif())
               u3 = dble(runif())
            end do
            if (u1.lt.0.986655477086949d0) then
               sndvec(i) = dsqrt(asq-2.d0*dlog(u3))
            else
               sndvec(i) = -dsqrt(asq-2.d0*dlog(u3))
            end if
            
         else if (u1.ge.0.958720824790463d0) then
            u2 = dble(runif())
            u3 = dble(runif())
            t = a - (0.630834801921960d0 * dmin1(u2,u3))
            g = .3989422804014327d0 * dexp(t * t * 0.5d0)
     &           - 0.180025191068563d0 * (a - dabs(t))
            do while ((dmax1(u2,u3).gt.0.755591531667601d0).and.
     &            (0.034240503750111d0*dabs(u2-u3).gt.g))
               u2 = dble(runif())
               u3 = dble(runif())
               t = a - (0.630834801921960d0 * dmin1(u2,u3))
               g = .3989422804014327d0 * dexp(t * t * 0.5d0)
     &               - 0.180025191068563d0 * (a - dabs(t))
            end do

            if (u2.lt.u3) then
               sndvec(i)=t
            else
               sndvec(i)=-t
            end if

         else if (u1.ge.0.911312780288703d0) then
            u2 = dble(runif())
            u3 = dble(runif())
            t = .479727404222441d0 +(1.105473661022070d0 * dmin1(u2,u3))
            g = .3989422804014327d0 * dexp(t * t * 0.5d0)
     &           - 0.180025191068563d0 * (a - dabs(t))
            do while ((dmax1(u2,u3).gt.0.872834976671790d0).and.
     &            (0.049264496373128d0*dabs(u2-u3).gt.g))
               u2 = dble(runif())
               u3 = dble(runif())
               t = .479727404222441d0+(1.105473661022070d0*dmin1(u2,u3))
               g = .3989422804014327d0 * dexp(t * t * 0.5d0)
     &               - 0.180025191068563d0 * (a - dabs(t))
            end do

            if (u2.lt.u3) then
               sndvec(i) = t
            else
               sndvec(i) = -t
            end if

         else
            u2 = dble(runif())
            u3 = dble(runif())
            do while (dmax1(u2,u3).gt.0.805577924423817d0)
               u2 = dble(runif())
               u3 = dble(runif())
            end do
            t = 0.47972740422241d0 - (0.595507138015940d0*dmin1(u2,u3))

            if (u2.lt.u3) then
               sndvec(i) = t
            else
               sndvec(i) = -t
            end if
         end if
      end do
      return
      end

    
c======================================================================      
      real*8 function rtnorm(mu,sd,a,b,ainf,binf)
c=======================================================================            
c     generate truncated(a,b) Normal(mu,sd**2) using the Geweke's 
c     algorithm.
c     mu is the mean of TN distribution
c     sd is standard deviation of TN distribution
c     a,b  = end points of interval (ainf = binf = .false.)   
c     ainf = true, if left endpoint is infinite; otherwise = false
c     binf = true, if right endpoint is infinite; otherwise = false      
c     A.J.V., 2005
      implicit none
      real*8 mu,sd,a,b,a1,b1,rtsnorm
      logical ainf,binf
      a1=(a-mu)/sd
      b1=(b-mu)/sd
      rtnorm=mu+sd*rtsnorm(a1,b1,ainf,binf)
      return
      end      
      
c======================================================================            
      real*8 function rtsnorm(a,b,la,lb)
c======================================================================            
c     generates a N(0,1) random variable
c     subject to the constraint that it be in an interval
c     (a,b), where the endpoints may be finite or infinite.
c     a, b    endpoints of interval; a < b if la = lb = .false.
c     la      .true. if left endpoint is - infinity; in this
c             case A is ignored.
c     lb      .true. if right endpoint is + infinity; in this
c             case B is ignored.
c     A.J.V., 2005
      implicit real*8 (a-h,o-z)
      logical la,lb,lflip
      real runif
      real*8 dexpone,rnorm
      data eps,t1,t2,t3,t4/2.0d0,.375d0,2.18d0,.725d0,.45d0/

      if(la.and.lb)go to 160
      lflip=.false.
      if(la.or.lb)go to 100
      if(b.le.a)go to 170
c ******* Finite interval
      c1=a
      c2=b
      if((c1*c2).gt.0.0d0)go to 30
c ++++ (A,B) includes 0
      if((c1.gt.-t1).and.(c2.lt.t1))go to 20
c -- F(A) or F(B) small: full normal with rejection
   10 x=rnorm(0.d0,1.d0)
      if(x.lt.c1)go to 10
      if(x.gt.c2)go to 10
      GO TO 150
c -- F(A) and F(B) large: uniform importance sampling
   20 cdel=c2-c1
   25 x=c1+cdel*dble(runif())
      if(dble(runif()).gt.dexpone(x))go to 25
      go to 150
c ++++ (A,B) excludes 0
c -- Transform to both positive
   30 if(c1.gt.0.0d0)go to 40
      c=c1
      c1=-c2
      c2=-c
      lflip=.true.
   40 f1=dexpone(c1)
      f2=dexpone(c2)
      if(f2.lt.eps)go to 60
      if((f1/f2).gt.t2)go to 60
c  -- F(A)/F(B) not large: uniform importance sampling
   50 cdel=c2-c1
   55 x=c1+cdel*runif()
      if(dble(runif()).gt.(dexpone(x)/f1))go to 55
      go to 140
   60 if(c1.gt.t3)go to 80
c -- P(X>A) and F(A)/F(B) large: half-normal with rejection
   70 x=abs(rnorm(0.d0,1.d0))
      if(x.lt.c1)go to 70
      if(x.gt.c2)go to 70
      go to 140
c -- P(X>A) small, F(A)/F(B) large: exponential importance
c    sampling with rejection
   80 c=c2-c1
   90 z=-log(runif())/c1
      if(z.gt.c)go to 90
      if(dble(runif()).gt.dexpone(z))GO TO 90
      x=c1+z
      go to 140
c ****** Half-line interval
  100 c1=a
c -- Transform to bound from below if A = -infinity
      if(lb)go to 110
      c1=-b
      lflip=.true.
  110 if(c1.gt.t4)go to 130
c -- A not large: full normal with rejection
  120 x=rnorm(0.d0,1.d0)
      if(x.lt.c1)go to 120
      go to 140
c -- A small: exponential importance sampling
  130 z=-log(runif())/c1
      if(dble(runif()).gt.dexpone(z))go to 130
      x=c1+z
  140 if(lflip)x=-x
  150 rtsnorm=X
      return
c ****** Whole interval
  160 rtsnorm=rnorm(0.d0,1.d0)
      return
c  ***** Singleton
  170 rtsnorm=A
      return
      end

c=======================================================================                        
      real*8 function dexpone(x)
c=======================================================================                       
c     evaluate a exponential function
c     A.J.V., 2005
      implicit none
      real*8 x,expin
      expin=-.5d0*x**2
      if (expin .le. -50.0d0) then
        dexpone=0.0d0
      else
        dexpone=dexp(expin)
      end if
      return
      end
      
c=======================================================================            
      subroutine rmvnorm(n,mean,sigma,work1,work2,y)
c=======================================================================      
c     Subroutine to generate vector of N normal variates with 
c     mean = MEAN and variance = SIGMA
c     A.J.V., 2005
      implicit none
      integer n
      real*8 mean(n),sigma(n,n),work1(n*(n+1)/2),work2(n),y(n)

      call cholesky(n,sigma,work1)
      call mvnchol(n,work1,mean,work2,y)
      return
      end      
      
c=======================================================================      
      subroutine mvnchol(n,l,mean,work,y)
c=======================================================================      
c     Subroutine to generate vector of N normal variates with 
c     mean = MEAN and variance = LL', i.e., L is the Cholesky
c     decomposition of the desired variance-covariance structure.
c     WORK is a double precision work vector of at least N elements.
c     The subroutine calls NORMAL for a vector of N iid normal(0,1)
c     deviates (stored in work).  The new variables are calculated as 
c     MEAN + L*WORK.
c     A.J.V., 2005
      implicit none
      integer n,i,j,jj
      real*8 l(n*(n+1)/2),mean(n),work(n),y(n)
      
      call normalvec(n,work)
      
      do i=1,n
         y(i)=mean(i)
      end do
      
      jj = 0
      
      do i = 1,n
         do j = i,n
            jj = jj + 1
            y(j) = y(j) + l(jj)*work(i)
         end do
      end do
      return
      end
 
c=======================================================================      
      subroutine normalvec(n,work)
c=======================================================================
c     generates a vector of normal variables
c     A.J.V., 2005
      implicit none
      integer i,n
      real*8 work(n),rnorm
      
      do i=1,n
         work(i)=rnorm(0.d0,1.d0)
      end do
      return
      end

c=======================================================================                        
      subroutine rbinom(n,p,evali)
c=======================================================================                  
c     This subroutine generates a Binomial(n,p) random variable 
c     A.J.V., 2005
      implicit none 
      integer n,evali,i
      real*8 uni,p,c,pr,f
      real runif
      
      if(p.lt.0.d0)p=0.d0
      if(p.gt.1.d0)p=1.d0
      
      uni=dble(runif())
      c=exp(log(p)-log(1-p))
      i=0
      pr=exp(n*log(1-p))
      f=pr
      do while(uni.ge.f.and.i.le.n)
         pr=pr*(c*(n-i)/(i+1))
         f=f+pr
         i=i+1
      end do
      evali=i
      return
      end         


c======================================================================      
      real*8 function rtlnorm(mu,sd,a,b,ainf,binf)
c=======================================================================            
c     generate truncated(a,b) LogNormal(mu,sd**2) using the Geweke's 
c     algorithm.
c     mu is the mean of log(variable)
c     sd is standard deviation of log(variable)
c     a,b  = end points of interval (ainf = binf = .false.)   
c     ainf = true, if left endpoint is infinite; otherwise = false
c     binf = true, if right endpoint is infinite; otherwise = false      
c     A.J.V., 2005
      implicit none
      real*8 mu,sd,a,b,a1,b1,rtsnorm,x
      logical ainf,binf
      
      if(ainf)then
         a1=0.d0
        else
         a1=(log(a)-mu)/sd  
      end if  

      if(binf)then
         b1=0.d0
        else
         b1=(log(b)-mu)/sd  
      end if  

      x=mu+sd*rtsnorm(a1,b1,ainf,binf)
      rtlnorm=exp(x)
      
      return
      end   
      

c=======================================================================                  
      double precision function rchisq(nu)
c=======================================================================                  
c     This function generates a random chi2 value.
c     A.J.V., 2005 
      implicit none 
      real*8 nu,rgamma
      rchisq=rgamma(0.5d0*nu,0.5d0)
      return
      end            
      
      
c=======================================================================      
      subroutine riwishart(n,nu0,K0,workm1,workm2,workv1,workhs1,
     &                     workhs2,iflag) 
c=======================================================================      
c     generate a (co)variance matrix using a inverted wishart 
c     distribution
c
c     Arguments:
c              n=dimension
c              K0= QUADF + T-1 
c              nu0=tb+nsubject 
c
c     Parametrization of the prior should be:
c              E(SIGMA)=1/(tb-n-1) * T-1
c
c     Output:
c              K0     = sigma  
c              workm1 = sigmainv
c
c     A.J.V., 2005
      implicit none
      integer i,j,n,nu0,iflag(n),ihmssf
      real*8 K0(n,n),workm1(n,n),workm2(n,n),workv1(n)
      real*8 workhs1(n*(n+1)/2),workhs2(n*(n+1)/2),detlog

      do i=1,n
         do j=1,n
            workm1(i,j)=k0(i,j)
         end do 
      end do

      call invdet(workm1,n,workm2,detlog,iflag,workv1)

      call cholesky(n,workm2,workhs1)

      call wishrt(workhs1,workhs2,n,n,nu0+1)

      do i=1,n
         do j=1,n
            workm1(i,j)=workhs2(ihmssf(i,j,n))
            workm2(i,j)=workhs2(ihmssf(i,j,n))
         end do 
      end do
      
      call invdet(workm2,n,k0,detlog,iflag,workv1)
      
      return
      end


c=======================================================================      
      subroutine wishrt(l,s,maxn,nv,n)
c=======================================================================      
c     Subroutine to generate a matrix from a Wishart distribution, i.e.
c     a matrix of sums of squares and cross products.  This is done for
c     a specific variance covariance matrix and number of traits.  The
c     algorithm is based on:
c
c     Odell, P. L., and A. H. Feiveson.  1966.  A numerical procedure to
c        to generate a sample covariance matrix.  Journal of the
c        American Statistical Association (JASA) 61:199-203.
c
c     The subroutine was derived from one written by Karin Meyer, which
c     was based on the algorithm adapted in Kennedy & Gentle in
c     'Statistical Computing,' Marcel Dekker 1980, p.231-232.
c
c     Variables passed:
c
c     Input
c        l  - lower triangular matrix, halfstored, resulting from the 
c             cholesky decomposition of the variance/covariance matrix in 
c             the population
c
c        nv - number of variables (= p in paper)
c
c        n  - sample size (no. of records) 
c====>        n is not the degrees of freedom!!
c
c      Output
c        s -  matrix of SS/CP sampled from the Wishart distr. This matrix
c             is halfstored (upper triangle)

      implicit none
      integer*4 maxvar,nv,n,nv1,i,j,k,ndf,jk,cnt,maxn
      parameter(maxvar=15)
      real*8 l,s,z,v,b,acc,chi(1)
      dimension l(maxn*(maxn+1)/2),s(maxn*(maxn+1)/2)
      dimension z(maxvar,maxvar),v(maxvar),b(maxvar,maxvar)
      real*8 rnorm,rchisq	
	  
      if(nv.gt.maxvar)then
           call intpr('"wishrt": rotuine dimension exceeded!',-1,0,0)
           call intpr('current maximum is',-1,maxvar,0)            
           call intpr('required size is',-1,nv,0)            
           return      
      end if

      nv1=nv+1
            
c     -----------------------------------------------------------------
c     set up vectors (z-alpha in paper) of normal deviates
c     -----------------------------------------------------------------
      do j=2,nv
	     do k=1,nv
		    v(k)=rnorm(0.d0,1.d0)
		 end do
         do i=1,j-1
            z(i,j)=v(i)
         end do
      end do

c     -----------------------------------------------------------------
c     set up vector of chi-square deviations (v-i in paper)
c     -----------------------------------------------------------------
      do i=1,nv
         ndf=n-i
         chi(1)=rchisq(dble(ndf))
         v(i)=chi(1)
      end do

c     -----------------------------------------------------------------
c     calculate diagonals of B
c     -----------------------------------------------------------------
      b(1,1) = v(1)
      v(1) = dsqrt(v(1))
      do j = 2,nv
         acc = 0.d0
         do i = 1,j-1
            acc = acc + z(i,j)*z(i,j)
         end do
         b(j,j) = v(j) + acc
         v(j) = dsqrt(v(j))
      end do

c     -----------------------------------------------------------------
c     calculate off diagonals of B
c     -----------------------------------------------------------------
      do j=2,nv
         b(1,j) = z(1,j)*v(1)
         b(j,1) = b(1,j)
         do i=2,j-1
            acc = 0.d0
            do k = 1,i-1
               acc = acc + z(k,i)*z(k,j)
            end do
            b(i,j) = z(i,j)*v(i) + acc
            b(j,i) = b(i,j)
         end do
      end do

c     -----------------------------------------------------------------
c     calculate L*B*L'
c     -----------------------------------------------------------------

c     FIRST CALCULATE L*B AND STORE IN B - REMEMBER THAT L IS LOWER
c     TRIANGULAR

      nv1 = nv + 1
      do i = 1,nv
         do j = 1,nv
            jk = j - nv
            acc = 0.d0
            do k = 1,j
               jk = jk + nv1 - k
               acc = acc +  l(jk) * b(k,i)
            end do
            v(j) = acc
         end do
         do j = 1,nv
            b(j,i) = v(j)
         end do
      end do

c     NEXT, CALCULATE LB*L' WHERE LB IS STORED IN B, NOW L' IS UPPER 
c     TRIANGULAR, AND THAT LBL' IS SYMMETRIC

      cnt = 0
      do i = 1,nv
         do j = i,nv
            cnt = cnt + 1
            acc = 0.d0
            jk = j - nv
            do k = 1,j
               jk = jk + nv1 - k
               acc = acc + b(i,k) * l(jk)
            end do
            s(cnt) = acc
         end do
      end do
      
      return
      end
      

c=======================================================================      
      integer function ihmssf(i,j,n)
c=======================================================================      
      integer i,j,n,i1,j1
      if(i.le.j)then
         i1=i-1
         ihmssf=n*i1-i*i1/2+j
       else
         j1=j-1
         ihmssf=n*j1-j*j1/2+i
      end if
      return
      end      
      
      
c=======================================================================                  
      subroutine samalph(alpha,aa0,ab0,ndis,k)
c=======================================================================            
c     This routine samples another alpha value using the technique
c     in Escobar and West (TR 533). The prior distribution for alpha
c     is a gamma(aa0, ab0).
c
c     ndis: the number of clusters.
c     k   : Sample size
c     A.J.V., 2005
      integer ndis,k
      real*8 alpha, aa0, ab0, xalp,s,e,rgamma
      real runif

      xalp=rgamma(1.d0+alpha,1.d0)
      xalp=xalp/(xalp+rgamma(dble(k),1.d0))
      s=ab0-dlog(xalp)
      e=aa0+dble(ndis)-1.d0
      if(dble(runif()).lt.e/(e+dble(k)*s))then
         e=e+1.d0
      end if
      alpha=rgamma(e,s)
      return
      end

