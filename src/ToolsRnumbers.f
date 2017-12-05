
c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES AND FUNCTIONS FOR RN GENERATION
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
      real function runif()
c=======================================================================                  
c     This function generates a uniform random variable
c     A.J.V., 2006
      real ranf
      runif=ranf()
      return
      end

c=======================================================================                        
      double precision function rexpo(lambda)
c=======================================================================                  
c     This function generates a exp(lambda) random variable  
c     A.J.V., 2005
      implicit none
      double precision lambda
      real runif
      rexpo=-log(1.d0-dble(runif()))/lambda
      return
      end     

c=======================================================================                        
      subroutine rdisc(imin,imax,evali)
c=======================================================================                  
c     This subroutine generates a discrite uniform random variable in 
c     {imin,...,imax}
c     A.J.V., 2006
      implicit none 
      integer imin,imax,evali,ignuin

      evali=ignuin(imin,imax)
      if(evali.lt.imin)evali=imin
      if(evali.gt.imax)evali=imax
      return
      end    

c=======================================================================                  
c      double precision function rgamma(alpha,beta)
c=======================================================================                  
c     This function generates a random gamma value.
c     The parametrization is such that E(X)=alpha/beta
c     A.J.V., 2006 
c      implicit none 
c      double precision beta,alpha
c      real a,r,gengam
c      a=beta
c      r=alpha
c      rgamma = gengam(a,r)
c      return
c      end      


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
      double precision beta,alpha,gamdv
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
      double precision a,aa,ea,u0,u1,u2,c1,c2,c3,c4,c5,w
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
      double precision function rbeta(a0,b0)
c=======================================================================                  
c     This function generates a beta random variable  
c     A.J.V., 2006
      implicit none
      double precision a0,b0
      real genbet,a,b
      a=a0
      b=b0
      rbeta=genbet(a,b)
      return
      end         

c======================================================================      
      double precision function rtbeta(alpha,beta,a,b,ainf,binf)
c=======================================================================            
c     generate truncated(a,b) Beta(alpha,beta)
c     a,b  = end points of interval (ainf = binf = .false.)   
c     ainf = true, if left endpoint is 0; otherwise = false
c     binf = true, if right endpoint is 1; otherwise = false      
c     A.J.V., 2006
      implicit none
      double precision alpha,beta,a,b
      double precision rbeta,invcdfbetas,cdfbetas
      double precision uni,tmp,tmp1,tmp2
      logical ainf,binf
      real runif

      uni=dble(runif())
      rtbeta=0.d0
      if(ainf.and.binf) go to 110
      if(ainf.or.binf) go to 100
      
      if(a.gt.b)then
        call rexit("error in limits rtbeta")
        rtbeta=a
        return
      end if  

      tmp1=cdfbetas(a,alpha,beta,1,0)
      tmp2=cdfbetas(b,alpha,beta,1,0)
      tmp=tmp1+uni*(tmp2-tmp1)
     
      rtbeta=invcdfbetas(tmp,alpha,beta,1,0)
      go to 120

100   if(ainf)then

         tmp2=cdfbetas(b,alpha,beta,1,0)
         tmp=uni*tmp2
         rtbeta=invcdfbetas(tmp,alpha,beta,1,0)
         go to 120
      end if
      if(binf)then
         tmp1=cdfbetas(a,alpha,beta,1,0)      
         tmp=uni+(1.d0-uni)*tmp1
         rtbeta=invcdfbetas(tmp,alpha,beta,1,0)
         go to 120
      end if

110   rtbeta=rbeta(alpha,beta)

120   continue

      return
      end      
      
c======================================================================      
      double precision function rtbeta2(alpha,beta,a,b)
c=======================================================================            
c     generate truncated(a,b) Beta(alpha,beta) using a AR
c     a,b  = end points of interval
c     A.J.V., 2006
      implicit none
      integer exit
      double precision alpha,beta,a,b
      double precision rbeta

      exit=0
      do while(exit.eq.0)
         rtbeta2=rbeta(alpha,beta)
         if(rtbeta2.gt.a.and.rtbeta2.le.b)exit=1
      end do

      return
      end      


c=======================================================================                        
      subroutine dirichlet(alpha,kreal,k,x)
c=======================================================================                  
c     This subroutine generates a dirichlet random vector x  
c     A.J.V., 2005
      implicit none
      integer i,k,kreal
      double precision alpha(kreal),x(kreal),tmp,rgamma,a0

      tmp=0.d0

      do i=1,k
         x(i)=0.d0 
         a0=alpha(i)
         if(a0.gt.0.d0)then
            x(i)=rgamma(a0,1.d0) 
            tmp=tmp+x(i)
         end if   
      end do
      
      do i=1,k
         x(i)=x(i)/tmp
      end do   
      return
      end         

c=======================================================================
      subroutine simdisc(prob,n,m,val)
c=======================================================================
c     generates a sample from a discrete distribution
c     n= real dimension
c     m= used dimension      
c     A.J.V., 2006
      implicit none 
      integer n,m,val,i1,ok
      double precision prob(n),temp1,u,total
      real runif
      
      total=0.d0
      do i1=1,m
         total=total+prob(i1)
      end do

      if(total.eq.0.d0)then
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
      subroutine simdiscint(prob,n,imin,imax,val)
c=======================================================================
c     generates a sample from a discrete distribution
c     n= real dimension
c     m= used dimension      
c     A.J.V., 2006
      implicit none 
      integer n,val,i1,ok,imin,imax
      double precision prob(n),temp1,u,total
      real runif
      
      total=0.d0
      do i1=imin,imax
         total=total+prob(i1)
      end do

      if(total.eq.0.d0)then
        call rdisc(imin,imax,val) 
        return
      end if  

c++++ Generating the rn
      temp1=0.d0
      
      u=dble(runif())
      i1=imin
      ok=1
      do while(ok.eq.1.and.i1.le.imax)
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
      double precision uni,eta,cdfslogistic,invcdfslogistic
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
      double precision uni,cdfslogistic,invcdfslogistic,a,b
      real runif
      logical ainf,binf
      
      uni=dble(runif())
      rtslogistic2=0.d0
      if(ainf.and.binf) go to 110
      if(ainf.or.binf) go to 100
      
      if(a.gt.b)then
        call rexit("error in limits rtslogistic2")
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
         rtslogistic2=invcdfslogistic(uni+(1.d0-uni)*cdfslogistic(a))
         go to 120
      end if

110   rtslogistic2=invcdfslogistic(uni)

120   continue
      return
      end      
      

c=======================================================================            
      double precision function rnorm(mu,sd)
c=======================================================================            
c     This function generates a N(mu,sd^2) random values.
c     A.J.V., 2006
      implicit none
      double precision mu,sd
      real gennor,av0,sd0
      
      av0=mu
      sd0=sd
      rnorm = gennor(av0,sd0)
      return
      end

c======================================================================      
      double precision function rtnorm(mu,sd,a,b,ainf,binf)
c=======================================================================            
c     generate truncated(a,b) Normal(mu,sd**2) using the Geweke's 
c     algorithm.
c     mu is the mean of TN distribution
c     sd is standard deviation of TN distribution
c     a,b  = end points of interval (ainf = binf = .false.)   
c     ainf = true, if left endpoint is infinite; otherwise = false
c     binf = true, if right endpoint is infinite; otherwise = false      
c     A.J.V., 2006
      implicit none
      double precision mu,sd,a,b,a1,b1,rtsnorm
      logical ainf,binf
      a1=(a-mu)/sd
      b1=(b-mu)/sd
      rtnorm=mu+sd*rtsnorm(a1,b1,ainf,binf)
      return
      end      
      
c======================================================================            
      double precision function rtsnorm(a,b,la,lb)
c======================================================================            
c     generates a N(0,1) random variable
c     subject to the constraint that it be in an interval
c     (a,b), where the endpoints may be finite or infinite.
c     a, b    endpoints of interval; a < b if la = lb = .false.
c     la      .true. if left endpoint is - infinity; in this
c             case A is ignored.
c     lb      .true. if right endpoint is + infinity; in this
c             case B is ignored.
c     A.J.V., 2006
      implicit double precision (a-h,o-z)
      logical la,lb,lflip
      real runif
      double precision dexpone,rnorm
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
      cdel=c2-c1
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
      double precision function dexpone(x)
c=======================================================================                       
c     evaluate a exponential function
c     A.J.V., 2006
      implicit none
      double precision x,expin
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
c     A.J.V., 2006
      implicit none
      integer n
      double precision mean(n),sigma(n,n),work1(n*(n+1)/2),work2(n),y(n)

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
c     A.J.V., 2006
      implicit none
      integer n,i,j,jj
      double precision l(n*(n+1)/2),mean(n),work(n),y(n)
      
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
      subroutine rmvnorm2(nr,n,mean,sigma,work1,work2,y)
c=======================================================================      
c     Subroutine to generate vector of N normal variates with 
c     mean = MEAN and variance = SIGMA
c     A.J.V., 2007
      implicit none
      integer nr,n
      double precision mean(nr),sigma(nr,nr),work1(nr*(nr+1)/2),
     1  work2(nr),y(nr)

      call cholesky2(nr,n,sigma,work1)
      call mvnchol2(nr,n,work1,mean,work2,y)
      return
      end      

c=======================================================================      
      subroutine mvnchol2(nr,n,l,mean,work,y)
c=======================================================================      
c     Subroutine to generate vector of N normal variates with 
c     mean = MEAN and variance = LL', i.e., L is the Cholesky
c     decomposition of the desired variance-covariance structure.
c     WORK is a double precision work vector of at least N elements.
c     The subroutine calls NORMAL for a vector of N iid normal(0,1)
c     deviates (stored in work).  The new variables are calculated as 
c     MEAN + L*WORK.
c     A.J.V., 2007
      implicit none
      integer nr,n,i,j,jj
      double precision l(nr*(nr+1)/2),mean(nr),work(nr),y(nr)
      double precision rnorm
      
      do i=1,n
         work(i)=rnorm(0.d0,1.d0)
      end do
      
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
c     A.J.V., 2006
      implicit none
      integer i,n
      double precision work(n),rnorm
      
      do i=1,n
         work(i)=rnorm(0.d0,1.d0)
      end do
      return
      end


c=======================================================================                        
      subroutine rbinom(n,p,evali)
c=======================================================================                  
c     This subroutine generates a Binomial(n,p) random variable 
c     A.J.V., 2006
      implicit none 
      integer n,evali,ignbin
      double precision p
      real pp
      
      pp=p

      evali=ignbin(n,pp) 

      return
      end        

c======================================================================      
      double precision function rtlnorm(mu,sd,a,b,ainf,binf)
c=======================================================================            
c     generate truncated(a,b) LogNormal(mu,sd**2) using the Geweke's 
c     algorithm.
c     mu is the mean of log(variable)
c     sd is standard deviation of log(variable)
c     a,b  = end points of interval (ainf = binf = .false.)   
c     ainf = true, if left endpoint is infinite; otherwise = false
c     binf = true, if right endpoint is infinite; otherwise = false      
c     A.J.V., 2006
      implicit none
      double precision mu,sd,a,b,a1,b1,rtsnorm,x
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
c     A.J.V., 2006
      implicit none 
      double precision nu
      double precision rgamma
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
      double precision K0(n,n),workm1(n,n),workm2(n,n),workv1(n)
      double precision workhs1(n*(n+1)/2),workhs2(n*(n+1)/2),detlog

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
      integer maxvar,nv,n,nv1,i,j,k,ndf,jk,cnt,maxn
      parameter(maxvar=70)
      double precision l,s,z,v,b,acc,chi(1)
      dimension l(maxn*(maxn+1)/2),s(maxn*(maxn+1)/2)
      dimension z(maxvar,maxvar),v(maxvar),b(maxvar,maxvar)
      double precision rnorm,rchisq

      if(nv.gt.maxvar)then
           call rexit("wishrt: routine dimension exceeded")
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
c     A.J.V., 2006
      integer ndis,k
      double precision alpha, aa0, ab0, xalp,s,e,rgamma
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


c=======================================================================            	        
      double precision function rtcauchy(mu,sd,a,b,ainf,binf)
c=======================================================================            	  
c     This function gerenates a truncated cauchy(alpha=0,beta=1) 
c     random variable. The truncation region is (a,b) 
c     a,b  = end points of interval (ainf = binf = .false.)   
c     ainf = true, if left endpoint is infinite; otherwise = false
c     binf = true, if right endpoint is infinite; otherwise = false      
c     A.J.V., 2006
      implicit none 
      double precision uni,mu,sd,cdfcauchy,invcdfcauchy,a,b
      real runif
      logical ainf,binf
      
      uni=dble(runif())
      rtcauchy=0.d0
      if(ainf.and.binf) go to 110
      if(ainf.or.binf) go to 100
      
      if(a.gt.b)then
        call rexit("error in limits rtcauchy")      
        rtcauchy=a
        return
      end if  

      rtcauchy=invcdfcauchy(cdfcauchy(a,mu,sd,1,0)+
     &         uni*(cdfcauchy(b,mu,sd,1,0)-cdfcauchy(a,mu,sd,1,0)),
     &         mu,sd,1,0)
      go to 120

100   if(ainf)then
         rtcauchy=invcdfcauchy(uni*cdfcauchy(b,mu,sd,1,0),mu,sd,1,0)
         go to 120
      end if
      if(binf)then
         rtcauchy=invcdfcauchy(uni+(1.d0-uni)*cdfcauchy(a,mu,sd,1,0),
     &                         mu,sd,1,0)
         go to 120
      end if

110   rtcauchy=invcdfcauchy(uni,mu,sd,1,0)

120   continue
      return
      end


c=======================================================================      
      integer function rpois(mu)
c=======================================================================      
c     This function gerenates a Poisson(mu) random variable. 
c     A.J.V., 2006
      integer ignpoi
      double precision  mu
      real mean

      mean=mu
      rpois=ignpoi(mean)
      return
      end


c=======================================================================      
      subroutine rmult(n,p,ncat,x)
c=======================================================================      
c     This function gerenates a Multinomial random vector. 
c     A.J.V., 2007
      integer n,ncat,x(ncat)
      real p(ncat-1)

      call genmul(n,p,ncat,x)

      return
      end


c======================================================================      
      double precision function rtchisq(nu,a,b,ainf,binf)
c=======================================================================            
c     generate truncated(a,b) Chisq(nu)
c     a,b  = end points of interval (ainf = binf = .false.)   
c     ainf = true, if left endpoint is 0; otherwise = false
c     binf = true, if right endpoint is 1; otherwise = false      
c     A.J.V., 2007
      implicit none
      double precision nu,a,b
      double precision rchisq,invcdfchisq,cdfchisq
      double precision uni,tmp
      logical ainf,binf
      real runif

      uni=dble(runif())

      rtchisq=0.d0
      if(ainf.and.binf) go to 110
      if(ainf.or.binf) go to 100
      
      if(a.gt.b)then
        call rexit("error in limits rtchisq")
        rtchisq=a
        return
      end if  

      tmp=cdfchisq(a,nu,1,0)+ 
     &    uni*(cdfchisq(b,nu,1,0)-cdfchisq(a,nu,1,0))
      rtchisq=invcdfchisq(tmp,nu,1,0)
      go to 120

100   if(ainf)then
         tmp=uni*cdfchisq(b,nu,1,0)
         rtchisq=invcdfchisq(tmp,nu,1,0)
         go to 120
      end if
      if(binf)then
         tmp=uni+(1.d0-uni)*cdfchisq(a,nu,1,0)
         rtchisq=invcdfchisq(tmp,nu,1,0)
         go to 120
      end if

110   rtchisq=rchisq(nu)

120   continue

      return
      end      


c======================================================================      
      subroutine rperm(nl,n,p)
c=======================================================================            
c     generate a random permutation of the first n integers
c     A.J.V., 2007
      implicit none
      integer n,nl,p(nl),i,j,k,ipj,itemp,m
      double precision u(100)
      real runif

      do i=1,n
         p(i)=i
      end do
      
      do i=1,n,100
         m=min(n-i+1,100)
         do j=1,100
            u(j)=dble(runif())
         end do
         
         do j=1,m
            ipj=i+j-1
            k=int(u(j)*(n-ipj+1))+ipj
            itemp=p(ipj)
            p(ipj)=p(k)
            p(k)=itemp
         end do
      end do
      return
      end


c=======================================================================
      subroutine rtmvn2(n,mu,sigma,lower,upper,l,typeint,y)
c=======================================================================
c     simulate multivariate truncated normal variables using a 
c     recursive approach. The samples can only be used in MH
c     steps. They are not proper samples from a TMVN distribution.
c
c     A.J.V. 2007
c=======================================================================
      implicit none

c+++++Input 
c     typeint=1 (left)
c     typeint=2 (interval)
c     typeint=3 (right)
c     typeint=4 (unconstrained support)    

      integer n,typeint(n)
      double precision mu(n),sigma(n,n),lower(n),upper(n),l(n*(n+1)/2)

c+++++Output  
      double precision y(n)

c+++++Working      
      integer i,j,maxn,ihmssf
      parameter(maxn=50)
      double precision a(maxn,maxn),slow,supp,tmp1,z(maxn)
      double precision rtnorm
      logical ainf,binf      

c+++++Check dimenions
      if(n.gt.maxn)then
         call rexit("increase dimension maxn in subroutine rtmvn2")
      end if

c+++++Algorithm 

      do i=1,n
         do j=1,n
            a(i,j)=0.d0
         end do 
      end do
 
      call cholesky(n,sigma,l)

      do i=1,n
         do j=1,i
            a(i,j)=l(ihmssf(i,j,n))
         end do 
      end do

      do i=1,n
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
      
         if(i.eq.1)then
            slow=(lower(i)-mu(i))/a(i,i)
            supp=(upper(i)-mu(i))/a(i,i)
          else
            tmp1=0.d0
            do j=1,i-1
               tmp1=tmp1+a(i,j)*z(j)
            end do
            slow=(lower(i)-mu(i)-tmp1)/a(i,i)
            supp=(upper(i)-mu(i)-tmp1)/a(i,i)
         end if
         
         z(i)=rtnorm(0.d0,1.d0,slow,supp,ainf,binf)
      end do
      
      do i=1,n
         tmp1=0.d0
         do j=1,n
            tmp1=tmp1+a(i,j)*z(j)   
         end do
         y(i)=mu(i)+tmp1
      end do

      return
      end

c=======================================================================                        
      integer function rpoiss(mu)
c=======================================================================                  
c     This function generates a poisson random variable  
c     A.J.V., 2006
      implicit none
      double precision mu
      integer ignpoi
      real mur
      
      mur=mu
      rpoiss=ignpoi(mur)
      return
      end         

c=======================================================================
      subroutine rtmvn(n,mu,sigma,lower,upper,work1,work2,typeint,
     &                 y,logcgko)
c=======================================================================
c     simulate multivariate truncated normal variables using a 
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

c+++++Output
      double precision y(n),logcgko

c+++++Working
      integer i,j,k,maxn
      parameter(maxn=50)
      double precision dnrm,cdfnorm
      double precision muv(maxn)
      double precision muc,sigmac
      double precision rtnorm
      double precision tmp1
      double precision slow,supp
      double precision yv(maxn)
      logical ainf,binf      

c+++++Algorithm

      if(maxn.lt.n)then
         call rexit("Increase 'maxn' in 'rtmvn'")
      end if   

      logcgko=0.d0 

      do i=1,n
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
               tmp1=tmp1+work2(j,k)*y(k) 
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
            y(i)=rtnorm(muc,sigmac,slow,supp,ainf,binf)
            logcgko=logcgko+dnrm(y(i),muc,sigmac,1)

            if(typeint(i).eq.1)then
               logcgko=logcgko-cdfnorm(supp,muc,sigmac,1,1)

             else if(typeint(i).eq.2)then
               logcgko=logcgko-
     &             log(               
     &                  cdfnorm(supp,muc,sigmac,1,0)-
     &                  cdfnorm(slow,muc,sigmac,1,0)
     &                ) 
             else if(typeint(i).eq.3)then
               logcgko=logcgko-cdfnorm(slow,muc,sigmac,0,1)
            end if

         end if  

      end do
      return
      end


c======================================================================      
      integer function rpoiss2(mu)
c======================================================================            
c     generate a Poisson(mu) random number
c     A.J.V., 2007
c======================================================================
      implicit none
      integer i
      double precision mu,uni,p,f
      real runif

      i=0
      p=dexp(-mu)
      f=p
      uni=dble(runif())

      rpoiss2=i
      do while(uni.ge.f)
         p=mu*p/dble(i+1)
         f=f+p
         i=i+1
         rpoiss2=i
      end do

      return
      end      


c======================================================================      
      subroutine rtbetas(alpha,beta,a,b,ainf,binf,val)
c=======================================================================            
c     generate truncated(a,b) Beta(alpha,beta)
c     a,b  = end points of interval (ainf = binf = .false.)   
c     ainf = true, if left endpoint is 0; otherwise = false
c     binf = true, if right endpoint is 1; otherwise = false      
c     A.J.V., 2006
      implicit none
      double precision alpha,beta,a,b
      double precision bound      
      double precision rbeta
      double precision val
      double precision uni,tmp,tmp1,tmp2,tmp3,tmp4
      logical ainf,binf
      real runif
      integer status

      uni=dble(runif())
      val=0.d0
      if(ainf.and.binf) go to 110
      if(ainf.or.binf) go to 100
      
      if(a.gt.b)then
        call rexit("error in limits rtbetas")
        val=a
        return
      end if  


      call cdfbet(1,tmp1,tmp3,a,1.d0-a,alpha,beta,status,bound)
      call cdfbet(1,tmp2,tmp3,b,1.d0-b,alpha,beta,status,bound)
      tmp=tmp1+uni*(tmp2-tmp1)
      
      tmp2=1.d0-tmp
      call cdfbet(2,tmp,tmp2,val,tmp4,alpha,beta,status,bound)

      go to 120

100   if(ainf)then

         call cdfbet(1,tmp2,tmp3,b,1.d0-b,alpha,beta,status,bound)
         tmp=uni*tmp2
         tmp2=1.d0-tmp
         call cdfbet(2,tmp,tmp2,val,tmp4,alpha,beta,status,bound)
         go to 120
      end if
      if(binf)then

         call cdfbet(1,tmp1,tmp3,a,1.d0-a,alpha,beta,status,bound)
         tmp=uni+(1.d0-uni)*tmp1
      
         tmp2=1.d0-tmp
         call cdfbet(2,tmp,tmp2,val,tmp4,alpha,beta,status,bound)
         go to 120
      end if

110   val=rbeta(alpha,beta)

120   continue

      return
      end   



c======================================================================            
      subroutine ginoe(n,variance,rmat)
c======================================================================      
c     Author: Valerio Cappellini                          March 2007
c     Version:: 1.0.0
c     Modified: Alejandro Jara (to use RNG) 
c
c     Generator of REAL N x N non symmetric matrices drawn according 
c     to the Ginibre ensemble called <GinOE> \subset GL(N,R). 
c     The distribution of matrix Z and its matrix 
c     elements are given by
c                                        z_{ij}^2
c                        1            - ----------
c     P(z_{ij}) = ----------------  e     2 s^2
c                  SQRT(2\pi s^2) 
c
c      and
c
c                       1 
c     P(Z) = -----------------------  exp[ - Tr Z^2 / (2 s^2) ]
c             (2\pi)^(N^2/2) s^(N^2)
c
c     where s denotes the input variable <VARIANCE> and Z the output 
c     variable <RMAT>.
c
c======================================================================      
      SAVE  S, T, A, B, R1, R2
      INTEGER N, emme, mu
      REAL U(2), S, T, A, B, R1, R2
      double precision RMAT(N,N)
      REAL V, X, Y, Q, DEVIAT, VARIANCE
      real ranf
      DATA  S, T, A, B / 0.449871, -0.386595, 0.19600, 0.25472/
      DATA  R1, R2 / 0.27597, 0.27846/
C         generate pair of uniform deviates

      DO 200 emme = 1, N
          DO 200 mu = 1, N
   50 u(1)=ranf()
      u(2)=ranf()
     
      V = 1.7156 * (U(2) - 0.5)
      X = U(1) - S
      Y = ABS(V) - T
      Q = X**2 + Y*(A*Y - B*X)

C           accept P if inside inner ellipse
      IF (Q .LT. R1)  GO TO 100
C           reject P if outside outer ellipse
      IF (Q .GT. R2)  GO TO 50
C           reject P if outside acceptance region
      IF (V**2 .GT. -4.0 *ALOG(U(1)) *U(1)**2)  GO TO 50
C           ratio of P's coordinates is normal deviate
  100 DEVIAT = V/U(1)*VARIANCE
  200 RMAT(emme,mu) = DEVIAT

      RETURN
      END


c======================================================================            
      subroutine rhaar(n,a,q)
c======================================================================      
c     Original: SUBROUTINE Haar_O_N(n,Q,ROUTIN)
c     Author: Valerio Cappellini                          March 2007
c     Version:: 1.0.0
c     Modified: Alejandro Jara (to use RNG) 
c
c     This Program produce N x N REAL Random Orthogonal Matrices 
c     distributed according
c     to the Haar measure, making use of the QR decomposition of N x N 
c     REAL Random (non symmetric) Matrices from the Ginibre Ensemble. 
c     The QR decomposition is performed
c     by means of (N-1) Householder reflections. 
c
c     Algorithm:
c     ^^^^^^^^^
c      + we start producing an N x N matrix A(0) = a(0)_{ij} from the 
c        Ginibre Ensemble.
c      + we fix Q(0) = I_N ( the N x N identity )
c
c         + we perform an iterated procedure on the index <k> running 
c           from 1 to N-1
c     
c             at each step a matrix H_k is produced and the running 
c             matrices A(k-1) 
c             and Q(k-1) are upgrated to A(k-1) -> A(k) = H_k * A(k-1) 
c             and Q(k-1) -> Q(k) = Q(k-1) * H_k
c
c         + end of iterated procedures
c
c      + If A(N-1)_{NN} is negative we change sign to the last 
c        column of Q(N-1)
c
c     Finally Q := Q(N-1) is the REAL Random Orthogonal Matrix 
c     distributed according to the Haar measure and 
c
c          t
c     R = Q  * A(0) (here and in the following the superscript t 
c     denotes the usual transposition) is an upper triangular matrix 
c     with positive diagonal entries so that this upper positive 
c     triangularity could be used as a test for the subroutine,
c     toghether with the orthogonality of Q.
c     
c     At each step the algorithm make use of the N-k+1 dimensional 
c     basis vector 
c
c                                    t
c     e_k : (1 , 0 , 0 , ... , 0 , 0)         , of the N-k+1 
c     dimensional vector
c
c                                                           t
c     v := (a_{k,k} , a_{k+1,k} , ... , a_{N-1,k} , a_{N,k})   ,
c
c     its Euclidean norm || v ||  ,  its first entry v_1 = a_{kk}   ,
c
c     the sign of the latter sgn(v_1).
c     
c     Then we construct the vector   u := v + sgn(v_1) * || v || * e_k 
c     and its 
c
c     positive rescaled-squared-norm c_k := || u ||^2 / 2  =
c     || v ||*(|| v || + | v_1 |)
c     so that finally
c
c
c                /                   |                   \
c               |       I_{k-1}      |         0          |
c               |                    |                    |
c      H_k :=   | ___________________|___________________ |
c               |                    |                    |
c               |                    |                    |
c               |          0         |        M_k         |           
c                \                   |                   /
c     
c     with I_{k-1} being the (k-1) x (k-1) identity , and M_k the 
c     (N-k+1) x (N-k+1)
c
c                                             /             u * u   \
c     dimensional matrix  M_k := -sgn(v_1) * | I_{N-k+1} - --------  |
c                                             \              c_k    /
c
c======================================================================      
      INTEGER n
      REAL dk,ck
      double precision a(n,n),Q(n,n),sum
      REAL VARIANCE,sigma,tau
      INTEGER i,j,k

      VARIANCE=1./SQRT(2.)
      call ginoe(n,VARIANCE,a)

      Do i=1,n
          Do j=1,n
              Q(i,j)=0.d0
          Enddo
          Q(i,i)=1.d0
      Enddo  

      sum=0.
      do 17 k=1,n-1
        do 11 i=k,n
          sum=max(sum,abs(a(i,k)))
11      continue
        if(sum.ne.0.)then         
          sum=0.
          do 13 i=k,n
            sum=sum+a(i,k)**2
13        continue                      
          sigma=sign(sqrt(sum),a(k,k))  
          a(k,k)=a(k,k)+sigma           
          ck=sigma*a(k,k)               
                                        
          dk=-sign(1.,sigma)            

          do 16 j=k+1,n                 
            sum=0.                      
            do 14 i=k,n                 
              sum=sum+a(i,k)*a(i,j)     
14          continue                    
            tau=sum/ck                  
            do 15 i=k,n                     
              a(i,j)=dk*(a(i,j)-tau*a(i,k)) 
15          continue                        
16        continue                           

          do 26 i=1,n                   
            sum=0.                      
            do 24 j=k,n                 
              sum=sum+Q(i,j)*a(j,k)     
24          continue                    
            tau=sum/ck                  
            do 25 j=k,n                     
              Q(i,j)=dk*(Q(i,j)-tau*a(j,k)) 
25          continue                        
26        continue

        endif
17    continue
      if(a(n,n).lt.0.) then  
          Do i=1,n           
              Q(i,n)=-Q(i,n) 
          Enddo              
      Endif                  
      
      return
      END


c======================================================================            
      subroutine rhaar2(xwork,x,n,q)
c======================================================================      
c     subroutine to generate a orthogonal random matrix with
c     haar distribution given a matrix of normally distributed
c     elements. This can be used for random walks in x.
c
c     A.J.V., 2008
c======================================================================      
      implicit none
      integer n
      double precision x(n,n),q(n,n)
      double precision xwork(n,n)
  
      integer i,j,k
      double precision ck,dk,scale,sigma,sums,tau   

      do i=1,n
         do j=1,n
            xwork(i,j)=x(i,j)
            q(i,j)=0.d0
         end do
         q(i,i)=1.d0
      end do

      do k=1,n-1
         scale=0.d0
         do i=1,k
            scale=max(scale,abs(xwork(i,k)))
         end do
         if(scale.ne.0.d0)then
            sums=0.d0  
            do i=k,n
               sums=sums+xwork(i,k)**2
            end do  
            sigma=sign(sqrt(sums),xwork(k,k)) 
            xwork(k,k)=xwork(k,k)+sigma
            ck=sigma*xwork(k,k) 

            dk=-sign(1.d0,sigma)

            do j=k+1,n
               sums=0.d0
               do i=k,n
                  sums=sums+xwork(i,k)*xwork(i,j) 
               end do
               tau=sums/ck
               do i=k,n
                  xwork(i,j)=dk*(xwork(i,j)-tau*xwork(i,k))
               end do
            end do

            do i=1,n
               sums=0.d0
               do j=k,n
                  sums=sums+q(i,j)*xwork(j,k)
               end do
               tau=sums/ck
               do j=k,n
                  q(i,j)=dk*(q(i,j)-tau*xwork(j,k)) 
               end do
            end do 

         end if  
      end do

      if(xwork(n,n).lt.0.d0)then
         do i=1,n
            q(i,n)=-q(i,n)
         end do
      end if  

      return
      end 



