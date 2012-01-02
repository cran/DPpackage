      SUBROUTINE advnst(k)
C**********************************************************************
C
C     SUBROUTINE ADVNST(K)
C               ADV-a-N-ce ST-ate
C
C     Advances the state  of  the current  generator  by 2^K values  and
C     resets the initial seed to that value.
C
C     This is  a  transcription from   Pascal to  Fortran    of  routine
C     Advance_State from the paper
C
C     L'Ecuyer, P. and  Cote, S. "Implementing  a  Random Number Package
C     with  Splitting   Facilities."  ACM  Transactions  on Mathematical
C     Software, 17:98-111 (1991)
C
C
C                              Arguments
C
C
C     K -> The generator is advanced by2^K values
C                                   INTEGER K
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalar Arguments ..
      INTEGER k
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER g,i,ib1,ib2
C     ..
C     .. External Functions ..
      INTEGER mltmod
      LOGICAL qrgnin
      EXTERNAL mltmod,qrgnin
C     ..
C     .. External Subroutines ..
      EXTERNAL getcgn,setsd
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/
C     ..
C     .. Executable Statements ..
C     Abort unless random number generator initialized
      IF (qrgnin()) GO TO 10
      call rexit("ADVNST called before RNG initialized")

   10 CALL getcgn(g)
C
      ib1 = a1
      ib2 = a2
      DO 20,i = 1,k
          ib1 = mltmod(ib1,ib1,m1)
          ib2 = mltmod(ib2,ib2,m2)
   20 CONTINUE
      CALL setsd(mltmod(ib1,cg1(g),m1),mltmod(ib2,cg2(g),m2))
C
C     NOW, IB1 = A1**K AND IB2 = A2**K
C
      RETURN

      END

      REAL FUNCTION genbet(aa,bb)
C**********************************************************************
C
C     REAL FUNCTION GENBET( A, B )
C               GeNerate BETa random deviate
C
C
C                              Function
C
C
C     Returns a single random deviate from the beta distribution with
C     parameters A and B.  The density of the beta is
C               x^(a-1) * (1-x)^(b-1) / B(a,b) for 0 < x < 1
C
C
C                              Arguments
C
C
C     A --> First parameter of the beta distribution
C                         REAL A
C     JJV                 (A > 1.0E-37)
C
C     B --> Second parameter of the beta distribution
C                         REAL B
C     JJV                 (B > 1.0E-37)
C
C
C                              Method
C
C
C     R. C. H. Cheng
C     Generating Beta Variates with Nonintegral Shape Parameters
C     Communications of the ACM, 21:317-322  (1978)
C     (Algorithms BB and BC)
C
C**********************************************************************
C     .. Parameters ..
C     Close to the largest number that can be exponentiated
      REAL expmax
C     JJV changed this - 89 was too high, and LOG(1.0E38) = 87.49823
      PARAMETER (expmax=87.49823)
C     Close to the largest representable single precision number
      REAL infnty
      PARAMETER (infnty=1.0E38)
C     JJV added the parameter minlog
C     Close to the smallest number of which a LOG can be taken.
      REAL minlog
      PARAMETER (minlog=1.0E-37)
C     ..
C     .. Scalar Arguments ..
      REAL aa,bb
C     ..
C     .. Local Scalars ..
      REAL a,alpha,b,beta,delta,gamma,k1,k2,olda,oldb,r,s,t,u1,u2,v,w,y,
     +     z
      LOGICAL qsame
C     ..
C     .. External Functions ..
      REAL ranf
      EXTERNAL ranf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC exp,log,max,min,sqrt
C     ..
C     .. Save statement ..
C     JJV added a,b
      SAVE olda,oldb,alpha,beta,gamma,k1,k2,a,b
C     ..
C     .. Data statements ..
C     JJV changed these to ridiculous values
      DATA olda,oldb/-1.0E37,-1.0E37/
C     ..
C     .. Executable Statements ..
      qsame = (olda.EQ.aa) .AND. (oldb.EQ.bb)
      IF (qsame) GO TO 20
C     JJV added small minimum for small log problem in calc of W
      IF (.NOT. (aa.LT.minlog.OR.bb.LT.minlog)) GO TO 10
      call rexit("AA or BB too small in GENBET - Abort")

   10 olda = aa
      oldb = bb
   20 IF (.NOT. (min(aa,bb).GT.1.0)) GO TO 100


C     Alborithm BB

C
C     Initialize
C
      IF (qsame) GO TO 30
      a = min(aa,bb)
      b = max(aa,bb)
      alpha = a + b
      beta = sqrt((alpha-2.0)/ (2.0*a*b-alpha))
      gamma = a + 1.0/beta
   30 CONTINUE
   40 u1 = ranf()
C
C     Step 1
C
      u2 = ranf()
      v = beta*log(u1/ (1.0-u1))
C     JJV altered this
      IF (v.GT.expmax) GO TO 55
C     JJV added checker to see if a*exp(v) will overflow
C     JJV 50 _was_ w = a*exp(v); also note here a > 1.0
      w = exp(v)
      IF (w.GT.infnty/a) GO TO 55
      w = a*w
      GO TO 60
 55   w = infnty

   60 z = u1**2*u2
      r = gamma*v - 1.3862944
      s = a + r - w
C
C     Step 2
C
      IF ((s+2.609438).GE. (5.0*z)) GO TO 70
C
C     Step 3
C
      t = log(z)
      IF (s.GT.t) GO TO 70
C
C     Step 4
C
C     JJV added checker to see if log(alpha/(b+w)) will 
C     JJV overflow.  If so, we count the log as -INF, and
C     JJV consequently evaluate conditional as true, i.e.
C     JJV the algorithm rejects the trial and starts over
C     JJV May not need this here since ALPHA > 2.0
      IF (alpha/(b+w).LT.minlog) GO TO 40

      IF ((r+alpha*log(alpha/ (b+w))).LT.t) GO TO 40
C
C     Step 5
C
   70 IF (.NOT. (aa.EQ.a)) GO TO 80
      genbet = w/ (b+w)
      GO TO 90

   80 genbet = b/ (b+w)
   90 GO TO 230


C     Algorithm BC

C
C     Initialize
C
  100 IF (qsame) GO TO 110
      a = max(aa,bb)
      b = min(aa,bb)
      alpha = a + b
      beta = 1.0/b
      delta = 1.0 + a - b
      k1 = delta* (0.0138889+0.0416667*b)/ (a*beta-0.777778)
      k2 = 0.25 + (0.5+0.25/delta)*b
  110 CONTINUE
  120 u1 = ranf()
C
C     Step 1
C
      u2 = ranf()
      IF (u1.GE.0.5) GO TO 130
C
C     Step 2
C
      y = u1*u2
      z = u1*y
      IF ((0.25*u2+z-y).GE.k1) GO TO 120
      GO TO 170
C
C     Step 3
C
  130 z = u1**2*u2
      IF (.NOT. (z.LE.0.25)) GO TO 160
      v = beta*log(u1/ (1.0-u1))

C     JJV instead of checking v > expmax at top, I will check
C     JJV if a < 1, then check the appropriate values

      IF (a.GT.1.0) GO TO 135
C     JJV A < 1 so it can help out if EXP(V) would overflow
      IF (v.GT.expmax) GO TO 132
      w = a*exp(v)
      GO TO 200
 132  w = v + log(a)
      IF (w.GT.expmax) GO TO 140
      w = exp(w)
      GO TO 200

C     JJV in this case A > 1
 135  IF (v.GT.expmax) GO TO 140
      w = exp(v)
      IF (w.GT.infnty/a) GO TO 140
      w = a*w
      GO TO 200
 140  w = infnty
      GO TO 200

  160 IF (z.GE.k2) GO TO 120
C
C     Step 4
C
C
C     Step 5
C
  170 v = beta*log(u1/ (1.0-u1))

C     JJV same kind of checking as above
      IF (a.GT.1.0) GO TO 175
C     JJV A < 1 so it can help out if EXP(V) would overflow
      IF (v.GT.expmax) GO TO 172
      w = a*exp(v)
      GO TO 190
 172  w = v + log(a)
      IF (w.GT.expmax) GO TO 180
      w = exp(w)
      GO TO 190

C     JJV in this case A > 1
 175  IF (v.GT.expmax) GO TO 180
      w = exp(v)
      IF (w.GT.infnty/a) GO TO 180
      w = a*w
      GO TO 190

  180 w = infnty

C     JJV here we also check to see if log overlows; if so, we treat it
C     JJV as -INF, which means condition is true, i.e. restart
  190 IF (alpha/(b+w).LT.minlog) GO TO 120
      IF ((alpha* (log(alpha/ (b+w))+v)-1.3862944).LT.log(z)) GO TO 120
C
C     Step 6
C
  200 IF (.NOT. (a.EQ.aa)) GO TO 210
      genbet = w/ (b+w)
      GO TO 220

  210 genbet = b/ (b+w)
  220 CONTINUE
  230 RETURN

      END
      REAL FUNCTION genchi(df)
C**********************************************************************
C
C     REAL FUNCTION GENCHI( DF )
C                Generate random value of CHIsquare variable
C
C
C                              Function
C
C
C     Generates random deviate from the distribution of a chisquare
C     with DF degrees of freedom random variable.
C
C
C                              Arguments
C
C
C     DF --> Degrees of freedom of the chisquare
C            (Must be positive)
C                         REAL DF
C
C
C                              Method
C
C
C     Uses relation between chisquare and gamma.
C
C**********************************************************************
C     .. Scalar Arguments ..
      REAL df
C     ..
C     .. External Functions ..
C      REAL gengam
C      EXTERNAL gengam
      REAL sgamma
      EXTERNAL sgamma
C     ..
C     .. Executable Statements ..
      IF (.NOT. (df.LE.0.0)) GO TO 10
      call rexit("DF <= 0 in GENCHI - ABORT")

C     JJV changed this to call sgamma directly
C   10 genchi = 2.0*gengam(1.0,df/2.0)
 10   genchi = 2.0*sgamma(df/2.0)
      RETURN

      END
      REAL FUNCTION genexp(av)

C**********************************************************************
C
C     REAL FUNCTION GENEXP( AV )
C
C                    GENerate EXPonential random deviate
C
C
C                              Function
C
C
C     Generates a single random deviate from an exponential
C     distribution with mean AV.
C
C
C                              Arguments
C
C
C     AV --> The mean of the exponential distribution from which
C            a random deviate is to be generated.
C                              REAL AV
C     JJV                      (AV >= 0)
C
C     GENEXP <-- The random deviate.
C                              REAL GENEXP
C
C
C                              Method
C
C
C     Renames SEXPO from TOMS as slightly modified by BWB to use RANF
C     instead of SUNIF.
C
C     For details see:
C
C               Ahrens, J.H. and Dieter, U.
C               Computer Methods for Sampling From the
C               Exponential and Normal Distributions.
C               Comm. ACM, 15,10 (Oct. 1972), 873 - 882.
C
C**********************************************************************
C     .. Scalar Arguments ..
      REAL av
C     ..
C     .. External Functions ..
      REAL sexpo
      EXTERNAL sexpo
C     ..
C     .. Executable Statements ..
C     JJV added check to ensure AV >= 0.0 
      IF (av.GE.0.0) GO TO 10
      call rexit("AV < 0.0 in GENEXP - ABORT")

 10   genexp = sexpo()*av
      RETURN

      END
      REAL FUNCTION genf(dfn,dfd)
C**********************************************************************
C
C     REAL FUNCTION GENF( DFN, DFD )
C                GENerate random deviate from the F distribution
C
C
C                              Function
C
C
C     Generates a random deviate from the F (variance ratio)
C     distribution with DFN degrees of freedom in the numerator
C     and DFD degrees of freedom in the denominator.
C
C
C                              Arguments
C
C
C     DFN --> Numerator degrees of freedom
C             (Must be positive)
C                              REAL DFN
C      DFD --> Denominator degrees of freedom
C             (Must be positive)
C                              REAL DFD
C
C
C                              Method
C
C
C     Directly generates ratio of chisquare variates
C
C**********************************************************************
C     .. Scalar Arguments ..
      REAL dfd,dfn
C     ..
C     .. Local Scalars ..
      REAL xden,xnum
C     ..
C     JJV changed this code to call sgamma directly
C     .. External Functions ..
C      REAL genchi
C      EXTERNAL genchi
      REAL sgamma
      EXTERNAL sgamma
C     ..
C     .. Executable Statements ..
      IF (.NOT. (dfn.LE.0.0.OR.dfd.LE.0.0)) GO TO 10
      call rexit("Degrees of freedom nonpositive in GENF - abort")

 10   xnum = 2.0*sgamma(dfn/2.0)/dfn

C      GENF = ( GENCHI( DFN ) / DFN ) / ( GENCHI( DFD ) / DFD )
      xden = 2.0*sgamma(dfd/2.0)/dfd
C     JJV changed constant so that it will not underflow at compile time
C     JJV while not slowing generator by using double precision or logs.
C      IF (.NOT. (xden.LE. (1.0E-38*xnum))) GO TO 20
      IF (.NOT. (xden.LE. (1.0E-37*xnum))) GO TO 20

      call rwarn("GENF - generated numbers would cause overflow")
      call rwarn("GGENF returning 1.0E37")
      genf = 1.0E37
      GO TO 30

   20 genf = xnum/xden
   30 RETURN

      END
      REAL FUNCTION gengam(a,r)
C**********************************************************************
C
C     REAL FUNCTION GENGAM( A, R )
C           GENerates random deviates from GAMma distribution
C
C
C                              Function
C
C
C     Generates random deviates from the gamma distribution whose
C     density is
C          (A**R)/Gamma(R) * X**(R-1) * Exp(-A*X)
C
C
C                              Arguments
C
C
C     JJV added the argument ranges supported
C     A --> Location parameter of Gamma distribution
C                              REAL A ( A > 0 )
C
C     R --> Shape parameter of Gamma distribution
C                              REAL R ( R > 0 )
C
C
C                              Method
C
C
C     Renames SGAMMA from TOMS as slightly modified by BWB to use RANF
C     instead of SUNIF.
C
C     For details see:
C               (Case R >= 1.0)
C               Ahrens, J.H. and Dieter, U.
C               Generating Gamma Variates by a
C               Modified Rejection Technique.
C               Comm. ACM, 25,1 (Jan. 1982), 47 - 54.
C     Algorithm GD
C
C     JJV altered the following to reflect sgamma argument ranges
C               (Case 0.0 < R < 1.0)
C               Ahrens, J.H. and Dieter, U.
C               Computer Methods for Sampling from Gamma,
C               Beta, Poisson and Binomial Distributions.
C               Computing, 12 (1974), 223-246/
C     Adapted algorithm GS.
C
C**********************************************************************
C     .. Scalar Arguments ..
      REAL a,r
C     ..
C     .. External Functions ..
      REAL sgamma
      EXTERNAL sgamma
C     ..
C     .. Executable Statements ..

C     JJV added argument value checker
      IF ( a.GT.0.0 .AND. r.GT.0.0 ) GO TO 10
      call rexit("Location or shape param out of range in GENGAM")

C     JJV end addition

 10   gengam = sgamma(r)/a
C      gengam = gengam/a
      RETURN

      END
      SUBROUTINE genmn(parm,x,work)
C**********************************************************************
C
C     SUBROUTINE GENMN(PARM,X,WORK)
C              GENerate Multivariate Normal random deviate
C
C
C                              Arguments
C
C
C     PARM --> Parameters needed to generate multivariate normal
C               deviates (MEANV and Cholesky decomposition of
C               COVM). Set by a previous call to SETGMN.
C               1 : 1                - size of deviate, P
C               2 : P + 1            - mean vector
C               P+2 : P*(P+3)/2 + 1  - upper half of cholesky
C                                       decomposition of cov matrix
C                                             REAL PARM(*)
C
C     X    <-- Vector deviate generated.
C                                             REAL X(P)
C
C     WORK <--> Scratch array
C                                             REAL WORK(P)
C
C
C                              Method
C
C
C     1) Generate P independent standard normal deviates - Ei ~ N(0,1)
C
C     2) Using Cholesky decomposition find A s.t. trans(A)*A = COVM
C
C     3) trans(A)E + MEANV ~ N(MEANV,COVM)
C
C**********************************************************************
C     .. Array Arguments ..
      REAL parm(*),work(*),x(*)
C     ..
C     .. Local Scalars ..
      REAL ae
      INTEGER i,icount,j,p
C     ..
C     .. External Functions ..
      REAL snorm
      EXTERNAL snorm
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC int
C     ..
C     .. Executable Statements ..
      p = int(parm(1))
C
C     Generate P independent normal deviates - WORK ~ N(0,1)
C
      DO 10,i = 1,p
          work(i) = snorm()
   10 CONTINUE
      DO 30,i = 1,p
C
C     PARM (P+2 : P*(P+3)/2 + 1) contains A, the Cholesky
C      decomposition of the desired covariance matrix.
C          trans(A)(1,1) = PARM(P+2)
C          trans(A)(2,1) = PARM(P+3)
C          trans(A)(2,2) = PARM(P+2+P)
C          trans(A)(3,1) = PARM(P+4)
C          trans(A)(3,2) = PARM(P+3+P)
C          trans(A)(3,3) = PARM(P+2-1+2P)  ...
C
C     trans(A)*WORK + MEANV ~ N(MEANV,COVM)
C
          icount = 0
          ae = 0.0
          DO 20,j = 1,i
              icount = icount + j - 1
              ae = ae + parm(i+ (j-1)*p-icount+p+1)*work(j)
   20     CONTINUE
          x(i) = ae + parm(i+1)
   30 CONTINUE
      RETURN
C
      END
      SUBROUTINE genmul(n,p,ncat,ix)
C**********************************************************************
C
C            SUBROUTINE GENMUL( N, P, NCAT, IX )
C     GENerate an observation from the MULtinomial distribution
C
C
C                              Arguments
C
C
C     N --> Number of events that will be classified into one of
C           the categories 1..NCAT
C                         INTEGER N
C
C     P --> Vector of probabilities.  P(i) is the probability that
C           an event will be classified into category i.  Thus, P(i)
C           must be [0,1]. Only the first NCAT-1 P(i) must be defined
C           since P(NCAT) is 1.0 minus the sum of the first
C           NCAT-1 P(i).
C                         REAL P(NCAT-1)
C
C     NCAT --> Number of categories.  Length of P and IX.
C                         INTEGER NCAT
C
C     IX <-- Observation from multinomial distribution.  All IX(i)
C            will be nonnegative and their sum will be N.
C                         INTEGER IX(NCAT)
C
C
C                              Method
C
C
C     Algorithm from page 559 of
C
C     Devroye, Luc
C
C     Non-Uniform Random Variate Generation.  Springer-Verlag,
C     New York, 1986.
C
C**********************************************************************
C     .. Scalar Arguments ..
      INTEGER n,ncat
C     ..
C     .. Array Arguments ..
      REAL p(*)
      INTEGER ix(*)
C     ..
C     .. Local Scalars ..
      REAL prob,ptot,sum
      INTEGER i,icat,ntot
C     ..
C     .. External Functions ..
      INTEGER ignbin
      EXTERNAL ignbin
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs
C     ..
C     .. Executable Statements ..

C     Check Arguments
      IF (n.LT.0) call rexit("N < 0 in GENMUL")
      IF (ncat.LE.1) call rexit("NCAT <= 1 in GENMUL")
      ptot = 0.0
      DO 10,i = 1,ncat - 1
          IF (p(i).LT.0.0) call rexit("Some P(i) < 0 in GENMUL")
          IF (p(i).GT.1.0) call rexit("Some P(i) > 1 in GENMUL")
          ptot = ptot + p(i)
   10 CONTINUE
      IF (ptot.GT.1.0) call rexit("Sum of P(i) > 1 in GENMUL")
C     Initialize variables
      ntot = n
      sum = 1.0
      DO 20,i = 1,ncat
          ix(i) = 0
   20 CONTINUE

C     Generate the observation
      DO 30,icat = 1,ncat - 1
          prob = p(icat)/sum
          ix(icat) = ignbin(ntot,prob)
          ntot = ntot - ix(icat)
          IF (ntot.LE.0) RETURN
          sum = sum - p(icat)
   30 CONTINUE
      ix(ncat) = ntot

C     Finished
      RETURN

      END
      REAL FUNCTION gennch(df,xnonc)
C**********************************************************************
C
C     REAL FUNCTION GENNCH( DF, XNONC )
C           Generate random value of Noncentral CHIsquare variable
C
C
C                              Function
C
C

C     Generates random deviate  from the  distribution  of a  noncentral
C     chisquare with DF degrees  of freedom and noncentrality  parameter
C     XNONC.
C
C
C                              Arguments
C
C
C     DF --> Degrees of freedom of the chisquare
C            (Must be >= 1.0)
C                         REAL DF
C
C     XNONC --> Noncentrality parameter of the chisquare
C               (Must be >= 0.0)
C                         REAL XNONC
C
C
C                              Method
C
C
C     Uses fact that  noncentral chisquare  is  the  sum of a  chisquare
C     deviate with DF-1  degrees of freedom plus the  square of a normal
C     deviate with mean sqrt(XNONC) and standard deviation 1.
C
C**********************************************************************
C     .. Scalar Arguments ..
      REAL df,xnonc
C     ..
C     .. External Functions ..
C     JJV changed these to call SGAMMA and SNORM directly
C      REAL genchi,gennor
C      EXTERNAL genchi,gennor
      REAL sgamma,snorm
      EXTERNAL sgamma,snorm
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC sqrt
C     ..
C     JJV changed abort to df < 1, and added case: df = 1 
C     .. Executable Statements ..
      IF (.NOT. (df.LT.1.0.OR.xnonc.LT.0.0)) GO TO 10
      call rexit("DF < 1 or XNONC < 0 in GENNCH - ABORT")

C     JJV changed this to call SGAMMA and SNORM directly
C      gennch = genchi(df-1.0) + gennor(sqrt(xnonc),1.0)**2

 10   IF (df.GE.1.000001) GO TO 20
C     JJV case DF = 1.0
      gennch = (snorm() + sqrt(xnonc))**2
      GO TO 30

C     JJV case DF > 1.0
 20   gennch = 2.0*sgamma((df-1.0)/2.0) + (snorm() + sqrt(xnonc))**2
 30   RETURN
      
      END
      REAL FUNCTION gennf(dfn,dfd,xnonc)

C**********************************************************************
C
C     REAL FUNCTION GENNF( DFN, DFD, XNONC )
C           GENerate random deviate from the Noncentral F distribution
C
C
C                              Function
C
C
C     Generates a random deviate from the  noncentral F (variance ratio)
C     distribution with DFN degrees of freedom in the numerator, and DFD
C     degrees of freedom in the denominator, and noncentrality parameter
C     XNONC.
C
C
C                              Arguments
C
C
C     DFN --> Numerator degrees of freedom
C             (Must be >= 1.0)
C                              REAL DFN
C      DFD --> Denominator degrees of freedom
C             (Must be positive)
C                              REAL DFD
C
C     XNONC --> Noncentrality parameter
C               (Must be nonnegative)
C                              REAL XNONC
C
C
C                              Method
C
C
C     Directly generates ratio of noncentral numerator chisquare variate
C     to central denominator chisquare variate.
C
C**********************************************************************
C     .. Scalar Arguments ..
      REAL dfd,dfn,xnonc
C     ..
C     .. Local Scalars ..
      REAL xden,xnum
      LOGICAL qcond
C     ..
C     .. External Functions ..
C     JJV changed the code to call SGAMMA and SNORM directly
C      REAL genchi,gennch
C      EXTERNAL genchi,gennch
      REAL sgamma,snorm
      EXTERNAL sgamma,snorm
C     ..
C     .. Executable Statements ..
C     JJV changed the argument checker to allow DFN = 1.0
C     JJV in the same way as GENNCH was changed.
      qcond = dfn .LT. 1.0 .OR. dfd .LE. 0.0 .OR. xnonc .LT. 0.0
      IF (.NOT. (qcond)) GO TO 10
      call rexit("DF or noncent param out of range in GENNF")

C      GENNF = ( GENNCH( DFN, XNONC ) / DFN ) / ( GENCHI( DFD ) / DFD )
C     JJV changed this to call SGAMMA and SNORM directly
C     xnum = gennch(dfn,xnonc)/dfn
 10   IF (dfn.GE.1.000001) GO TO 20
C     JJV case dfn = 1.0 - here I am treating dfn as exactly 1.0
      xnum = (snorm() + sqrt(xnonc))**2
      GO TO 30

C     JJV case dfn > 1.0
 20   xnum = (2.0*sgamma((dfn-1.0)/2.0) + (snorm()+sqrt(xnonc))**2)/dfn

C     xden = genchi(dfd)/dfd
 30   xden = 2.0*sgamma(dfd/2.0)/dfd
      
C     JJV changed constant so that it will not underflow at compile time
C     JJV while not slowing generator by using double precision or logs.
C      IF (.NOT. (xden.LE. (1.0E-38*xnum))) GO TO 40
      IF (.NOT. (xden.LE. (1.0E-37*xnum))) GO TO 40

      call rwarn("GENNF - generated numbers would cause overflow")
      call rwarn("GENNF returning 1.0E37")
      gennf = 1.0E37
      GO TO 50

   40 gennf = xnum/xden
   50 RETURN

      END
      REAL FUNCTION gennor(av,sd)
C**********************************************************************
C
C     REAL FUNCTION GENNOR( AV, SD )
C
C         GENerate random deviate from a NORmal distribution
C
C
C                              Function
C
C
C     Generates a single random deviate from a normal distribution
C     with mean, AV, and standard deviation, SD.
C
C
C                              Arguments
C
C
C     AV --> Mean of the normal distribution.
C                              REAL AV
C
C     SD --> Standard deviation of the normal distribution.
C                              REAL SD
C     JJV                      (SD >= 0)
C
C     GENNOR <-- Generated normal deviate.
C                              REAL GENNOR
C
C
C                              Method
C
C
C     Renames SNORM from TOMS as slightly modified by BWB to use RANF
C     instead of SUNIF.
C
C     For details see:
C               Ahrens, J.H. and Dieter, U.
C               Extensions of Forsythe's Method for Random
C               Sampling from the Normal Distribution.
C               Math. Comput., 27,124 (Oct. 1973), 927 - 937.
C
C
C**********************************************************************
C     .. Scalar Arguments ..
      REAL av,sd
C     ..
C     .. External Functions ..
      REAL snorm
      EXTERNAL snorm
C     ..
C     .. Executable Statements ..
C     JJV added check to ensure SD >= 0.0 
      IF (sd.GE.0.0) GO TO 10
      call rexit("SD < 0.0 in GENNOR - ABORT")

 10   gennor = sd*snorm() + av
      RETURN

      END
      SUBROUTINE genprm(iarray,larray)
C**********************************************************************
C
C    SUBROUTINE GENPRM( IARRAY, LARRAY )
C               GENerate random PeRMutation of iarray
C
C
C                              Arguments
C
C
C     IARRAY <--> On output IARRAY is a random permutation of its
C                 value on input
C                         INTEGER IARRAY( LARRAY )
C
C     LARRAY <--> Length of IARRAY
C                         INTEGER LARRAY
C
C**********************************************************************
C     .. Scalar Arguments ..
      INTEGER larray
C     ..
C     .. Array Arguments ..
      INTEGER iarray(larray)
C     ..
C     .. Local Scalars ..
      INTEGER i,itmp,iwhich
C     ..
C     .. External Functions ..
      INTEGER ignuin
      EXTERNAL ignuin
C     ..
C     .. Executable Statements ..
      DO 10,i = 1,larray
          iwhich = ignuin(i,larray)
          itmp = iarray(iwhich)
          iarray(iwhich) = iarray(i)
          iarray(i) = itmp
   10 CONTINUE
      RETURN

      END
      REAL FUNCTION genunf(low,high)
C**********************************************************************
C
C     REAL FUNCTION GENUNF( LOW, HIGH )
C
C               GeNerate Uniform Real between LOW and HIGH
C
C
C                              Function
C
C
C     Generates a real uniformly distributed between LOW and HIGH.
C
C
C                              Arguments
C
C
C     LOW --> Low bound (exclusive) on real value to be generated
C                         REAL LOW
C
C     HIGH --> High bound (exclusive) on real value to be generated
C                         REAL HIGH
C
C**********************************************************************
C     .. Scalar Arguments ..
      REAL high,low
C     ..
C     .. External Functions ..
      REAL ranf
      EXTERNAL ranf
C     ..
C     .. Executable Statements ..
      IF (.NOT. (low.GT.high)) GO TO 10
      call rexit("LOW > High in GENUNF - Abort")

   10 genunf = low + (high-low)*ranf()

      RETURN

      END
      SUBROUTINE getcgn(g)
      INTEGER g
C**********************************************************************
C
C      SUBROUTINE GETCGN(G)
C                         Get GeNerator
C
C     Returns in G the number of the current random number generator
C
C
C                              Arguments
C
C
C     G <-- Number of the current random number generator (1..32)
C                    INTEGER G
C
C**********************************************************************
C
      INTEGER curntg,numg
      SAVE curntg
      PARAMETER (numg=32)
      DATA curntg/1/
C
      g = curntg
      RETURN

      ENTRY setcgn(g)
C**********************************************************************
C
C     SUBROUTINE SETCGN( G )
C                      Set GeNerator
C
C     Sets  the  current  generator to G.    All references to a generat
C     are to the current generator.
C
C
C                              Arguments
C
C
C     G --> Number of the current random number generator (1..32)
C                    INTEGER G
C
C**********************************************************************
C
C     Abort if generator number out of range
C
      IF (.NOT. (g.LT.0.OR.g.GT.numg)) GO TO 10
      call rexit("Generator number out of range in SETCGN")

   10 curntg = g
      RETURN

      END
      SUBROUTINE getsd(iseed1,iseed2)
C**********************************************************************
C
C     SUBROUTINE GETSD(G,ISEED1,ISEED2)
C               GET SeeD
C
C     Returns the value of two integer seeds of the current generator
C
C     This  is   a  transcription from  Pascal   to  Fortran  of routine
C     Get_State from the paper
C
C     L'Ecuyer, P. and  Cote,  S. "Implementing a Random Number  Package
C     with   Splitting Facilities."  ACM  Transactions   on Mathematical
C     Software, 17:98-111 (1991)
C
C
C                              Arguments
C
C
C
C     ISEED1 <- First integer seed of generator G
C                                   INTEGER ISEED1
C
C     ISEED2 <- Second integer seed of generator G
C                                   INTEGER ISEED1
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalar Arguments ..
      INTEGER iseed1,iseed2
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER g
C     ..
C     .. External Functions ..
      LOGICAL qrgnin
      EXTERNAL qrgnin
C     ..
C     .. External Subroutines ..
      EXTERNAL getcgn
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/
C     ..
C     .. Executable Statements ..
C     Abort unless random number generator initialized
      IF (qrgnin()) GO TO 10
      call rexit("GETSD called before RNG initialized")

   10 CALL getcgn(g)
      iseed1 = cg1(g)
      iseed2 = cg2(g)
      RETURN

      END
      INTEGER FUNCTION ignbin(n,pp)
C**********************************************************************
C
C     INTEGER FUNCTION IGNBIN( N, PP )
C
C                    GENerate BINomial random deviate
C
C
C                              Function
C
C
C     Generates a single random deviate from a binomial
C     distribution whose number of trials is N and whose
C     probability of an event in each trial is P.
C
C
C                              Arguments
C
C
C     N  --> The number of trials in the binomial distribution
C            from which a random deviate is to be generated.
C                              INTEGER N
C     JJV                      (N >= 0)
C
C     PP --> The probability of an event in each trial of the
C            binomial distribution from which a random deviate
C            is to be generated.
C                              REAL PP
C     JJV                      (0.0 <= pp <= 1.0)
C
C     IGNBIN <-- A random deviate yielding the number of events
C                from N independent trials, each of which has
C                a probability of event P.
C                              INTEGER IGNBIN
C
C
C                              Note
C
C
C     Uses RANF so the value of the seeds, ISEED1 and ISEED2 must be set
C     by a call similar to the following
C          DUM = RANSET( ISEED1, ISEED2 )
C
C
C                              Method
C
C
C     This is algorithm BTPE from:
C
C         Kachitvichyanukul, V. and Schmeiser, B. W.
C
C         Binomial Random Variate Generation.
C         Communications of the ACM, 31, 2
C         (February, 1988) 216.
C
C**********************************************************************
C     SUBROUTINE BTPEC(N,PP,ISEED,JX)
C
C     BINOMIAL RANDOM VARIATE GENERATOR
C     MEAN .LT. 30 -- INVERSE CDF
C       MEAN .GE. 30 -- ALGORITHM BTPE:  ACCEPTANCE-REJECTION VIA
C       FOUR REGION COMPOSITION.  THE FOUR REGIONS ARE A TRIANGLE
C       (SYMMETRIC IN THE CENTER), A PAIR OF PARALLELOGRAMS (ABOVE
C       THE TRIANGLE), AND EXPONENTIAL LEFT AND RIGHT TAILS.
C
C     BTPE REFERS TO BINOMIAL-TRIANGLE-PARALLELOGRAM-EXPONENTIAL.
C     BTPEC REFERS TO BTPE AND "COMBINED."  THUS BTPE IS THE
C       RESEARCH AND BTPEC IS THE IMPLEMENTATION OF A COMPLETE
C       USABLE ALGORITHM.
C     REFERENCE:  VORATAS KACHITVICHYANUKUL AND BRUCE SCHMEISER,
C       "BINOMIAL RANDOM VARIATE GENERATION,"
C       COMMUNICATIONS OF THE ACM, FORTHCOMING
C     WRITTEN:  SEPTEMBER 1980.
C       LAST REVISED:  MAY 1985, JULY 1987
C     REQUIRED SUBPROGRAM:  RAND() -- A UNIFORM (0,1) RANDOM NUMBER
C                           GENERATOR
C     ARGUMENTS
C
C       N : NUMBER OF BERNOULLI TRIALS            (INPUT)
C       PP : PROBABILITY OF SUCCESS IN EACH TRIAL (INPUT)
C       ISEED:  RANDOM NUMBER SEED                (INPUT AND OUTPUT)
C       JX:  RANDOMLY GENERATED OBSERVATION       (OUTPUT)
C
C     VARIABLES
C       PSAVE: VALUE OF PP FROM THE LAST CALL TO BTPEC
C       NSAVE: VALUE OF N FROM THE LAST CALL TO BTPEC
C       XNP:  VALUE OF THE MEAN FROM THE LAST CALL TO BTPEC
C
C       P: PROBABILITY USED IN THE GENERATION PHASE OF BTPEC
C       FFM: TEMPORARY VARIABLE EQUAL TO XNP + P
C       M:  INTEGER VALUE OF THE CURRENT MODE
C       FM:  FLOATING POINT VALUE OF THE CURRENT MODE
C       XNPQ: TEMPORARY VARIABLE USED IN SETUP AND SQUEEZING STEPS
C       P1:  AREA OF THE TRIANGLE
C       C:  HEIGHT OF THE PARALLELOGRAMS
C       XM:  CENTER OF THE TRIANGLE
C       XL:  LEFT END OF THE TRIANGLE
C       XR:  RIGHT END OF THE TRIANGLE
C       AL:  TEMPORARY VARIABLE
C       XLL:  RATE FOR THE LEFT EXPONENTIAL TAIL
C       XLR:  RATE FOR THE RIGHT EXPONENTIAL TAIL
C       P2:  AREA OF THE PARALLELOGRAMS
C       P3:  AREA OF THE LEFT EXPONENTIAL TAIL
C       P4:  AREA OF THE RIGHT EXPONENTIAL TAIL
C       U:  A U(0,P4) RANDOM VARIATE USED FIRST TO SELECT ONE OF THE
C           FOUR REGIONS AND THEN CONDITIONALLY TO GENERATE A VALUE
C           FROM THE REGION
C       V:  A U(0,1) RANDOM NUMBER USED TO GENERATE THE RANDOM VALUE
C           (REGION 1) OR TRANSFORMED INTO THE VARIATE TO ACCEPT OR
C           REJECT THE CANDIDATE VALUE
C       IX:  INTEGER CANDIDATE VALUE
C       X:  PRELIMINARY CONTINUOUS CANDIDATE VALUE IN REGION 2 LOGIC
C           AND A FLOATING POINT IX IN THE ACCEPT/REJECT LOGIC
C       K:  ABSOLUTE VALUE OF (IX-M)
C       F:  THE HEIGHT OF THE SCALED DENSITY FUNCTION USED IN THE
C           ACCEPT/REJECT DECISION WHEN BOTH M AND IX ARE SMALL
C           ALSO USED IN THE INVERSE TRANSFORMATION
C       R: THE RATIO P/Q
C       G: CONSTANT USED IN CALCULATION OF PROBABILITY
C       MP:  MODE PLUS ONE, THE LOWER INDEX FOR EXPLICIT CALCULATION
C            OF F WHEN IX IS GREATER THAN M
C       IX1:  CANDIDATE VALUE PLUS ONE, THE LOWER INDEX FOR EXPLICIT
C             CALCULATION OF F WHEN IX IS LESS THAN M
C       I:  INDEX FOR EXPLICIT CALCULATION OF F FOR BTPE
C       AMAXP: MAXIMUM ERROR OF THE LOGARITHM OF NORMAL BOUND
C       YNORM: LOGARITHM OF NORMAL BOUND
C       ALV:  NATURAL LOGARITHM OF THE ACCEPT/REJECT VARIATE V
C
C       X1,F1,Z,W,Z2,X2,F2, AND W2 ARE TEMPORARY VARIABLES TO BE
C       USED IN THE FINAL ACCEPT/REJECT TEST
C
C       QN: PROBABILITY OF NO SUCCESS IN N TRIALS
C
C     REMARK
C       IX AND JX COULD LOGICALLY BE THE SAME VARIABLE, WHICH WOULD
C       SAVE A MEMORY POSITION AND A LINE OF CODE.  HOWEVER, SOME
C       COMPILERS (E.G.,CDC MNF) OPTIMIZE BETTER WHEN THE ARGUMENTS
C       ARE NOT INVOLVED.
C
C     ISEED NEEDS TO BE DOUBLE PRECISION IF THE IMSL ROUTINE
C     GGUBFS IS USED TO GENERATE UNIFORM RANDOM NUMBER, OTHERWISE
C     TYPE OF ISEED SHOULD BE DICTATED BY THE UNIFORM GENERATOR
C
C**********************************************************************

C
C
C
C*****DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY
C
C     ..
C     .. Scalar Arguments ..
      REAL pp
      INTEGER n
C     ..
C     .. Local Scalars ..
      REAL al,alv,amaxp,c,f,f1,f2,ffm,fm,g,p,p1,p2,p3,p4,psave,q,qn,r,u,
     +     v,w,w2,x,x1,x2,xl,xll,xlr,xm,xnp,xnpq,xr,ynorm,z,z2
      INTEGER i,ix,ix1,k,m,mp,nsave
C     ..
C     .. External Functions ..
      REAL ranf
      EXTERNAL ranf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,alog,amin1,iabs,int,sqrt
C     JJV ..
C     JJV .. Save statement ..
      SAVE p,q,m,fm,xnp,xnpq,p1,xm,xl,xr,c,xll,xlr,p2,p3,p4,qn,r,g,
     +     psave,nsave
C     JJV I am including the variables in data statements
C     ..
C     .. Data statements ..
C     JJV made these ridiculous starting values - the hope is that
C     JJV no one will call this the first time with them as args
      DATA psave,nsave/-1.0E37,-214748365/
C     ..
C     .. Executable Statements ..
      IF (pp.NE.psave) GO TO 10
      IF (n.NE.nsave) GO TO 20
      IF (xnp-30.) 150,30,30
C
C*****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE
C

C     JJV added the argument checker - involved only renaming 10
C     JJV and 20 to the checkers and adding checkers
C     JJV Only remaining problem - if called initially with the
C     JJV initial values of psave and nsave, it will hang
 10   IF (pp.LT.0.0) then
        call rexit("PP < 0.0 in IGNBIN - ABORT")
      end if

      IF (pp.GT.1.0)then
        call rexit("PP > 1.0 in IGNBIN - ABORT")
      end if

      psave = pp
      p = amin1(psave,1.-psave)
      q = 1. - p
 20   IF (n.LT.0) then
        call rexit("N < 0 in IGNBIN - ABORT")
      end if

      xnp = n*p
      nsave = n
      IF (xnp.LT.30.) GO TO 140
      ffm = xnp + p
      m = ffm
      fm = m
      xnpq = xnp*q
      p1 = int(2.195*sqrt(xnpq)-4.6*q) + 0.5
      xm = fm + 0.5
      xl = xm - p1
      xr = xm + p1
      c = 0.134 + 20.5/ (15.3+fm)
      al = (ffm-xl)/ (ffm-xl*p)
      xll = al* (1.+.5*al)
      al = (xr-ffm)/ (xr*q)
      xlr = al* (1.+.5*al)
      p2 = p1* (1.+c+c)
      p3 = p2 + c/xll
      p4 = p3 + c/xlr

   30 u = ranf()*p4
      v = ranf()
C
C     TRIANGULAR REGION
C
      IF (u.GT.p1) GO TO 40
      ix = xm - p1*v + u
      GO TO 170
C
C     PARALLELOGRAM REGION
C
   40 IF (u.GT.p2) GO TO 50
      x = xl + (u-p1)/c
      v = v*c + 1. - abs(xm-x)/p1
      IF (v.GT.1. .OR. v.LE.0.) GO TO 30
      ix = x
      GO TO 70
C
C     LEFT TAIL
C
   50 IF (u.GT.p3) GO TO 60
      ix = xl + alog(v)/xll
      IF (ix.LT.0) GO TO 30
      v = v* (u-p2)*xll
      GO TO 70
C
C     RIGHT TAIL
C
   60 ix = xr - alog(v)/xlr
      IF (ix.GT.n) GO TO 30
      v = v* (u-p3)*xlr
C
C*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST
C
   70 k = iabs(ix-m)
      IF (k.GT.20 .AND. k.LT.xnpq/2-1) GO TO 130
C
C     EXPLICIT EVALUATION
C
      f = 1.0
      r = p/q
      g = (n+1)*r
      IF (m-ix) 80,120,100
   80 mp = m + 1
      DO 90 i = mp,ix
          f = f* (g/i-r)
   90 CONTINUE
      GO TO 120

  100 ix1 = ix + 1
      DO 110 i = ix1,m
          f = f/ (g/i-r)
  110 CONTINUE
  120 IF (v-f) 170,170,30
C
C     SQUEEZING USING UPPER AND LOWER BOUNDS ON ALOG(F(X))
C
  130 amaxp = (k/xnpq)* ((k* (k/3.+.625)+.1666666666666)/xnpq+.5)
      ynorm = -k*k/ (2.*xnpq)
      alv = alog(v)
      IF (alv.LT.ynorm-amaxp) GO TO 170
      IF (alv.GT.ynorm+amaxp) GO TO 30
C
C     STIRLING'S FORMULA TO MACHINE ACCURACY FOR
C     THE FINAL ACCEPTANCE/REJECTION TEST
C
      x1 = ix + 1
      f1 = fm + 1.
      z = n + 1 - fm
      w = n - ix + 1.
      z2 = z*z
      x2 = x1*x1
      f2 = f1*f1
      w2 = w*w
      IF (alv- (xm*alog(f1/x1)+ (n-m+.5)*alog(z/w)+ (ix-
     +    m)*alog(w*p/ (x1*q))+ (13860.- (462.- (132.- (99.-
     +    140./f2)/f2)/f2)/f2)/f1/166320.+ (13860.- (462.- (132.- (99.-
     +    140./z2)/z2)/z2)/z2)/z/166320.+ (13860.- (462.- (132.- (99.-
     +    140./x2)/x2)/x2)/x2)/x1/166320.+ (13860.- (462.- (132.- (99.-
     +    140./w2)/w2)/w2)/w2)/w/166320.)) 170,170,30
C
C     INVERSE CDF LOGIC FOR MEAN LESS THAN 30
C
  140 qn = q**n
      r = p/q
      g = r* (n+1)
  150 ix = 0
      f = qn
      u = ranf()
  160 IF (u.LT.f) GO TO 170
      IF (ix.GT.110) GO TO 150
      u = u - f
      ix = ix + 1
      f = f* (g/ix-r)
      GO TO 160

  170 IF (psave.GT.0.5) ix = n - ix
      ignbin = ix
      RETURN

      END
      INTEGER FUNCTION ignlgi()
C**********************************************************************
C
C     INTEGER FUNCTION IGNLGI()
C               GeNerate LarGe Integer
C
C     Returns a random integer following a uniform distribution over
C     (1, 2147483562) using the current generator.
C
C     This is a transcription from Pascal to Fortran of routine
C     Random from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER curntg,k,s1,s2,z
      LOGICAL qqssd
C     ..
C     .. External Functions ..
      LOGICAL qrgnin
      EXTERNAL qrgnin
C     ..
C     .. External Subroutines ..
      EXTERNAL getcgn,inrgcm,rgnqsd,setall
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/
C     ..
C     .. Executable Statements ..
C
C     IF THE RANDOM NUMBER PACKAGE HAS NOT BEEN INITIALIZED YET, DO SO.
C     IT CAN BE INITIALIZED IN ONE OF TWO WAYS : 1) THE FIRST CALL TO
C     THIS ROUTINE  2) A CALL TO SETALL.
C
      IF (.NOT. (qrgnin())) CALL inrgcm()
      CALL rgnqsd(qqssd)
      IF (.NOT. (qqssd)) CALL setall(1234567890,123456789)
C
C     Get Current Generator
C
      CALL getcgn(curntg)
      s1 = cg1(curntg)
      s2 = cg2(curntg)
      k = s1/53668
      s1 = a1* (s1-k*53668) - k*12211
      IF (s1.LT.0) s1 = s1 + m1
      k = s2/52774
      s2 = a2* (s2-k*52774) - k*3791
      IF (s2.LT.0) s2 = s2 + m2
      cg1(curntg) = s1
      cg2(curntg) = s2
      z = s1 - s2
      IF (z.LT.1) z = z + m1 - 1
      IF (qanti(curntg)) z = m1 - z
      ignlgi = z
      RETURN

      END
      INTEGER FUNCTION ignnbn(n,p)
C**********************************************************************
C
C     INTEGER FUNCTION IGNNBN( N, P )
C
C                GENerate Negative BiNomial random deviate
C
C
C                              Function
C
C
C     Generates a single random deviate from a negative binomial
C     distribution.
C
C
C                              Arguments
C
C
C     N  --> Required number of events.
C                              INTEGER N
C     JJV                      (N > 0)
C
C     P  --> The probability of an event during a Bernoulli trial.
C                              REAL P
C     JJV                      (0.0 < P < 1.0)
C
C
C
C                              Method
C
C
C     Algorithm from page 480 of
C
C     Devroye, Luc
C
C     Non-Uniform Random Variate Generation.  Springer-Verlag,
C     New York, 1986.
C
C**********************************************************************
C     ..
C     .. Scalar Arguments ..
      REAL p
      INTEGER n
C     ..
C     .. Local Scalars ..
      REAL y,a,r
C     ..
C     .. External Functions ..
C     JJV changed to call SGAMMA directly
C     REAL gengam
      REAL sgamma
      INTEGER ignpoi
C      EXTERNAL gengam,ignpoi
      EXTERNAL sgamma,ignpoi
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC real
C     ..
C     .. Executable Statements ..
C     Check Arguments
C     JJV changed argumnet checker to abort if N <= 0
      IF (n.LE.0)then
        call rexit("N <= 0 in IGNNBN")
      end if

      IF (p.LE.0.0)then
        call rexit("P <= 0.0 in IGNNBN")
      end if

      IF (p.GE.1.0)then
        call rexit("P >= 1.0 in IGNNBN")
      end if


C     Generate Y, a random gamma (n,(1-p)/p) variable
C     JJV Note: the above parametrization is consistent with Devroye,
C     JJV       but gamma (p/(1-p),n) is the equivalent in our code
      r = real(n)
      a = p/ (1.0-p)
C      y = gengam(a,r)
      y = sgamma(r)/a

C     Generate a random Poisson(y) variable
      ignnbn = ignpoi(y)
      RETURN

      END
      INTEGER FUNCTION ignpoi(mu)
C**********************************************************************
C
C     INTEGER FUNCTION IGNPOI( MU )
C
C                    GENerate POIsson random deviate
C
C
C                              Function
C
C
C     Generates a single random deviate from a Poisson
C     distribution with mean MU.
C
C
C                              Arguments
C
C
C     MU --> The mean of the Poisson distribution from which
C            a random deviate is to be generated.
C                              REAL MU
C     JJV                    (MU >= 0.0)
C
C     IGNPOI <-- The random deviate.
C                              INTEGER IGNPOI (non-negative)
C
C
C                              Method
C
C
C     Renames KPOIS from TOMS as slightly modified by BWB to use RANF
C     instead of SUNIF.
C
C     For details see:
C
C               Ahrens, J.H. and Dieter, U.
C               Computer Generation of Poisson Deviates
C               From Modified Normal Distributions.
C               ACM Trans. Math. Software, 8, 2
C               (June 1982),163-179
C
C**********************************************************************
C**********************************************************************C
C**********************************************************************C
C                                                                      C
C                                                                      C
C     P O I S S O N  DISTRIBUTION                                      C
C                                                                      C
C                                                                      C
C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     FOR DETAILS SEE:                                                 C
C                                                                      C
C               AHRENS, J.H. AND DIETER, U.                            C
C               COMPUTER GENERATION OF POISSON DEVIATES                C
C               FROM MODIFIED NORMAL DISTRIBUTIONS.                    C
C               ACM TRANS. MATH. SOFTWARE, 8,2 (JUNE 1982), 163 - 179. C
C                                                                      C
C     (SLIGHTLY MODIFIED VERSION OF THE PROGRAM IN THE ABOVE ARTICLE)  C
C                                                                      C
C**********************************************************************C
C
C      INTEGER FUNCTION IGNPOI(IR,MU)
C
C     INPUT:  IR=CURRENT STATE OF BASIC RANDOM NUMBER GENERATOR
C             MU=MEAN MU OF THE POISSON DISTRIBUTION
C     OUTPUT: IGNPOI=SAMPLE FROM THE POISSON-(MU)-DISTRIBUTION
C
C
C
C     MUPREV=PREVIOUS MU, MUOLD=MU AT LAST EXECUTION OF STEP P OR CASE B
C     TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
C     COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL
C
C
C
C     SEPARATION OF CASES A AND B
C
C     .. Scalar Arguments ..
      REAL mu
C     ..
C     .. Local Scalars ..
      REAL a0,a1,a2,a3,a4,a5,a6,a7,b1,b2,c,c0,c1,c2,c3,d,del,difmuk,e,
     +     fk,fx,fy,g,muold,muprev,omega,p,p0,px,py,q,s,t,u,v,x,xx
C     JJV I added a variable 'll' here - it is the 'l' for CASE A
      INTEGER j,k,kflag,l,ll,m
C     ..
C     .. Local Arrays ..
      REAL fact(10),pp(35)
C     ..
C     .. External Functions ..
      REAL ranf,sexpo,snorm
      EXTERNAL ranf,sexpo,snorm
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,alog,exp,float,ifix,max0,min0,sign,sqrt
C     ..
C     JJV added this for case: mu unchanged
C     .. Save statement ..
      SAVE s, d, l, ll, omega, c3, c2, c1, c0, c, m, p, q, p0,
     +     a0, a1, a2, a3, a4, a5, a6, a7, fact, muprev, muold
C     ..
C     JJV end addition - I am including vars in Data statements
C     .. Data statements ..
C     JJV changed initial values of MUPREV and MUOLD to -1.0E37
C     JJV if no one calls IGNPOI with MU = -1.0E37 the first time,
C     JJV the code shouldn't break
      DATA muprev,muold/-1.0E37,-1.0E37/
      DATA a0,a1,a2,a3,a4,a5,a6,a7/-.5,.3333333,-.2500068,.2000118,
     +     -.1661269,.1421878,-.1384794,.1250060/
      DATA fact/1.,1.,2.,6.,24.,120.,720.,5040.,40320.,362880./
C     ..
C     .. Executable Statements ..

      ignpoi=0
      difmuk=0.0
      e=0.0
      fk=0.0
      u=0.0
      
      IF (mu.EQ.muprev) GO TO 10
      IF (mu.LT.10.0) GO TO 120
C
C     C A S E  A. (RECALCULATION OF S,D,LL IF MU HAS CHANGED)
C
C     JJV This is the case where I changed 'l' to 'll'
C     JJV Here 'll' is set once and used in a comparison once

      muprev = mu
      s = sqrt(mu)
      d = 6.0*mu*mu
C
C             THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
C             PROBABILITIES FK WHENEVER K >= M(MU). LL=IFIX(MU-1.1484)
C             IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .
C
      ll = ifix(mu-1.1484)
C
C     STEP N. NORMAL SAMPLE - SNORM(IR) FOR STANDARD NORMAL DEVIATE
C
   10 g = mu + s*snorm()
      IF (g.LT.0.0) GO TO 20
      ignpoi = ifix(g)
C
C     STEP I. IMMEDIATE ACCEPTANCE IF IGNPOI IS LARGE ENOUGH
C
      IF (ignpoi.GE.ll) RETURN
C
C     STEP S. SQUEEZE ACCEPTANCE - SUNIF(IR) FOR (0,1)-SAMPLE U
C
      fk = float(ignpoi)
      difmuk = mu - fk
      u = ranf()
      IF (d*u.GE.difmuk*difmuk*difmuk) RETURN
C
C     STEP P. PREPARATIONS FOR STEPS Q AND H.
C             (RECALCULATIONS OF PARAMETERS IF NECESSARY)
C             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
C             THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
C             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
C             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.
C
   20 IF (mu.EQ.muold) GO TO 30
      muold = mu
      omega = .3989423/s
      b1 = .4166667E-1/mu
      b2 = .3*b1*b1
      c3 = .1428571*b1*b2
      c2 = b2 - 15.*c3
      c1 = b1 - 6.*b2 + 45.*c3
      c0 = 1. - b1 + 3.*b2 - 15.*c3
      c = .1069/mu
   30 IF (g.LT.0.0) GO TO 50
C
C             'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)
C
      kflag = 0
      GO TO 70
C
C     STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)
C
   40 IF (fy-u*fy.LE.py*exp(px-fx)) RETURN
C
C     STEP E. EXPONENTIAL SAMPLE - SEXPO(IR) FOR STANDARD EXPONENTIAL
C             DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
C             (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)
C
   50 e = sexpo()
      u = ranf()
      u = u + u - 1.0
      t = 1.8 + sign(e,u)
      IF (t.LE. (-.6744)) GO TO 50
      ignpoi = ifix(mu+s*t)
      fk = float(ignpoi)
      difmuk = mu - fk
C
C             'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)
C
      kflag = 1
      GO TO 70
C
C     STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)
C
   60 IF (c*abs(u).GT.py*exp(px+e)-fy*exp(fx+e)) GO TO 50
      RETURN
C
C     STEP F. 'SUBROUTINE' F. CALCULATION OF PX,PY,FX,FY.
C             CASE IGNPOI .LT. 10 USES FACTORIALS FROM TABLE FACT
C
   70 IF (ignpoi.GE.10) GO TO 80
      px = -mu
      py = mu**ignpoi/fact(ignpoi+1)
      GO TO 110
C
C             CASE IGNPOI .GE. 10 USES POLYNOMIAL APPROXIMATION
C             A0-A7 FOR ACCURACY WHEN ADVISABLE
C             .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)
C
   80 del = .8333333E-1/fk
      del = del - 4.8*del*del*del
      v = difmuk/fk
      IF (abs(v).LE.0.25) GO TO 90
      px = fk*alog(1.0+v) - difmuk - del
      GO TO 100

   90 px = fk*v*v* (((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0) -
     +     del
  100 py = .3989423/sqrt(fk)
  110 x = (0.5-difmuk)/s
      xx = x*x
      fx = -0.5*xx
      fy = omega* (((c3*xx+c2)*xx+c1)*xx+c0)
      IF (kflag) 40,40,60
C
C     C A S E  B. (START NEW TABLE AND CALCULATE P0 IF NECESSARY)
C
C     JJV changed MUPREV assignment from 0.0 to initial value
  120 muprev = -1.0E37
      IF (mu.EQ.muold) GO TO 130
C     JJV added argument checker here
      IF (mu.GE.0.0) GO TO 125
      call rexit("MU < 0 in IGNPOI - ABORT")

C     JJV added line label here
 125  muold = mu
      m = max0(1,ifix(mu))
      l = 0
      p = exp(-mu)
      q = p
      p0 = p
C
C     STEP U. UNIFORM SAMPLE FOR INVERSION METHOD
C
  130 u = ranf()
      ignpoi = 0
      IF (u.LE.p0) RETURN
C
C     STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
C             PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
C             (0.458=PP(9) FOR MU=10)
C
      IF (l.EQ.0) GO TO 150
      j = 1
      IF (u.GT.0.458) j = min0(l,m)
      DO 140 k = j,l
          IF (u.LE.pp(k)) GO TO 180
  140 CONTINUE
      IF (l.EQ.35) GO TO 130
C
C     STEP C. CREATION OF NEW POISSON PROBABILITIES P
C             AND THEIR CUMULATIVES Q=PP(K)
C
  150 l = l + 1
      DO 160 k = l,35
          p = p*mu/float(k)
          q = q + p
          pp(k) = q
          IF (u.LE.q) GO TO 170
  160 CONTINUE
      l = 35
      GO TO 130

  170 l = k
  180 ignpoi = k
      RETURN

      END
      INTEGER FUNCTION ignuin(low,high)
C**********************************************************************
C
C     INTEGER FUNCTION IGNUIN( LOW, HIGH )
C
C               GeNerate Uniform INteger
C
C
C                              Function
C
C
C     Generates an integer uniformly distributed between LOW and HIGH.
C
C
C                              Arguments
C
C
C     LOW --> Low bound (inclusive) on integer value to be generated
C                         INTEGER LOW
C
C     HIGH --> High bound (inclusive) on integer value to be generated
C                         INTEGER HIGH
C
C
C                              Note
C
C
C     If (HIGH-LOW) > 2,147,483,561 prints error message on * unit and
C     stops the program.
C
C**********************************************************************

C     IGNLGI generates integers between 1 and 2147483562
C     MAXNUM is 1 less than maximum generable value
C     .. Parameters ..
      INTEGER maxnum
      PARAMETER (maxnum=2147483561)
      CHARACTER*(*) err1,err2
      PARAMETER (err1='LOW > HIGH in IGNUIN',
     +          err2=' ( HIGH - LOW ) > 2,147,483,561 in IGNUIN')
C     ..
C     .. Scalar Arguments ..
      INTEGER high,low
C     ..
C     .. Local Scalars ..
      INTEGER err,ign,maxnow,range,ranp1
C     ..
C     .. External Functions ..
      INTEGER ignlgi
      EXTERNAL ignlgi
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC mod
C     ..
C     .. Executable Statements ..
      IF (.NOT. (low.GT.high)) GO TO 10
      err = 1
C      ABORT-PROGRAM
      GO TO 80

   10 range = high - low
      IF (.NOT. (range.GT.maxnum)) GO TO 20
      err = 2
C      ABORT-PROGRAM
      GO TO 80

   20 IF (.NOT. (low.EQ.high)) GO TO 30
      ignuin = low
      RETURN

C     Number to be generated should be in range 0..RANGE
C     Set MAXNOW so that the number of integers in 0..MAXNOW is an
C     integral multiple of the number in 0..RANGE

   30 ranp1 = range + 1
      maxnow = (maxnum/ranp1)*ranp1
   40 ign = ignlgi() - 1
      IF (.NOT. (ign.LE.maxnow)) GO TO 40
      ignuin = low + mod(ign,ranp1)
      RETURN

   80 IF (.NOT. (err.EQ.1)) GO TO 90
      call rexit("LOW > HIGH in IGNUIN")
      GO TO 100

C     TO ABORT-PROGRAM
   90 continue
  100 continue
      IF (.NOT. (err.EQ.1)) GO TO 110
      call rexit("LOW > HIGH in IGNUIN")

  110 call rexit("( HIGH - LOW ) > 2,147,483,561 in IGNUIN")

      END
      SUBROUTINE initgn(isdtyp)
C**********************************************************************
C
C     SUBROUTINE INITGN(ISDTYP)
C          INIT-ialize current G-e-N-erator
C
C     Reinitializes the state of the current generator
C
C     This is a transcription from Pascal to Fortran of routine
C     Init_Generator from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C
C                              Arguments
C
C
C     ISDTYP -> The state to which the generator is to be set
C
C          ISDTYP = -1  => sets the seeds to their initial value
C          ISDTYP =  0  => sets the seeds to the first value of
C                          the current block
C          ISDTYP =  1  => sets the seeds to the first value of
C                          the next block
C
C                                   INTEGER ISDTYP
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalar Arguments ..
      INTEGER isdtyp
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER g
C     ..
C     .. External Functions ..
      LOGICAL qrgnin
      INTEGER mltmod
      EXTERNAL qrgnin,mltmod
C     ..
C     .. External Subroutines ..
      EXTERNAL getcgn
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/
C     ..
C     .. Executable Statements ..
C     Abort unless random number generator initialized
      IF (qrgnin()) GO TO 10
      call rexit("INITGN called before RNG initialized")

   10 CALL getcgn(g)
      IF ((-1).NE. (isdtyp)) GO TO 20
      lg1(g) = ig1(g)
      lg2(g) = ig2(g)
      GO TO 50

   20 IF ((0).NE. (isdtyp)) GO TO 30
      CONTINUE
      GO TO 50
C     do nothing
   30 IF ((1).NE. (isdtyp)) GO TO 40
      lg1(g) = mltmod(a1w,lg1(g),m1)
      lg2(g) = mltmod(a2w,lg2(g),m2)
      GO TO 50

   40 call rexit("ISDTYP NOT IN RANGE")

   50 cg1(g) = lg1(g)
      cg2(g) = lg2(g)
      RETURN

      END
      SUBROUTINE inrgcm()
C**********************************************************************
C
C     SUBROUTINE INRGCM()
C          INitialize Random number Generator CoMmon
C
C
C                              Function
C
C
C     Initializes common area  for random number  generator.  This saves
C     the  nuisance  of  a  BLOCK DATA  routine  and the  difficulty  of
C     assuring that the routine is loaded with the other routines.
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER i
      LOGICAL qdum
C     ..
C     .. External Functions ..
      LOGICAL qrgnsn
      EXTERNAL qrgnsn
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/
C     ..
C     .. Executable Statements ..
C     V=20;                            W=30;
C
C     A1W = MOD(A1**(2**W),M1)         A2W = MOD(A2**(2**W),M2)
C     A1VW = MOD(A1**(2**(V+W)),M1)    A2VW = MOD(A2**(2**(V+W)),M2)
C
C   If V or W is changed A1W, A2W, A1VW, and A2VW need to be recomputed.
C    An efficient way to precompute a**(2*j) MOD m is to start with
C    a and square it j times modulo m using the function MLTMOD.
C
      m1 = 2147483563
      m2 = 2147483399
      a1 = 40014
      a2 = 40692
      a1w = 1033780774
      a2w = 1494757890
      a1vw = 2082007225
      a2vw = 784306273
      DO 10,i = 1,numg
          qanti(i) = .FALSE.
   10 CONTINUE
C
C     Tell the world that common has been initialized
C
      qdum = qrgnsn(.TRUE.)
      RETURN

      END
      INTEGER FUNCTION lennob(string)
      IMPLICIT INTEGER (a-p,r-z),LOGICAL (q)
C**********************************************************************
C
C     INTEGER FUNCTION LENNOB( STRING )
C                LENgth NOt counting trailing Blanks
C
C
C                              Function
C
C
C     Returns the length of STRING up to and including the last
C     non-blank character.
C
C
C                              Arguments
C
C
C     STRING --> String whose length not counting trailing blanks
C                is returned.
C
C**********************************************************************
      CHARACTER*(*) string

      end = len(string)
      DO 20,i = end,1,-1
          IF (.NOT. (string(i:i).NE.' ')) GO TO 10
          lennob = i
          RETURN

   10     CONTINUE
   20 CONTINUE
      lennob = 0
      RETURN

      END
      INTEGER FUNCTION mltmod(a,s,m)
C**********************************************************************
C
C     INTEGER FUNCTION MLTMOD(A,S,M)
C
C                    Returns (A*S) MOD M
C
C     This is a transcription from Pascal to Fortran of routine
C     MULtMod_Decompos from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C
C                              Arguments
C
C
C     A, S, M  -->
C                         INTEGER A,S,M
C
C**********************************************************************
C     .. Parameters ..
      INTEGER h
      PARAMETER (h=32768)
C     ..
C     .. Scalar Arguments ..
      INTEGER a,m,s
C     ..
C     .. Local Scalars ..
      INTEGER a0,a1,k,p,q,qh,rh
C     ..
C     .. Executable Statements ..
C
C     H = 2**((b-2)/2) where b = 32 because we are using a 32 bit
C      machine. On a different machine recompute H
C
      IF (.NOT. (a.LE.0.OR.a.GE.m.OR.s.LE.0.OR.s.GE.m)) GO TO 10
      call rexit("A, M, S out of order in MLTMOD - ABORT")

   10 IF (.NOT. (a.LT.h)) GO TO 20
      a0 = a
      p = 0
      GO TO 120

   20 a1 = a/h
      a0 = a - h*a1
      qh = m/h
      rh = m - h*qh
      IF (.NOT. (a1.GE.h)) GO TO 50
      a1 = a1 - h
      k = s/qh
      p = h* (s-k*qh) - k*rh
   30 IF (.NOT. (p.LT.0)) GO TO 40
      p = p + m
      GO TO 30

   40 GO TO 60

   50 p = 0
C
C     P = (A2*S*H)MOD M
C
   60 IF (.NOT. (a1.NE.0)) GO TO 90
      q = m/a1
      k = s/q
      p = p - k* (m-a1*q)
      IF (p.GT.0) p = p - m
      p = p + a1* (s-k*q)
   70 IF (.NOT. (p.LT.0)) GO TO 80
      p = p + m
      GO TO 70

   80 CONTINUE
   90 k = p/qh
C
C     P = ((A2*H + A1)*S)MOD M
C
      p = h* (p-k*qh) - k*rh
  100 IF (.NOT. (p.LT.0)) GO TO 110
      p = p + m
      GO TO 100

  110 CONTINUE
  120 IF (.NOT. (a0.NE.0)) GO TO 150
C
C     P = ((A2*H + A1)*H*S)MOD M
C
      q = m/a0
      k = s/q
      p = p - k* (m-a0*q)
      IF (p.GT.0) p = p - m
      p = p + a0* (s-k*q)
  130 IF (.NOT. (p.LT.0)) GO TO 140
      p = p + m
      GO TO 130

  140 CONTINUE
  150 mltmod = p
C
      RETURN

      END
      SUBROUTINE phrtsd(phrase,seed1,seed2)
C**********************************************************************
C
C     SUBROUTINE PHRTSD( PHRASE, SEED1, SEED2 )
C               PHRase To SeeDs
C
C
C                              Function
C
C
C     Uses a phrase (character string) to generate two seeds for the RGN
C     random number generator.
C
C
C                              Arguments
C
C
C     PHRASE --> Phrase to be used for random number generation
C                         CHARACTER*(*) PHRASE
C
C     SEED1 <-- First seed for RGN generator
C                         INTEGER SEED1
C
C     SEED2 <-- Second seed for RGN generator
C                         INTEGER SEED2
C
C
C                              Note
C
C
C     Trailing blanks are eliminated before the seeds are generated.
C
C     Generated seed values will fall in the range 1..2^30
C     (1..1,073,741,824)
C
C**********************************************************************
C     .. Parameters ..
      CHARACTER*(*) table
      PARAMETER (table='abcdefghijklmnopqrstuvwxyz'//
     +          'ABCDEFGHIJKLMNOPQRSTUVWXYZ'//'0123456789'//
     +          '!@#$%^&*()_+[];:''"<>?,./')
      INTEGER twop30
      PARAMETER (twop30=1073741824)
C     ..
C     .. Scalar Arguments ..
      INTEGER seed1,seed2
      CHARACTER phrase* (*)
C     ..
C     .. Local Scalars ..
      INTEGER i,ichr,j,lphr
C     ..
C     .. Local Arrays ..
      INTEGER shift(0:4),values(5)
C     ..
C     .. External Functions ..
      INTEGER lennob
      EXTERNAL lennob
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC index,mod
C     ..
C     JJV added Save statement for variable in Data statement 
C     .. Save statements ..
      SAVE shift
C     JJV end addition 
C     .. 
C     .. Data statements ..
      DATA shift/1,64,4096,262144,16777216/
C     ..
C     .. Executable Statements ..
      seed1 = 1234567890
      seed2 = 123456789
      lphr = lennob(phrase)
      IF (lphr.LT.1) RETURN
      DO 30,i = 1,lphr
          ichr = mod(index(table,phrase(i:i)),64)
          IF (ichr.EQ.0) ichr = 63
          DO 10,j = 1,5
              values(j) = ichr - j
              IF (values(j).LT.1) values(j) = values(j) + 63
   10     CONTINUE
          DO 20,j = 1,5
              seed1 = mod(seed1+shift(j-1)*values(j),twop30)
              seed2 = mod(seed2+shift(j-1)*values(6-j),twop30)
   20     CONTINUE
   30 CONTINUE
      RETURN

      END
      LOGICAL FUNCTION qrgnin()
C**********************************************************************
C
C     LOGICAL FUNCTION QRGNIN()
C               Q Random GeNerators INitialized?
C
C     A trivial routine to determine whether or not the random
C     number generator has been initialized.  Returns .TRUE. if
C     it has, else .FALSE.
C
C**********************************************************************
C     .. Scalar Arguments ..
      LOGICAL qvalue
C     ..
C     .. Local Scalars ..
      LOGICAL qinit
C     ..
C     .. Entry Points ..
      LOGICAL qrgnsn
C     ..
C     .. Save statement ..
      SAVE qinit
C     ..
C     .. Data statements ..
      DATA qinit/.FALSE./
C     ..
C     .. Executable Statements ..
      qrgnin = qinit
      RETURN

      ENTRY qrgnsn(qvalue)
C**********************************************************************
C
C     LOGICAL FUNCTION QRGNSN( QVALUE )
C               Q Random GeNerators Set whether iNitialized
C
C     Sets state of whether random number generator is initialized
C     to QVALUE.
C
C     This routine is actually an entry in QRGNIN, hence it is a
C     logical function.  It returns the (meaningless) value .TRUE.
C
C**********************************************************************
      qinit = qvalue
      qrgnsn = .TRUE.
      RETURN

      END
      REAL FUNCTION ranf()
C**********************************************************************
C
C     REAL FUNCTION RANF()
C                RANDom number generator as a Function
C
C     Returns a random floating point number from a uniform distribution
C     over 0 - 1 (endpoints of this interval are not returned) using the
C     current generator
C
C     This is a transcription from Pascal to Fortran of routine
C     Uniform_01 from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C**********************************************************************
C     .. External Functions ..
      INTEGER ignlgi
      EXTERNAL ignlgi
C     ..
C     .. Executable Statements ..
C
C     4.656613057E-10 is 1/M1  M1 is set in a data statement in IGNLGI
C      and is currently 2147483563. If M1 changes, change this also.
C
      ranf = ignlgi()*4.656613057E-10
      RETURN

      END
      REAL FUNCTION sdot(n,sx,incx,sy,incy)
      REAL sx(1),sy(1),stemp
      INTEGER i,incx,incy,ix,iy,m,mp1,n

      stemp = 0.0E0
      sdot = 0.0E0
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) GO TO 20
      ix = 1
      iy = 1
      IF (incx.LT.0) ix = (-n+1)*incx + 1
      IF (incy.LT.0) iy = (-n+1)*incy + 1
      DO 10 i = 1,n
          stemp = stemp + sx(ix)*sy(iy)
          ix = ix + incx
          iy = iy + incy
   10 CONTINUE
      sdot = stemp
      RETURN

   20 m = mod(n,5)
      IF (m.EQ.0) GO TO 40
      DO 30 i = 1,m
          stemp = stemp + sx(i)*sy(i)
   30 CONTINUE
      IF (n.LT.5) GO TO 60
   40 mp1 = m + 1
      DO 50 i = mp1,n,5
          stemp = stemp + sx(i)*sy(i) + sx(i+1)*sy(i+1) +
     +            sx(i+2)*sy(i+2) + sx(i+3)*sy(i+3) + sx(i+4)*sy(i+4)
   50 CONTINUE
   60 sdot = stemp
      RETURN

      END
      SUBROUTINE setall(iseed1,iseed2)
C**********************************************************************
C
C      SUBROUTINE SETALL(ISEED1,ISEED2)
C               SET ALL random number generators
C
C     Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
C     initial seeds of the other generators are set accordingly, and
C     all generators states are set to these seeds.
C
C     This is a transcription from Pascal to Fortran of routine
C     Set_Initial_Seed from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C
C                              Arguments
C
C
C     ISEED1 -> First of two integer seeds
C                                   INTEGER ISEED1
C
C     ISEED2 -> Second of two integer seeds
C                                   INTEGER ISEED1
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalar Arguments ..
      INTEGER iseed1,iseed2
      LOGICAL qssd
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER g,ocgn
      LOGICAL qqssd
C     ..
C     .. External Functions ..
      INTEGER mltmod
      LOGICAL qrgnin
      EXTERNAL mltmod,qrgnin
C     ..
C     .. External Subroutines ..
      EXTERNAL getcgn,initgn,inrgcm,setcgn
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/,qqssd
C     ..
C     .. Data statements ..
      DATA qqssd/.FALSE./
C     ..
C     .. Executable Statements ..
C
C     TELL IGNLGI, THE ACTUAL NUMBER GENERATOR, THAT THIS ROUTINE
C      HAS BEEN CALLED.
C
      qqssd = .TRUE.
      CALL getcgn(ocgn)
C
C     Initialize Common Block if Necessary
C
      IF (.NOT. (qrgnin())) CALL inrgcm()
      ig1(1) = iseed1
      ig2(1) = iseed2
      CALL initgn(-1)
      DO 10,g = 2,numg
          ig1(g) = mltmod(a1vw,ig1(g-1),m1)
          ig2(g) = mltmod(a2vw,ig2(g-1),m2)
          CALL setcgn(g)
          CALL initgn(-1)
   10 CONTINUE
      CALL setcgn(ocgn)
      RETURN

      ENTRY rgnqsd(qssd)
C**********************************************************************
C
C     SUBROUTINE RGNQSD
C                    Random Number Generator Query SeeD set?
C
C     Returns (LOGICAL) QSSD as .TRUE. if SETALL has been invoked,
C     otherwise returns .FALSE.
C
C**********************************************************************
      qssd = qqssd
      RETURN

      END
      SUBROUTINE setant(qvalue)
C**********************************************************************
C
C      SUBROUTINE SETANT(QVALUE)
C               SET ANTithetic
C
C     Sets whether the current generator produces antithetic values.  If
C     X   is  the value  normally returned  from  a uniform [0,1] random
C     number generator then 1  - X is the antithetic  value. If X is the
C     value  normally  returned  from a   uniform  [0,N]  random  number
C     generator then N - 1 - X is the antithetic value.
C
C     All generators are initialized to NOT generate antithetic values.
C
C     This is a transcription from Pascal to Fortran of routine
C     Set_Antithetic from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C
C                              Arguments
C
C
C     QVALUE -> .TRUE. if generator G is to generating antithetic
C                    values, otherwise .FALSE.
C                                   LOGICAL QVALUE
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalar Arguments ..
      LOGICAL qvalue
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER g
C     ..
C     .. External Functions ..
      LOGICAL qrgnin
      EXTERNAL qrgnin
C     ..
C     .. External Subroutines ..
      EXTERNAL getcgn
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/
C     ..
C     .. Executable Statements ..
C     Abort unless random number generator initialized
      IF (qrgnin()) GO TO 10
      call rexit("SETANT called before RNG initialized")

   10 CALL getcgn(g)
      qanti(g) = qvalue
      RETURN

      END
      SUBROUTINE setgmn(meanv,covm,ldcovm,p,parm)
C      SUBROUTINE setgmn(meanv,covm,p,parm)
C     JJV changed this routine to take leading dimension of COVM
C     JJV argument and pass it to SPOFA, making it easier to use
C     JJV if the COVM which is used is contained in a larger matrix
C     JJV and to make the routine more consistent with LINPACK.
C     JJV Changes are in comments, declarations, and the call to SPOFA.
C**********************************************************************
C
C     SUBROUTINE SETGMN( MEANV, COVM, LDCOVM, P, PARM)
C            SET Generate Multivariate Normal random deviate
C
C
C                              Function
C
C
C      Places P, MEANV, and the Cholesky factoriztion of COVM
C      in PARM for GENMN.
C
C
C                              Arguments
C
C
C     MEANV --> Mean vector of multivariate normal distribution.
C                                        REAL MEANV(P)
C
C     COVM   <--> (Input) Covariance   matrix    of  the  multivariate
C                 normal distribution.  This routine uses only the
C                 (1:P,1:P) slice of COVM, but needs to know LDCOVM.
C
C                 (Output) Destroyed on output
C                                        REAL COVM(LDCOVM,P)
C
C     LDCOVM --> Leading actual dimension of COVM.
C                                        INTEGER LDCOVM
C
C     P     --> Dimension of the normal, or length of MEANV.
C                                        INTEGER P
C
C     PARM <-- Array of parameters needed to generate multivariate
C                normal deviates (P, MEANV and Cholesky decomposition
C                of COVM).
C                1 : 1                - P
C                2 : P + 1            - MEANV
C                P+2 : P*(P+3)/2 + 1  - Cholesky decomposition of COVM
C                                             REAL PARM(P*(P+3)/2 + 1)
C
C**********************************************************************
C     .. Scalar Arguments ..
C      INTEGER p
      INTEGER p, ldcovm
C     ..
C     .. Array Arguments ..
C      REAL covm(p,p),meanv(p),parm(p* (p+3)/2+1)
      REAL covm(ldcovm,p),meanv(p),parm(p* (p+3)/2+1)
C     ..
C     .. Local Scalars ..
      INTEGER i,icount,info,j
C     ..
C     .. External Subroutines ..
      EXTERNAL spofa
C     ..
C     .. Executable Statements ..
C
C
C     TEST THE INPUT
C
      IF (.NOT. (p.LE.0)) GO TO 10
      call rexit("P nonpositive in SETGMN")

   10 parm(1) = p
C
C     PUT P AND MEANV INTO PARM
C
      DO 20,i = 2,p + 1
          parm(i) = meanv(i-1)
   20 CONTINUE
C
C      Cholesky decomposition to find A s.t. trans(A)*(A) = COVM
C
C      CALL spofa(covm,p,p,info)
      CALL spofa(covm,ldcovm,p,info)
      IF (.NOT. (info.NE.0)) GO TO 30
      call rexit("COVM not positive definite in SETGMN")

   30 icount = p + 1
C
C     PUT UPPER HALF OF A, WHICH IS NOW THE CHOLESKY FACTOR, INTO PARM
C          COVM(1,1) = PARM(P+2)
C          COVM(1,2) = PARM(P+3)
C                    :
C          COVM(1,P) = PARM(2P+1)
C          COVM(2,2) = PARM(2P+2)  ...
C
      DO 50,i = 1,p
          DO 40,j = i,p
              icount = icount + 1
              parm(icount) = covm(i,j)
   40     CONTINUE
   50 CONTINUE
      RETURN
C
      END
      SUBROUTINE setsd(iseed1,iseed2)
C**********************************************************************
C
C     SUBROUTINE SETSD(ISEED1,ISEED2)
C               SET S-ee-D of current generator
C
C     Resets the initial  seed of  the current  generator to  ISEED1 and
C     ISEED2. The seeds of the other generators remain unchanged.
C
C     This is a transcription from Pascal to Fortran of routine
C     Set_Seed from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C
C                              Arguments
C
C
C     ISEED1 -> First integer seed
C                                   INTEGER ISEED1
C
C     ISEED2 -> Second integer seed
C                                   INTEGER ISEED1
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalar Arguments ..
      INTEGER iseed1,iseed2
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER g
C     ..
C     .. External Functions ..
      LOGICAL qrgnin
      EXTERNAL qrgnin
C     ..
C     .. External Subroutines ..
      EXTERNAL getcgn,initgn
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/
C     ..
C     .. Executable Statements ..
C     Abort unless random number generator initialized
      IF (qrgnin()) GO TO 10
      call rexit("SETSD called before RNG initialized")

   10 CALL getcgn(g)
      ig1(g) = iseed1
      ig2(g) = iseed2
      CALL initgn(-1)
      RETURN

      END
      REAL FUNCTION sexpo()
C**********************************************************************C
C                                                                      C
C                                                                      C
C     (STANDARD-)  E X P O N E N T I A L   DISTRIBUTION                C
C                                                                      C
C                                                                      C
C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     FOR DETAILS SEE:                                                 C
C                                                                      C
C               AHRENS, J.H. AND DIETER, U.                            C
C               COMPUTER METHODS FOR SAMPLING FROM THE                 C
C               EXPONENTIAL AND NORMAL DISTRIBUTIONS.                  C
C               COMM. ACM, 15,10 (OCT. 1972), 873 - 882.               C
C                                                                      C
C     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM       C
C     'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)       C
C                                                                      C
C     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   C
C     SUNIF.  The argument IR thus goes away.                          C
C                                                                      C
C**********************************************************************C
C
C
C     Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
C     (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
C
C     JJV added a Save statement for q (in Data statement)
C     .. Local Scalars ..
      REAL a,q1,u,umin,ustar
      INTEGER i
C     ..
C     .. Local Arrays ..
      REAL q(8)
C     ..
C     .. External Functions ..
      REAL ranf
      EXTERNAL ranf
C     ..
C     .. Equivalences ..
      EQUIVALENCE (q(1),q1)
C     ..
C     .. Save statement ..
      SAVE q
C     ..
C     .. Data statements ..
      DATA q/.6931472,.9333737,.9888778,.9984959,.9998293,.9999833,
     +     .9999986,.9999999/
C     ..
C
      a = 0.0
      u = ranf()
      GO TO 30

   20 a = a + q1
   30 u = u + u
C     JJV changed the following to reflect the true algorithm and
C     JJV prevent unpredictable behavior if U is initially 0.5.
C      IF (u.LE.1.0) GO TO 20
      IF (u.LT.1.0) GO TO 20
      u = u - 1.0
      IF (u.GT.q1) GO TO 60
      sexpo = a + u
      RETURN

   60 i = 1
      ustar = ranf()
      umin = ustar
   70 ustar = ranf()
      IF (ustar.LT.umin) umin = ustar
      i = i + 1
      IF (u.GT.q(i)) GO TO 70
      sexpo = a + umin*q1
      RETURN

      END
      REAL FUNCTION sgamma(a)
C**********************************************************************C
C                                                                      C
C                                                                      C
C     (STANDARD-)  G A M M A  DISTRIBUTION                             C
C                                                                      C
C                                                                      C
C**********************************************************************C
C**********************************************************************C
C                                                                      C
C               PARAMETER  A >= 1.0  !                                 C
C                                                                      C
C**********************************************************************C
C                                                                      C
C     FOR DETAILS SEE:                                                 C
C                                                                      C
C               AHRENS, J.H. AND DIETER, U.                            C
C               GENERATING GAMMA VARIATES BY A                         C
C               MODIFIED REJECTION TECHNIQUE.                          C
C               COMM. ACM, 25,1 (JAN. 1982), 47 - 54.                  C
C                                                                      C
C     STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER     C
C                                 (STRAIGHTFORWARD IMPLEMENTATION)     C
C                                                                      C
C     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   C
C     SUNIF.  The argument IR thus goes away.                          C
C                                                                      C
C**********************************************************************C
C                                                                      C
C               PARAMETER  0.0 < A < 1.0  !                            C
C                                                                      C
C**********************************************************************C
C                                                                      C
C     FOR DETAILS SEE:                                                 C
C                                                                      C
C               AHRENS, J.H. AND DIETER, U.                            C
C               COMPUTER METHODS FOR SAMPLING FROM GAMMA,              C
C               BETA, POISSON AND BINOMIAL DISTRIBUTIONS.              C
C               COMPUTING, 12 (1974), 223 - 246.                       C
C                                                                      C
C     (ADAPTED IMPLEMENTATION OF ALGORITHM 'GS' IN THE ABOVE PAPER)    C
C                                                                      C
C**********************************************************************C
C
C
C     INPUT: A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION
C     OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION
C
C     COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
C     COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
C     COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
C
C     .. Scalar Arguments ..
      REAL a
C     ..
C     .. Local Scalars .. (JJV added B0 to fix rare and subtle bug)
      REAL a1,a2,a3,a4,a5,a6,a7,aa,aaa,b,b0,c,d,e,e1,e2,e3,e4,e5,p,q,q0,
     +     q1,q2,q3,q4,q5,q6,q7,r,s,s2,si,sqrt32,t,u,v,w,x
C     ..
C     .. External Functions ..
      REAL ranf,sexpo,snorm
      EXTERNAL ranf,sexpo,snorm
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,alog,exp,sign,sqrt
C     ..
C     .. Save statement ..
C     JJV added Save statement for vars in Data satatements
      SAVE aa,aaa,s2,s,d,q0,b,si,c,q1,q2,q3,q4,q5,q6,q7,a1,a2,a3,a4,a5,
     +     a6,a7,e1,e2,e3,e4,e5,sqrt32
C     ..
C     .. Data statements ..
C
C     PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
C     SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
C
      DATA q1,q2,q3,q4,q5,q6,q7/.04166669,.02083148,.00801191,.00144121,
     +     -.00007388,.00024511,.00024240/
      DATA a1,a2,a3,a4,a5,a6,a7/.3333333,-.2500030,.2000062,-.1662921,
     +     .1423657,-.1367177,.1233795/
      DATA e1,e2,e3,e4,e5/1.,.4999897,.1668290,.0407753,.0102930/
      DATA aa/0.0/,aaa/0.0/,sqrt32/5.656854/
C     ..
C     .. Executable Statements ..
C
      IF (a.EQ.aa) GO TO 10
      IF (a.LT.1.0) GO TO 130
C
C     STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
C
      aa = a
      s2 = a - 0.5
      s = sqrt(s2)
      d = sqrt32 - 12.0*s
C
C     STEP  2:  T=STANDARD NORMAL DEVIATE,
C               X=(S,1/2)-NORMAL DEVIATE.
C               IMMEDIATE ACCEPTANCE (I)
C
   10 t = snorm()
      x = s + 0.5*t
      sgamma = x*x
      IF (t.GE.0.0) RETURN
C
C     STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
C
      u = ranf()
      IF (d*u.LE.t*t*t) RETURN
C
C     STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
C
      IF (a.EQ.aaa) GO TO 40
      aaa = a
      r = 1.0/a
      q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r
C
C               APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
C               THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
C               C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
C
      IF (a.LE.3.686) GO TO 30
      IF (a.LE.13.022) GO TO 20
C
C               CASE 3:  A .GT. 13.022
C
      b = 1.77
      si = .75
      c = .1515/s
      GO TO 40
C
C               CASE 2:  3.686 .LT. A .LE. 13.022
C
   20 b = 1.654 + .0076*s2
      si = 1.68/s + .275
      c = .062/s + .024
      GO TO 40
C
C               CASE 1:  A .LE. 3.686
C
   30 b = .463 + s + .178*s2
      si = 1.235
      c = .195/s - .079 + .16*s
C
C     STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
C
   40 IF (x.LE.0.0) GO TO 70
C
C     STEP  6:  CALCULATION OF V AND QUOTIENT Q
C
      v = t/ (s+s)
      IF (abs(v).LE.0.25) GO TO 50
      q = q0 - s*t + 0.25*t*t + (s2+s2)*alog(1.0+v)
      GO TO 60

   50 q = q0 + 0.5*t*t* ((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
C
C     STEP  7:  QUOTIENT ACCEPTANCE (Q)
C
   60 IF (alog(1.0-u).LE.q) RETURN
C
C     STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
C               U= 0,1 -UNIFORM DEVIATE
C               T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
C
   70 e = sexpo()
      u = ranf()
      u = u + u - 1.0
      t = b + sign(si*e,u)
C
C     STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
C
      IF (t.LT. (-.7187449)) GO TO 70
C
C     STEP 10:  CALCULATION OF V AND QUOTIENT Q
C
      v = t/ (s+s)
      IF (abs(v).LE.0.25) GO TO 90
      q = q0 - s*t + 0.25*t*t + (s2+s2)*alog(1.0+v)
      GO TO 100

   90 q = q0 + 0.5*t*t* ((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
C
C     STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
C
  100 IF (q.LE.0.0) GO TO 70
      IF (q.LE.0.5) GO TO 110
C
C     JJV modified the code through line 125 to handle large Q case
C
      IF (q.LT.15.0) GO TO 105
C
C     JJV Here Q is large enough that Q = log(exp(Q) - 1.0) (for real Q)
C     JJV so reformulate test at 120 in terms of one EXP, if not too big
C     JJV 87.49823 is close to the largest real which can be
C     JJV exponentiated (87.49823 = log(1.0E38))
C
      IF ((q+e-0.5*t*t).GT.87.49823) GO TO 125
      IF (c*abs(u).GT.exp(q+e-0.5*t*t)) GO TO 70
      GO TO 125

 105  w = exp(q) - 1.0
      GO TO 120

  110 w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q
C
C               IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
C
  120 IF (c*abs(u).GT.w*exp(e-0.5*t*t)) GO TO 70
 125  x = s + 0.5*t
      sgamma = x*x
      RETURN
C
C     ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))
C
C     JJV changed B to B0 (which was added to declarations for this)
C     JJV in 130 to END to fix rare and subtle bug.
C     JJV Line: '130 aa = 0.0' was removed (unnecessary, wasteful).
C     JJV Reasons: the state of AA only serves to tell the A .GE. 1.0
C     JJV case if certain A-dependant constants need to be recalculated.
C     JJV The A .LT. 1.0 case (here) no longer changes any of these, and
C     JJV the recalculation of B (which used to change with an
C     JJV A .LT. 1.0 call) is governed by the state of AAA anyway.
C
 130  b0 = 1.0 + .3678794*a
  140 p = b0*ranf()
      IF (p.GE.1.0) GO TO 150
      sgamma = exp(alog(p)/a)
      IF (sexpo().LT.sgamma) GO TO 140
      RETURN

  150 sgamma = -alog((b0-p)/a)
      IF (sexpo().LT. (1.0-a)*alog(sgamma)) GO TO 140
      RETURN

      END
      REAL FUNCTION snorm()
C**********************************************************************C
C                                                                      C
C                                                                      C
C     (STANDARD-)  N O R M A L  DISTRIBUTION                           C
C                                                                      C
C                                                                      C
C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     FOR DETAILS SEE:                                                 C
C                                                                      C
C               AHRENS, J.H. AND DIETER, U.                            C
C               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             C
C               SAMPLING FROM THE NORMAL DISTRIBUTION.                 C
C               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          C
C                                                                      C
C     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  C
C     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  C
C                                                                      C
C     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   C
C     SUNIF.  The argument IR thus goes away.                          C
C                                                                      C
C**********************************************************************C
C
C
C     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
C     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
C
C     .. Local Scalars ..
      REAL aa,s,tt,u,ustar,w,y
      INTEGER i
C     ..
C     .. Local Arrays ..
      REAL a(32),d(31),h(31),t(31)
C     ..
C     .. External Functions ..
      REAL ranf
      EXTERNAL ranf
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC float,int
C     ..
C     .. Save statement ..
C     JJV added a Save statement for arrays initialized in Data statmts
      SAVE a,d,t,h
C     ..
C     .. Data statements ..
      DATA a/0.0,.3917609E-1,.7841241E-1,.1177699,.1573107,.1970991,
     +     .2372021,.2776904,.3186394,.3601299,.4022501,.4450965,
     +     .4887764,.5334097,.5791322,.6260990,.6744898,.7245144,
     +     .7764218,.8305109,.8871466,.9467818,1.009990,1.077516,
     +     1.150349,1.229859,1.318011,1.417797,1.534121,1.675940,
     +     1.862732,2.153875/
      DATA d/5*0.0,.2636843,.2425085,.2255674,.2116342,.1999243,
     +     .1899108,.1812252,.1736014,.1668419,.1607967,.1553497,
     +     .1504094,.1459026,.1417700,.1379632,.1344418,.1311722,
     +     .1281260,.1252791,.1226109,.1201036,.1177417,.1155119,
     +     .1134023,.1114027,.1095039/
      DATA t/.7673828E-3,.2306870E-2,.3860618E-2,.5438454E-2,
     +     .7050699E-2,.8708396E-2,.1042357E-1,.1220953E-1,.1408125E-1,
     +     .1605579E-1,.1815290E-1,.2039573E-1,.2281177E-1,.2543407E-1,
     +     .2830296E-1,.3146822E-1,.3499233E-1,.3895483E-1,.4345878E-1,
     +     .4864035E-1,.5468334E-1,.6184222E-1,.7047983E-1,.8113195E-1,
     +     .9462444E-1,.1123001,.1364980,.1716886,.2276241,.3304980,
     +     .5847031/
      DATA h/.3920617E-1,.3932705E-1,.3950999E-1,.3975703E-1,
     +     .4007093E-1,.4045533E-1,.4091481E-1,.4145507E-1,.4208311E-1,
     +     .4280748E-1,.4363863E-1,.4458932E-1,.4567523E-1,.4691571E-1,
     +     .4833487E-1,.4996298E-1,.5183859E-1,.5401138E-1,.5654656E-1,
     +     .5953130E-1,.6308489E-1,.6737503E-1,.7264544E-1,.7926471E-1,
     +     .8781922E-1,.9930398E-1,.1155599,.1404344,.1836142,.2790016,
     +     .7010474/
C     ..
C     .. Executable Statements ..
C
      u = ranf()
      s = 0.0
      IF (u.GT.0.5) s = 1.0
      u = u + u - s
      u = 32.0*u
      i = int(u)
      IF (i.EQ.32) i = 31
      IF (i.EQ.0) GO TO 100
C
C                                START CENTER
C
      ustar = u - float(i)
      aa = a(i)
   40 IF (ustar.LE.t(i)) GO TO 60
      w = (ustar-t(i))*h(i)
C
C                                EXIT   (BOTH CASES)
C
   50 y = aa + w
      snorm = y
      IF (s.EQ.1.0) snorm = -y
      RETURN
C
C                                CENTER CONTINUED
C
   60 u = ranf()
      w = u* (a(i+1)-aa)
      tt = (0.5*w+aa)*w
      GO TO 80

   70 tt = u
      ustar = ranf()
   80 IF (ustar.GT.tt) GO TO 50
      u = ranf()
      IF (ustar.GE.u) GO TO 70
      ustar = ranf()
      GO TO 40
C
C                                START TAIL
C
  100 i = 6
      aa = a(32)
      GO TO 120

  110 aa = aa + d(i)
      i = i + 1
  120 u = u + u
      IF (u.LT.1.0) GO TO 110
      u = u - 1.0
  140 w = u*d(i)
      tt = (0.5*w+aa)*w
      GO TO 160

  150 tt = u
  160 ustar = ranf()
      IF (ustar.GT.tt) GO TO 50
      u = ranf()
      IF (ustar.GE.u) GO TO 150
      u = ranf()
      GO TO 140

      END
*DECK SPOFA
      SUBROUTINE spofa(a,lda,n,info)
      INTEGER lda,n,info
      REAL a(lda,1)
C
C     SPOFA FACTORS A REAL SYMMETRIC POSITIVE DEFINITE MATRIX.
C
C     SPOFA IS USUALLY CALLED BY SPOCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR SPOCO) = (1 + 18/N)*(TIME FOR SPOFA) .
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE SYMMETRIC MATRIX TO BE FACTORED.  ONLY THE
C                DIAGONAL AND UPPER TRIANGLE ARE USED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = TRANS(R)*R
C                WHERE  TRANS(R)  IS THE TRANSPOSE.
C                THE STRICT LOWER TRIANGLE IS UNALTERED.
C                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
C
C        INFO    INTEGER
C                = 0  FOR NORMAL RETURN.
C                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
C                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SDOT
C     FORTRAN SQRT
C
C     INTERNAL VARIABLES
C
      REAL sdot,t
      REAL s
      INTEGER j,jm1,k
C     BEGIN BLOCK WITH ...EXITS TO 40
C
C
      DO 30 j = 1,n
          info = j
          s = 0.0E0
          jm1 = j - 1
          IF (jm1.LT.1) GO TO 20
          DO 10 k = 1,jm1
              t = a(k,j) - sdot(k-1,a(1,k),1,a(1,j),1)
              t = t/a(k,k)
              a(k,j) = t
              s = s + t*t
   10     CONTINUE
   20     CONTINUE
          s = a(j,j) - s
C     ......EXIT
          IF (s.LE.0.0E0) GO TO 40
          a(j,j) = sqrt(s)
   30 CONTINUE
      info = 0
   40 CONTINUE
      RETURN

      END
