static const char css_c_sccs_id[] = "@(#)$Workfile: rand.c$ $Revision: 5$"; 
 
/*+ 
 * rand.c - Random Variate Generators 
 * 
 * $Revision: 5$ 
 * 
 * Description: 
 *		Holds a number of functions concerned with random variate 
 *		generation. 
 * 
 * Modification History: 
 * <unknown>	V1	Peter Mueller 
 *		Created. 
 * 
 * 25-JUN-1999	V2		Stephen Morris (Tessella Support services plc) 
 *		1) Removed nested comments (The C slash-star comment introducer within 
 *		a comment block).  The Gnu C compiler issued warnings on these. 
 *		2) Added "stdlib.h" for a definition of "exit()". 
 *		3) Include "mathconst.h" to define M_PI for Win32 platforms. 
 *		4) Include "unixrep.h" to define replacement functions for Unix. 
 *		5) Include "fortif.h" to handle C/Fortran interfacing issues.  (As part 
 *		of this, the call to cdfnor_ has been replaced with a call to cdfnor, 
 *		and the type of the "status" argument in this call corrected.) 
 *		6) Removed unreferenced local variables. 
 *		7) In "beta_pars", replaced (r < var*(1+r)^2) by (r < var*(1+r)*(1+r)). 
 *		"^" is the bit-wise exclusive-OR operator, which does not seem 
 *		appropriate in this case. 
 *		8) The functions ld_mixture_normal and rcond_mixture_normal_I do not 
 *		return a value, so declare them as "void". 
 * 
 * 15-JUL-1999	V4		Stephen Morris (Tessella Support Services plc) 
 *		1) Within Response\01, tests showed that a lot of time is being spent 
 *		in mvnS_rand() allocating and deallocating memory.  Other tests showed 
 *		that many calls within the program use the same size elements.  To 
 *		speed up the routine, the memory is cached, not being freed on routine 
 *		exit.  If a subsequent call is for the same size memory block, the 
 *		memory is reused.  An exit handler is added to ensure that allocated 
 *		memory is freed on exit. 
 *		2) Modified header to use "Versions" revision numbering. 
 * 
 * 12-JAN-2000	V5    P. J. L. Heesterman (Tessella Support Services plc) 
 *		    Constants that were assigned to float made explicit float to 
 *		    suppress contstan double to float convertion warning. 
 *		    Warnings on non-constant assignments of double to float are 
 *		    suppressed by disabling the warning. 
 * 
-*/ 
  
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
 
#include "ToolsRanlib.h" 
#include "ToolsCMatrix.h" 
#include "ToolsVector.h" 
#include "ToolsNrutil.h" 
#include "ToolsRand.h" 
#include "ToolsMess.h" 
#include "ToolsMathconst.h" 
#include "ToolsFortif.h" 
#include "ToolsUnixrep.h" 
#include "ToolsCCdflib.h"
 
 
static long is1=123456789, is2=981963; 
static int set=0; 
 
/* Memory cache for mvnS_rand, all prefixed by mvnsr_ */ 
 
static double *mvnsr_z = NULL; 
static double **mvnsr_S_cdt = NULL; 
static double **mvnsr_S_cd = NULL; 
static int mvnsr_p = -1;			/* Saved value of "p" */ 
 
#if defined(WIN32) 
#pragma warning(disable : 4244) 
#endif 
 
/******************************************** 
 *         Setseed 
 ********************************************/ 
/* call setall(is1,is2) */ 
void setseed(long is1, long is2) 
{ 
  set=1;   
  setall(is1,is2); 
} 
 
void getseed(long *is1, long *is2) 
{ 
  set=1;  
  getsd(is1,is2); 
} 
 
/* =====================================================================  
   Univariate RV generation & pdfs 
   ===================================================================== */  
/******************************************** 
 *         uniform 
 ********************************************/ 
/* returns uniform r.v. of size n in x			  
   calls IMSL routine 
   d_uniform: double 
   uniform:   float 
  ------------------------------------------*/ 
 
void duniform(int n, double *x) 
/* double */ 
{ 	 
  int i; 
 
  if (set==0){ 
    setall(is1,is2); 
    set=1; 
  } 
  for(i=0;i<n;i++) 
    x[i]=genunf(0.0,1.0); 
} 
 
void	uniform(int n, float *x) 
/* float */ 
{ int  
    i; 
 
  if (set==0){ 
    setall(is1,is2); 
    set=1; 
  } 
  for(i=0;i<n;i++) 
    x[i]=genunf(0.0,1.0); 
} 
 
double runif() 
{ 
  double  
    x; 
 
  if (set==0){ 
    setall(is1,is2); 
    set=1; 
  } 
   
  /* assign to double x for conversion */ 
  x = genunf(0.0,1.0); 
  return(x); 
} 
 
 
/******************************************* 
         stdnormal 
 ********************************************/ 
/* returns standard normal r.vector of size n in x  
   calls IMSL routine 
  ------------------------------------------*/ 
 
void	dstdnormal(int n,double *x) 
{	 
  int  
    i; 
 
  for(i=0;i<n;i++) 
    x[i] = gennor(0.0,1.0); 
} 
 
 
void	stdnormal(int n,float *x) 
/* float */ 
{	 
  int  
    i; 
 
  for(i=0;i<n;i++) 
    x[i] = gennor(0.0,1.0); 
} 
 
double norm_rand() 
{ 
 
  return gennor(0.0,1.0); 
} 
 
/* Poisson Distribution */ 
long l_poisson(double mu)
{
    return ignpoi(mu);
}

/******************************************** 
         normal_dens 
 ********************************************/ 
double d_normal_dens(double y, double m, double s) 
     /* returns log pdf of normal y ~ Normal(mu,s^2) */ 
{	 
  double log(); 
 
  return(-0.5*log(2.0*M_PI) - log(1.0*s) -  
	 (y-m)*(y-m)/(2.0*s*s)); 
} 
double pdfnorm(double y, double m, double s) 
{ 
  return(exp(d_normal_dens(y,m,s))); 
} 
 
/******************************************* 
         log_normal_dens 
 ********************************************/ 
double	d_lognormal_dens(y,m,s) 
/* returns pdf of lognormal y = exp(x), 
   x ~ Normal(mu,s^2) */ 
double y, m, s; 
{	 
  double x, f, log(), d_normal_dens(); 
 
  if (y > 0){ 
    x = log(y); 
    f = d_normal_dens(x,m,s) - x; 
  } 
  else f = -10e10; 
   
  return(f); 
} 
 
 
/* *************************************************' 
   normal cdf and inv cdf 
 ************************************************* */ 
double	cdfnorm(double y, double m, double s) 
/* returns cdf of normal N(m,s^2) at x */ 
{	 
  double  /* used to all be float... */
    p, q, mean,sd,bound,x,z; 
  double  
    cdf; 
  int  
    status,which; 
 
  /* primitive type conversion */ 
  x = y;  
  mean = m; 
  sd = s; 
  which = 1; 
  z = (x-mean)/sd; 
 
  if (z < -5.0) 
    p = 2.86e-7F; 
  else if (z > 5.0) 
    p = 0.9999997F; 
  else 
    cdfnor(&which,&p, &q,     &x,     &mean,  &sd,    &status,    &bound); 
   
  cdf = p; /* another primitive type conversion */ 
  return cdf; 
} 
 
double	invcdfn (double cdf, double m, double s) 
/* returns inv cdf of normal N(m,s^2) at p */ 
{	 
  double /* used to all be floats.. */
    p,q,mean,sd,bound,x; 
  double  
    y; 
  int  
    status,which; 
 
  if( (cdf < 0.0) | (cdf > 1.0) ){ 
    error("invcdfn",  
	  "Tried inverse cdf with p<0 or p>1", 
	  1); 
  }  
  /* par check */ 
  if (cdf <= 2.86e-07) 
    y = -5.0*s+m; 
  else if (cdf >= 0.9999997) 
    y = 5.0*s+m; 
  else{ 
    /* primitive type conversion */ 
    p = cdf; 
    q = 1.0-p;
    mean = m; 
    sd = s; 
    which = 2; 
     
    cdfnor(&which,&p,&q, &x,&mean,&sd,&status,&bound); 
   
    y = x; /* another primitive type conversion */ 
  } 
  return y; 
} 
 
 
/* ************************************************* 
   trunc normal 
 ************************************************* */ 
double norm_trunc_left(double m, double s,  
		       double cut) 
/* returns y ~ N(m,s^2), with cut < y */ 
{ 
  double u, pcut, y; 
 
  if (s==0.0){ /* degenerate s=0 */ 
    return (m>cut) ? m : cut; 
  } 
 
  pcut = cdfnorm(cut, m, s); 
  if (pcut>0.9999) 
    return cut+0.0001*s; 
  duniform(1, &u); 
  u = pcut + (1-pcut)*u; 
  y = invcdfn(u,m,s); 
 
  return y; 
} 
double norm_trunc_right(double m, double s,  
		       double cut) 
/* returns y ~ N(m,s^2), with y < cut */ 
{ 
  double  
    u, pcut, y; 
 
  if (s==0.0){ /* degenerate s=0 */ 
    return (m<cut) ? m : cut; 
  } 
  pcut = cdfnorm(cut, m, s); 
  if (pcut < 0.0001) 
    return cut-0.0001*s; 
  duniform(1, &u); 
  u = u*pcut; 
  y = invcdfn(u,m,s); 
 
  return y; 
} 
double norm_trunc(double m, double s,  
	       double left, double right) 
/* returns y ~ N(m,s^2), with left < y < right */ 
{ 
  double u, pleft, pright, y; 
 
  pleft = cdfnorm(left, m, s); 
  pright = cdfnorm(right, m, s); 
  if (pleft>0.9999) 
    return left+0.0001*max_xy(right-left,s); 
  if (pright<0.0001) 
    return right-0.0001*max_xy(right-left,s);  
  duniform(1, &u); 
  u = pleft + (pright-pleft)*u; 
  y = invcdfn(u,m,s); 
 
  return y; 
} 

double norm_trunc2(double m, double s,  
		   double left, double right,
		   double *lpdf)
/* returns y ~ N(m,s^2), with left < y < right */ 
{ 
  double u, pleft, pright, y; 
 
  pleft = cdfnorm(left, m, s); 
  pright = cdfnorm(right, m, s); 
  if (pleft>0.99) 
    return left; 
  if (pright<0.01) 
    return right; 
  duniform(1, &u); 
  u = pleft + (pright-pleft)*u; 
  y = invcdfn(u,m,s); 

  *lpdf = (-0.5*log(2.0*M_PI) - log(1.0*s) -  
	  (y-m)*(y-m)/(2.0*s*s))
    - log(pright-pleft);

  return y; 
} 

double pdf_norm_trunc(double y, double m, double s,  
		      double left, double right)
     /* evaluate turncated normal pdf */
{ 
  double u, pleft, pright, lpdf;
 
  pleft = cdfnorm(left, m, s); 
  pright = cdfnorm(right, m, s); 

  lpdf = (-0.5*log(2.0*M_PI) - log(1.0*s) -  
	  (y-m)*(y-m)/(2.0*s*s))
    - log(pright-pleft);

  return y; 
} 
 
 
double m_norm_trunc_left(double m, double s, double cut)  
/* returns E(y) ~ N(m,s^2), with cut < y  
   E(y) = m + w(2*pi*w)^{-1/2} exp{-(m-cut)^2/(2w) / Phi{(m-c)/s} 
 
   w = s^2 (from O'Hagan, p73) */ 
 
{ 
  double  
    w, pcut, E, 
    sqrt2pi = 2.5066282746310002; 
 
  if (s==0.0){ /* degenerate s=0 */ 
    E = (m>cut) ? m : cut; 
  } else { 
    pcut = cdfnorm((m-cut)/s, 0.0, 1.0); 
    w = s*s; 
    E = m + s/sqrt2pi*exp(-0.5*(m-cut)*(m-cut)/w)/pcut; 
  } 
  return(E);  
} 
 
/* *************************************************  
`t 
 ************************************************* */  
 
double t_logpdf(double x, double m, double sd, int nu) 
/* without normalizing constant */ 
{ 
  double y; 
 
  y = -0.5*(nu+1.0)*log(1.0+sqr((x-m)/sd)/nu); 
  return y; 
} 
 
double tc_logpdf(double sd, int nu) 
/* log normalizing constant */ 
{ 
  double y; 
 
  y = lgamma(0.5*(nu+1.0))-lgamma(0.5*nu)-0.5*log(nu*M_PI)-log(sd); 
  return y; 
} 
 
/* *********************************************************** 
   beta 
   *********************************************************** */ 
 
double betadev(double alpha, double beta) 
{ 
      double gamdev(); 
      double x, y; 
       
      x = gamdev(alpha);               /* X ~ gamma(alpha) */ 
      y = gamdev(beta);                /* Y ~ gamma(beta)  */ 
 
      return (x/(x+y));     /* X/(X+y) ~ ebta(alpha, beta) */ 
} 
 
void beta_pars(double mu, double var, double *a, double *b) 
{ 
  /* returns beta par's a,b as a function of mu & var */ 
  double r; 
  double r1; 
 
  r = mu/(1-mu); 
  r1 = r + 1.0; 
  if (r < var*r1*r1){ 
    printf("\n impossible beta pars! \n"); 
    exit(1); 
  } 
   * b = (r - var * r1 * r1) / (var * r1 * r1 * r1); 
  *a = *b * r; 
} 
 
 
double beta_logpdf(double x, double a, double b) 
{ 
  	double y; 

    	if(x <=0 || x >=1)
        	return -30.0;     //log of zero.

  	y = (a-1)*log(x)+(b-1)*log(1.0-x); 
  
  	//normalization constant: added by Raj -- 08/22/2003	
  	y += (lgamma(a+b) - lgamma(a) -lgamma(b));
  
	return y; 
} 


/* *************************************************' 
   beta cdf and inv cdf 
 ************************************************* */ 
double	cdfbeta(double y, double a, double b) 
/* returns cdf of Beta(a,b) at x */ 
{	 
  double /* float   */
    p, q, af,bf,x, x1, bound;
  double  
    cdf; 
  int  
    status,which; 
 
  if( (y < 0.0) | (y > 1.0) ){ 
    error("cdfbeta",  
	  "Tried cdf beta with y<0 or y>1 (it's okay, but weird)", 
	  1); 
  }  
  if( (a < 1E-30) | (a > 1E30) ){ 
    error("cdfbeta",  
	  "Tried cdf beta with a outside [1E-30,1E30]",
	  1); 
  }  
  if( (b < 1E-30) | (b > 1E30) ){ 
    error("cdfbeta",  
	  "Tried cdf beta with b outside [1E-30,1E30]",
	  1); 
  }  

  /* primitive type conversion */ 
  x = y;  
  x1 = 1.0-y;
  af = a; bf = b;
  which = 1; 
 
  cdfbet(&which, &p, &q, &x,&y, &af,&bf,&status,&bound); 
  if(status < 0){ 
    error("cdfbeta",  
	  "Input parameter out of range",
	  1); 
  }  
  if(status == 1)
    p = 0.0;
  else if (status == 2)
    p = 1.0;

  cdf = p; /* another primitive type conversion */ 
  return cdf; 
} 
 
double	invcdfbeta (double cdf, double a, double b) 
/* returns inv cdf of beta Be(a,b) at p */ 
{	 
  double /* float   */
    p,q,af,bf,bound,x,x1; 
  double  
    y; 
  int  
    status,which; 
 
  /* par check */ 
  if( (cdf < 0.0) | (cdf > 1.0) ){ 
    error("invcdfbeta",  
	  "Tried invcdf beta with y<0 or y>1",
	  1); 
  }  
  if( (a < 1E-30) | (a > 1E30) ){ 
    error("invcdfbeta",  
	  "Tried invcdf beta with a outside [1E-30,1E30]",
	  1); 
  }  
  if( (b < 1E-30) | (b > 1E30) ){ 
    error("invcdfbeta",  
	  "Tried invcdf beta with b outside [1E-30,1E30]",
	  1); 
  }  
  /* primitive type conversion */ 
  p = cdf; 
  q = 1.0-p;
  af = a; bf = b;
  which = 2; 
     
  cdfbet(&which,&p,&q,&x,&x1,&af,&bf,&status,&bound); 
  if(status < 0){ 
    error("cdfbeta",  
	  "Input parameter out of range",
	  1); 
  }  
  if(status == 1)
    x = 0.0;
  else if (status == 2)
    x = 1.0;
   
  y = x; /* another primitive type conversion */ 
  return y; 
} 
 
void dirichlet(double *alpha, double *w, int p) 
/* ****!! not tested *****/ 
{ 
  double s, a,b, W; 
  int j; 
 
  for(s=0,j=0;j<p;j++) s += alpha[j]; 
  for(j=0,W=1.0,b=s;j<p-1;j++){ 
    a = alpha[j]; 
    b -= alpha[j]; 
    w[j] = betadev(a,b)*W; 
    W -= w[j]; 
  } 
  w[p-1] = W; 
  if (W <= 0) 
    printf("\n\n **** non-pos dirich gen!! **\n"); 
} 
   
 
/* *********************************************************** 
   gamma 
   *********************************************************** */ 
double gamdev(double alpha) 
{ 
  double  
    value; 
  float 
    a; 
  float  
    gengam(float,float); 
   
  a = alpha; /* type conversion */ 
  value = gengam(1.0,a); 
  return value; 
} 

double cdfgam_ab(double y, double a, double b)
/* untested */
{ // returns cdf of a G(a,b), E(X)=a/b, V(x)=a/b^2
  double  /* used to all be float ???? */
    p, q, shape, scale,bound,x,z; 
  double  
    cdf; 
  int  
    status,which; 
 
  /* primitive type conversion */ 
  x = y;  
  scale=b;
  shape=a;
  which = 1; 
  z = y*b;       // gamma(a) 

  if (z > 10.0*sqrt(a))  // to avoid infty, V(z)=a
    p = 0.9999997; 
  else 
    cdfgam(&which,&p, &q, &x, &shape, &scale, &status, &bound); 
   
  cdf = p; /* another primitive type conversion */ 
  return cdf; 
} 
 
double	invcdfgam (double cdf, double a, double b) 
/* untested */
/* returns inv cdf of gamma G(a,b), with E(x)=a/b */
{	 
  double /* used to all be floats.. */
    p,q,shape,scale,bound,x; 
  double  
    y; 
  int  
    status,which; 
 
  if( (cdf < 0.0) | (cdf > 1.0) ){ 
    error("invcdfn",  
	  "Tried inverse cdf with p<0 or p>1", 
	  1); 
  }  
  /* par check */ 
  if (cdf <= 0)
    y = 0;
  else if (cdf >= 0.9999997) 
     y = 5.0*sqrt(a)/b;  // 5 sd's out...
  else{ 
    /* primitive type conversion */ 
    p = cdf; 
    q = 1.0-p;
    scale=b;
    shape = a; 
    which = 2; 
     
    cdfgam(&which,&p,&q, &x,&shape,&scale,&status,&bound); 
   
    y = x; /* another primitive type conversion */ 
  } 
  return y; 
} 
 
 
/* ************************************************* 
   trunc normal 
 ************************************************* */ 
double gam_trunc_left(double a, double b,  
		       double cut) 
/* untested */
/* returns y ~ Ga(a,b), with cut < y */ 
{ 
  double u, pcut, y; 
 
  pcut = cdfgam_ab(cut, a, b);
  if (pcut>0.99) 
    return cut; 
  duniform(1, &u); 
  u = pcut + (1-pcut)*u; 
  y = invcdfgam(u,a,b);
 
  return y; 
} 

double gam_trunc_right(double a, double b, double cut) 
/* untested */
/* returns y ~ Ga(a,b), with cut > y */ 
{ 
  double u, pcut, y; 
 
  pcut = cdfgam_ab(cut, a, b);
  if (pcut < 1.0E-10)
    return 0.5*cut; 
  duniform(1, &u); 
  u = pcut*u; 
  y = invcdfgam(u,a,b);
 
  return y; 
} 

double gam_trunc(double a, double b, double left, double right) 
/* untested */
/* returns y ~ Ga(a,b), with cut > y */ 
{ 
  double u, pleft, pright, y; 
 
  pleft = cdfgam_ab(left, a, b);
  pright = cdfgam_ab(right, a, b);
  duniform(1, &u); 
  u = pleft + u*(pright-pleft);
  y = invcdfgam(u,a,b);
 
  return y; 
} 

/* ************************************************* */

double expdev(double beta)
{
  double value;	
  float 
    b, genexp(float);
  
  b=beta; // type conversion 
  value = genexp(beta);
  return(value);
}

double expdev_truncleft(double beta, double x0)
{
  double x;
  x = genexp(beta)+x0;
  return(x);
}
  
  
  
 
/* *********************************************************** 
   gamma_ab 
   *********************************************************** */ 
/* returns a Gamma(alpha,beta) r.v. */ 
/* f(x; a,b) = b*(bx)^(a-1)*exp(-bx) / Gamma(a) 
   E(x) = a/b, V(x) = a/b^2 */ 
double gamdev_ab(double alpha, double beta) 
{ 
  double value; 
   
  value = gamdev(alpha)/beta; 
  return value; 
} 
 
 
double gamma_logpdf(double x,double a, double b) 
{ 
  return a*log(b) + log(x)*(a-1.0) - b*x - lgamma(a); 
} 
 
/******************************************** 
 *         Weibull 
 ********************************************/ 
/* generates  
   X = Weibull(a,b) 
       ba(bx)^(a-1) exp[-(bx)^a] 
  ------------------------------------------   
   Args:  
   x:      random variate (OUTPUT) 
   a,b:    parameters of Weibul 
  ------------------------------------------   
   Reference: Bratley, Fox and Schrage, "Simulation" 
 ********************************************/ 
 
double d_weibull(double *x, double a, double b) 
{ 
  double u; 
  
  duniform(1,&u); 
  *x = pow(-log(1-u), 1.0/a) / b; 
  return(*x); 
} 
 
/********************************************  
 *         binomial 
 ********************************************/ 
/* returns binomial r.v. Bin(n,p)  
  ------------------------------------------*/ 
 
int	binomial(int n, double p) 
{ 	 
  int i,x; 
  for(i=0,x=0;i<n;i++) 
    x += (runif() < p) ? 1 : 0; 
  return x; 
} 
 
int binomial2(int n, double p) 
{ 	 
  int  
    x; 
  float 
    pf; 
  long 
    nl, xl; 
 
  nl = n; 
  pf = p; 
  xl = ignbin(nl,pf); 
  x = xl; 
  return(x); 
} 
 
double logp_binomial(int y, int n, double p) 
{ 
  double lp; 
  lp = lgamma(n+1.0)-lgamma(n-y+1.0)-lgamma(y+1.0) + 
    y*log(p) + (n-y)*log(1.0-p); 
  return(lp); 
} 
 
int rand_bebinomial(int n, double a, double b) 
{ 
  int i,x; 
  double th; 
 
  th = betadev(a,b); 
  for(i=0,x=0;i<n;i++) 
    x += (runif() < th) ? 1 : 0; 
  return x; 
} 
   
double logpdf_bebinomial(int y, int n, double a, double b, int c) 
/* evaluates beta-binom log pdf. If const=0 without const's in a,b */ 
{ 
  double lp1,lp2,lp3; 
 
  if (y>n) 
    return(-10.0e10); /* outside support */ 
  lp1 = lgamma(n+1.0)-lgamma(y+1.0)-lgamma(n-y+1.0); 
  lp2 = lgamma(a+y)+lgamma(n+b-y)-lgamma(a+b+n); 
  if (c==1) 
    lp3 = lgamma(a+b)-lgamma(a)-lgamma(b); 
  else  
    lp3 = 0; 
  return(lp1+lp2+lp3); 
} 
   
 
double lg(double x) /* for debugging */ 
{ 
  printf("%5.3f\n",lgamma(x)); 
  return(lgamma(x)); 
} 
 
/******************************************** 
 *         multinomial 
 ********************************************/ 
/* returns draw from multinomial with cell prob pr  
  ------------------------------------------   
  Value: 
    x:   vector of indices indicating draws 
         x in [0,..n_cell-1] 
  ------------------------------------------   
  Args: 
    n_draw: number of draws 
    n_cell: number of cells 
    pr:     n_cell vector of cell probs 
            (not necessarily standardized) 
    x:      n_draw vector of indices 
            (OUTPUT) 
 ********************************************/ 
void multinomial(int n_draw, int n_cell, double *pr, int *x) 
{ 
  double *cum_p, uj; 
  int i,j; 
 
  cum_p = dvector(0,n_cell); 
 
  for(i=1,cum_p[0]=pr[0];i<n_cell;i++)     
    cum_p[i] = cum_p[i-1]+pr[i]; 
  for(j=0;j<n_draw;j++){ 
    uj = runif()*cum_p[n_cell-1]; 
    for(i=0; ((uj > cum_p[i]) & (i<n_cell)); i++); 
    x[j] = i; 
  } 
 
  free_dvector(cum_p,0,n_cell); 
} 
 
void multinomial_I(int n_draw, int n_cell, int *pr, int *x) 
     /* same as multinomial(), but with int prob vector */ 
{ 
  double *cum_p, uj; 
  int i,j; 
 
  cum_p = dvector(0,n_cell); 
 
  for(i=1,cum_p[0]=pr[0];i<n_cell;i++)     
    cum_p[i] = cum_p[i-1]+pr[i]; 
  for(j=0;j<n_draw;j++){ 
    uj = runif()*cum_p[n_cell-1]; 
    for(i=0; ((uj > cum_p[i]) & (i<n_cell)); i++); 
    x[j] = i; 
  } 
 
  free_dvector(cum_p,0,n_cell); 
} 
 
void multinomial_table(int n_draw, int n_cell, double *pr, int *ct) 
{ /* same as multinomial, but returns counts instead of indices */ 
  double *cum_p, uj; 
  int i,j; 
 
  cum_p = dvector(0,n_cell); 
 
  ct[0]=0; 
  for(i=1,cum_p[0]=pr[0];i<n_cell;i++){ 
    cum_p[i] = cum_p[i-1]+pr[i]; 
    ct[i] = 0; 
  } 
  for(j=0;j<n_draw;j++){ 
    uj = runif()*cum_p[n_cell-1]; 
    for(i=0; ((uj > cum_p[i]) & (i<n_cell)); i++); 
    ct[i] += 1; 
  } 
 
  free_dvector(cum_p,0,n_cell); 
} 
 
/* =====================================================================  
   multivariate rv generation & pdfs 
   ===================================================================== */  
 
/* *********************************************************** 
   normal 
   *********************************************************** */ 
void mvn_rand(int p, double *mean, double **L, double *vector)    
       /* LL' = Sigma cov matrix */ 
{ 
 
     double *z; 
 
     z = dvector(0,p-1); 
 
     dstdnormal(p,z);  
     Ax_plus_y(L,z,mean,vector,p); 
     free_dvector(z,0,p-1); 
     return; 
}  
 
 
 
/*+ 
 * mvnS_rand_free - Free Memory Allocated By mvnS_rand 
 * 
 * Description: 
 *		Frees up any memory allocated by mvnS_rand.  mvnS_rand registers this 
 *		as an exit handler the first time it is called. 
 * 
 * Declaration: 
 *		void mvnS_rand_free() 
 * 
 * Arguments: 
 *		None. 
 * 
 * Side Effects: 
 *		Sets the pointer variables to NULL, and invalidates the saved value of 
 *		"p" (sets it to -1). 
 *		 
-*/ 
 
static void mvnS_rand_free(void) 
{ 
	/* Free memory */ 
 
	if (mvnsr_z != NULL) { 
		free_dvector(mvnsr_z, 0, (mvnsr_p - 1)); 
		mvnsr_z = NULL; 
	} 
 
	if (mvnsr_S_cd != NULL) { 
		free_dmatrix(mvnsr_S_cd, 0, (mvnsr_p - 1), 0, (mvnsr_p - 1)); 
		mvnsr_S_cd = NULL; 
	} 
 
	if (mvnsr_S_cdt != NULL) { 
		free_dmatrix(mvnsr_S_cdt, 0, (mvnsr_p - 1), 0, (mvnsr_p - 1)); 
		mvnsr_S_cdt = NULL; 
	} 
 
	/* Invalidate saved size of memory allocated */ 
 
	mvnsr_p = -1; 
} 
 
 
/* like normal, just with cov-matrix instead of L */   
double mvnS_rand(int p, double *mean, double **S, double *vector) 
{ 
	static int first_time = 1;		/* First-time flag */ 
	double
	  zsq, log_coef, log_value,
	  det_half;

	if (p != mvnsr_p) { 
 
		/* 
		 * Memory required different to that cached, so free memory and 
		 * allocate more. 
		 */ 
 
		mvnS_rand_free(); 
 
		mvnsr_S_cdt = dmatrix(0, (p - 1), 0, (p - 1)); 
		mvnsr_S_cd  = dmatrix(0, (p - 1), 0, (p - 1)); 
		mvnsr_z = dvector(0, (p - 1)); 
 
		mvnsr_p = p;		/* Save indicator of amount stored */ 
 
		/* 
		 * If this is the first time through the routine, register the exit 
		 * handler to deallocate the memory on program exit. 
		 */ 
 
		if (first_time) { 
			atexit(mvnS_rand_free); 
			first_time = 0; 
		} 
	} 
 
 
	d_chol_decomp(p, S, mvnsr_S_cdt); 
	At(mvnsr_S_cdt, mvnsr_S_cd, p); 
	dstdnormal(p, mvnsr_z); 
	Ax_plus_y(mvnsr_S_cd, mvnsr_z, mean, vector, p); 
 
	/* Don't free memory.  This is freed at program termination */ 
 
	xtx(mvnsr_z,&zsq,p); 
	det_half = fabs(prod_diag(mvnsr_S_cdt,p)); 
	if (det_half==0) /* no determinant for improper mvn */ 
	  log_coef = -0.5*p*log(2.0*M_PI); 
	else  
	  log_coef = -log(det_half)-0.5*p*log(2.0*M_PI); 
	log_value = log_coef-0.5*zsq; 
	return (log_value); 
}    
 
double mvn_cdf(int p, double *mu, double **V, double **L, 
	       double *left, double *right, int I) 
/* compute p(left<x<right) by importance sampling using 
   h(x) = 1/2*N(mu,VV')+1/2*N(0,VV')  
   taylor made for left or right approx= 0 */ 
{ 
  int i,j; 
  double W,w,*eps,*x,*z,u,f,s,*Lmu,h1,h2; 
 
  eps = dvector(0,p); 
  Lmu = dvector(0,p); 
  x = dvector(0,p); 
  z = dvector(0,p); 
 
  Ax(L,mu,Lmu,p); 
  for(W=0,s=0,i=0;i<I;i++){ 
    dstdnormal(p,eps); 
    u = runif(); 
    if (u<0.5){ 
      Ax(V,eps,x,p); 
      xtx(eps,&h1,p); 
      x_min_y(eps,Lmu,z,p); 
      xtx(z,&h2,p); 
    } 
    else{ 
      Ax_plus_y(V,eps,mu,x,p); 
      xtx(eps,&h2,p); 
      x_plus_y(eps,Lmu,z,p); 
      xtx(z,&h1,p); 
    } 
    w = exp(-0.5*h2)/(exp(-0.5*h1)+exp(-0.5*h2)); 
    W += w; 
    for(j=0,f=1;j<p;j++) 
      if ((x[j]<left[j])|(x[j]>right[j])) 
	f=0; 
    s += w*f; 
  } 
 
  free_dvector(eps,0,p); 
  free_dvector(Lmu,0,p); 
  free_dvector(x,0,p); 
  free_dvector(z,0,p); 
 
  return(s/W); 
} 
 
void mvntrunc_rand(int p, double *mu, double **L,  
		 double *left, double *right,  
		 double *x, double *z, 
		 int niter)    
/* LL' = Sigma cov matrix  
  carefull -  univ marginal of a mv trunc normal is not nec 
  trunc normal!!  
  niter steps of a gibbs sampler in the truncated mv normal 
  in standardized parametrization  
  (1) compute left & right cut for z[0] given z[1..p-1] 
  (2) generate z[0] 
  (3) repeat for z[1..p-1] 
*/ 
{ 
  double lo,hi,loj,hij; 
  int 
    i,j,iter; 
   
  if ((L[0][1] != 0.0)|(prod_diag(L,p)==0)){ 
    printf("non lower triangular L or 0 on diag in mvntrunc!\n"); 
    exit(1); 
  } 
 
  for(iter=0;iter<niter;iter++){ 
    for(i=0;i<p;i++){ 
      inequ(left[i],right[i],mu[i],L[i],z,i,&lo,&hi,i); 
      for(j=i+1;j<p;j++){ /*compute bounds */ 
	if (L[j][i]==0) continue; 
	inequ(left[j],right[j],mu[j],L[j],z,i,&loj,&hij,j+1); 
	lo = max_xy(lo,loj); 
	hi = min_xy(hi,hij); 
      } 
      z[i] = norm_trunc(0.0,1.0,lo,hi); 
    } 
  } 
  Ax_plus_y(L,z,mu,x,p); 
} 
 
 
void mvntrunc_rand_init(int p, double *mean, double **L,  
  double *left, double *right, double *x, double *z) 
/* aux for mvn_trunc */ 
{ 
  double **L_inv; 
  int i; 
 
  L_inv = dmatrix(0,p,0,p); 
 
  d_inv_triang_lo(p,L,L_inv); 
  for(i=0;i<p;i++){ 
    if ((mean[i] <= right[i]) & (left[i]<=mean[i])) 
      x[i] = mean[i]; 
    else if (right[i] < mean[i]) 
      x[i] = right[i]; 
    else if (left[i] > mean[i]) 
      x[i] = left[i]; 
  } 
  A_xminusy(L_inv,x,mean,z,p); 
  free_dmatrix(L_inv, 0,p,0,p); 
} 
 
void inequ(double left, double right, double m, double *L,  
	   double *x, int i, 
	   double *lo, double *hi, int p) 
/* auxilary to mvn_trunc; computes lower & upper bound 
   on x[i] in an inequality: 
   left < m+L'x < right */ 
{ 
  double Lx; 
  int ii; 
  for(ii=0,Lx=0;ii<p;ii++){ 
    if (ii==i) continue; 
    Lx += L[ii]*x[ii]; 
  } 
  left = (left-m-Lx)/L[i]; 
  right = (right-m-Lx)/L[i]; 
 
  *lo = (L[i]>0) ? left : right; 
  *hi = (L[i]>0) ? right : left; 
} 
 
/* *********************************************************** 
   MV-Normal pdf 
   *********************************************************** */ 
/* evaluates MV-Normal pdf with mean=m, CovMat = S, dim = p at 
   point x. 
   L_inv is the inverse of the chol-decomp of S.  
   LL' = S and L_inv = 1/L.  
   don't include coeff for improper mvn dist */ 
 
double mvn_logpdf(double *x, double *m, double **L_inv, int p) 
{ 
  double *z, zsq, det_minushalf,  
    log_coef, log_value; 
 
  /* Allocate memory for aux vector */ 
  z = dvector(0,p-1); 
 
  /* compute z=L-inv*(x-m) */ 
  A_xminusy(L_inv,x,m,z,p); 
 
  /* log-pdf=log(det(L-inv))-p/2(log(2*pi))-z'z/2 */ 
  xtx(z,&zsq,p); 
  det_minushalf = fabs(prod_diag(L_inv,p)); 
  if (det_minushalf==0) /* no determinant for improper mvn */ 
    log_coef = -0.5*p*log(2.0*M_PI); 
  else  
    log_coef = log(det_minushalf)-0.5*p*log(2.0*M_PI); 
  log_value = log_coef-0.5*zsq; 
 
  /* free z */ 
  free_dvector(z,0,p-1); 
 
  return(log_value); 
} 
 
double mvn_pdf(double *x, double *m, double **L_inv, int p) 
{ 
  return exp(mvn_logpdf(x,m,L_inv,p)); 
} 
 
double mvnS_logpdf(double *x, double *m, double **S, int p) 
{ 
  double 
    **L, **V,lpdf; 
 
  L = dmatrix(0,p,0,p); 
  V = dmatrix(0,p,0,p); 
 
  /* S = V'V */ 
  d_chol_decomp(p,S,V); 
  /* L = 1/V */ 
  d_inv_triang(p,V,L); 
  /* LSL' = I */ 
  At_self(L,p); 
  lpdf = mvn_logpdf(x,m,L,p); 
 
  free_dmatrix(L,0,p,0,p); 
  free_dmatrix(V,0,p,0,p); 
  return lpdf; 
} 
 
/* ************************************************* 
   mvn cond moments 
 ************************************************* */ 
void cond_nmoments(int k, int p, double *x, 
		   double *mu, double **S, double **S22_inv, 
		   double *s, double *m) 
/* returns cond moments of x[k] | x[-k]  
   in s and m (variance and mean) */ 
{ 
  int  
    i,j,i2,j2; 
   
  /* A_xminusy(S22_inv,x[-k],mu[-k] */ 
  *m = mu[k]; 
  *s = S[k][k]; 
  for(i=0,i2=0;i<p;i++){ 
    if (i==k) continue; 
    for(j=j2=0; j<p; j++){ 
      if (j==k) continue; 
      *m += S[i][k]*S22_inv[i2][j2]*(x[j]-mu[j]); 
      *s -= S[i][k]*S22_inv[i2][j2]*S[k][j]; 
      j2 += 1; 
    } 
    i2 += 1; 
  } 
} 
 
void cond_n2moments(int k, int p, double y, double *mu, double **V, 
		    double *s, double *m) 
{ 
  int j; 
 
  j = 1-k; /* i.e. 0 if k=1, 1 if k=0 */ 
  *m = mu[k] + V[0][1]/V[j][j]*(y-mu[j]); 
  *s = V[k][k] - V[0][1]/V[j][j]*V[1][0]; 
} 
   
void cond_mvnmoments(int p, int q, double *x, double *mu, double **S, 
		     double **S22_inv, int *qi, int *qi_c, double **V, 
		     double *m) 
/* returns cond moments of  
   x[qi] | x[-qi] of a N(mu,S) */ 
/* q=dim of x[qi], p=dim of mu 
   qi = indices of x[qi] 
   qi_c = indices of x[-qi] 
   S22_inv = S[qi_c,qi_c]^-1 
   m, V: moments of conditional m.v.n. dist */ 
{ 
  int  
    i,i2,j,j2,k,k2,l,l2; 
 
  /* m = mu1 + S12*S22_inv*(x2-mu2) */ 
  for(k2=0;k2<q;k2++){ 
    k = qi[k2]; 
    m[k2] = mu[k]; 
    for(i2=0;i2<p-q;i2++){ 
      i = qi_c[i2]; 
      for(j2=0; j2<p-q; j2++){ 
	j = qi_c[j2]; 
	m[k2] += S[k][i]*S22_inv[i2][j2]*(x[j]-mu[j]); 
      } 
    } 
  } 
  /* S = S11 - S12*S22_inv*S12 */ 
  for(l2=0;l2<q;l2++){ 
    l = qi[l2]; 
    for(k2=0;k2<q;k2++){ 
      k = qi[k2]; 
      V[k2][l2] = S[k][l]; 
      for(i2=0;i2<p-q;i2++){ 
	i = qi_c[i2]; 
	for(j2=0; j2<p-q; j2++){ 
	  j = qi_c[j2]; 
	  V[k2][l2] -= S[k][i]*S22_inv[i2][j2]*S[j][l]; 
	} 
      } 
    } 
  } 
} 
 
void cond_mvnmean(int p, int q, double *x, double *mu, double **S, 
                  double **S22_inv, int *qi, int *qi_c, double *m) 
/* same as above, but mean only */ 
{ 
  int  
    i,i2,j,j2,k,k2; 
 
  /* m = mu1 + S12*S22_inv*(x2-mu2) */ 
  for(k2=0;k2<q;k2++){ 
    k = qi[k2]; 
    m[k2] = mu[k]; 
    for(i2=0;i2<p-q;i2++){ 
      i = qi_c[i2]; 
      for(j2=0; j2<p-q; j2++){ 
	j = qi_c[j2]; 
	m[k2] += S[k][i]*S22_inv[i2][j2]*(x[j]-mu[j]); 
      } 
    } 
  } 
} 
 
/*============ untested ====================*/ 
/********************************************  
 *         normal-logdens 
 ********************************************/ 
/* log MV normal density                    
   For x[1,..p] ~ MVN(mu, HH') 
  ------------------------------------------   
   Args: 
   p:     dimensionality 
   n:     number of rows in x 
   x:     n x p matrix of n r.v. vectors 
   mu:    mean (p-vector) 
   H:     chol-decomp of cov matrix (HH'=S) 
   f:     vector of log densities (OUTPUT) 
  ------------------------------------------*/ 
void normal_logdens(int p,int n, double *x, double *mu, double **H, 
		    double *f) 
{ 
  int i,j; 
  double *xi, *z,log_f,two_pi, det_H, *y; 
 
  z = dvector(0,p); 
  y = dvector(0,p); 
 
  two_pi = -0.5*p*log(2*M_PI); 
 
  for(j=0,det_H=0;j<p;j++) 
    det_H += log(H[j][j]); 
  for(i=0;i<n;i++){ 
    xi = &x[p*i]; 
    x_min_y(xi,mu,y,p); 
    Ax(H,y,z,p); 
    xtx(z,&log_f,p); 
    log_f *= -0.5; 
    f[i] = det_H+log_f+two_pi; 
  } 
   
  free_dvector(z,0,p); 
  free_dvector(y,0,p); 
} 
 
 
/* *********************************************************** 
   MV-t pdf 
   *********************************************************** */ 
/* evaluates MV-t pdf with location = m, scale = S, df = nu, 
   dim = p, at point x.  
   L_inv is the inverse of the chol-decomposition of S.  
   LL' = S, and L_inv = 1/L */ 
 
double mvt_pdf(double *x, double *m, double **L_inv, int nu, 
	       int p) 
{ 
  double *z, zsq, det_minushalf,  
      value, log_coef, log_value; 
 
 
  /* allocate memory for auxilary array */ 
  z = dvector(0,p-1); 
 
  /* compute z = L-inv*(x-m) */ 
  A_xminusy(L_inv,x,m,z,p); 
 
  /* log-pdf= det(L-inv) + lgamma([nu+p]/2) - 
     p/2*log(pi*nu) - lgamma(p/2) + 
     - (nu+p)/2*log(1+z'z/nu) */ 
  xtx(z,&zsq,p); 
  det_minushalf = fabs( prod_diag(L_inv,p) ); 
  log_coef = lgamma( (nu+p)*0.5 ) + log(det_minushalf) + 
    log(nu*M_PI)*(-p*0.5) - lgamma(nu*0.5); 
  log_value = log_coef + log(1.0+zsq/nu)*(-0.5*(nu+p)); 
  value = exp(log_value); 
 
  /* free z */ 
  free_dvector(z,0,p-1); 
 
  return(value); 
} 
 
 
double mvt_logpdf(double *x, double *m, double **L_inv, int nu, 
	       int p) 
{ 
  double *z, zsq, det_minushalf,  
      log_coef, log_value; 
 
 
  /* allocate memory for auxilary array */ 
  z = dvector(0,p-1); 
 
  /* compute z = L-inv*(x-m) */ 
  A_xminusy(L_inv,x,m,z,p); 
 
  /* log-pdf= det(L-inv) + lgamma([nu+p]/2) - 
     p/2*log(pi*nu) - lgamma(p/2) + 
     - (nu+p)/2*log(1+z'z/nu) */ 
  xtx(z,&zsq,p); 
  det_minushalf = fabs( prod_diag(L_inv,p) ); 
  log_coef = lgamma( (nu+p)*0.5 ) + log(det_minushalf) + 
    log(nu*M_PI)*(-p*0.5) - lgamma(nu*0.5); 
  log_value = log_coef + log(1.0+zsq/nu)*(-0.5*(nu+p)); 
 
  /* free z */ 
  free_dvector(z,0,p-1); 
 
  return(log_value); 
} 
 
double mvtS_logpdf(double *x, double *m, double **S, int nu, int p) 
{ 
  double 
    **L, **V,lpdf; 
 
  L = dmatrix(0,p,0,p); 
  V = dmatrix(0,p,0,p); 
 
  /* S = V'V */ 
  d_chol_decomp(p,S,V); 
  /* L = 1/V */ 
  d_inv_triang(p,V,L); 
  /* LSL' = I */ 
  At_self(L,p); 
  lpdf = log(mvt_logpdf(x,m,L,nu,p)); 
 
  free_dmatrix(L,0,p,0,p); 
  free_dmatrix(V,0,p,0,p); 
  return lpdf; 
} 
 
 
/* *********************************************************** 
   MV-Normal cond-pdf 
   *********************************************************** 
   evaluates univariate conditional in a MV normal density. 
   p(y | x), where x = z[xi], y = z[yi], and 
   z is N(m,S) */ 
 
/* ************* untested ******************************* */ 
double mvn_cond_pdf(double *z, double *m, double **S, int p,  
		    int px, int* xi, int yi) 
{ 
  double **Sxx, **Sxx_inv, *mx, *e, *Le, *x, *Sxy, *h, y; 
  double dm, Syy, sy, ds, my, pdf, logpdf; 
   
  /* Allocate aux memory */ 
  Sxx = dmatrix(0,px,0,px); 
  Sxx_inv = dmatrix(0,px,0,px); 
  mx = dvector(0,px); 
  e = dvector(0,px); 
  Le = dvector(0,px); 
  x = dvector(0,px); 
  Sxy = dvector(0,px); 
  h = dvector(0,px); 
 
  /* compute submatrices */ 
  dsubmatrix(S,p,p,Sxx,px,px,xi, xi); 
  dsubvector(S[yi],p,Sxy,px, xi); 
  dsubvector(m,p,mx,px,xi); 
  dsubvector(z,p,x,px,xi); 
  Syy = S[yi][yi]; 
  d_invert(px,Sxx,Sxx_inv); 
 
  /* E(y|x) = my + Sxy/Sxx(x-mx) */ 
  x_min_y(mx,x,e,px); 
  Ax(Sxx_inv,e,Le,px); 
  xy(Le,Sxy,&dm,px); 
  my = m[yi]+dm; 
 
  /* Var(y|x) = Syy - Sxy*Sxx_inv*Sxy */ 
  Ax(Sxx_inv,Sxy,h,px); 
  xy(h,Sxy,&ds,px); 
  sy = Syy - ds; 
 
  /* free memory */ 
  free_dmatrix(Sxx,0,px,0,px); 
  free_dmatrix(Sxx_inv,0,px,0,px); 
  free_dvector(mx ,0,px); 
  free_dvector(e  ,0,px); 
  free_dvector(Le ,0,px); 
  free_dvector(x  ,0,px); 
  free_dvector(Sxy,0,px); 
  free_dvector(h  ,0,px); 
 
  /* compute normal pdf */ 
  y = z[yi]; 
  logpdf = d_normal_dens(y,my,sy); 
  pdf = exp(logpdf); 
  return(pdf); 
} 
 
 
/********************************************  
 *         normal_marg_logdens 
 ********************************************/ 
/* log density of x[1,..,i-1,i+1,..p], 
   if x[1,..,p] ~ N(mu,HH')                  
  ------------------------------------------    
   Args: 
   p:      dimensionality 
   i_marg: i in above description 
   n:      number of rows in x 
   x:      n x p matrix of n vectors x[i,] 
   H:      chol-decomp of S, i.e. HH' = S 
   f:      vector of log marg densities 
           (OUTPUT)                       
  ----------------------------------------- */ 
void mvn_marglogdens(int p, int i_marg, int n, double *x,  
		     double *mu, double **H, double *f) 
{ 
  int i,j; 
  double *xi, *z,log_f,two_pi, det_H, *y; 
 
  z = dvector(0,p); 
  y = dvector(0,p); 
 
  two_pi = -0.5*(p-1)*log(2.0*M_PI); 
 
  for(j=0,det_H=0;j<i_marg;j++) 
    det_H += log(H[j][j]); 
  for(j=i_marg+1;j<p;j++) 
    det_H += log(H[j][j]); 
  for(i=0;i<n;i++){ 
    xi = &x[ i * p]; 
    x_min_y(xi,mu,y,p); /* check this !!!!! */ 
    Ax_j(H,xi,z,i_marg,p); 
    xx_j(z,&log_f,i_marg,p); 
    log_f *= -0.5; 
    f[i] = det_H+log_f+two_pi; 
  } 
 
  free_dvector(z,0,p); 
  free_dvector(y,0,p); 
} 
 
/*------------------------------------------*  
 * same as normal_marg_logdens, 
 * but assumes S = I 
 *----------------------------------------- */ 
 
void mvnI_marglogdens (int p, int i_marg, int n, double *x,  
		       double *mu, double *f) 
{ 
  int i; 
  double *xi, log_f,two_pi; 
  double *z; 
 
  z = dvector(0,p); 
 
  two_pi = -0.5*(p-1)*log(2*M_PI); 
 
  for(i=0;i<n;i++){ 
    xi = &x[ i * p]; 
    x_min_y(xi,mu,z,p); 
    xx_j(z,&log_f,i_marg,p); 
    log_f *= -0.5; 
    f[i] = log_f+two_pi; 
  } 
 
  free_dvector(z,0,p); 
 
} 
/* =====================================================================  
   Mixture normals 
   ===================================================================== */  
/******************************************** 
 *         ld_mixture_normal 
 ********************************************/ 
/* computes log density for mixture of normals: 
   
      p(y) =  
      a/(n+a)*N(m0,tau*S) + \sum_{i=1}^k ai[i]/(n+a)*N(mu[i],S) 
     
   where S = V'V and L = 1/V 
   theta is a m x p matrix of parameter vectors theta[j,] 
  ------------------------------------------   
   VALUE: 
      list with log density values 
  ------------------------------------------   
   ARGS: 
   m:   number of parameter vectors 
   k:   number of terms in the mixture 
   a,ai:       coefficients of mixture terms 
   theta:      m x p matrix of parameter vectors 
   m0:  first term in mixture 
   mu:  k x p matrix with means of mixture terms 
   tau: scalar coefficient for cov matrix of first mixture term 
   L:   L = V^(-1), where V'V = R, correlation matrix for mixture terms 
   sd:  together with R = (1/L)'(1/L) defines cov matrix S 
   logdens:      returns logdensity values for the m parameter vectors 
 ********************************************/ 
 
void ld_mixture_normal (int m,double **theta,int p, 
		       double *m0, double tau, int k,int a, double *ai, 
		       /* note: !! check ai as double array !! */ 
		       /******** !!!!!!!! check int a */ 
		       double **mu, double **L, double *sd, 
		       double *logdens) 
{ 
  int i,j; 
  double *logd_j, *pr_j, logd0, dens0, expt, 
  *thi, *x, *y, *z, *prior, *muj; 
   
  int dbg = 0; 
 
  logd_j = dvector(0,m); 
  pr_j = dvector(0,m); 
  x = dvector(0,p); 
  y = dvector(0,p); 
  z = dvector(0,p); 
  prior = dvector(0,p); 
 
  /* loop over all m par vectors */ 
  for(i=0;i<m;i++){ 
    thi = theta[i]; 
 
    /* loop over all k mixture terms */ 
    for(prior[i]=0,j=0; j<k; j++){ 
      muj = mu[j]; 
 
      /* compute z = L*diag(1/sd)*(thi-muj) */ 
      x_min_y(thi,muj,x,p); 
      x_div_y(x,sd,y,p); 
      xA(y,L,z,p); 
       
      /* evaluate log posterior z'z/2 */ 
      xtx(z,&expt,p); 
      logd_j[j] = -expt/2; 
 
      /* add up prior for i-th par vector */ 
      pr_j[j] = exp(logd_j[j]) * ai[j]; 
      if (dbg){ 
	messtabdouble("logd[0] ", logd_j[j]); 
	messtabdouble("prj[0] ", pr_j[j]);  
      } 
      prior[i] += pr_j[j]; 
    } 
 
    /* now add pr(thi) if thi comes from G0 */ 
    /* first compute again z = L*diag(1/sd)*(thi-m0) */ 
    x_min_y(thi,m0,x,p); 
    x_div_y(x,sd,y,p); 
    xA(y,L,z,p); 
 
    /* compute pr0 and add to prior[i] */ 
    xtx(z,&expt,p); 
    logd0 = -expt*0.5/(1+tau) - p*0.5*(log(1.0+tau)); 
    dens0  = a*exp(logd0); 
    prior[i] += dens0; 
    if (dbg) messdouble(" m0: ", dens0);  
    if (prior[i] < 1e-10) 
      logdens[i] = -1e10; 
    else  
      logdens[i] = log(prior[i]); 
  } 
 
  free_dvector(logd_j,0,m); 
  free_dvector(pr_j,0,m); 
  free_dvector(x,0,p); 
  free_dvector(y,0,p); 
  free_dvector(z,0,p); 
  free_dvector(prior,0,p); 
} 
 
/******************************************** 
 *         mdp_logpdf 
 ********************************************/ 
/* computes log density for mixture of normals: 
   
      p(y) =  
      a/(n+a)*N(m0,tau*S) + \sum_{i=1}^k ai[i]/(n+a)*N(mu[i],S) 
     
   where S = VV' and L = 1/V 
   theta is a m x p matrix of parameter vectors theta[j,] 
  ------------------------------------------   
   VALUE: 
      list with log density values 
  ------------------------------------------   
   ARGS: 
   m:   number of parameter vectors 
   k:   number of terms in the mixture 
   a,ai:       coefficients of mixture terms 
   theta:      m x p matrix of parameter vectors 
   m0:  first term in mixture 
   mu:  k x p matrix with means of mixture terms 
   tau: scalar coefficient for cov matrix of first mixture term 
   L:   L = V^(-1), where VV' = R, correlation matrix for mixture terms 
   sd:  together with R = (1/L)(1/L)' defines cov matrix S 
   logdens:      returns logdensity values for the m parameter vectors 
 ********************************************/ 
 
double mdp_logpdf(double *theta,int p, 
	       double *m0, double **B_cd_inv, int k, double a, double *ai, 
	       double **mu, double **L) 
{ 
  int j; 
  double *logd_j, *pr_j, logd0, dens0, expt, 
  *x, *z, *muj, dens, logdens, 
  det_minushalf, log_coef, 
  det_minushalf0, log_coef0; 
   
  logd_j = dvector(0,k); 
  pr_j = dvector(0,k); 
  x = dvector(0,p); 
  z = dvector(0,p); 
   
  /* loop over all k mixture terms */ 
  det_minushalf = fabs(prod_diag(L,p)); 
  log_coef = log(det_minushalf)-0.5*p*log(2.0*M_PI); 
 
  for(dens = 0.0,j=0; j<k; j++){ 
    muj = mu[j]; 
     
    /* compute z = L*diag(1/sd)*(thi-muj) */ 
    x_min_y(theta,muj,x,p); 
    Ax(L,x,z,p); 
     
    /* evaluate log posterior z'z/2 */ 
    xtx(z,&expt,p); 
    logd_j[j] = -0.5*expt; /* + log_coef */ 
     
    /* add up prior for i-th par vector */ 
    pr_j[j] = exp(logd_j[j]) * ai[j]; 
    dens += pr_j[j]; 
  } 
   
  /* now add pr(thi) if thi comes from G0 */ 
  /* first compute again z = L*diag(1/sd)*(thi-m0) */ 
  det_minushalf0 = fabs(prod_diag(B_cd_inv,p)); 
  log_coef0 = log(det_minushalf0)-0.5*p*log(2.0*M_PI); 
  x_min_y(theta,m0,x,p); 
  Ax(B_cd_inv,x,z,p); 
   
  /* compute pr0 and add to prior[i] */ 
  xtx(z,&expt,p); 
  logd0 = -expt*0.5 + log_coef0 - log_coef; /*scaled by log_coef */ 
  dens0  = a*exp(logd0); 
  dens += dens0; 
  if (dens <= 0) 
    logdens = -1e10; 
  else  
    logdens = log(dens); 
   
   
  free_dvector(logd_j,0,k); 
  free_dvector(pr_j,0,k); 
  free_dvector(x,0,p); 
  free_dvector(z,0,p); 
 
  return(logdens); 
} 
 
/******************************************** 
 *         mvnmixI_logdens 
 ********************************************/ 
/* computes log density for mixture of mvnormals with diag cov matrix 
   
      p(x) =  
       sum N(x; mu[i],diag(sd)) = 
       sum N(z; nu[i],I) 
   where z[j] = x[j]/sd[j], and nu[j] = mu[j]/sd[j] */ 
/********************************************/ 
 
double mvnmixI_logdens(double *x, int k, double **nu, double *sd, int p) 
/* nu[i] are the means of the normal kernels/sd[i] 
   sd    is the diagonal of the cov matrix (same for all terms) 
   p     is the dimension  
   x     is the location at which the log dens is evaluated */ 
{ 
  int j; 
  double  
    *logd,expt,*z,*delta,dens,logdens; 
   
  /* alloc mem */ 
  logd = dvector(0,k); 
  z    = dvector(0,p); 
  delta= dvector(0,p); 
 
  /* tra sform x to z */ 
  x_div_y(x,sd,z,p); 
 
    /* loop over all k mixture terms */ 
    for(j=0,dens=0.0; j<k; j++){ 
 
      /* compute z = L*diag(1/sd)*(thi-muj) */ 
      x_min_y(z,nu[j],delta,p); 
       
      /* evaluate log posterior delta'delta/2 */ 
      xtx(delta,&expt,p); 
      logd[j] = -expt*0.5; 
 
      /* add up prior for i-th par vector */ 
      dens += exp(logd[j]); 
    } 
 
 
  logdens = log(dens); 
 
  /* free mem */ 
  free_dvector(logd,0,k); 
  free_dvector(delta,0,p); 
  free_dvector(z,0,p); 
   
  return logdens; 
} 
 
/********************************************  
 *         rcond_mixture_normal 
 ********************************************/ 
/* simulates a draw from a complete conditional 
   of a mixture of normals. 
   
      p(y) =  
      a/(n+a)*N(m0,tau*S) + \sum_{i=1}^k ai[i]/(n+a)*N(mu[i],S) 
     
   where S = V'V and L = 1/V 
   theta is a m x p matrix of parameter vectors theta[j,] 
  ------------------------------------------   
   VALUE: 
     none 
  ------------------------------------------   
   ARGS: 
   m:   number of parameter vectors 
   k:   number of terms in the mixture 
   a,ai:       coefficients of mixture terms 
   theta:      m x p matrix of parameter vectors 
   m0:  first term in mixture 
   mu:  k x p matrix with means of mixture terms 
   tau: scalar coefficient for cov matrix of first mixture term 
   L:   L = V^(-1), where V'V = R, correlation matrix for mixture terms 
   sd:  together with R = (1/L)'(1/L) defines cov matrix S 
   logdens:      returns logdensity values for the m parameter vectors 
 ********************************************  
 * version of rcond_mixture_normal with S=I 
   assume m0 = (0,..0) 
   only one draw  
  ------------------------------------------*/ 
void rcond_mixture_normal_I(int p, double *x_cond, int i_cond, 
			    int k, double *aj, double **mu, double tau, 
			    /* note: aj is a double array !!! */ 
			    double a, 
			    double *x, double *pdf_old, double *pdf_new) 
{ 
  int i,j, j_chosen; 
  double  
    *pr, *logf_marg, max_logf=0,z,*z_cond, mn, 
    x_old, *dummy; 
 
  pr = dvector(0,k+1); 
  logf_marg = dvector(0,k+1); 
  z_cond = dvector(0,p); 
  dummy = dvector(0,p); 
 
/*------------------------------------------*/ 
/* compute marginal prob of each mixture term */ 
/* first the k terms: */ 
  for(i=0, max_logf = -10e10;i<k;i++){ 
    mvnI_marglogdens(p,i_cond,1,x_cond, mu[i], &logf_marg[i]); 
    if (logf_marg[i] > max_logf) max_logf = logf_marg[i]; 
  } 
 
 
/* then the m0 term: */ 
  for(j=0;j<p;j++) dummy[j] = 0.0; 
  y_assgn_x(z_cond,x_cond,p); 
  x_div_r(z_cond,sqrt(tau),z_cond,p); 
  mvnI_marglogdens(p,i_cond,1,z_cond, dummy, &logf_marg[k]); 
  logf_marg[k] -= 0.5*(p-1)*log(tau); 
 
  if (logf_marg[i] > max_logf) max_logf = logf_marg[i]; 
  for(i=0;i<k;i++) 
    pr[i] = aj[i]*exp(logf_marg[i]-max_logf); 
  pr[k] = a*exp(logf_marg[k]-max_logf); 
   
 
/*------------------------------------------*/ 
/* choose mixture term and sample */ 
  multinomial(1,k+1,pr,&j_chosen); 
  dstdnormal(1,&z); 
  if (j_chosen == k) 
    *x =  z*sqrt(tau); 
  else 
    *x = mu[j_chosen][i_cond] + z; 
 
/*------------------------------------------*/ 
/* compute pdf of old point x_cond[i_cond] and new point x */ 
/* first: the k mixture terms: */ 
  x_old = x_cond[i_cond]; 
  for(i=0,*pdf_old = *pdf_new = 0.0;i<k;i++){ 
    mn = mu[i][i_cond]; 
    *pdf_old += pr[i]*exp(-(x_old-mn)*(x_old-mn)/2); 
    *pdf_new += pr[i]*exp(-(*x-mn)*(*x-mn)/2); 
  } 
 
/* then: the m0 term: */ 
  *pdf_old += pr[k]*exp( -x_old*x_old/(2*tau))/sqrt(tau); 
  *pdf_new += pr[k]*exp( -(*x)*(*x)/(2*tau))/sqrt(tau); 
 
  free_dvector(pr,0,k+1); 
  free_dvector(logf_marg,0,k+1); 
  free_dvector(z_cond,0,p); 
  free_dvector(dummy,0,p); 
 
} 
 
/* =====================================================================  
   Random Matrices 
   ===================================================================== */  
 
/* *********************************************************** 
   wishart 
   *********************************************************** */ 
 
void wishart(int nu,             /* degrees of freedom */ 
	int p,              /* dimension  */ 
	double **L,         /* LL' = Sigma */ 
	double **W)         /* Wishart variate (output) */ 
/* 
                  (nu-p-1)/2               -1 
               |A|            exp(-0.5 tr(S  A) 
p(A; S,nu,p)= ------------- 
	          nu/2 
	       |S| 
E(A) = nu*S, nu >= p,  
       
*/ 
{ 
  double *zeta, **eps, **V, **LV; 
  double sumeps2, sumepsij; 
  int i, j, k, nu_zeta; 
   
  zeta = dvector(0,p-1); 
  eps = dmatrix(0,p-1,0,p-1); 
  V = dmatrix(0,p-1,0,p-1); 
  LV = dmatrix(0,p-1,0,p-1); 
 
  /* generate zeta[i] from Chi-Sq(nu-i) */ 
  for(i=0;i<p;i++){ 
    nu_zeta = nu-i; 
    zeta[i] = gamdev_ab(nu_zeta*0.5,0.5); 
 
    /* eps is a right upper triang array of normals */ 
    dstdnormal(p-i,&eps[i][i]); 
  } 
   
  /* V[i][i] = zeta[i] + sum{k<i} eps[k,i]^2 */ 
  for(i=0;i<p;i++){ 
    for(k=0,sumeps2=0; k<i; k++) 
      sumeps2 += eps[k][i]*eps[k][i]; 
    V[i][i] = zeta[i] + sumeps2; 
 
    /* For i<j: 
       V[i][j] = eps[i,j]*sqrt(zeta[i]) + sum{k<i} eps[k,i]*eps[k,j] */ 
    for(j=i+1; j<p; j++){ 
      for(k=0,sumepsij=0; k<i; k++) 
	sumepsij += eps[k][i]*eps[k][j]; 
      V[i][j] = V[j][i] = eps[i][j]*sqrt(zeta[i]) + sumepsij; 
    } 
  } 
   
  AB(L,p,p,V,p,p,LV); 
  ABt(LV,p,p,L,p,p,W); 
 
  free_dvector(zeta,0,p-1); 
  free_dmatrix(eps,0,p-1,0,p-1); 
  free_dmatrix(V,0,p-1,0,p-1); 
  free_dmatrix(LV,0,p-1,0,p-1); 
 
} 
 
double log_dwishart(int nu,int p,double **V_inv,double **B) 
/* Sigma = VV', A=BB' is the Wishart covariate */ 
{ 
  double  
    f,tr,vb; 
  int  
    i,j,k; 
 
  /* tr = tr S^-1A = tr (V_inv*B)*(V_inv*B)' = 
     = sum_i (sum_j (sum_k V_inv[i][k]*B[k][j])^2) */ 
 
  for(i=0,tr=0; i<p; i++){ 
    for(j=0;j<p;j++){ 
      for(k=0,vb=0; k<p; k++) 
	vb += V_inv[i][k]*B[k][j]; 
      tr += vb*vb; 
    } 
  } 
   
  f = (nu-p-1.0)*log(prod_diag(B,p)) + 
    nu*log(prod_diag(V_inv,p)) - 
      0.5*tr; 
  return(f); 
} 
     
 
   
 
/* ***********************************************************  
   empirical dist 
   *********************************************************** */ 
void moments(double *mn, double *var, double *x, int n) 
/* returns mean & variance of x[0..n] */ 
{ 
  int i; 
  double m1,m2; 
 
  for(i=0,m1=m2=0;i<n;i++){ 
    m1 += x[i]; 
    m2 += x[i]*x[i]; 
  } 
  *mn = m1/(1.0*n); 
  *var = (m2/(1.0*n) - m1*m1/(1.0*n*n)); 
} 
 
void std(double *x, double *z, int n) 
/* returns in z the standardized version of x */ 
{ 
  double m,s,v; 
  int i; 
 
  moments(&m,&v,x,n); 
  s = sqrt(v); 
  for(i=0;i<n;i++) 
    z[i] = (x[i]-m)/s; 
} 
