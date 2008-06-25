#if !defined(RAND_H_INCLUDED)
#define RAND_H_INCLUDED

/*+
 * rand.h - Random Number Routines
 *
 * Description:
 *		Declaration of various random number routines in the library.
 *
 * Modification History:
 * <unknown>	1.1		Peter Mueller
 *		Created.
 *
 * 29-JUN-1999	1.2		Stephen Morris (Tessella Support Services plc)
 *		Remove embedded comments (slash-star in a comment block).
 *		Made ld_mixture_normal and rcond_mixture_normal_I of type void.
 *		Added header and #if...#endif sentinels.
-*/

void setseed(long, long);
void getseed(long *is1, long *is2);
/* =========================================================== *
 *		 UNIVARIATE
 * =========================================================== */
void	duniform(int , double *);
double  runif();
void	uniform(int , float *);
void	dstdnormal(int ,double *);
void	stdnormal(int , float *);
long l_poisson(double mu);
double  norm_rand();

double	cdfnorm(double , double , double );
double	invcdfn (double , double , double );
double norm_trunc_left(double , double , double);
double norm_trunc_right(double , double , double);
double norm_trunc(double , double , double, double);
double norm_trunc2(double m, double s,  
		   double left, double right,
		   double *lpdf);
double m_norm_trunc_left(double m, double s, double cut) ;

double logpdfnorm(double , double , double );
double d_normal_dens(double , double , double );
double pdfnorm(double , double , double );
double	pdflognormal(double, double, double);
double t_logpdf(double x, double m, double sd, int nu);
double tc_logpdf(double sd, int nu);
double betadev(double , double );
double beta_logpdf(double x, double a, double b);
double	cdfbeta(double y, double a, double b);
double	invcdfbeta (double cdf, double a, double b);

void beta_pars(double mu, double var, double *a, double *b);
void dirichlet(double *alpha, double *w, int p);
double gamdev(double );
double gamdev_ab(double , double );
double gamma_logpdf(double x,double a, double b);
double gam_trunc_left(double a, double b,  double cut) ;
double gam_trunc_right(double a, double b, double cut) ;
double gam_trunc(double a, double b, double left, double right) ;
double invcdfgam (double cdf, double a, double b) ;
double cdfgam_ab(double y, double a, double b);
double expdev(double beta);
double expdev_truncleft(double beta, double x0);

double d_weibull(double *, double, double);
int	binomial(int , double );
int binomial2(int n, double p);
double logp_binomial(int y, int n, double p);
int rand_bebinomial(int n, double a, double b);
double logpdf_bebinomial(int y, int n, double a, double b, int c);

void multinomial(int, int, double *, int *);
void multinomial_I(int, int, int *, int *);
void multinomial_table(int n_draw, int n_cell, double *pr, int *ct);

/* =========================================================== *
 *		 MULTIVARIATE
 * =========================================================== */
void mvn_rand(int,double *,double **,double *);
void mvntrunc_rand_init(int p, double *mean, double **L, 
  double *left, double *right, double *x, double *z);
void mvntrunc_rand(int p, double *mu, double **L, 
		 double *left, double *right, 
		 double *x, double *z,
		 int niter)   ;
void inequ(double left, double right, double m, double *L, 
	   double *x, int i,
	   double *lo, double *hi, int p);
double mvn_cdf(int p, double *mu, double **V, double **L,
	       double *left, double *right, int I);
double mvnS_rand(int,double *,double **,double *);   // added return value
double mvn_pdf(double *, double *, double **, int );
double mvn_logpdf(double *, double *, double **, int );
double mvnS_logpdf(double *, double *, double **, int );
void mvn_marglogdens(int , int , int , double *,double *, double **, double *);
void mvnI_marglogdens(int, int, int, double *, double *, double *);
double mvn_cond_pdf(double *, double *, double **, int , int , int* , int );
double mvt_pdf(double *, double *, double **, int, int );
double mvt_logpdf(double *x, double *m, double **L_inv, int nu,int p);
double mvtS_logpdf(double *x, double *m, double **S, int nu, int p);
void cond_nmoments(int k, int p, double *x,double *mu, double **S, 
		   double **S22_inv, double *s, double *m);
void cond_n2moments(int k, int p, double y, double *mu, double **V,
		    double *s, double *m);
void cond_mvnmoments(int p, int q, double *x, double *mu, double **S,
		     double **S22_inv, int *qi, int *qi_c, double **V,
		     double *m);
void ld_mixture_normal (int ,double **, int ,
		       double *, double , int ,int , double *,
		       double **, double **, double *,
		       double *);
double mdp_logpdf(double *theta,int p,
	       double *m0, double **B_cd_inv, int k,
		  double a, double *ai, double **mu, double **L);
void rcond_mixture_normal_I(int , double *, int ,
			    int , double *, double **, double ,
			    double ,
			    double *, double *, double *);
double mvnmixI_logdens(double *x, int k, double **nu, double *sd, int p);

/* =========================================================== *
 *		 MATRIX VARIATE
 * =========================================================== */
void wishart(int   ,        /* degrees of freedom */
	int  ,              /* dimension  */
	double ** ,         /* LL' = Sigma */
	double ** );        /* Wishart variate (output) */
double log_dwishart(int nu,int p,double **V_inv,double **B);

void moments(double *mn, double *var, double *x, int n);
void std(double *x, double *z, int n);

#endif

