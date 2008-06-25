#if !defined(REGR_H_INCLUDED)
#define REGR_H_INCLUDED

/*+
 * regr.h - Regression
 *
 * Description:
 *		Declaration of various regression routines in the library.
 *
 * Modification History:
 * <unknown>	1.1		Peter Mueller
 *		Created.
 *
 * 29-JUN-1999	1.2		Stephen Morris (Tessella Support Services plc)
 *		Remove embedded comments (slash-star in a comment block).
 *		Added header and #if...#endif sentinels.
-*/

/* **********************************************************************
 * Multiple regression 
 * ********************************************************************** */
void regr(double *mu, double **H_inv_cd, double **H_cdt,
	  double *y, double *h, double **x, int n, int p,
	  int var_known);
void regr1(double *mu, double *sd, double *prec,
	  double *y, double *h, double *x, int n, int var_known);
void bayes_regr(double *mu2, double **V, double **L,
	  double *y, double *h, double **x, 
	  double *mu1, double **S1_inv, int n, int p,
	  int var_known);
void bayes_regr1(double *mu, double *sd, double *prec,
	  double *y, double *h, double *x,
	  double mu1, double s1, int n);
void logit_setup(double *mu, double **V, double **L,
		 double *y, double **x, double *offset, int N, int n,
		 int p, double *sig); 

double lp_logistic(double *beta, double *y, double **x, 
		   double *offset, int N, int n, int p, double *sig);

void regr_logistic(double *beta, 
		   double *y, double **x, double *offset, int N,
		   int n, int p, double *sig, int n_iter, int
		   *n_jump);

void bayes_regr_rand(double *th,
	  double *y, double *h, double **x, 
	  double *mu1, double **S1_inv, int n, int p,
	  int var_known);
void bayes_regr_old(double *mu2, double **V, double **L,
	  double *y, double *h, double **x, 
	  double *mu1, double **S1_inv, int n, int p,
	  int var_known);
void bayes_logit_setup(double *mu, double **V, double **L,
		       double *m1, double **S1_inv,
		       double *y, double **x, double *offset, int N,
		       int n, int p, double *sig);
void bayes_logit_setup2(double *mu, double **V, double **L,
		       double *m1, double **S1_inv,
		       double *y, double **x, double *offset, int *N,
		       int n, int p, double *sig, int normal);
double lp_logistic2(double *beta, double *y, double **x, 
		    double *offset, int *N,int n, int p, double *sig,
		    int normal);
void bayes_regr_logistic(double *beta, 
   double *m1, double **S1_inv,		       
   double *y, double **x, double *offset, int N,
   int n, int p, double *sig, int n_iter, int *n_jump); 
void bayes_regr_logistic2(double *beta, 
   double *m1, double **S1_inv,		       
   double *y, double **x, double *offset, int *N,
   int n, int p, double *sig, int n_iter, int *n_jump, int normal);

#endif

