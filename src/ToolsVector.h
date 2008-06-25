#if !defined(VECTOR_H_INCLUDED)
#define VECTOR_H_INCLUDED

/*+
 * vector.h - Matrix Algebra
 *
 * Description:
 *		Declaration of various matrix routines in the library.
 *
 * Modification History:
 * <unknown>	1.1		Peter Mueller
 *		Created.
 *
 * 29-JUN-1999	1.2		Stephen Morris (Tessella Support Services plc)
 *		Remove embedded comments (slash-star in a comment block).
 *		Added header and #if...#endif sentinels.
-*/

/* -----------------------------------------------------------
   matrix algebra
   ----------------------------------------------------------- */

double sqr(double x);
int ipow(int r, int s);
int quad_eq(double a, double b, double c, 
	     double *x1, double *x2);
int ifloor(double x);
int iround(double x);

void r_swap_s(double *r, double *s);
void i_swap_j(int *i, int *j);
double mod(double x, double y);
int int_div(double x, double y);
double max_xy(double x, double y);
double min_xy(double x, double y);

void x_zero(double *, int);
double x_max(double *x, int *i, int p);
double i_min(int *x, int *i, int p);
double i_max(int *x, int *i, int p);
double x_min(double *x, int *i, int p);
void i_zero(int *x, int p);
void y_assgn_x(double *, double *, int);
void j_assgn_i(int *y, int *x, int p);
void x_assgn_r(double *, double, int );
void i_assgn_r(int *x, int r, int p);
void x_swap_y(double *x, double *y, int p);
void xtx (double *, double *value, int);
double norm (double *x, int p) ;
void xxt (double *, double **A, int);
void prod (double *x, double *value, int p);
void sumx (double *x, double *value, int p);
void sumi (int *x, int *value, int p);
void xy (double *, double *, double *value, int);
void x_plus_y(double *, double *, double *, int);
void x_plus_r(double *x, double r, double *z, int p);
double x_dist_y(double *x, double *y, int p);
void rx(double, double *, int);
void rx_plus_m(double, double *, double, double *, int);
void sort (double *x, int *ind, int p) ;
void sorti (int *x, int *ind, int p, int up) ;
void sort2 (double *x, int p) ;
void grid (double x0, double xn, int n, double *x);

void x_min_y(double *, double *, double *, int);
void x_div_y(double *, double *, double *, int);
void x_mult_y(double *x, double *y, double *z, int p);
void x_div_r(double *, double, double *, int);
void x_plus_ry(double *, double, double *, double *, int);
void x_plus_ryminz(double *x, double r, double *y, double *z, 
		   double *w, int p) ;
void a_plus_rxy(double, double, double *, double *, double *, int);
void rx_plus_sy (double, double *, double, double *, double *, int);

/* -----------------------------------------------------------
   matrix only
   ----------------------------------------------------------- */

void A_zero(double **, int);
double R_max(double **x, int *i, int *j, int p, int q);
double R_min(double **x, int *i, int *j, int p, int q);
void I_zero(int **A, int p,int q);
void R_zero(double **, int, int);
void A_unit(double **, int);
void A_swap_B(double **A, double **B, int p, int q);
void Aj_swap_Ak(double **A, int j, int k, int p, int q);
void At(double **, double **, int);
void Rt(double **, double **, int, int);
void At_self(double **, int);
void B_assgn_A(double **, double **, int);
void A_assgn_r(double **A, double r, int p);
void R_assgn_r(double **A, double r, int p, int q);
void y_assgn_Aj(double *y, double **A, int j, int p);
void R_assgn_T(double **, double **, int, int);
void  A_plus_B(double **, double **, double **, int, int);
void  I_plus_J(int **A, int **B, int **C, int p, int q);
void A_minus_B(double **A, double **B, double **C, int p, int q);
void AplusBr(double **, double **, double, double **, int);
void AplusBi(double **A, double **B, int i, double **C, int p);
void rA_plus_sB(double, double **, double, double **, double **, int);
void  AB(double **, int, int, double **, int, int, double **);
void  AtB(double **, int, int, double **, int, int, double **);
void  ABt(double **, int, int, double **, int, int, double **);
void AtBA(double **A, int p, double **B, double **C);
void RtBR(double **A, int p, int q, double **B, double **C) ;
void rSTU (double r, 
	   double **S, int pS, int qS, 
	   double **T, int pT, int qT, 
	   double **U, int pU, int qU,
	   double **V);
void rR_plus_sSTU (double r, double **R, double s,
		  double **S, int pS, int qS, 
		  double **T, int pT, int qT, 
		  double **U, int pU, int qU,
		  double **V);
void rx_plus_sSTy (double r, double *x, 
		   double s,
		  double **S, int pS, int qS, 
		  double **T, int pT, int qT, 
		  double *y, int pU, 
		  double *z);
void rR_plus_sSTSt (double r, double **R, double s,
		   double **S, int pS, int qS, 
		   double **T, int pT, int qT, 
		   double **V);
void sSTSt (double s,
		   double **S, int pS, int qS, 
		   double **T, int pT, int qT, 
		   double **V);
double prod_diag(double **, int);
double trace(double **A, int p);

/* ------------------------------------------------- 
   diag matrix
   ------------------------------------------------- */ 

void L_assgn_x(double **L, int p, double *x);
void L_assgn_recx(double **L, int p, double *x);
void  AtLB(double **, int, int, double *, 
	   double **, int, int, double **);
void  ALBt(double **, int, int, double *, 
	   double **, int, int, double **);
void  LminusHalfAB(double *,double **, int, int, 
		     double **, int, int, double **);
void  LminusHalfAt(double *,double **, int, int, double **);
void  AdivDiagB(double **, double *, int, int);

/* -----------------------------------------------------------
   matrix/vector
   ----------------------------------------------------------- */

void  xA(double *,double **,double *, int);	
void  Ax(double **,double *,double *, int);	
void  Rx(double **, int, int,double *, double *);	
void  Ri(double **R, int p, int q, int *x, double *z)	 ;
void  Rtx(double **, int, int,double *, double *);	
void  rA(double,double **, double **, int, int);	
void Ax_plus_y(double **, double *, double *, double *, int);
void Rx_plus_y(double **, double *, double *, double *, int, int);
void rAx_plus_sBy(double , double **, double *, 
		  double , double **, double *, double *, int);
void x_plus_rAy(double *, double, double **, double *, 
		double *, int, int);
void rAx_plus_y(double, double **, double *, double *, 
		double *, int);
void A_plus_rxxt(double **, double, double *, double **, int);
void A_xminusy(double **, double *, double *, double *, int);
void A_xminusy_plus_z(double **A, double *x, double *y, double *z,
		      double *w, int p);
double xtAy (double *x, double *y, double **A, int p);
void ABAt(double **A, double **B, double **C, int p);
void ABAt_plus_C(double **A,double **B,double **C,double **D,int p);
void RSRt_plus_C(double **A,double **B,double **C,double **D,int p,
		 int q);
void C_min_RSRt(double **A,double **B,double **C,double **D,int p,
		 int q);

/* -----------------------------------------------------------
   submatrix/vector
   ----------------------------------------------------------- */

void  Ax_j(double ** ,double *,double *, int, int);
void  xx_j(double * , double *, int, int); 
void  dsubmatrix(double **, int, int, 
		 double **, int, int, int *, int *);
void  dsubvector(double *, int, double *, int, int *);
void dsubmatrix_j(double **a, int rc, double **m, int p);

void permutate (int *y, int n);

#endif

