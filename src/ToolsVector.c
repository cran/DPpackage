static const char vector_c_sccs_id[] = "%W%"; 
 
/*+ 
 * vector.c - Vector Routines 
 * 
 * Description: 
 *		Miscellaneous vector routines. 
 * 
 * Modification History: 
 * <unknown>	1.1		Peter Mueller 
 *		Created.  Modified on 21.6.96 to add norm(). 
 * 
 * 29-JUN-1999	1.2		Stephen Morris (Tessella Support Services plc) 
 *		Removed embedded comments (slash-star in a comment block). 
 *		Removed unused variables. 
 *		Include "mess.h" to define error(). 
 *		Add a few casts to "int" to remove warnings from MSVC compiler. 
-*/ 
 
#include <math.h> 
 
#include "ToolsVector.h" 
#include "ToolsMess.h" 
#include "ToolsNrutil.h" 
#include "ToolsRand.h" 
 
/* ===================================================================== 
   Matrix algebra 
   ===================================================================== */ 
 
int	_i, _j, _k, _l; 
double *_x; 
 
/* ------------------------------------------------- 
   scalars 
 ------------------------------------------------- */ 
double sqr(double x) 
{ 
  return x*x; 
} 
 
int ipow(int r, int s) 
{ _i = (int) floor( exp(log(1.0*r)*s*1.0)+0.5 ); 
  return(_i); 
} 
 
int quad_eq(double a, double b, double c, 
	     double *x1, double *x2) 
{ /* solves quadratic equation 
     a*x^2+b*x+c = 0 
     solution in x1,x2 
     returns 
     -1 if not solvable */ 
  double 
    p,q,p4q; 
 
  p = b/a; 
  q = c/a; 
  p4q = p*p*0.25-q; 
  if (p4q < 0) 
    return(-1); 
  *x1 = -0.5*p-sqrt(p4q); 
  *x2 = -0.5*p+sqrt(p4q); 
  return(0); 
} 
 
 
int ifloor(double x) 
{ _i = (int) floor(x+0.5); /*!!!??? */ 
  return(_i); 
} 
 
int iround(double x) 
{ _i = (int) floor(x+0.5); 
  return(_i); 
} 
 
 
double mod(double x, double y) 
{ 
  return(x-y*floor(x/y)); 
} 
 
int int_div(double x, double y) 
{ 
  _i = (int) floor(x/y); 
  return(_i); 
} 
void r_swap_s(double *r, double *s) 
{ 
  double 
    x; 
  x = *r; 
  *r = *s; 
  *s = x; 
} 
 
void i_swap_j(int *i, int *j) 
{ 
  int 
    x; 
  x = *i; 
  *i = *j; 
  *j = x; 
} 
 
double max_xy(double x, double y) 
{ 
  return (x>y) ? x : y; 
} 
 
double min_xy(double x, double y) 
{ 
  return (x<y) ? x : y; 
} 
 
 
 
/* ----------------------------------------------------------- 
   vectors only 
   ----------------------------------------------------------- */ 
void x_zero(double *x, int p) 
{ 
  for(_i=0;_i<p;_i++) 
    x[_i] = 0.0; 
} 
 
 
double x_max(double *x, int *i, int p) 
{ 
  double mx; 
  for(mx=x[0],*i=0,_i=1; _i<p; _i++) 
    if (x[_i]>mx){ 
      mx=x[_i]; 
      *i = _i; 
    } 
  return(mx); 
} 
 
double x_min(double *x, int *i, int p) 
{ 
  double mx; 
  for(mx=x[0],*i=0,_i=1; _i<p; _i++) 
    if (x[_i]<mx){ 
      *i=_i; 
      mx=x[_i]; 
    } 
  return(mx); 
} 
 
double i_min(int *x, int *i, int p) 
{ 
  int mx; 
  for(mx=x[0],*i=0,_i=1; _i<p; _i++) 
    if (x[_i]<mx){ 
      *i=_i; 
      mx=x[_i]; 
    } 
  return(mx); 
} 
 
double i_max(int *x, int *i, int p) 
{ 
  int mx; 
  for(mx=x[0],*i=0,_i=1; _i<p; _i++) 
    if (x[_i]>mx){ 
      *i=_i; 
      mx=x[_i]; 
    } 
  return(mx); 
} 
 
void i_zero(int *x, int p) 
{ 
  for(_i=0;_i<p;_i++) 
    x[_i] = 0; 
} 
 
void x_assgn_r(double *x, double r, int p) 
{ 
  for(_i=0;_i<p;_i++) 
    x[_i] = r; 
} 
 
void   y_assgn_x(double *y, double *x, int p) 
{ 
  for(_i=0;_i<p;_i++) y[_i]=x[_i]; 
} 
 
void   j_assgn_i(int *y, int *x, int p) 
{ 
  for(_i=0;_i<p;_i++) y[_i]=x[_i]; 
} 
 
void   i_assgn_r(int *y, int r, int p) 
{ 
  for(_i=0;_i<p;_i++) y[_i]=r; 
} 
 
void x_swap_y (double *x, double *y, int p) 
{ double h; 
  for(_i=0;_i<p;_i++){ 
    h = x[_i]; 
    x[_i] = y[_i]; 
    y[_i] = h; 
  } 
} 
 
void   rx(double r, double *x, int p) 
{ 
  
  for(_i=0;_i<p;_i++) x[_i] *= r; 
} 
 
void xtx (double *x, double *value, int p) 
{ 
  
  for(_i=0,*value=0;_i<p;_i++) *value += x[_i]*x[_i]; 
} 
 
double norm (double *x, int p) 
{ 
  double value; 
  
  for(_i=0,value=0;_i<p;_i++) value += x[_i]*x[_i]; 
  return(value); 
} 
 
void sumx (double *x, double *value, int p) 
{ 
  for(_i=0,*value=0;_i<p;_i++) *value += x[_i]; 
} 
 
void sumi (int *x, int *value, int p) 
{ 
  for(_i=0,*value=0;_i<p;_i++) *value += x[_i]; 
} 
 
void prod (double *x, double *value, int p) 
{ 
  for(_i=0,*value=1;_i<p;_i++) *value *= x[_i]; 
} 
 
void xxt (double *x, double **A, int p) 
{ 
  for(_i=0;_i<p;_i++) 
    for(_j=0;_j<p;_j++) 
      A[_i][_j] = x[_i]*x[_j]; 
} 
 
void xy (double *x, double *y, double *value, int p) 
{ 
  for(_i=0,*value=0;_i<p;_i++) *value += x[_i]*y[_i]; 
} 
 
void x_plus_y(double *x, double *y, double *z, int p) 
{ 
  for(_i=0;_i<p;_i++) z[_i] = x[_i]+y[_i]; 
} 
 
void x_plus_r(double *x, double r, double *z, int p) 
{ 
  for(_i=0;_i<p;_i++) z[_i] = x[_i]+r; 
} 
 
double x_dist_y(double *x, double *y, int p) 
{ 
  double d; 
 
  for(_i=0,d=0;_i<p;_i++) 
    d += (x[_i]-y[_i])*(x[_i]-y[_i]); 
  return(d); 
} 
 
void rx_plus_m(double r, double *x, double m, double *z, int p) 
{ 
  for(_i=0;_i<p;_i++) z[_i] = r*x[_i]+m; 
} 
 
void x_min_y(double *x, double *y, double *z, int p) 
{ 
  for(_i=0;_i<p;_i++) z[_i] = x[_i]-y[_i]; 
} 
void x_div_y(double *x, double *y, double *z, int p) 
{ 
  for(_i=0;_i<p;_i++) z[_i] = x[_i]/y[_i]; 
} 
 
void x_mult_y(double *x, double *y, double *z, int p) 
{ 
  for(_i=0;_i<p;_i++) z[_i] = x[_i]*y[_i]; 
} 
 
 
void x_div_r(double *x, double r, double *z, int p) 
{ 
  for(_i=0;_i<p;_i++) z[_i] = x[_i]/r; 
} 
void x_plus_ry(double *x, double r, double *y, double *z, int p) 
{ 
  for(_i=0;_i<p;_i++) z[_i] = x[_i]+r*y[_i]; 
} 
void x_plus_ryminz(double *x, double r, double *y, double *z, 
		   double *w, int p) 
{ 
  for(_i=0;_i<p;_i++) w[_i] = x[_i]+r*(y[_i]-z[_i]); 
} 
void a_plus_rxy(double a, double r, double *x, double *y, double *b, int p) 
{ 
  for(_i=0, *b=a; _i<p; _i++) *b += r*x[_i]*y[_i]; 
} 
void rx_plus_sy (double r, double *x, double s, double *y, 
		 double *z, int p) 
{ 
   
  for(_i=0;_i<p;_i++) z[_i] = r*x[_i]+s*y[_i]; 
} 
 
void sort (double *x, int *ind, int p) 
{ /* naive sort in ascending order 
     ind gives the ascending indices */ 
  int i,j; 
 
  for(i=0;i<p;i++) ind[i]=i; 
  for(i=0;i<p-1;i++) 
    for(j=i+1;j<p;j++) 
      if(x[ind[j]]<x[ind[i]]){ 
	i_swap_j(&ind[i],&ind[j]); 
      } 
} 
 
void sorti (int *x, int *ind, int p, int up) 
{ /* naive sort in ascending order 
     ind gives the ascending indices (up=1) 
     else in descending order */ 
  int i,j; 
 
  for(i=0;i<p;i++) ind[i]=i; 
  if (up==1){ /* ascending order */ 
  for(i=0;i<p-1;i++) 
    for(j=i+1;j<p;j++) 
      if (x[ind[j]]<x[ind[i]]) 
	i_swap_j(&ind[i],&ind[j]); 
  } else { /* descending order */ 
  for(i=0;i<p-1;i++) 
    for(j=i+1;j<p;j++) 
      if (x[ind[j]]>x[ind[i]]) 
	i_swap_j(&ind[i],&ind[j]); 
  } 
} 
 
void sort2 (double *x, int p) 
{ /* sort x ascending 
     and return sorted vector instead of index */ 
  int i,j,*ind; 
  double *y; 
 
  ind = ivector(0,p); 
  y = dvector(0,p); 
  for(i=0;i<p;i++) ind[i]=i; 
  for(i=0;i<p-1;i++) 
    for(j=i+1;j<p;j++) 
      if(x[ind[j]]<x[ind[i]]){ 
	i_swap_j(&ind[i],&ind[j]); 
      } 
  y_assgn_x(y,x,p); 
  for(i=0;i<p;i++) 
    x[i] = y[ ind[i] ]; 
 
  free_ivector(ind,0,p); 
  free_dvector(y,0,p); 
} 
 
void grid (double x0, double xn, int n, double *x) 
{ 
  int i; 
  double dx,xi; 
  dx = (xn-x0)/(n-1.0); 
  for(i=0,xi=x0;i<n;i++,xi+=dx) 
    x[i]=xi; 
} 
 
/* **********************************************************************/ 
/* 		macros for matrix manipulations     			*/ 
/* **********************************************************************/ 
void A_zero(double **A, int p) 
{ 
  for(_i=0;_i<p;_i++) 
    for(_j=0;_j<p;_j++) 
      A[_i][_j] = 0.0; 
} 

double R_max(double **x, int *i, int *j, int p, int q) 
{ 
  double mx; 
  for(mx=x[0][0],*i=*j=0,_i=0; _i<p; _i++){
    for(_j=0; _j<q;_j++){
      if (x[_i][_j]>mx){ 
	mx=x[_i][_j]; 
	*i = _i;  *j=_j;
      }// if >mx
    } //_j
  }// _i
  return(mx); 
} 

double R_min(double **x, int *i, int *j, int p, int q) 
{ 
  double mx; 
  for(mx=x[0][0],*i=*j=0,_i=0; _i<p; _i++){
    for(_j=0; _j<q;_j++){
      if (x[_i][_j]<mx){ 
	mx=x[_i][_j]; 
	*i = _i;  *j=_j;
      }
    } 
  }
  return(mx); 
} 
 
void A_diag_x(double **A, double *x, int p) 
{ 
  for(_i=0;_i<p;_i++){ 
    for(_j=0;_j<p;_j++) 
      A[_i][_j] = 0.0; 
    A[_i][_i] = x[_i]; 
  } 
} 
 
 
void R_zero(double **A, int p,int q) 
{ 
  for(_i=0;_i<p;_i++) 
    for(_j=0;_j<q;_j++) 
      A[_i][_j] = 0.0; 
} 
 
void I_zero(int **A, int p,int q) 
{ 
  for(_i=0;_i<p;_i++) 
    for(_j=0;_j<q;_j++) 
      A[_i][_j] = 0.0; 
} 
 
void A_unit(double **A, int p) 
{ 
  for(_i=0;_i<p;_i++){ 
    for(_j=0;_j<p;_j++) 
      A[_i][_j] = 0.0; 
    A[_i][_i] = 1.0; 
  } 
} 
 
void A_swap_B(double **A, double **B, int p, int q) 
{ 
  double h; 
  for(_i=0;_i<p;_i++) 
    for(_j=0; _j<q;_j++){ 
      h = A[_i][_j]; 
      A[_i][_j] = B[_i][_j]; 
      B[_i][_j] = h; 
    } 
} 
 
 
void Aj_swap_Ak(double **A, int j, int k, int p, int q) 
{ 
  double h; 
  for(_i=0;_i<p;_i++){ 
    h = A[_i][j]; 
    A[_i][j] = A[_i][k]; 
    A[_i][k] = h; 
  } 
} 
 
void At(double **A, double **B, int p) 
{ 
  for(_i=0;_i<p;_i++) 
    for(_j=0;_j<p;_j++) 
      B[_i][_j] = A[_j][_i]; 
} 
 
void Rt(double **A, double **B, int p, int q) 
{ 
  for(_i=0;_i<p;_i++) 
    for(_j=0;_j<q;_j++) 
      B[_i][_j] = A[_j][_i]; 
} 
 
void At_self(double **A, int p) 
{ 
  double x; 
  for(_i=0;_i<p;_i++) 
    for(_j=0;_j<_i;_j++){ 
      x = A[_j][_i]; 
      A[_j][_i] = A[_i][_j]; 
      A[_i][_j] = x; 
    } 
} 
	 
void B_assgn_A(double **B, double **A, int p) 
{ 
  
  for(_i=0;_i<p;_i++) 
    for(_j=0;_j<p;_j++) 
      B[_i][_j] = A[_i][_j]; 
} 
 
void A_assgn_r(double **A, double r, int p) 
{ 
  
  for(_i=0;_i<p;_i++) 
    for(_j=0;_j<p;_j++) 
      A[_i][_j] = r; 
} 
 
void R_assgn_r(double **A, double r, int p, int q) 
{ 
  
  for(_i=0;_i<p;_i++) 
    for(_j=0;_j<q;_j++) 
      A[_i][_j] = r; 
} 
 
void y_assgn_Aj(double *y, double **A, int j, int p) 
{ 
  for(_i=0;_i<p;_i++) 
    y[_i] = A[_i][j]; 
} 
 
void R_assgn_T(double **B, double **A, int p, int q) 
{ 
  
  for(_i=0;_i<p;_i++) 
    for(_j=0;_j<q;_j++) 
      B[_i][_j] = A[_i][_j]; 
} 
void  A_plus_B(double **A, double **B, double **C, int p, int q) 
{ 
  for(_i=0;_i<(p);_i++)  
    for(_j=0;_j<(q);_j++)			 
      C[_i][_j]=A[_i][_j]+B[_i][_j];  
} 
 
void  I_plus_J(int **A, int **B, int **C, int p, int q) 
{ 
  for(_i=0;_i<(p);_i++)  
    for(_j=0;_j<(q);_j++)			 
      C[_i][_j]=A[_i][_j]+B[_i][_j];  
} 
 
void A_minus_B(double **A, double **B, double **C, int p, int q) 
{ 
  for(_i=0;_i<(p);_i++)  
    for(_j=0;_j<(q);_j++)			 
      C[_i][_j]=A[_i][_j]-B[_i][_j];  
} 
 
void AplusBr(double **A, double **B, double r, double **C, int p) 
{ 
  for(_i=0;_i<p;_i++) 
    for(_j=0; _j<p; _j++) 
      C[_i][_j] = A[_i][_j]+r*B[_i][_j]; 
} 
void AplusBi(double **A, double **B, int i, double **C, int p) 
{ 
  for(_i=0;_i<p;_i++) 
    for(_j=0; _j<p; _j++) 
      C[_i][_j] = A[_i][_j]+i*B[_i][_j]; 
} 
void  rA_plus_sB(double r, double **A, double s, double **B, 
		 double **C, int p) 
{ 
  for(_i=0;_i<(p);_i++)  
    for(_j=0;_j<(p);_j++)			 
      C[_i][_j]= r*A[_i][_j]+s*B[_i][_j];  
} 
void  AB(double **A, int p, int q, double **B, int r, int s, double **C) 
{ 
  if (q != r) 
    error("AB","dimensions don't match",1); 
  for(_i=0;_i<(p);_i++)				 
    for(_j=0;_j<(s);_j++)			 
      for(C[_i][_j]=0,_k=0;_k<(q);_k++)	 
	C[_i][_j]+=A[_i][_k]*B[_k][_j]; 
} 
 
void  AtB(double **A, int p, int q, double **B, int r, int s, double **C) 
{ 
  if (p != r) 
    error("AtB","dimensions don't match",1); 
  for(_i=0;_i<(q);_i++)			 
    for(_j=0;_j<(s);_j++)			 
      for(C[_i][_j]=0,_k=0;_k<(p);_k++)	 
	C[_i][_j]+=A[_k][_i]*B[_k][_j]; 
} 
 
void AtBA(double **A, int p, double **B, double **C) 
{ 
  double AtBik; 
 
  for(_i=0;_i<p;_i++) 
    for(_j=0;_j<p;_j++){ 
      C[_i][_j] = 0.0; 
      for(_k=0;_k<p;_k++){ 
	for(_l=0, AtBik=0.0; _l<p; _l++) 
	  AtBik += A[_l][_i]*B[_l][_k]; 
	C[_i][_j] += AtBik*A[_k][_j]; 
      } 
    } 
} 
 
void RtBR(double **A, int p, int q, double **B, double **C) 
{ 
  double AtBik; 
 
  for(_i=0;_i<q;_i++) 
    for(_j=0;_j<q;_j++){ 
      C[_i][_j] = 0.0; 
      for(_k=0;_k<p;_k++){ 
	for(_l=0, AtBik=0.0; _l<p; _l++) 
	  AtBik += A[_l][_i]*B[_l][_k]; 
	C[_i][_j] += AtBik*A[_k][_j]; 
      } 
    } 
} 
 
void  ABt(double **A, int p, int q, double **B, int r, int s, double **C) 
{ 
  if (q != s) 
    error("AtB","dimensions don't match",1); 
  for(_i=0;_i<(p);_i++)			 
    for(_j=0;_j<(r);_j++)			 
      for(C[_i][_j]=0,_k=0;_k<(q);_k++)	 
	C[_i][_j]+=A[_i][_k]*B[_j][_k]; 
} 
 
void rSTU (double r, 
	   double **S, int pS, int qS, 
	   double **T, int pT, int qT, 
	   double **U, int pU, int qU, 
	   double **V) 
{ 
  int i,j,k,l; 
 
  R_zero(V,pS,qU); 
  for(i=0;i<pS;i++) 
    for(j=0;j<qU;j++) 
      for(k=0;k<pT;k++) 
	for(l=0;l<qT;l++) 
	  V[i][j] += r*S[i][k]*T[k][l]*U[l][j]; 
} 
 
void rR_plus_sSTU (double r, double **R, 
		   double s, 
		  double **S, int pS, int qS, 
		  double **T, int pT, int qT, 
		  double **U, int pU, int qU, 
		  double **V) 
{ 
  int i,j,k,l; 
 
  if (r==0) 
    R_zero(V,pS,qU); 
  else 
    rA(r,R,V,pS,qU); 
  for(i=0;i<pS;i++) 
    for(j=0;j<qU;j++) 
      for(k=0;k<pT;k++) 
	for(l=0;l<qT;l++) 
	  V[i][j] += s*S[i][k]*T[k][l]*U[l][j]; 
} 
 
void rx_plus_sSTy (double r, double *x, 
		   double s, 
		  double **S, int pS, int qS, 
		  double **T, int pT, int qT, 
		  double *y, int pU, 
		  double *z) 
{ 
  int i,k,l; 
 
  if (r==0) 
    x_zero(z,pS); 
  else 
    for(i=0;i<pS;i++) 
      z[i] = r*x[i]; 
  for(i=0;i<pS;i++) 
    for(k=0;k<pT;k++) 
      for(l=0;l<qT;l++) 
	z[i] += s*S[i][k]*T[k][l]*y[l]; 
} 
 
void rR_plus_sSTSt (double r, double **R, double s, 
		   double **S, int pS, int qS, 
		   double **T, int pT, int qT, 
		   double **V) 
{ 
  int i,j,k,l; 
 
  if (r==0) 
    R_zero(V,pS,pS); 
  else 
    rA(r,R,V,pS,pS); 
  for(i=0;i<pS;i++) 
    for(j=0;j<pS;j++) 
      for(k=0;k<pT;k++) 
	for(l=0;l<qT;l++) 
	  V[i][j] += s*S[i][k]*T[k][l]*S[j][l]; 
} 
 
void sSTSt (double s, 
		   double **S, int pS, int qS, 
		   double **T, int pT, int qT, 
		   double **V) 
{ 
  int i,j,k,l; 
 
  R_zero(V,pS,pS); 
  for(i=0;i<pS;i++) 
    for(j=0;j<pS;j++) 
      for(k=0;k<pT;k++) 
	for(l=0;l<qT;l++) 
	  V[i][j] += s*S[i][k]*T[k][l]*S[j][l]; 
} 
 
 
double prod_diag(double **A, int p) 
{ 
  double value; 
  for(value=1.0,_i=0;_i<p;_i++) 
    value *= A[_i][_i]; 
  return value; 
} 
 
double trace(double **A, int p) 
{ 
  double value; 
  for(value=1.0,_i=0;_i<p;_i++) 
    value += A[_i][_i]; 
  return value; 
} 
 
/* ----------------------------------------------------------- 
   diag matrix 
   ----------------------------------------------------------- */ 
 
void L_assgn_x(double **L, int p, double *x) 
{ 
  A_zero(L,p); 
  for(_i=0;_i<(p);_i++) 
    L[_i][_i] = x[_i]; 
} 
 
void L_assgn_recx(double **L, int p, double *x) 
{ 
  A_zero(L,p); 
  for(_i=0;_i<(p);_i++) 
    L[_i][_i] = 1.0/x[_i]; 
} 
 
void  AtLB(double **A, int p, int q, double *L, 
	   double **B, int r, int s, double **C) 
{ 
  for(_i=0;_i<(p);_i++)			 
    for(_j=0;_j<(s);_j++)			 
      for(C[_i][_j]=0,_k=0;_k<(q);_k++)	 
	C[_i][_j]+=A[_k][_i]*L[_k]*B[_k][_j]; 
} 
void  ALBt(double **A, int p, int q, double *L, 
	   double **B, int r, int s, double **C) 
{ 
  for(_i=0;_i<(p);_i++)			 
    for(_j=0;_j<(s);_j++)			 
      for(C[_i][_j]=0,_k=0;_k<(q);_k++)	 
	C[_i][_j]+=A[_i][_k]*L[_k]*B[_j][_k]; 
} 
void  LminusHalfAB(double *L,double **A, int p, int q, 
		     double **B, int r, int s, double **C) 
{ 
  for(_i=0;_i<(p);_i++)		 
    for(_j=0;_j<(s);_j++)			 
      for(C[_i][_j]=0,_k=0;_k<(q);_k++)	 
	C[_i][_j]+=A[_i][_k]/sqrt(L[_i])*B[_k][_j]; 
} 
void  LminusHalfAt(double *L,double **A, int p, int q, 
		     double **C) 
{ 
  for(_i=0;_i<(p);_i++){ 
    if (L[_i] <= 0.0) error("LminusHalfAt","non-pos sqrt",1); 
    for(_j=0;_j<(q);_j++)			 
      C[_i][_j] = A[_j][_i]/sqrt(L[_i]); 
  } 
} 
void  AmltDiagB(double **A, double *B, int p, int m) 
{ 
  for(_i=0;_i<(m);_i++) 
    for(_k=0;_k<(p);_k++)				 
      A[_i][_k] *= B[_k];			 
} 
void  AdivDiagB(double **A, double *B, int p, int m) 
{ 
  for(_i=0;_i<(m);_i++)			 
    for(_k=0;_k<(p);_k++)				 
      A[_i][_k] /= B[_k];			 
} 
/* ----------------------------------------------------------- 
   vector & matrix 
   ----------------------------------------------------------- */ 
 
void  xA(double *x,double **A,double *z, int p)	 
{ 
  for(_i=0;_i<(p);_i++){				 
    for(z[_i]=0,_j=0; _j<(p); _j++)		 
      z[_i]+=A[_j][_i]*x[_j];	 
  } 
} 
void  Ax(double **A,double *x,double *z, int p)	 
{ 
  for(_i=0;_i<(p);_i++){				 
    for(z[_i]=0,_j=0; _j<(p); _j++)		 
      z[_i]+=A[_i][_j]*x[_j];	 
  } 
} 
void  Rx(double **R, int p, int q,double *x, double *z)	 
{ 
  for(_i=0;_i<(p);_i++){				 
    for(z[_i]=0,_j=0; _j<(q); _j++)		 
      z[_i]+=R[_i][_j]*x[_j];	 
  } 
} 
void  Ri(double **R, int p, int q, int *x, double *z)	 
{ 
  for(_i=0;_i<(p);_i++){				 
    for(z[_i]=0,_j=0; _j<(q); _j++)		 
      z[_i]+=R[_i][_j]*x[_j];	 
  } 
} 
void  Rtx(double **R, int p, int q,double *x, double *z)	 
{ 
  for(_i=0;_i<(q);_i++){				 
    for(z[_i]=0,_j=0; _j<(p); _j++)		 
      z[_i]+=R[_j][_i]*x[_j];	 
  } 
} 
void  rA(double r,double **A, double **B, int p, int q)	 
{ 
  for(_i=0; _i<(p); _i++){			 
    for(_j=0; _j<(q); _j++)			 
      B[_i][_j] = r* A[_i][_j];	 
  } 
} 
 
void Ax_plus_y(double **A, double *x, double *y, double *z, int p) 
{ 
  
  for(_i=0;_i<p;_i++) 
    for(z[_i]=y[_i],_j=0; _j<p; _j++) 
      z[_i] += A[_i][_j]*x[_j]; 
} 
 
void Rx_plus_y(double **A, double *x, double *y, double *z, int p, int q) 
{ 
  for(_i=0;_i<p;_i++) 
    for(z[_i]=y[_i],_j=0; _j<q; _j++) 
      z[_i] += A[_i][_j]*x[_j]; 
} 
 
void rAx_plus_sBy(double r, double **A, double *x, 
		  double s, double **B, double *y, 
		  double *z, int p) 
{ 
  for(_i=0;_i<p;_i++) 
    for(z[_i]=0,_j=0; _j<p; _j++) 
      z[_i] += r*A[_i][_j]*x[_j] + s*B[_i][_j]*y[_j]; 
} 
double xtAy (double *x, double *y, double **A, int p) 
{ 
  double z; 
 
  for(z=0,_i=0;_i<p;_i++) 
    for(_j=0; _j<p; _j++) 
      z += A[_i][_j]*x[_j]*y[_i]; 
  return(z); 
} 
void ABAt(double **A, double **B, double **C, int p) 
{ 
  int _i, _j; 
  for(_i=0;_i<p;_i++) 
    for(_j=0; _j<p; _j++) 
      C[_i][_j] = xtAy(A[_i],A[_j],B,p); 
} 
void ABAt_plus_C(double **A, double **B, double **C, double **D, int p) 
{ 
  int _i,_j; 
 
  for(_i=0;_i<p;_i++) 
    for(_j=0; _j<p; _j++) 
      { 
	D[_i][_j] = C[_i][_j]+xtAy(A[_i],A[_j],B,p); 
      } 
} 
void RSRt_plus_C(double **A, double **B, double **C, double **D, int 
		 p, int q) 
{ 
  int _i,_j; 
 
  for(_i=0;_i<p;_i++) 
    for(_j=0; _j<p; _j++) 
      { 
	D[_i][_j] = C[_i][_j]+xtAy(A[_i],A[_j],B,q); 
      } 
} 
void C_min_RSRt(double **A, double **B, double **C, double **D, int 
		 p, int q) 
{ 
  int _i,_j; 
 
  for(_i=0;_i<p;_i++) 
    for(_j=0; _j<p; _j++) 
      { 
	D[_i][_j] = C[_i][_j] - xtAy(A[_i],A[_j],B,q); 
      } 
} 
void A_plus_rxxt (double **A, double r, double *x, double **C, int p) 
{ 
  for(_i=0;_i<p;_i++) 
    for(_j=0;_j<p;_j++) 
      C[_i][_j] = A[_i][_j]+r*x[_i]*x[_j]; 
} 
 
void A_xminusy(double **A, double *x, double *y, double *z, int p) 
{ 
  
  for(_i=0;_i<p;_i++) 
    for(z[_i]=0,_j=0; _j<p; _j++) 
      z[_i] += A[_i][_j]*(x[_j]-y[_j]); 
} 
 
void A_xminusy_plus_z(double **A, double *x, double *y, double *z, 
		      double *w, int p) 
{ 
  for(_i=0;_i<p;_i++) 
    for(w[_i]=z[_i],_j=0; _j<p; _j++) 
      w[_i] += A[_i][_j]*(x[_j]-y[_j]); 
} 
 
 
/* ----------------------------------------------------------- 
   submatrices/vecyors 
   ----------------------------------------------------------- */ 
void  Ax_j(double **A,double *x,double *z, int skip, int p) 
{ 
  for(_i=0;_i<skip;_i++){ 
    for(z[_i]=0,_j=0; _j<(skip); _j++) 
      z[_i]+=A[_i][_j]*x[_j];       
    for(_j=(skip+1); _j<(p); _j++)     
      z[_i]+=A[_i][_j]*x[_j];}      
  for(_i=skip+1;_i<p;_i++){          
    for(z[_i]=0,_j=0; _j<(skip); _j++) 
      z[_i]+=A[_i][_j]*x[_j];       
    for(_j = (skip+1); _j<(p); _j++)   
      z[_i]+=A[_i][_j]*x[_j];}  
}  
void  xx_j(double *x, double *r, int skip, int p) 
{ 
  for(_i=0, *r=0;_i<skip;_i++)  
    *r += x[_i]*x[_i];         
  for(_i=skip+1; _i<p; _i++)              
    *r += x[_i]*x[_i]; 
}    
 
void dsubmatrix(double **a, int r, int c, 
		   double **m, int rnew, int cnew, int *inew, int *jnew) 
/* returns a submatrix with rnew rows inew[0..rnew-1] and 
   cnew cols jnew[0..cnew-1] */ 
{ 
  for(_i=0; _i<rnew; _i++) 
    for(_j=0; _j<cnew; _j++) 
      m[_i][_j] = a[inew[_i]][jnew[_j]]; 
} 
 
void dsubmatrix_j(double **a, int rc, double **m, int p) 
/* returns matrix a without rc-th row and column */ 
{ 
  int i,i2,j,j2; 
 
  for(i=0,i2=0;i<p;i++){ 
    if (i==rc) continue; 
    for(j=j2=0;j<p;j++){ 
      if (j==rc) continue; 
      m[i2][j2] = a[i][j]; 
      j2 += 1; 
    } 
    i2 += 1; 
  } 
} 
 
 
void dsubvector(double *a, int r, 
		   double *m, int rnew, int *inew) 
/* returns a submatrix with rnew rows inew[0..rnew-1] and 
   cnew cols jnew[0..cnew-1] */ 
{ 
  for(_i=0; _i<rnew; _i++) 
      m[_i] = a[inew[_i]]; 
} 
 
void permutate (int *y, int n) 
{ 
  int N,i; 
 
  for(i=0;i<n;i++) 
    y[i] = i; 
 
  for(N=n; N>1; N--){ 
    i = floor(runif()*N); 	/* random (not yet used) index */ 
    i_swap_j(&y[i],&y[N-1]); 	/* swap with last one */ 
  } 
}   
 
   
