static const char matric_c_sccs_id[] = "%W%";

/*+
 * matrix.c - Matrix Routines
 *
 * Version:	%I%
 *
 * Description:
 *		Holds a number of functions needed to manipulate matrices.
 *
 *		N.B. all these routines assume that matrices are passed as
 *		arrays of arrays.
 *
 *		(Using lapack routines are compiled but not yet tested !!!
 *		Test with test-matrix)
 *
 * Modification History:
 * <unknown>	1.1		Peter Muller
 *		Created.
 *
 * 25-JUN-1999	1.2		Stephen Morris (Tessella Support services plc)
 *		1) Removed nested comments (The C slash-star comment introducer
 *		within a comment block).  The Gnu C compiler issued warnings on
 *		these.
 *		2) Added "fortif.h" for a definition of the Lapack FORTRAN routines;
 *		inclusion of this file eliminates the need to suffix the names
 *		of FORTRAN routines with an underscore.
 *		3) Add "stdlib.h" for a definition of exit().
 *		4) cdinv and invcd are now "void", as they do not return a value.
 *		   The result of invtr is now passed back to the caller of
 *		   d_inv_triang, and d_inv_triang_lo.
-*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ToolsMess.h"
#include "ToolsRand.h"
#include "ToolsNrutil.h"
#include "ToolsVector.h"
#include "ToolsCMatrix.h"
#include "ToolsFortif.h"

/* here is some dummmy */
int seed;


/********************************************
 *         eigen
 ********************************************/
int	d_eigen(int n,double **a_in, double *lambda, double **e_out)
/* returns vector of eigenvalues, lambda
   and unit eigenvectors as rows of E */
/*!!!!! WARNING: might have a bug, compare with Splus eigenvalues */
{ 
  int 
    info, i, j, lwork;
  double 
    *a, *work;
  char
    jobz, uplo;

  a = dvector(0,n*n);
  work = dvector(0,6*n);
  lwork = 6*n;

  /* copy 2-dim matrix to 1-dim vector a 
     by row */
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      a[i*n+j] = a_in[i][j];

  /* set up pars for and do IMSL call */
  jobz = 'V';
  uplo = 'U';
  dsyev(&jobz,&uplo,&n,a,&n,lambda,work,&lwork,&info); 
  /* devcsf(&n,a,&lda,lambda,e,&lde); */

  /* dlftds returns matrix of eigenvectors, but by 
     columns!! */
  if (info > 0){
    lp_error("d_eigen","compute eigen values","dsyev",info);
  } 
  
  /* copy e into 2-dim matrix e_out */
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      e_out[i][j] = a[i*n+j];
  
  /* release memory */
  free_dvector(a, 0, n*n);
  free_dvector(work,0,6*n);

  return(info);
}

/********************************************
 *         chol-decomp
 ********************************************/
int	d_chol_decomp(int n, double **a_in, double **v_out)
/* returns in v the chol-decomp of a, i.e. V'V=A
   and V is upper triangular */
{ 

  /* dummy prototype */

  int 
    info, i, j;
  double 
    *a;
  char
    uplo[2];

  a = dvector(0,n*n);

  /* copy 2-dim matrix to 1-dim vector a 
     by row */
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      a[i*n+j] = a_in[i][j];
  uplo[0] = 'U';
  uplo[1] = 'P';

  /* Set up pars and do lapack call */
  dpotrf(uplo, &n, a, &n, &info);

  /* dlftds returns upper triangular, but by 
     columns!! (FORTRAN) */
  for(i=0;i<n;i++){
    v_out[i][i] = a[i*n+i];
    for(j=0;j<i;j++){
      v_out[j][i] = a[i*n+j];
      v_out[i][j] = 0.0;
    }
  }
  
  free_dvector(a, 0,n*n);

  if (info != 0){
    printf("\n****** spurious chol decomp \n");
  }
  if (info != 0){
    lp_error("d_chol_decomp","compute choleski decomp","dpotrf",info);
  } 
  return(info);
}

void cd(double **A, double **V, int p)
/* returns lower triangular cd, i.e. VV'=A, V lower triang */
{
  int i,j;
  double x;

  d_chol_decomp(p,A,V); /* returns V st V'V = A, V upper triang */
  for(i=0;i<p;i++)
    for(j=(i+1);j<p;j++){
      x = V[i][j];
      V[i][j] = V[j][i];
      V[j][i] = x;
    }
}
  
/* ***********************************************************
   Chol-decomp-inv
 * *********************************************************** */
/* Returns L, if 1/A = LL', i.e. L=1/V where A=V'V
 * i.e. L = cd(inv(A)) = A-inv-cd */

void cdinv(int p, double **A, double **L)
{
  double **V;

  V = dmatrix(0,p-1,0,p-1);
  d_chol_decomp(p,A,V); /* A = V'V */
  d_inv_triang(p,V,L);  /* L = 1/V */
  free_dmatrix(V,0,p-1,0,p-1);
}
/* ***********************************************************
   Inv-Chol-decomp
   *********************************************************** */
/* Returns L=1/V, if A=VV', i.e. 1/A = L'L
   i.e. L = inv(cd(A)) = A-cd-inv */

void invcd(int p, double **A, double **L)
{
  double **A_inv;

  A_inv = dmatrix(0,p-1,0,p-1);
  d_invert(p,A,A_inv);
  d_chol_decomp(p,A_inv,L);

  free_dmatrix(A_inv,0,p-1,0,p-1);
}

/******************************************** 
 *         invert
 ********************************************/
int	d_invert(int n, double **a_in, double **v_out)
/* returns in v the inverse of a, i.e. v = a^-1 
   for a real symmetric pos def */
{ 
  int 
    infoc, info, i, j;
  double 
    *a;
  char
    uplo='U';

  a = dvector(0,n*n+1);

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      a[i*n+j] = a_in[i][j];

  dpotrf(&uplo, &n, a, &n, &infoc); /* chol decomp */
  dpotri(&uplo,&n,a,&n,&info);
  /*   dlinds(&n,a,&lda,v,&ldv); */

  if (infoc != 0)
    lp_error("d_invert","invert a matrix","dpotrf",infoc);
  if (info != 0)
    lp_error("d_invert","invert a matrix","dpotri",info);
      
  for(i=0;i<n;i++){
    for(j=0;j<=i;j++){
      v_out[i][j] = a[i*n+j];
      v_out[j][i] = a[i*n+j];
    }
  }

  free_dvector(a, 0, n*n+1);
  if ((info != 0)|(infoc != 0))
    printf("\n Warning: spurious matrix inversion.\n");
  return(info);
}

/********************************************
 *         inv_triang
 ********************************************/
int	d_inv_triang(int n, double **a_in, double **v_out)
/* returns in v the inverse of the upper t
   triangular matrix a */
{ 
  return invtr(n,a_in,v_out,'L');
}

int d_inv_triang_lo(int n, double **a_in, double **v_out)
/* returns in v the inverse of the LOWER triang matrix a */
{
  return invtr(n,a_in,v_out,'U');
}

int invtr(int n, double **a_in, double **v_out, char uplo)
/* returns in v the inverse of the triangular matrix a 
   note: uplo='L' means high, uplo='H' means low 'cos of
   Fortran/C conversion */
{
  int 
    info, i, j;
  double 
    *a;
  char
    diag='N';

  a = dvector(0,n*n);

  /* copy 2-dim matrix to 1-dim vector a 
     by row */
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      a[i*n+j] = a_in[i][j];

  /* set up pars for and do IMSL call */
  /* dlinrt(&n,a,&lda,&ipath,v,&ldv); */
  dtrtri(&uplo,&diag,&n,a,&n,&info);
  if (info != 0)
    lp_error("d_invtr","invert a triang matrix","dpotri",info);

  /* convert vector v into 2-dim matrix v_out */
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      v_out[i][j] = a[i*n+j];

  free_dvector(a, 0, n*n);
  return(info);
}

/* *************************************************
   corr
 ************************************************* */

void corr(double **S, double **R, int p, int prt)
{
  int 
    i,j;
  double 
    *s;

  s = dvector(0,p);

  for(i=0;i<p;i++) 
    s[i] = sqrt(S[i][i]);
  for(i=0;i<p;i++)
    for(j=0;j<p;j++)
      R[i][j] = S[i][j]/(s[i]*s[j]);
  
  if (prt==1){
    messdoublematrix2("R ",R,p,p);
    messdoublevec("s ",s,p);
  }
  free_dvector(s,0,p);
}
  

/* *************************************************
   error
 ************************************************* */
int _mt_exit=1;

void matrix_noexit()
{
  /* prints warning only, without exiting.. */
  _mt_exit=0;
}

void lp_error(char *proc, char *act, char *lpfct, int info)
{	
  /* error in a call to a lapack function */

  fprintf(stderr, "\n ** Error ");
  if (proc[0]!='\0') /* not empty */
    fprintf(stderr, " in function '%s', ", proc);
  if (act[0]!='\0') /* not empty */
    fprintf(stderr, " trying to  %s", act);
  if (lpfct[0]!='\0') /* not empty */
    fprintf(stderr, "\n ** using the LAPACK function '%s' ", lpfct);
  else
    fprintf(stderr, ", ");
  fprintf(stderr, "(info=%d),", info);
  fprintf(stderr, "\n ** .. exiting program");
  fprintf(stderr, " (from a function in 'matrix.c').\n");
  if (_mt_exit==1)
    exit(1);
}
