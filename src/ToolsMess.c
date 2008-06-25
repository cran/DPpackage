static const char mess_c_sccs_id[] = "%W%";

/*+
 * mess.c - Random Auxillary C Functions
 *
 * Version:	%I%
 *
 * Description:
 *		Set of miscellaneous functions concerned with printing messages.
 *
 * Modification History:
 * 20-NOV-1998	1.1		Peter Muller
 *		Created.
 *
 * 25-JUN-1999	1.2		Stephen Morris (Tessella Support services plc)
 *		1) Removed nested comments (The C slash-star comment introducer
 *		within a comment block).  The Gnu C compiler issued warnings on
 *		these.
 *		2) Added "stdlib.h" for a definition of "exit()".
 *		3) Changed return type of functions not returning a value to void.
 *		4) Remove unreferenced variables.
 *		5) Include matrix.h and vector.h to define some referenced functions.
-*/

#include <stdio.h>
#include <stdlib.h>

#include "ToolsCMatrix.h"
#include "ToolsMess.h"
#include "ToolsNrutil.h"
#include "ToolsVector.h"

/* static variables */

static int fatal = 0;

/********************************************
 *         warning
 ********************************************/

void warning(char *module, char *mess, int nr)              
{
   	printf("\n*** WARNING # %d in %s***\n %s \n",
	       nr,module,mess);
}

/********************************************
 *         error
 ********************************************/   

void error(char *module, char *mess, int nr)              
{
   	printf("\n *** ERROR # %d in %s***\n %s",nr,module, mess);
	printf(  " exiting program \n");
	exit(1);
}

void error3(char *mess)
{
  printf("\n *** ERROR %s \n",mess);
  exit(1);
}

void err_msg(char *fct, char *txt, int n1, int n2, int n3)
{ /* print error message */
  printf("\n\n *** Error in %s \n", fct);
  printf(txt,n1,n2,n3); /* n1,n2 allows to include
			numbers in txt */
  printf("\n");
  exit(1);
}

/* ************************************* err     ************** */     
int	err()
{
	return(fatal);
}

/*
 * ********* various functions to print debug messages ********
 * to turn off these debug messages, replace mess.. by dummies.
 */

/********************************************
 *         message
 ********************************************/

void messblk(char *txt) 	/* prints message on stderr (for debugging) */           
{
   	printf(" %s ",txt);
}

void message(char *txt) 	/* prints message on stderr (for debugging) */           
{
   	printf("%s\n",txt);
}

void messtab(float x)              
{
   	printf("%6.4f\t",x);
}

void header(char *txt)	/* prints header with text */
{
	printf
	("\n******************************************************************\n");
	printf
	("*               %-40s         *\n", txt);
	printf
	("******************************************************************\n");
}

void header2(char *txt)	/* prints header with text */
{
	printf
	("\n================= ");
	printf
	(" %-30s ", txt);
	printf
	("===================\n");
}

void header3()	/* prints header with text */
{
	printf
	("\n=============================================\n ");
}

/********************************************
 *         messint
 ********************************************/

void messintblk(char *txt, int i)              
{
   	printf("%s %d ",txt,i);
}

void messinttab(char *txt, int i)              
{
   	printf("%s %d\t",txt,i);
}

void messint(char *txt, int i)              
{
   	printf("%s %d\n",txt,i);
}

void messintvec(char *txt, int *i, int p)
{
  int j;
  printf("\n %s ", txt);
  for(j=0;j<p;j++){
    printf("%4d ",i[j]);
  }
  printf("\n");
}

void messintmatrix2(char *txt, int **x, int p, int q)              
{
	int 	i,j;

	printf("%s\n",txt);
	for(i=0;i<p;i++){
		for(j=0;j<q;j++)
			printf("%d ",x[i][j]);
		printf("\n");
	}
}

/*********************************************
 *         messFloat
 *********************************************/

void messfloatblk(char *txt, float x)              
{
   	printf("%s %8.4e ",txt,x);
}

void messfloat(char *txt, float x)              
{
  printf("%s %8.4e\n",txt,x);
}

void messfloatvec(char *txt, float *x,  int p)              
{
  int 	i;
  
  printf("%s\n",txt);
  for(i=0;i<p;i++)
    printf("%6.3f ",x[i]);
  printf("\n");
}

void messfloattab(char *txt, float x)              
{
  printf("%s %8.4e\t",txt,x);
}

void messmatrix(char *txt, float *x, int p, int q)              
{
	int 	i,j;

	printf("%s\n",txt);
	for(i=0;i<p;i++){
		for(j=0;j<q;j++)
			printf("%8.4e ",x[i*q+j]);
		printf("\n");
	}
}

/********************************************
 *         messDouble
 ********************************************/

void messdouble(char *txt, double x)              
{
  printf("%s %6.2f\n",txt,x);
}

void messtabdouble(char *txt, double x)              
{
  printf("%s %6.4f\t ",txt,x);
}

void messblkdouble(char *txt, double x)              
{
  printf("%s %6.4f ",txt,x);
}

void messdoublematrix(char *txt, double *x, int p, int q)              
{
	int 	i,j;

	printf("%s\n",txt);
	for(i=0;i<p;i++){
		for(j=0;j<q;j++)
			printf("%6.4f ",x[i*q+j]);
		printf("\n");
	}
}

void messdoublematrix2(char *txt, double **x, int p, int q)              
{
	int 	i,j;

	printf("%s\n",txt);
	for(i=0;i<p;i++){
		for(j=0;j<q;j++)
			printf("%6.4f ",x[i][j]);
		printf("\n");
	}
}

void messdoublematrix3(char *txt, double ***x, int r, int p, int q)              
{
	int 	k,i,j;

	printf("%s\n",txt);
	for(k=0;k<r;k++){
	  printf("\n ------- %2d ---------- \n", k);
	  for(i=0;i<p;i++){
	    for(j=0;j<q;j++)
	      printf("%5.2f ",x[k][i][j]);
	    printf("\n");
	  }
	}
}

void messdoublevec(char *txt, double *x,  int p)
{
  int 	i;
  
  printf("%s: ",txt);
  for(i=0;i<p;i++)
    printf(" %5.2f ",x[i]);
  printf("\n");
}

void messdoublevectab(char *txt, double *x,  int p)              
{
  int 	i;
  
  printf("%s: ",txt);
  for(i=0;i<p;i++)
    printf(" %5.2f ",x[i]);
  printf("\t");
}

/* -----------------------------------------------------------
   for debugging
   ---------------------------------------------------------- */

void dbl(double x)
{
  messdouble(" ",x);
}

void dvec(double *x, int p)
{
  messdoublevec(" ",x,p);
}

void dvecl(double *x, int p)
{
  int
    i;
  for(i=0;i<p;i++)
    printf("%7.5f ",x[i]);
  printf("\n");
}

void ivec(int *x, int p)
{
  messintvec(" ",x,p);
}

void dmat(double *x, int p, int q)
{
  messdoublematrix(" ",x,p,q);
}

void dmat2(double **x, int p, int q)
{
  messdoublematrix2(" ",x,p,q);
}

void dmat2l(double **x, int p, int q)
{
  int
    i,j;
  for(i=0;i<p;i++){
    for(j=0;j<q;j++)
      printf("%8.6f ",x[i][j]);
    printf("\n");
  }
}

void dmat2s(double **x, int p, int q)
{
  int
    i,j;
  for(i=0;i<p;i++){
    for(j=0;j<q;j++)
      printf("% 4.2f ",x[i][j]);
    printf("\n");
  }
}

void imat2(int **x, int p, int q)
{
  int
    i,j;
  for(i=0;i<p;i++){
    for(j=0;j<q;j++)
      printf("% 4d ",x[i][j]);
    printf("\n");
  }
}

void dmat3(double ***x, int r,int p, int q)
{
  messdoublematrix3(" ",x,r,p,q);
}

void dR(double **V, int p)
{
  double **S;
  int i;

  S = dmatrix(0,p,0,p);
  
  ABt(V,p,p,V,p,p,S);
  for(i=0;i<p;i++)
    printf("% 4.2f ",S[i][i]);
  printf("\n");
  corr(S,S,p,1);
  
  free_dmatrix(S,0,p,0,p);
}

void dAB(double **A, double **B, int p)
{
  double **C;

  C = dmatrix(0,p,0,p);
  AB(A,p,p,B,p,p,C);
  dmat2(C,p,p);

  free_dmatrix(C,0,p,0,p);
}

void dAtB(double **A, double **B, int p)
{
  double **C;

  C = dmatrix(0,p,0,p);
  AtB(A,p,p,B,p,p,C);
  dmat2(C,p,p);

  free_dmatrix(C,0,p,0,p);
}

void dABt(double **A, double **B, int p)
{
  double **C;

  C = dmatrix(0,p,0,p);
  ABt(A,p,p,B,p,p,C);
  dmat2(C,p,p);

  free_dmatrix(C,0,p,0,p);
}

void dA_inv(double **A, int p)
{
  double **C;

  C = dmatrix(0,p,0,p);
  d_invert(p,A,C);
  dmat2(C,p,p);

  free_dmatrix(C,0,p,0,p);
}

void dVup_inv(double **A, int p)
{
  double **C;

  C = dmatrix(0,p,0,p);
  d_inv_triang(p,A,C);
  dmat2(C,p,p);

  free_dmatrix(C,0,p,0,p);
}

void dVlo_inv(double **A, int p)
{
  double **C;

  C = dmatrix(0,p,0,p);
  d_inv_triang_lo(p,A,C);
  dmat2(C,p,p);

  free_dmatrix(C,0,p,0,p);
}

void dA_invB(double **A, double **B, int p)
{
  double **C,**D;

  C = dmatrix(0,p,0,p);
  D = dmatrix(0,p,0,p);

  d_invert(p,A,C);
  AB(C,p,p,B,p,p,D);
  dmat2(D,p,p);

  free_dmatrix(C,0,p,0,p);
  free_dmatrix(D,0,p,0,p);
}

void dABAt(double **A, double **B, int p)
{
  double **C,**D;

  C = dmatrix(0,p,0,p);
  D = dmatrix(0,p,0,p);

  AB (A,p,p,B,p,p,C);
  ABt(C,p,p,A,p,p,D);
  dmat2(D,p,p);

  free_dmatrix(C,0,p,0,p);
  free_dmatrix(D,0,p,0,p);
}

void dA_invA_invt(double **A, int p)
{
  double **C,**D;

  C = dmatrix(0,p,0,p);
  D = dmatrix(0,p,0,p);

  d_invert(p,A,C);
  ABt(C,p,p,C,p,p,D);
  dmat2(D,p,p);

  free_dmatrix(C,0,p,0,p);
  free_dmatrix(D,0,p,0,p);
}


void drA(double r, double **A, int p)
{
  double **C;

  C = dmatrix(0,p,0,p);

  rA(r,A,C,p,p);
  dmat2(C,p,p);

  free_dmatrix(C,0,p,0,p);
}
