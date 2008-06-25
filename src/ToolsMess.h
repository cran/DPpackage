#if !defined(MESS_H_INCLUDED)
#define MESS_H_INCLUDED

/*+
 * mess.h - Message Routines
 *
 * Description:
 *		defines a number of routines for printing information to stdout.
 *
 * Modification History:
 * <unknown>	1.1		Peter Muller
 *		Created.
 *
 * 25-JUN-1999	1.2		Stephen Morris (Tessella Support Services plc)
 *		Included #if....#endif sentinels.
 *		Removed nested comments.
 *		Define functions that do not take a value to be void.
-*/


void error(char *,char *, int);
void warning(char *,char *, int);
void error2(char *, char *, int);
void error3(char *);
void err_msg(char *fct, char *txt, int n1, int n2, int n3);
void messblk(char *);
void message(char *); 	
void messtab(float);              
void header(char *);	/* prints header with text */
void header2(char *);	/* prints header with text */
void header3();	/* prints header with text */

/********************************************
 *         messint
 ********************************************/
void messintblk(char *, int);              
void messinttab(char *,  int);              
void messint(char *,  int);              
void messintvec(char *, int *, int);              
/********************************************
 *         messFloat
 ********************************************/
void messfloatblk(char *, float);              
void messfloat(char *, float);              
void messfloatvec(char *, float*,  int);              
void messfloattab(char *,  float);              
void messmatrix(char *, float*,  int, int);              
/********************************************
 *         messDouble
 ********************************************/
void messdouble(char *, double);              
void messtabdouble(char *, double);              
void messblkdouble(char *, double);              
void messdoublematrix(char *, double *,  int, int);              
void messdoublematrix2(char *, double **, int, int);
void messdoublevec(char *, double *,  int);              

void dbl(double );
void dvec(double *,int );
void ivec(int *,int );
void dmat(double *,int, int);
void dmat2(double **, int, int);
void dmat2s(double **, int, int);
void dmat3(double ***, int, int, int);

void dABAt(double **A, double **B, int p);

#endif





