#if !defined(INTERFACE_H_INCLUDED)
#define INTERFACE_H_INCLUDED

/*+
 * interface.h - Miscellaneous Useful Routines
 *
 * Description:
 *		Declaration of various routines in the library, primarily
 *		concerned with I/O.
 *
 * Modification History:
 * <unknown>	1.1		Peter Mueller
 *		Created.
 *
 * 29-JUN-1999	1.2		Stephen Morris (Tessella Support Services plc)
 *		Remove embedded comments (slash-star in a comment block).
 *		Added header and #if...#endif sentinels.
-*/

#include <stdio.h>



void swap(int *i, int *j);
void swapDoubleArray(double **x, double **y);
void swapDoubleMatrix(double ***x, double ***y);

void trim(char *in, char *out);
void trim2(char *in, char *out, char *excl, int nexcl);

/* ***********************************************************
   open/close files
   *********************************************************** */

void openInterface_toC();
void closeInterface_toC();
void openInterface_fromC();
void openInterface_from(char *);
FILE *openIn(char *);
void closeIn();
void closeOut();
FILE *openOut(char *);
FILE *openAppend(char *);
void openInterface_fromAppend(char *);
void closeInterface_fromC();

/* ***********************************************************
   read scalars
   *********************************************************** */

void scanFloat(char *, float *);
void scanDouble(char *, double *);
void scanInt(char *, int *);
void fscanDouble(FILE *,char *, double *);
void fscanInt(FILE *,char *, int *);
void scanLong(char *, long *);

/* ***********************************************************
   read arrays
   *********************************************************** */

void scanFloatArray(char *,float *, int);
void scanArray(char *,float *, int);
void scanDoubleArray(char *,double *, int);
void scanString(char *txt, char *s, int n);
void fscanString(FILE *,char *txt, char *s, int n);
void fscanDoubleArray(FILE *,double *, int);
void scanDoubleMatrix(char *, double **, int, int);
void fscanDoubleMatrix(FILE *ifile, double **x,int r,int c);
void scanIntArray(char *,int *, int);
void fscanIntArray(FILE *ifile, int *x, int n);

/* ***********************************************************
   write scalar
   *********************************************************** */

void writeInt(int);
void writeLong(long i);
void writeFloat(float);
void writeDouble(double);

/* ***********************************************************
   write array
   *********************************************************** */

void writeIntArray(int *,int, int);
void fwriteIntArray(FILE *, int *,int, int);
void fwriteIntMatrix(FILE *f, int **x, int rows, int cols);
void writeIntMatrix(int **x, int rows, int cols);
void writeDoubleArray(double *,int, int);
void writeDoubleMatrix2(double **, int , int);
void fwriteDoubleArray(FILE *, double *,int, int);
void fwriteDoubleMatrix2(FILE *, double **, int , int);
void writeDoubleMatrix(double **, int, int);
void writeFloatArray(float *, int, int);
void writeArray(float *, int, int); 

void fserror(char *proc, char *act, char *what);

#endif
