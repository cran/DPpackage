static const char interface_c_sccs_id[] = "%W%";

/*+
 * interface.c - Interface Routines
 *
 * Version:	%I%
 *
 * Description:
 *		Holds a number of functions needed to access files.
 *
 * Modification History:
 * <unknown>	1.1		Peter Muller
 *		Created.
 *
 * 25-JUN-1999	1.2		Stephen Morris (Tessella Support services plc)
 *		1) Removed nested comments (The C slash-star comment introducer
 *		within a comment block).  The Gnu C compiler issued warnings on
 *		these.
 *		2) Added "stdlib.h" for a definition of "exit()".
 *		3) Removed unreferenced local variables.
 *		4) In "writeLong", use "%ld" for output of a long variable.
 *		5) Initialize s1 in fwrite[Int|Double]Array, to avoid warnings about
 *		possible uninitialized use.
-*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ToolsNrutil.h"
#include "ToolsInterface.h"

FILE	*ifile, *fopen(), *ofile;

/* -------------------------------------------------
   swap
 ------------------------------------------------- */
void swap(int *i, int *j)
{
  int k;
  k= *i; *i= *j; *j= k;
}
void swapDoubleArray(double **x, double **y)
{
  double *z;
  z = *x; *x = *y; *y = z;
}
void swapDoubleMatrix(double ***x, double ***y)
{
  double **z;
  z = *x; *x = *y; *y = z;
}

/* -----------  trim comments from a file  ---------------------- */ 
void trim(char *in, char *out){  
  /* trims away comments from "init.upd" 
     discards comments from any "#" to "\n" */
  char c; 
  FILE *f,*g;

  f = openIn(in);
  g = openOut(out);
  for(;feof(f)==0;){
    fscanf(f,"%c",&c);
    if (c=='#') /* read input until newline and discard */
      for(;c!='\n';fscanf(f,"%c",&c));
    fprintf(g,"%c",c);
  }
  fclose(g);
  fclose(f); 
}

void trim2(char *in, char *out, char *excl, int nexcl){
  /* trims away comments from "init.upd" 
     discards comments from any "#" to "\n" 
     excludes chars in excl[1..nexcl] */
  char c; 
  int  skp, j;
  FILE *f,*g;

  f = openIn(in);
  g = openOut(out);
  for(;feof(f)==0;){
    fscanf(f,"%c",&c);
    if (c=='#') /* read input until newline and discard */
      for(;c!='\n';fscanf(f,"%c",&c));
    for(j=skp=0;j<nexcl;j++)
      if (c==excl[j])
	skp=1;
    if (skp==0)
      fprintf(g,"%c",c);
  }
  fclose(g);
  fclose(f); 
}
  
/* -----------------------------------------------------------
   open/close files
   ----------------------------------------------------------- */

void openInterface_toC()
{
	if ( (ifile=fopen("toC","r") ) == NULL) {
		fserror("","open file","toC");
	}
}

void closeInterface_toC()
{
	fclose(ifile);
}


void openInterface_fromC()
{
	if ( (ofile=fopen("fromC","w") ) == NULL) {
		fserror("","open file","fromC");
	}
}

void openInterface_from(char *name)
{
	if ( (ofile=fopen(name,"w") ) == NULL) {
		fserror("","open file",name); 
	}
}

FILE *openIn(char *name)
{
  if ((ifile=fopen(name,"r"))==NULL){
    fserror("openIn","open file",name);
  }
  return(ifile);
}

void closeIn()
{
  fclose(ifile);
}

void closeOut()
{
  fclose(ofile);
}

FILE *openOut(char *name)
{
  if ((ofile=fopen(name,"w"))==NULL){
    fserror("openOut","open file",name);
  }
  return ofile;
}

FILE *openAppend(char *name)
{
  if ((ofile=fopen(name,"a"))==NULL){
    fserror("openAppend","open file",name);
  }
  return ofile;
}


void openInterface_fromAppend(char *name)
{
  if ((ofile=fopen(name,"a"))==NULL){
    fserror("openInterface_fromAppend","open file",name);
  }
}

void closeInterface_fromC()
{
	fclose(ofile);
}

/* -----------------------------------------------------------
   read in
   ----------------------------------------------------------- */

void scanFloat(char *txt, float *f)
{
  fscanf(ifile,txt);
  if (fscanf(ifile," %f ",f) != 1) {
    fserror ("scanFloat","read float",txt);
  }
}

void scanDouble(char *txt,double *f)
{
  fscanf(ifile,txt);
  if (fscanf(ifile," %lf ",f) != 1) {
    fserror ("scanDouble","read double",txt);
  }
}

void fscanDouble(FILE *ifile, char *txt, double *f)
{
  fscanf(ifile,txt);
  if (fscanf(ifile," %lf ",f) != 1) {
    fserror ("fscanDouble","read double",txt);
  }
}

void scanInt(char *txt, int *n)
{
  int 	s;

  fscanf(ifile,txt);
  if ((s = fscanf(ifile," %d ",n)) != 1) {
    fserror ("scanInt","read int",txt);
  }
}

void fscanInt(FILE *ifile, char *txt, int *n)
{
	int 	s;

  fscanf(ifile,txt);
  if ((s = fscanf(ifile," %d ",n)) != 1) {
    fserror ("fscanInt","read int",txt);
  }
}


void scanLong(char *txt, long *n)
{
        int     s;

  fscanf(ifile,txt);
  if ((s = fscanf(ifile," %ld ",n)) != 1) {
    fserror ("scanLong","read long",txt);
  }
}

/* -----------------------------------------------------------
   read arrays
   ----------------------------------------------------------- */
void scanFloatArray(char *txt, float *x, int n)
{
	scanArray(txt,x,n);
}

void scanArray(char *txt, float *x, int n)
{
	int	i; 

  fscanf(ifile,txt);
  for(i=0;i<n;i++){
    if (fscanf(ifile," %f ",&x[i]) != 1) {
      fserror ("scanArray","read float array",txt);
    }
  }
}

void scanDoubleArray(char *txt, double *x, int n)
{
  int	i;

  fscanf(ifile,txt);
  for(i=0;i<n;i++){
    if (fscanf(ifile," %lg ",&x[i]) != 1) {
      fserror ("scanDoubleArray",
	       "read double array",txt);
    }
  }
}

void fscanDoubleArray(FILE *in, double *x, int n)
{
  int	i;

  for(i=0;i<n;i++){
    if (fscanf(in," %lg ",&x[i]) != 1) {
      /* printf("i=%d\n",i); */
      fserror("fscanDoubleArray","read double array","");
    }
  }
}

void scanString(char *txt, char *s, int n)
{
  fgets(s,n,ifile);
}
  

void fscanString(FILE *ifile, char *txt, char *s, int n)
{
  fgets(s,n,ifile);
}
  
void scanDoubleMatrix(char *txt,double **x,int r,int c)
{
  int	i,j;
  
  fscanf(ifile,txt);
  for(i=0;i<r;i++)
    for(j=0;j<c;j++){
      if (fscanf(ifile," %lg ",&x[i][j]) != 1) {
	fserror ("scanDoubleMatrix","read double matrix",txt);
      }
    }
}

void fscanDoubleMatrix(FILE *ifile, double **x,int r,int c)
{
  int	i,j;
  
  for(i=0;i<r;i++)
    for(j=0;j<c;j++){
      if (fscanf(ifile," %lg ",&x[i][j]) != 1) {
	exit(1);
      }
    }
}

void scanIntArray(char *txt, int *x, int n)
{
  int	i;

  fscanf(ifile,txt);
  for(i=0;i<n;i++){
    if (fscanf(ifile," %d ",&x[i]) != 1) {
      fserror ("scanIntArray","read int array",txt);
    }
  }
}

void fscanIntArray(FILE *ifile, int *x, int n)
{
  int	i;

  for(i=0;i<n;i++){
    if (fscanf(ifile," %d ",&x[i]) != 1) {
      fserror("fscanIntArray","read int array","");
    }
  }
}

/* -----------------------------------------------------------
   write scalars
   ----------------------------------------------------------- */
void writeInt(int i)
{
  int s;
  s=fprintf(ofile,"%d\n",i);
  if(s<0)
    fserror("writeInt","write int","");
}
void writeLong(long i)
{
  int s;
  s=fprintf(ofile,"%ld\n",i);
  if (s<0)
    fserror("writeLong","write long","");
    
}

void writeIntTab(int i)
{
  int s;
  s=fprintf(ofile,"%d\t",i);
  if (s<0)
    fserror("writeIntTab","write int and <tab>","");
}

void writeFloat(float x)
{
  int s;
  s=fprintf(ofile,"%f\n",x);
  if (s<0)
    fserror("writeFloat","write float","");
  
}

void writeDouble(double x)
{
  int s;
  s=fprintf(ofile,"%5.3e\n",x);
  if (s<0)
    fserror("writeDouble","write double","");
}

/* -----------------------------------------------------------
   write arrays
   ----------------------------------------------------------- */

void fwriteDoubleArray(FILE *f, double *x, int rows, int cols)
{
  int	i,j,s1,s2;
  
  s1 = 0;
  for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
      if(j%10 == 9)
	fprintf(f,"\n\t");
      s1=fprintf(f,"%5.3e ",x[i*cols+j]);
      if (s1<0) break;
    }
    s2=fprintf(f,"\n");
    if ((s2<0)|(s1<0))
      fserror("fwriteDoubleArray","write double array","");
  }
}

void fwriteIntArray(FILE *f, int *x, int rows, int cols)
{
  int	i,j,s1,s2;
  
  s1 = 0;
  for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
      if(j%10 == 9)
	fprintf(f,"\n\t");
      s1=fprintf(f,"%d\t",x[i*cols+j]);
      if (s1<0) break;
    }
    s2=fprintf(f,"\n");
    if ((s2<0)|(s1<0))
      fserror("fwriteIntArray","write int array","");
  }
}


void writeIntArray(int *x, int rows, int cols)
{
  fwriteIntArray(ofile,x,rows,cols);
}

void fwriteIntMatrix(FILE *f, int **x, int rows, int cols)
{
  int	i,j,s1;
  
  for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
      if(j%10 == 9)
	fprintf(f,"\n\t");
      s1=fprintf(f,"%d\t",x[i][j]);
      if (s1<0) 
	fserror("fwriteIntMatrix","write int matrix","");
    }
    fprintf(f,"\n");
  }
}

void writeIntMatrix(int **x, int rows, int cols)
{
  fwriteIntMatrix(ofile,x,rows,cols);
}

void writeDoubleArray(double *x,int rows,int cols)
{
  fwriteDoubleArray(ofile,x,rows,cols);
}

void fwriteDoubleMatrix2(FILE *f, double **x, int rows, int cols)
{
  int	i, j, s;
  
  for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
      if(j%10 == 9){
	fprintf(f,"\n\t");
      }
      s=fprintf(f,"%5.3e ",x[i][j]);
      if (s<0)
	fserror("fwriteDoubleMatrix2","write double matrix","");
    }
    fprintf(f,"\n");
  }
}

void writeDoubleMatrix2(double **x, int rows, int cols)
{
  fwriteDoubleMatrix2(ofile,x,rows,cols);
}

void writeDoubleMatrix(double **x, int rows, int cols)
{
  int	i,j, c,s;
  
  for(i=0;i<rows;i++){
    for(j=0, c=0;j<cols;j++){
      if(++c > 10){
	fprintf(ofile,"\n\t");
	c = 0;
      }
      s=fprintf(ofile,"%5.3e ",x[i][j]);
      if (s<0)
	fserror("fwriteDoubleMatrix","write double matrix","");
    }
    fprintf(ofile,"\n");
  }
}

void writeFloatArray(float *x,int rows,int cols)
{
 writeArray(x,rows,cols);
}

void writeArray(float *x,int rows,int cols)
{
  int	
    i,j, c,s;
  
  for(i=0;i<rows;i++){
    for(j=0,c=0;j<cols;j++){
      if (c++>9){
	fprintf(ofile,"\n\t");
	c = 0;
      }
      s=fprintf(ofile,"%5.3e ",x[i*cols+j]);
      if (s<0)
	fserror("fwriteDoubleMatrix","write double matrix","");
    }
    fprintf(ofile,"\n");
  }
}

/* ************************************************* 
   error
 ************************************************* */ 

void fserror(char *proc, char *act, char *what){
  /* writes error message and aborts */
  
  fprintf(stderr, "\n ** Error ");
  if (proc[0]!='\0') /* not empty */
    fprintf(stderr, " in function '%s', ", proc);
  if (act[0]!='\0') /* not empty */
    fprintf(stderr, " trying to %s", act);
  if (what[0]!='\0') /* not empty */
    fprintf(stderr, " '%s'", what);
  fprintf(stderr, "\n ** .. exiting program");
  fprintf(stderr, " (from a function in 'interface.c').\n");
  
  exit(1);
}
