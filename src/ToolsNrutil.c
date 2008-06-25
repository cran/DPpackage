static const char nrutil_c_sccs_id[] = "%W%"; 
 
/*+ 
 * nrutil.c - NR Utility Routines 
 * 
 * Description: 
 *		Various auxilary functions for NR programs (syntax is changed to fit 
 *		the compiler on our system). 
 * 
 * Modification History: 
 * <unknown>	1.1		Peter Mueller 
 *		Created. 
 * 
 * 29-JUN-1999	1.2		Stephen Morris (Tessella Support Services plc) 
 *		Include stdlib.h for the definition of "malloc()". 
-*/ 
 
#include <stdio.h> 
#include <stdlib.h> 
 
#include "ToolsNrutil.h" 
 
int nv = 0; 
 
void nvprt(){printf("\n nv=%d\n",nv);} 
int nvn(){return(nv);} 
 
void nrerror(char *proc, char *act, char *what) 
{ 
  void exit(); 
 
  fprintf(stderr, "\n ** Error "); 
  if (proc[0]!='\0') /* not empty */ 
    fprintf(stderr, " in function '%s', ", proc); 
  if (act[0]!='\0') /* not empty */ 
    fprintf(stderr, " trying to %s", act); 
  if (what[0]!='\0') /* not empty */ 
    fprintf(stderr, " '%s',", what); 
  else 
    fprintf(stderr, ", "); 
  fprintf(stderr, "\n ** .. exiting program.\n"); 
  fprintf(stderr, "\n ** (a function in 'nrutil.c').\n"); 
  exit(1); 
} 
 
 
 
float *vector(int nl,int nh) 
{ 
        float *v; 
 
        v=(float *)calloc((unsigned) (nh-nl+1),sizeof(float)); 
        if (!v) nrerror("vector","allocate a float vector",""); 
        return v-nl; 
} 
 
 
int  *ivector(int nl,int nh) 
{ 
        int  *v; 
 
	nv += (nh-nl+1); 
        v=(int  *)calloc((unsigned) (nh-nl+1),sizeof(int)); 
        if (!v) nrerror("ivector","allocate an int vector",""); 
        return v-nl; 
} 
 
char  *cvector(int nl,int nh) 
{ 
        char  *v; 
 
	nv += (nh-nl+1); 
        v=(char  *)calloc((unsigned) (nh-nl+1),sizeof(char)); 
        if (!v) nrerror("cvector","allocate a char vector",""); 
        return v-nl; 
} 
 
char  **cmatrix(int nrl,int nrh,int ncl,int nch) 
{ 
        int i; 
	char **m; 
 
	nv += (nrh-nrl+1)*(nch-ncl+1); 
        m=(char  **)calloc((unsigned) (nrh-nrl+1),sizeof(char  *)); 
        if (!m) nrerror("cmatrix","allocate a char matrix (1st dim)", 
			""); 
        m -= nrl; 
 
        for(i=nrl;i<=nrh;i++) { 
                m[i]=(char  *)calloc((unsigned) (nch-ncl+1),sizeof(char)); 
                if (!m[i]) nrerror("cmatrix",  
		     "allocate a char matrix (2nd dim)",""); 
                m[i] -= ncl; 
        } 
        return m; 
} 
 
char **carray_2(int lo, int hi) 
{ 
  char **m; 
 
  nv += (hi-lo+1); 
  m=(char **)calloc((unsigned) (hi-lo+1),sizeof(char *)); 
  m -= lo; 
  return m; 
} 
 
double  *dvector(int nl,int nh) 
{ 
        double  *v; 
 
	nv += (nh-nl+1); 
        v=(double  *)calloc((unsigned) (nh-nl+1),sizeof(double)); 
        if (!v)  
	  nrerror("dvector","allocate a double vector",""); 
        return v-nl; 
} 
 
 
 
float  **matrix(int nrl,int nrh,int ncl,int nch) 
{ 
        int i; 
        float  **m; 
 
	nv += (nrh-nrl+1)*(nch-ncl+1); 
        m=(float  **) calloc((unsigned) (nrh-nrl+1),sizeof(float  *)); 
        if (!m)  
	  nrerror("matrix","allocate a float matrix (1st dim)",""); 
        m -= nrl; 
 
        for(i=nrl;i<=nrh;i++) { 
                m[i]=(float  *)calloc((unsigned) (nch-ncl+1), 
					    sizeof(float)); 
                if (!m[i])  
		  nrerror("matrix","allocate a float matrix (2nd dim)",""); 
                m[i] -= ncl; 
        } 
        return m; 
} 
 
double  **dmatrix(int nrl,int nrh,int ncl,int nch) 
{ 
        int i; 
        double  **m; 
	 
	nv += (nrh-nrl+1)*(nch-ncl+1); 
        m=(double  **)calloc((unsigned) (nrh-nrl+1), 
				   sizeof(double  *)); 
        if (!m)  
	  nrerror("dmatrix","allocate a double matrix (1st dim)",""); 
        m -= nrl; 
 
        for(i=nrl;i<=nrh;i++) { 
                m[i]=(double  *)calloc((unsigned) (nch-ncl+1), 
					     sizeof(double)); 
                if (!m[i])  
	  	  nrerror("dmatrix","allocate a double matrix (2nd dim)",""); 
                m[i] -= ncl; 
        } 
        return m; 
} 
 
double ***darray3(int n, int p, int q) 
{ 
  double ***a; 
  int i; 
  a = darray_3(0,n); 
  for(i=0;i<n;i++) 
    a[i] = dmatrix(0,p,0,q); 
  return(a); 
} 
 
double ****darray4(int m, int n, int p, int q) 
{ 
  double ****a; 
  int i; 
  a = darray_4(0,m); 
  for(i=0;i<m;i++) 
    a[i] = darray3(n,p,q); 
  return(a); 
} 
 
double ****darray_4(int lo, int hi) 
{ 
  double ****m; 
 
  nv += (hi-lo+1); 
  m=(double ****)calloc((unsigned) (hi-lo+1),sizeof(double ***)); 
  if (!m)  
    nrerror("darray_4","allocate a 4dim double array ",""); 
  m -= lo; 
  return m; 
} 
 
double ***darray_3(int lo, int hi) 
{ 
  double ***m; 
 
  nv += (hi-lo+1); 
  m=(double ***)calloc((unsigned) (hi-lo+1),sizeof(double **)); 
  if (!m)  
    nrerror("darray_3","allocate a 3dim double array ",""); 
  m -= lo; 
  return m; 
} 
 
double **darray_2(int lo, int hi) 
{ 
  double **m; 
 
  nv += (hi-lo+1); 
  m=(double **)calloc((unsigned) (hi-lo+1),sizeof(double *)); 
  if (!m)  
    nrerror("darray_2","allocate a 2dim double array ",""); 
  m -= lo; 
  return m; 
} 
 
int  **imatrix(int nrl,int nrh,int ncl,int nch) 
{ 
        int i, **m; 
 
	nv += (nrh-nrl+1)*(nch-ncl+1); 
        m=(int  **)calloc((unsigned) (nrh-nrl+1),sizeof(int  *)); 
        if (!m)  
	  nrerror("imatrix","allocate a int matrix (1st dim).",""); 
        m -= nrl; 
 
        for(i=nrl;i<=nrh;i++) { 
                m[i]=(int  *)calloc((unsigned) (nch-ncl+1),sizeof(int)); 
                if (!m[i])  
	  	  nrerror("imatrix","allocate a int matrix (2nd dim).",""); 
                m[i] -= ncl; 
        } 
        return m; 
} 
 
 
int **iarray_2(int lo, int hi) 
{ 
  int **m; 
 
  nv += (hi-lo+1); 
  m=(int **)calloc((unsigned) (hi-lo+1),sizeof(int *)); 
  if (!m)  
    nrerror("iarray_2","allocate a 2dim int array ",""); 
  m -= lo; 
  return m; 
} 
 
int ****iarray_4(int lo, int hi) 
{ 
  int ****m; 
 
  nv += (hi-lo+1); 
  m=(int ****)calloc((unsigned) (hi-lo+1),sizeof(int ***)); 
  if (!m)  
    nrerror("iarray_4","allocate a 4dim int array ",""); 
  m -= lo; 
  return m; 
} 
 
int ***iarray_3(int lo, int hi) 
{ 
  int ***m; 
 
  nv += (hi-lo+1); 
  m=(int ***)calloc((unsigned) (hi-lo+1),sizeof(int **)); 
  if (!m)  
    nrerror("iarray_3","allocate a 3dim int array ",""); 
  m -= lo; 
  return m; 
} 
 
int ***iarray3(int p1, int p2, int p3) 
{ 
  int ***m, i; 
 
  m = iarray_3(0,p1); 
  for (i=0;i<p1;i++) 
    m[i] = imatrix(0,p2,0,p3); 
  return m; 
} 
 
int ****iarray4(int p1, int p2, int p3, int p4) 
{ 
  int ****m, i; 
 
  m = iarray_4(0,p1); 
  for(i=0;i<p1;i++) 
    m[i] = iarray3(p2,p3,p4); 
  return m; 
} 
 
 
 
float  **submatrix(float  **a,int oldrl,int oldrh,int oldcl, 
		      int oldch,int newrl,int newcl) 
{ 
  int i,j; 
  float  **m; 
 
  m=(float  **)calloc((unsigned) (oldrh-oldrl+1), 
		      sizeof(float  *)); 
  if (!m)  
    nrerror("submatrix","allocate a submatrix",""); 
  m -= newrl; 
 
  for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl; 
 
  return m; 
} 
 
 
void free_vector(float  *v,int nl,int nh) 
{ 
        if( (v+nl) != NULL ) free((char  *) (v+nl)); 
	nv -= (nh-nl+1); 
} 
 
void free_ivector(int  *v,int nl,int nh) 
{ 
        if( (v+nl) != NULL ) free((char  *) (v+nl)); 
	nv -= (nh-nl+1); 
} 
 
void free_dvector(double  *v,int nl,int nh) 
{ 
        if( (v+nl) != NULL ) free((char  *) (v+nl)); 
	nv -= (nh-nl+1); 
} 
 
 
 
void free_matrix(float  **m,int nrl,int nrh,int ncl,int nch) 
{ 
        int i; 
 
        for(i=nrh;i>=nrl;i--) {if( (m[i]+ncl) != NULL )  
				 free((char  *) (m[i]+ncl));} 
        if( (m+nrl) != NULL ) free((char  *) (m+nrl)); 
        nv -= (nch-ncl+1)*(nrh-nrl+1); 
} 
 
void free_dmatrix(double  **m,int nrl,int nrh,int ncl,int nch) 
{ 
        int i; 
 
        for(i=nrh;i>=nrl;i--) {if( (m[i]+ncl) != NULL )  
				 free((char  *) (m[i]+ncl));} 
        if( (m+nrl) != NULL ) free((char  *) (m+nrl)); 
        nv -= (nch-ncl+1)*(nrh-nrl+1); 
} 
 
void free_imatrix(int  **m,int nrl,int nrh,int ncl,int nch) 
{ 
        int i; 
 
        for(i=nrh;i>=nrl;i--) {if( (m[i]+ncl) != NULL )  
				 free((char  *) (m[i]+ncl));} 
        if( (m+nrl) != NULL ) free((char  *) (m+nrl)); 
        nv -= (nch-ncl+1)*(nrh-nrl+1); 
} 
 
 
 
void free_submatrix(float  **b,int nrl,int nrh,int ncl,int nch) 
{ 
        if( (b+nrl) != NULL ) free((char *) (b+nrl)); 
} 
 
 
 
 
 
 
 
 
 
 
 
