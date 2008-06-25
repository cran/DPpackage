void    nvprt();
int     nvn();
float   *vector(int,int);
float   **matrix(int,int,int,int);
float   **convert_matrix(float  *,int,int,int,int);
char  *cvector(int nl,int nh);
char **carray_2(int lo, int hi);
char  **cmatrix(int nrl,int nrh,int ncl,int nch);
double  *dvector(int,int);
double  **dmatrix(int,int,int,int);
double ***darray3(int n, int p, int q);
double ****darray_4(int lo, int hi);
double  ***darray_3(int, int);
double  **darray_2(int, int);
int     *ivector(int,int);
int     **imatrix(int,int,int,int);
int **iarray_2(int lo, int hi);
int ****iarray_4(int lo, int hi);
int ***iarray_3(int lo, int hi);
int ***iarray3(int p1, int p2, int p3);
int ****iarray4(int p1, int p2, int p3, int p4);
double ****darray4(int m, int n, int p, int q);

float   **submatrix(float  **,int,int,int,int,int,int);
void free_vector(float  *,int,int);
void free_dvector(double  *,int,int);
void free_ivector(int  *,int,int);
void free_matrix(float  **,int,int,int,int);
void free_dmatrix(double  **,int,int,int,int);
void free_imatrix(int  **,int,int,int,int);
void free_submatrix(float  **,int,int,int,int);
void free_convert_matrix(float  **,int,int,int,int);
void nrerror(char *proc, char *act, char *what);


