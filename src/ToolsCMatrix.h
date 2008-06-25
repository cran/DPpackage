#if !defined(MATRIX_H_INCLUDED)
#define MATRIX_H_INCLUDED

/*+
 * matrix.h - Matrix Routine Definitions
 *
 * Description:
 *		Declarations of matrix-manipulation routines.
 *
 * Modification History:
 * <unknown>	1.1		Peter Muller
 *		Created.
 *
 * 25-JUN-1999	1.2		Stephen Morris (Tessella Support services plc)
 *		Add declaration of cd().
 *		Make cdinv() and invcd() "void".
-*/

int	d_eigen(int, double **, double *, double **);
int	d_chol_decomp(int, double **, double **);
int	d_invert(int, double **, double **);
int	d_inv_triang(int, double **, double **);
void cd(double **A, double **V, int p);
void cdinv(int, double **, double **);
void invcd(int, double **, double **);
void corr(double **S, double **R, int p, int prt);

int	d_inv_triang_lo(int n, double **a_in, double **v_out);
int	invtr(int n, double **a_in, double **v_out, char uplo);

void lp_error(char *proc, char *act, char *lpfct, int info);
void matrix_noexit();

#endif
