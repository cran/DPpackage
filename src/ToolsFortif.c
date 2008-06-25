static const char fortif_c_sccs_id[] = "%W%";

/*+
 * fortif.c - Fortran Interface
 *
 * Version:	%I%
 *
 * Description:
 *		Handles the interface between FORTRAN and C on the two platforms, by
 *		providing a single C interface to the Fortran routine, and dispatching
 *		according to the platform.
 *
 *		This particular file only contains Win32 interface routines for FORTRAN
 *		subroutines that have a character arguments.  When a FORTRAN routine is
 *		compiled, for each character argument in the signature the compiler
 *		adds another unsigned long argument immediately following it, that
 *		holds the length of the char array being passed to it.
 *
 *		In Response\01, we have C routines calling FORTRAN ones, and the size
 *		of the char array is generally unknown (it may be passed in from the
 *		caller's caller).  However, for the FORTRAN routines we are interested
 *		in, the character argument is input only, not being changed by the
 *		routine.  Therefore, we can use the length of the string represented by
 *		the "char*" argument as a suitable length to pass to the FORTRAN code.
 *
 *		Note that this ".c" file does not contain interfaces for routines that
 *		do not pass in character arguments; these interfaces can be handled
 *		by simple definitions in the associated header file.
 *
 * Modification History:
 * 10-JUN-1999	1.1			Stephen Morris (Tessella Support Services plc)
 *		Created.
-*/

#include <string.h>

#define FORTIF_H_SCCS_ID
#include "ToolsFortif.h"

#if defined(_WIN32)

/*+
 * dpotrf - Lapack Routine
-*/

extern void __stdcall DPOTRF(char *uplo, unsigned int uplo_length, int *n,
	double *a, int *ncol_a, int *info);

void dpotrf(char *uplo, int *n, double *a, int *ncol_a, int *info)
{
	DPOTRF(uplo, (unsigned int) strlen(uplo), n, a, ncol_a, info);
}

/*+
 * dpotri - Lapack Routine
-*/

extern void __stdcall DPOTRI(char *uplo, unsigned int uplo_length, int *n,
	double *a, int *ncol_a,	int *info);

void dpotri(char *uplo, int *n, double *a, int *ncol_a, int *info)
{
	DPOTRI(uplo, (unsigned int) strlen(uplo), n, a, ncol_a, info);
}

/*+
 * dtrtri - Lapack Routine
-*/

extern void __stdcall DTRTRI(char *uplo, unsigned int uplo_length,
	char *diag, unsigned int diag_len, int *n, double *a, int *ncol_a, int *info);

void dtrtri(char *uplo, char *diag, int *n, double *a, int *ncol_a, int *info)
{
	DTRTRI(uplo, (unsigned int) strlen(uplo), diag, (unsigned int) strlen(diag),
		n, a, ncol_a, info);
}

/*+
 * dsyev - Lapack Routine
 */

extern void __stdcall DSYEV(char *jobz, unsigned int jobz_length,
	char *uplo, unsigned int uplo_length, int *n, double *a, int *lda,
	double *w, double *work, int *lwork, int *info);

void dsyev(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work,
	int *lwork, int *info)
{
	DSYEV(jobz, (unsigned int) strlen(jobz), uplo, (unsigned int) strlen(uplo),
		n, a, lda, w, work,	lwork, info);
}


#endif

