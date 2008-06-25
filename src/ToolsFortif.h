#if !defined(FORTIF_H_INCLUDED)
#define FORTIF_H_INCLUDED

/*+
 * fortif.h - Fortran Interface
 *
 * Description:
 *		Handles the interface between Fortran and C on multiple platforms.
 *		The following rules are used:
 *
 *			True64 Unix (or DEC Unix, whatever name it's curently known by).
 *				Fortran routine XYZ() is called from C as xyz_()
 *				A CHARACTER argument in a FORTRAN routine is handled by passing a char*
 *
 *			Win32
 *				Fortran routine XYZ() is called from C as XYZ()
 *				A CHARACTER argument in a FORTRAN routine is handled by passing a char*,
 *					and an extra unsigned int argument immediately after the char* giving
 *					the length of the string.
 *
 *		Handling this requires the following strategy:
 *
 *		a) From C, the routine XYZ() is called as xyz().
 *
 *		Then what we do depends on the operating system:
 *
 *		True64 Unix
 *			b) In the header file, we #define xyz as xyz_
 *
 *		Win32
 *			There are two cases to consider:
 *
 *				i) If the routine does not contain a CHARACTER argument
 *					a)	We #define xyz as XYZ
 *					b)	We declare XYZ using __stdcall calling sequence
 *
 *				ii) If the routine does contain a CHARACTER argument
 *					a)	We declare xyz() as routine in fortif.c
 *					b)	In fortif.c, we declare XYZ using the __stdcall calling sequence,
 *						and include the extra length argument.
 *					c)	In fortif.c, we put a routine xyz() that calls XYZ, passing
 *						(unsigned int) strlen(arg) as the length argument for the string.
 *
 * Modification History:
 * 10-JUN-1999	1.1			Stephen Morris (Tessella Support Services plc)
 *		Created.
-*/

#if defined(FORTIF_H_SCCS_ID)
static const char fortif_h_sccs_id[] = "%W%";
#endif


#if defined(_WIN32)

#define STDCALL __stdcall

#else
#define STDCALL
#define dpotrf	dpotrf_
#define dpotri	dpotri_
#define dtrtri	dtrtri_
#define dsyev	dsyev_

#endif


/*
 * The following routines include a char* argument, so are defined in fortif.c
 * in the Win32 environment, and are external (noting that xyz() has become
 * xyz_() because of the previous #define's) in the UNIX environment.
 */

void dpotrf(char *uplo, int *n, double *a, int *ncol_a, int *info);
void dpotri(char *uplo, int *n, double *a, int *ncol_a, int *info);
void dtrtri(char *uplo, char *diag, int *n, double *a, int *ncol_a, int *info);
void dsyev(char *jobz, char *uplo, int *n, double *a, int *lda, double *w,
	double *work, int *lwork, int *info);

#endif

