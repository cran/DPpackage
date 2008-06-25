#if !defined(UNIXREP_H_INCLUDED)
#define UNIXREP_H_INCLUDED

/*+
 * unixrep.h - Unix Replacement Functions
 *
 * Version:	%I%
 *
 * Description:
 *		This routine declares functions that are present under True64 UNIX,
 *		but which are not implemented in Visual C++.  Currently, these are:
 *
 *			lgamma		Return the natural logarithm of the Gamma function
 *			rint		Return nearest integer
 *
 *		Note that his header file declares nothing if compiled under UNIX.
 *
 * Modification History:
 * 11-JUN-1999	1.1			Stephen Morris (Tessella Support Services plc)
 *		Created.
-*/

#if defined(UNIXREP_H_SCCS_ID)
static const char unixrep_h_sccs_id[] = "%W%";
#endif

#if defined(_WIN32)
double lgamma(double x);
double rint(double x);
#endif

#endif

