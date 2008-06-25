#if !defined(MATHCONST_H_INCLUDED)
#define MATHCONST_H_INCLUDED

/*+
 * mathconst.h - Maths Constants
 *
 * Version:	%I%
 *
 * Description:
 *		Supplies definitions of maths constants that are not available on some machines.
 *
 * Modification History:
 * 10-JUN-1999	1.1			Stephen Morris (Tessella Support Services plc)
 *		Created.
-*/

#if defined(MATHCONST_H_SCCS_ID)
static const char mathconst_h_sccs_id[] = "%W%";
#endif

#if defined(_WIN32)

#define M_PI (3.1415926535897932385)

#else

/* On UNIX, most constants are in "math.h" */

#include <math.h>

#endif

#endif

