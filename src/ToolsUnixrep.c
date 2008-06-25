static const char unixrep_c_sccs_id[] = "%W%";

/*+
 * unixrep.c - UNIX Replacement Routines
 *
 * Version:	%I%
 *
 * Description:
 *		This routine defines functions that are present under True64 UNIX,
 *		but which are not implemented in Visual C++.  Currently, these are:
 *
 *			lgamma		Return the natural logarithm of the Gamma function
 *
 *		Note that this file defines nothing if compiled under UNIX.
 *
 * Modification History:
 * 11-JUN-1999	1.1		Stephen Morris (Tessella Support Services plc)
 *		Created.
 *
 * 02-JUL-1999	1.2		Stephen Morris (Tessella Support services plc)
 *		Change the reflection formula so that it only comes into effect for
 *		values less than 1 (rather than less than or equal to one).  As it stood,
 *		it could get into an infinite recursion loop (as lgamma(2 - x) = lgamma(x)
 *		when x was 1.0).
-*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ToolsMathconst.h"

#define UNIXREP_H_SCCS_ID
#include "ToolsUnixrep.h"

#if defined(_WIN32)

/*+
 * lgamma - Natural Logarithm of Gamma Function
 *
 * Description:
 *		This uses the approximation for the Gamma function derived by Lanczos.
 *
 *		The code below (and a description, but not a derivation), comes from
 *		the book "Numerical Recipes in C" by Press W, Flannery B,
 *		Teukolsky S, Vetterling W (Cambridge University Press, 1988).
 *
 *		The routine is valid for all values of the argument (xx) > 0.  The
 *		routine described in the book gives full accuracy for xx > 1.  For
 *		values of xx in	the range 0 to 1, we let y = (1 - xx), then use the
 *		fact that
 *
 *			Gamma(1 - y) = (pi * y) / (Gamma(1 + y) * sin(pi * y))
 *
 *		giving:
 *
 *			Gamma(xx) = (pi * (1 - xx)) / (Gamma(2 - xx) * sin(pi * (1 - xx)))
 *
 *		... and hence:
 *
 *			ln(Gamma(xx)) = ln(pi * (1 - xx)) - ln(Gamma(2 - xx) -
 *									ln(sin(pi * (1 - xx)))
 *
 *		Since 0 < xx < 1, then 1 < (2 - xx) < 2.  For this region, the Gamma
 *		function has values in the range of the normal "double" data type.
 *
 * Declaration:
 *		double lgamma(double xx)
 *
 * Arguments:
 *		xx (double, input)
 *			Function argument.  This should be > 0.  (The name "xx" is that
 *			used in the "Numerical Recipes" book).
 *
 * Returns:
 *		double
 *			Natural logarithm of the gamma function.
-*/

double lgamma(double xx)
{
	static double cof[6] = {76.18009173, -86.50532033, 24.01409822,
		-1.231739516, 0.120858003e-2, -0.536382e-5};
	double	x, tmp, result, ser;
	int		j;	


	if (xx <= 0.0) {
		printf("lgamma error - negative value of argument\n");
		exit(1);
	}
	else if (xx < 1.0) {

		/* Low accuracy, use reflection formula to increase it +*/

		result = log(M_PI * (1.0 - xx)) - lgamma(2.0 - xx) -
			log(sin(M_PI * (1.0 - xx)));
	}
	else {

		/* Maximum accuracy regime - code copied from Numerical Recipes */

		x = xx - 1.0;
		tmp = x + 5.5;
		tmp -= (x + 0.5) * log(tmp);
		ser = 1.0;
		for (j = 0; j <= 5; ++j) {
			x += 1.0;
			ser += cof[j] / x;
		}
		result = log(2.50662827465 * ser) - tmp;
	}

	return result;
}



/*+
 * rint - Rounding
 *
 * Description:
 *		Returns a double value, rounded upwards according to the modulus, i.e.
 *		rint(3.5) = 4.0, rint(-3.5) = -4.0;
 *
 * Declaration:
 *		double rint(double x)
 *
 * Arguments:
 *		x (double, input)
 *			Number to round.
 *
 * Returns:
 *		double
 *			Rounded number.
-*/

double rint(double x)
{
	double result;

	if (x >= 0.0) {
		result = floor(x + 0.5);
	}
	else {
		result = -floor(-x + 0.5);
	}

	return result;
}

#endif

