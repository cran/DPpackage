static const char bayes_c_sccs_id[] = "%W%";

/*+
 * bayes.c - Basic Conjugate Posteriors
 *
 * Description:
 *		Various basic conjugate posteriors
 *
 *		NOTE: all these routines assume that matrices are passed as arrays of
 *		arrays
 *
 * Modification History:
 * <unknown>	1.1		Peter Mueller
 *		Created.
 *
 * 29-JUN-1999	1.2		Stephen Morris (Tessella Support Services plc)
 *		Remove embedded comments (slash-star within a comment block)
 *		Added includes of "bayes.h" and "rand.h" to define some functions.
-*/

#include <math.h>
#include <stdio.h>

#include "ToolsBayes.h"
#include "ToolsNrutil.h"
#include "ToolsCMatrix.h"
#include "ToolsRand.h"
#include "ToolsVector.h"

/********************************************
 *         normal_normal
 ********************************************/
void nn_bayes(double r1, double **Spr_inv, double *mpr, 
		  double r2, double **Slik_inv, double *y,
		  double **Spo, double **Spo_inv, double *mpo, 
		  int p)
/* prior: N(x; mpr, r1*Spr)
   likl:  N(y; x, r2*Slik)
   returns:
   post:  N(x; mpo,Spo)
   Spo = (1/r1*Spr_inv + 1/r2*Slik_inv)^-1
   mpo = Spo*(1/r1*Spr*mpr + 1/r2*Slik*y)
*/
{ double *z;

  z = dvector(0,p-1);

  rA_plus_sB(1.0/r1, Spr_inv, 1.0/r2, Slik_inv, Spo_inv, p);
  d_invert(p, Spo_inv, Spo);
  rAx_plus_sBy(1.0/r1, Spr_inv, mpr, 1.0/r2, Slik_inv, y, z, p);
  Ax(Spo, z, mpo, p);
  
  free_dvector(z,0,p-1);
}

void nn_bayes2(double **Spr_inv, double *mpr, 
	       double **Slik_inv, double *y,
	       double **Spo, double **Spo_inv, double *mpo, 
		  int p)
/* same as nn_baeys, but assumes Sy = Slikl_inv*y */
{ double *z;

  z = dvector(0,p-1);

  rA_plus_sB(1.0, Spr_inv, 1.0, Slik_inv, Spo_inv, p);
  d_invert(p, Spo_inv, Spo);
  Ax_plus_y(Spr_inv, mpr, y, z, p);
  Ax(Spo, z, mpo, p);
  
  free_dvector(z,0,p-1);
}

void nn_bayes_rand(double r1, double **Spr_inv, double *mpr, 
		  double r2, double **Slik_inv, double *y,
		  double *theta, int p)
/* same as normal_normal, but returns a draw only 
   returns:
   theta: draw from the posterior N(theta; mpo, Spo):
*/
{ double *z, **S, **S_inv, *m;

  /* allocate memory */
  z = dvector(0,p-1);
  m = dvector(0,p-1);
  S = dmatrix(0,p-1,0,p-1);
  S_inv = dmatrix(0,p-1,0,p-1);

  rA_plus_sB(1.0/r1, Spr_inv, 1.0/r2, Slik_inv, S_inv, p);
  d_invert(p, S_inv, S);
  rAx_plus_sBy(1.0/r1, Spr_inv, mpr, 1.0/r2, Slik_inv, y, z, p);
  Ax(S, z, m, p);

  mvnS_rand(p,m,S,theta);

  free_dvector(z,0,p-1);
  free_dvector(m,0,p-1);
  free_dmatrix(S,0,p-1,0,p-1);
  free_dmatrix(S_inv,0,p-1,0,p-1);
}

void nn_bayes1(double r1, double Spr_inv, double mpr, 
		  double r2, double Slik_inv, double y,
		  double *Spo, double *mpo)
/* same as nn_bayes, just for 1-dim */
{ 
  double 
    sd;

  *Spo = 1.0/(Spr_inv/r1 + Slik_inv/r2);
  sd = sqrt(*Spo);
  *mpo = *Spo*( Spr_inv/r1*mpr + Slik_inv/r2*y);

}

void nn_bayes1_rand(double r1, double Spr_inv, double mpr, 
		  double r2, double Slik_inv, double y,
		  double *theta)
/* same as draw_normal_normal, but univariate 
   returns:
   theta: draw from the posterior N(theta; mpo, Spo):
*/
{ double S, m, sd, z;

  dstdnormal(1,&z);

  S = 1.0/(Spr_inv/r1 + Slik_inv/r2);
  sd = sqrt(S);
  m = S*( Spr_inv/r1*mpr + Slik_inv/r2*y);
  *theta = m+sd*z;
}
