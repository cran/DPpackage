static const char regr_c_sccs_id[] = "%W%";

/*+
 * regr.c - Regression
 *
 * Description:
 *
 * Modification History:
 * <unknown>	1.1		Peter Mueller
 *		Created.  Corrected on 20/6/96 with:
 *			logit_setup:        added normal error, 
 *			logit_bayes_setup   same,
 *			lp_logistic2()		corrected normal pdf!!! 
 *
 * 29-JUN-1999	1.3		Stephen Morris (Tessella Support Services plc)
 *		Added this header.
 *		Removed embedded comments (slash-star in a comment block)
 *		Include bayes.h to define nn_bayes, and stdlib.h for exit().
+*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ToolsBayes.h"
#include "ToolsInterface.h"
#include "ToolsVector.h"
#include "ToolsCMatrix.h"
#include "ToolsNrutil.h"
#include "ToolsRand.h"
#include "ToolsRegr.h"



/* *************************************************
 *
 *   multiple regression
 *
 * *************************************************
   computes mean and cov matrix for:
MODEL
  y    = x'beta + eps, with eps ~ N(0,h[t])

ARGUMENTS
 mu            mean vector of normal posterior
 V=H_inv_cd,
 L=H_cdt       cov matrix: S=VV', L=1/V
 y             dep vars
 h             varianes (or inverse weights)
 x             matrix of covariates
 n             n obs
 p             n covariates
*/

/* *************************************************
   multiple regressiuon
 ************************************************* */
void regr(double *mu, double **H_inv_cd, double **H_cdt,
	  double *y, double *h, double **x, int n, int p,
	  int var_known)
{
  double 
    *xy, W, sig, SS, yhat,
    **H, **H_inv;
  int 
    i,j,k;

  /* allocate memory */
  xy = dvector(0,p);
  H = dmatrix(0,p,0,p);
  H_inv = dmatrix(0,p,0,p);

  for(j=0;j<p;j++){
    xy[j] = 0.0;
    for(k=0;k<p;k++)
      H[j][k] = 0.0;
  }

  for(i=0;i<n;i++){
    W = (var_known==1) ? h[i] : 1.0;
    for(j=0;j<p;j++)
      for(k=0;k<p;k++)
	{
	  H[j][k] += x[i][j]*x[i][k]/W;
	  if (H[j][k] == 0)
	    {
	      j = j;
	    }
	}
    for(j=0;j<p;j++)      
      xy[j] += x[i][j]*y[i]/W;
  }

  /* compute V_cd */
  d_chol_decomp(p,H,H_cdt); 
  d_inv_triang(p,H_cdt,H_inv_cd);
  ABt(H_inv_cd,p,p,H_inv_cd,p,p,H_inv);
  Ax(H_inv,xy,mu,p);

  /* compute s2 if var is unkwnown */
  if (var_known==0){
    for(i=0,SS=0;i<n;i++){
      for(j=0,yhat=0.0;j<p;j++)
	yhat += x[i][j]*mu[j];
      SS += (yhat-y[i])*(yhat-y[i]);
    }
    sig = SS/(n-p*1.0);
  }
  else sig=1.0;

  rA(sqrt(sig),H_inv_cd,H_inv_cd,p,p);
  rA(1.0/sqrt(sig),H_cdt,H_cdt,p,p);

  /* free memory */
  free_dvector(xy,0,p);

  free_dmatrix(H, 0,p,0,p);
  free_dmatrix(H_inv, 0,p,0,p);
}

/* *************************************************
   simple regr
 ************************************************* */

void regr1(double *mu, double *sd, double *prec,
	  double *y, double *h, double *x, int n, 
	  int var_known)
{
  double 
    xy=0, W, sig, SS, yhat, xtx=0;
  int 
    i;

  /* compute xtx = x'x = sum(x*x), xy = x'y = sum(x*y) */
  for(i=0,xtx=0,xy=0;i<n;i++){
    W = (var_known==1) ? h[i] : 1.0;
    xtx += x[i]*x[i]/W;
    xy  += x[i]*y[i]/W;
  }

  /* mu = (x'x)^-1 x'y */
  *mu = xy/xtx;

  /* compute var if h[i]'s are unknown */
  if (var_known==0){
    for(i=0,SS=0;i<n;i++){
      yhat = x[i]*(*mu);
      SS += (yhat-y[i])*(yhat-y[i]);
    }
    sig = SS/(n-1.0);
  }
  else sig=1.0;

  /* poserior sd and prec */
  *sd = sqrt(sig/xtx);
  *prec = 1.0/(*sd);

}


/* *************************************************
   bayes multiple regressiuon
 ************************************************* */
void bayes_regr(double *mu2, double **V, double **L,
	  double *y, double *h, double **x, 
	  double *mu1, double **S1_inv, int n, int p,
	  int var_known)
{
  double 
    *b,**S2_inv,**S2;
  int 
    i;
  double 
    *xy, W, **H;

  /* var unknown not yet implemented ! */
  if (var_known==0){
    printf(
   "\n*** Error in bayes_regr: unknown var not yet impl\n");
    exit(1);
  }

  /* alloc mem */
  S2_inv = dmatrix(0,p,0,p);
  S2 = dmatrix(0,p,0,p);
  b = dvector(0,p);
  xy = dvector(0,p);
  H = dmatrix(0,p,0,p);

  /* get likelihood from regression */
  x_zero(xy,p);
  A_zero(H,p);

  for(i=0;i<n;i++){
    W = (var_known==1) ? h[i] : 1.0;
    A_plus_rxxt(H,1.0/W,x[i],H,p);
    x_plus_ry(xy,y[i]/W,x[i],xy,p);
  }

  /* S2_inv = S1_inv + H*sig^(-2) */
  A_plus_B(S1_inv,H,S2_inv,p,p);
  Ax_plus_y(S1_inv,mu1,xy,b,p);

  /* compute V_cd */
  d_chol_decomp(p,S2_inv,L);
  d_inv_triang(p,L,V);
  ABt(V,p,p,V,p,p,S2);
  Ax(S2,b,mu2,p);

  /* free memory */
  free_dvector(xy,0,p);
  free_dmatrix(H, 0,p,0,p);
  free_dmatrix(S2_inv ,0,p,0,p);
  free_dmatrix(S2 ,0,p,0,p);
  free_dvector(b,0,p);
}

void bayes_regr_rand(double *th,
	  double *y, double *h, double **x, 
	  double *mu1, double **S1_inv, int n, int p,
	  int var_known)
{
  double *m2,**V,**L;

  /*alloc mem */
  m2 = dvector(0,p);
  V = dmatrix(0,p,0,p);
  L = dmatrix(0,p,0,p);

  bayes_regr(m2,V,L,y,h,x,mu1,S1_inv,n,p,var_known);
  mvn_rand(p,m2,V,th);

  free_dvector(m2,0,p);
  free_dmatrix(V,0,p,0,p);
  free_dmatrix(L,0,p,0,p);
}


/* -------------------------------------------------
   old version
 ------------------------------------------------- */
void bayes_regr_old(double *mu2, double **V, double **L,
	  double *y, double *h, double **x, 
	  double *mu1, double **S1_inv, int n, int p,
	  int var_known)
{
  double 
    *mu,**Slikl_inv,**S2_inv,**S2;

  /* alloc mem */
  Slikl_inv = dmatrix(0,p,0,p);
  S2_inv = dmatrix(0,p,0,p);
  S2 = dmatrix(0,p,0,p);
  mu = dvector(0,p);

  /* get likelihood from regression */
  regr(mu,V,L,y,h,x,n,p,var_known);
  AtB(L,p,p,L,p,p,Slikl_inv);

  nn_bayes(1.0,S1_inv,mu1,1.0,Slikl_inv,mu,S2,S2_inv,mu2,p);
  d_chol_decomp(p,S2_inv,L);
  d_inv_triang(p,L,V);

  /* free mem */
  free_dmatrix(Slikl_inv ,0,p,0,p);
  free_dmatrix(S2_inv ,0,p,0,p);
  free_dmatrix(S2 ,0,p,0,p);
  free_dvector(mu,0,p);
}


/* *************************************************
   bayes simple regr
 ************************************************* */

void bayes_regr1(double *mu, double *sd, double *prec,
	  double *y, double *h, double *x,
	  double mu1, double s1, int n)
{
  double 
    xy=0, W, sig, xtx=0, mu2, s2, spost;
  int 
    i;

  /* compute xtx = x'x = sum(x*x), xy = x'y = sum(x*y) */
  for(i=0,xtx=0,xy=0;i<n;i++){
    W = h[i];
    xtx += x[i]*x[i]/W;
    xy  += x[i]*y[i]/W;
  }

  /* mu = (x'x)^-1 x'y */
  mu2 = xy/xtx;

   sig=1.0;

  /* poserior sd and prec */
  s2 = sig/xtx;
  
  spost = 1.0/(1.0/s1 + 1.0/s2);
  *sd = sqrt(spost);
  *mu = (mu1/s1 + mu2/s2)*spost;
  *prec = 1.0/(*sd);
}





/* *************************************************
 *
 *   indep chain for logistic regr
 * *************************************************

 V,L           cov matrix: S=VV', L=1/V
 p:            dim par vector
 N             binomial N (=0 for normal)
 n             n obs
 sig           var for normal error (ignored for binom)
*/

void logit_setup(double *mu, double **V, double **L,
		 double *y, double **x, double *offset, int N,
		 int n, int p, double *sig)
{
  double 
    e,D,W, eta0, mu0; /* auxilary */
  double 
    *dmu, **H, **H_cdt, **H_inv_cd, **H_inv, 
    *z, *h;
  int 
    i,j,k, normal /* indicator for normal */;

  /* allocate memory */
  h = dvector(0,p);
  dmu = dvector(0,p);
  z = dvector(0,n);

  H = dmatrix(0,p,0,p);
  H_cdt = L; 
  H_inv_cd = V;
  H_inv = dmatrix(0,p,0,p);

  normal = (N==0); /* normal error */
  if (normal)
    N=1; /* to avoid splitting cases.. */

  for(j=0;j<p;j++){
    h[j] = 0.0;
    for(k=0;k<p;k++)
      H[j][k] = 0.0;
  }

  for(i=0;i<n;i++){
    /* -------------------------------------------------
       (1) compute z
       ------------------------------------------------- */
    xy(mu,x[i],&eta0,p);
    eta0 += offset[i];
    e = exp(-eta0);
    mu0 = N*1.0/(1.0+e);
    if (mu0 > 0.999*N) mu0 = 0.999*N;
    if (mu0 < 0.001*N) mu0 = 0.001*N; 
    if (normal){
      D = 1.0/(mu0*(1-mu0));
      W = D*D*sig[i];
    }
    else {
      D = 1.0/(mu0*(1.0-mu0/N));
      W = D*D*mu0*(1.0-mu0/N); /* could simplify .. */
    }
    z[i] = D*(y[i]-mu0);

    /* -------------------------------------------------
       (2) regression
       ------------------------------------------------- */
    for(j=0;j<p;j++)
      for(k=0;k<p;k++)
	H[j][k] += x[i][j]*x[i][k]/W;
    for(j=0;j<p;j++)
      h[j] += x[i][j]*z[i]/W;
  }
  /* compute V_cd */
  d_chol_decomp(p,H,H_cdt); 
  d_inv_triang(p,H_cdt,H_inv_cd);
  ABt(H_inv_cd,p,p,H_inv_cd,p,p,H_inv);
  Ax(H_inv,h,dmu,p);

  for(j=0;j<p;j++)
    mu[j] += dmu[j];

  /* free memory */
  free_dvector(h,0,p);
  free_dvector(dmu,0,p);
  free_dvector(z,0,n);

  free_dmatrix(H, 0,p,0,p);
  free_dmatrix(H_inv, 0,p,0,p);
}

/* -------------------------------------------------
   bayes_logit_setup
 ------------------------------------------------- */
void bayes_logit_setup(double *mu, double **V, double **L,
		       double *m1, double **S1_inv,
		       double *y, double **x, double *offset, int N,
		       int n, int p, double *sig)
{ /* old version of bayes_logit_setup2  - kept for compatibility */
  int *Nvec, normal;
  Nvec = ivector(0,n);

  i_assgn_r(Nvec,N,n);
  /* normal samplg dist */
  normal = (N==0);
  bayes_logit_setup2(mu,V,L,m1,S1_inv,y,x,offset,Nvec,n,p,sig,normal);
  free_ivector(Nvec,0,n);
}

void bayes_logit_setup2(double *mu, double **V, double **L,
		       double *m1, double **S1_inv,
		       double *y, double **x, double *offset, int *N,
		       int n, int p, double *sig, int normal)
{
  double 
    e,D,W, eta0, mu0; /* auxilary */
  double 
    *dmu, **H,
    *z, *h, *dmmy,
    **S2, **S2_inv, *dm1, *b;
  int 
    i,j,k;

  /* allocate memory */
  h = dvector(0,p);
  dmu = dvector(0,p);
  dm1 = dvector(0,p);
  b = dvector(0,p);
  dmmy = dvector(0,100);
  z = dvector(0,n);

  H = dmatrix(0,p,0,p);
  S2 = dmatrix(0,p,0,p);
  S2_inv = dmatrix(0,p,0,p);

  /* had here: if (normal) i_assgn_r(N2,1,n); */

  /* initialize to zero */
  for(j=0;j<p;j++){
    h[j] = 0.0;
    for(k=0;k<p;k++)
      H[j][k] = 0.0;
  }

  /* prior mean for mu - mu0 */
  x_min_y(m1,mu,dm1,p);

  for(i=0;i<n;i++){
    if ((!normal)&(N[i]==0)) 
      continue; /* void observation - allowed to be included */
    /* -------------------------------------------------
       (1) compute z
       ------------------------------------------------- */
    xy(mu,x[i],&eta0,p);
    eta0 += offset[i];
    e = exp(-eta0);
    if (normal){
      mu0 = 1.0/(1.0+e);
      if (mu0 > 0.999) mu0 = 0.999;
      if (mu0 < 0.001) mu0 = 0.001;
      D = 1.0/(mu0*(1-mu0));
      W = D*D*sig[i];
    }
    else {
      mu0 = N[i]*1.0/(1.0+e);
      if (mu0 > 0.999*N[i]) mu0 = 0.999*N[i];
      if (mu0 < 0.001*N[i]) mu0 = 0.001*N[i]; 
      D = 1.0/(mu0*(1.0-mu0/N[i]));
      W = D*D*mu0*(1.0-mu0/N[i]); /* could simplify .. */
    }
    z[i] = D*(y[i]-mu0);

    /* -------------------------------------------------
       (2) regression
       ------------------------------------------------- */
    for(j=0;j<p;j++)
      for(k=0;k<p;k++)
	H[j][k] += x[i][j]*x[i][k]/W;
    for(j=0;j<p;j++)
      h[j] += x[i][j]*z[i]/W;
  }

  /* S2_inv = S1_inv + H*sig^(-2) */
  A_plus_B(S1_inv,H,S2_inv,p,p);
  Ax_plus_y(S1_inv,dm1,h,b,p);
  /* note: dm1 instead of mu1! */

  /* compute V_cd */
  d_chol_decomp(p,S2_inv,L);
  d_inv_triang(p,L,V);
  ABt(V,p,p,V,p,p,S2);
  Ax(S2,b,dmu,p);

  for(j=0;j<p;j++)
    mu[j] += dmu[j];

  /* free memory */
  free_dvector(dm1,0,p);
  free_dvector(b,0,p);
  free_dvector(h,0,p);
  free_dvector(dmu,0,p);
  free_dvector(z,0,n);

  free_dvector(dmmy,0,100);

  free_dmatrix(S2_inv, 0,p,0,p);
  free_dmatrix(H, 0,p,0,p);
  free_dmatrix(S2, 0,p,0,p);
}

/* -------------------------------------------------
   lp logistic
 ------------------------------------------------- */
/* returns logistic log pdf */
double lp_logistic(double *beta, double *y, double **x, 
		   double *offset, int N,int n, int p, double *sig)
{ /* old version of lp_logistic2 - kept for compatibility */
  int 
    *Nvec,normal;
  double f;

  Nvec = ivector(0,n);
  normal = (N==0) ? 1 : 0;

  i_assgn_r(Nvec,N,n);
  normal = (N==0) ? 1 : 0;
  f = lp_logistic2(beta, y, x, offset, Nvec, n, p, sig, normal);
  free_ivector(Nvec,0,n);
  return (f);
}

double lp_logistic2(double *beta, double *y, double **x, 
		    double *offset, int *N,int n, int p, double *sig,
		    int normal)
{
  int 
    i;
  double 
    lp,pi,eta;

  for(i=0,lp=0.0;i<n;i++){
    if ((!normal) & (N[i]==0))
      continue; /* void observation */
    xy(beta,x[i],&eta,p);
    eta += offset[i];
    pi = 1.0/(1.0+exp(-eta));
    if (normal)
      lp += -0.5*(y[i]-pi)*(y[i]-pi)/sig[i];
    /* CORRECTED from missing "-0.5*" !!!! */
    else{
      if ( (pi==0.0) | (pi == 1.0) ) lp = -10e10;
      else
      lp += y[i]*log(pi) + (N[i]-y[i])*(log(1.0-pi));
    }
  }
  return (lp);
}

/* *************************************************;
   regr_logistic
 ************************************************* 
PROGRAM FLOW
   Loop over n_iter iterations
   (1) if necesary, update approximation
   (2) draw bnew from N(mub_hat, sigb_hat^2)
   (3) compare p(bnew,..) with p(b,...) and accept with
       appropriate prob (independence chain)
ARGUMENTS
 n_iter        number its in indep chain
 beta:         initial guess
 p:            dim par vector
 n:            n obs
 y:            0/1 obs
 x             covariates
 n_jump        incremented by 1 if jump accepted, unchanged else
*/
void regr_logistic(double *beta, 
		   double *y, double **x, double *offset, int N,
		   int n, int p, double *sig,
		   int n_iter, int *n_jump)
{ 
  int 
    it,
    updated;    /* indicator whether envelope is up-to-date */

  double
    *e, *bnew,
    qz,q,lp,lpz,dpdq,u,probjump,
    **L, **V, *mu_hat;


  /* allocate memory */
  L = dmatrix(0,p,0,p);
  V = dmatrix(0,p,0,p);
  mu_hat = dvector(0,p);
  e = dvector(0,p);
  bnew = dvector(0,p);

  lp = lp_logistic(beta,y,x,offset,N,n,p,sig);
  y_assgn_x(mu_hat,beta,p);
  for(it=0,updated=0;it<n_iter;it++){
    /* -------------------------------------------------
       (1) update approx
       ------------------------------------------------- */
    if (updated==0){
      logit_setup(mu_hat,V,L,y,x,offset,N,
		  n,p,sig);
      updated = 1;
    }
    /* -------------------------------------------------
       (2) draw new b1..6
       ------------------------------------------------- */
    dstdnormal(p,e);
    xtx(e,&qz,p);
    qz *= -0.5;
    Ax_plus_y(V,e,mu_hat,bnew,p);
    
    /* -------------------------------------------------
       (3) compare
       ------------------------------------------------- */
    lpz = lp_logistic(bnew,y,x,offset,N,n,p,sig);
    A_xminusy(L,beta,mu_hat,e,p);
    xtx(e,&q,p);
    q *= -0.5;
    dpdq = lpz-lp + q - qz;
    probjump = (dpdq>0.0) ? 1.0 : exp(dpdq);
    duniform(1,&u);
    if (u<probjump){
      *n_jump += 1;
      y_assgn_x(beta,bnew,p);
      lp = lpz;
      updated = 0; /* approx is outdated */
    }
  }

  /* free memory */
  free_dvector(mu_hat,0,p);
  free_dvector(e,0,p);
  free_dvector(bnew,0,p);

  free_dmatrix(L, 0,p,0,p);
  free_dmatrix(V, 0,p,0,p);

}

/* -------------------------------------------------
   bayes_regr_logistic
 ------------------------------------------------- */
void bayes_regr_logistic(double *beta, 
   double *m1, double **S1_inv,		       
   double *y, double **x, double *offset, int N,
   int n, int p, double *sig, int n_iter, int *n_jump)
{ /* old version of bayes_regr_logistic2: kept for compatibility */
  int *Nvec, normal;
  Nvec = ivector(0,n);

  i_assgn_r(Nvec,N,n);
  normal = (N==0) ? 1 : 0;
  bayes_regr_logistic2(beta, m1, S1_inv, y, x, offset, Nvec,
		       n,  p, sig, n_iter, n_jump, normal);
  free_ivector(Nvec,0,n);
}

void bayes_regr_logistic2(double *beta, 
   double *m1, double **S1_inv,		       
   double *y, double **x, double *offset, int *N,
   int n, int p, double *sig, int n_iter, int *n_jump, int normal)
{
  int 
    it;

  double
    *e, *bnew,lpr,lprz,
    qz,q,lp,lpz,dpdq,u,probjump,
    **L, **Ltilde, **V, *mu_hat, *mu_hattilde, **L1;


  /* allocate memory */
  L = dmatrix(0,p,0,p);
  L1 = dmatrix(0,p,0,p);
  V = dmatrix(0,p,0,p);
  mu_hat = dvector(0,p);
  Ltilde = dmatrix(0,p,0,p);
  mu_hattilde = dvector(0,p);
  e = dvector(0,p);
  bnew = dvector(0,p);

  d_chol_decomp(p,S1_inv,L1);
  lp = lp_logistic2(beta,y,x,offset,N,n,p,sig,normal);
  lpr = mvn_logpdf(beta,m1,L1,p); /* log prior */

  y_assgn_x(mu_hat,beta,p);
  bayes_logit_setup2(mu_hat,V,L,m1,S1_inv,y,x,offset,N,
		     n,p,sig,normal); /* V is dummy only */
  for(it=0;it<n_iter;it++){
    /* -------------------------------------------------
       (2) draw new b1..6
       ------------------------------------------------- */
    dstdnormal(p,e);
    xtx(e,&qz,p);
    qz *= -0.5;
    Ax_plus_y(V,e,mu_hat,bnew,p);
    /* -------------------------------------------------
       (1) update approx
       ------------------------------------------------- */
    y_assgn_x(mu_hattilde,bnew,p);
    bayes_logit_setup2(mu_hattilde,V,Ltilde,m1,S1_inv,y,x,offset,N,
      n,p,sig,normal); /* V is dummy only */
    
    /* -------------------------------------------------
       (3) compare
       ------------------------------------------------- */
    lpz = lp_logistic2(bnew,y,x,offset,N,n,p,sig,normal);
    lprz = mvn_logpdf(bnew, m1, L1, p); /* log prior */
    A_xminusy(Ltilde,beta,mu_hattilde,e,p);
    xtx(e,&q,p);
    q *= -0.5;
    dpdq = lpz+lprz-lp-lpr + q - qz;
    probjump = (dpdq>0.0) ? 1.0 : exp(dpdq);
    duniform(1,&u);
    if (u<probjump){
      *n_jump += 1;
      y_assgn_x(beta,bnew,p);
      lp = lpz;
      y_assgn_x(mu_hat,mu_hattilde,p);
      B_assgn_A(L,Ltilde,p);
    }
  }

  /* free memory */
  free_dvector(mu_hat,0,p);
  free_dvector(mu_hattilde,0,p);
  free_dvector(e,0,p);
  free_dvector(bnew,0,p);

  free_dmatrix(L1, 0,p,0,p);
  free_dmatrix(L, 0,p,0,p);
  free_dmatrix(Ltilde, 0,p,0,p);
  free_dmatrix(V, 0,p,0,p);

}

void dummy()
{
  printf("\n i'm the dummy \n");
}
