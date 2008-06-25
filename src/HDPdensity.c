/* =========================================================== 
							       
  meta.c                 				       
							       
 =========================================================== */

/* MODIFIED: for generic normal-normal eps-DP mixture,
   w/o longitudinal model */
   
/* HAVE commented out "sample_sig, mdp_samplealpha & sample_S" !! 
   for meta-analysis (ENAR ..) */

/* Meta analysis over related non-par Bayes models:


   Obs:
      th[ji]|mu[ji] ~ N(th[ji]; mu[ji],S)
      mu[ji]|x[ji] ~ 
        phi(x[ji]; mu[ji],S)*{(1-eps)*F0(mu) + eps*Fj(mu)} / m(x[ji]|F)
      where phi(x;m,S) = N(x;m,S) and
                       /			   /
      m(x|F) = (1-eps) | phi(x;mu,S) dF0(mu) + eps | phi(x;mu,S) dFj(mu)
                       /			   /
      m(X|F) = Prod{j,i} m(x[ji]|F)

   Population distr:
      (F0,F1...FJ) ~ DP(F0|alpha G0)*Prod{j} DP(Fj|alpha G0) * m(X|F)

   Hyperprior: ....
*/
				   
				   
#include "math.h"
#include "stdio.h"
#include "ToolsInterface.h"
#include "ToolsVector.h"
#include "ToolsCMatrix.h"
#include "ToolsRand.h"
#include "HDPdensity.h"
#include "ToolsNrutil.h"
#include "ToolsBayes.h"
#include "ToolsRegr.h"
#include "ToolsMess.h"
#include <stdlib.h>



/* =========================================================== */
/*
/*		 main 
/*
/* =========================================================== */


/* *************************************************
   global vars
 ************************************************* */
struct MDP_PARS mdp;
struct MCMC g;
struct Model mdl;
struct Data dta;
static int line=1; /* counts number of writes to
		      output files */

/********************************************/
/*         main */
/********************************************/

void hdpmn()
{
  int j;

  init();
  mdp_init();
  for(g.iter=0;g.iter<g.niter;g.iter++){ 
    mdp_gibbs();
    if (g.sampleeps==1) sample_eps();
    

    if (g.verbose>0){
      printf(" %d:%d=",g.iter,mdp.n_class);
      for(j=0;j<=mdp.J;j++)
	if (j==mdp.J) printf("%d",mdp.K[j]);
	else printf("%d+",mdp.K[j]);
      printf("  ");
      /* printf("%2.0f,%2.0f",100.0*nb_jump/(dta.npa*g.nbatch*1.0),
	 100.0*nt_jump/(dta.npa*g.nbatch*1.0)); */
      if (g.iter % 5 == 0) printf("\n");
    }
  }
  finish();
}

  
/* =========================================================== 

		 initialization & finish

 ===========================================================  */
/* ******************************************
         init
 *******************************************/

void init()
{
  int	i,j,k,s,currpa,pa,ni,npa;
  
/* -------------------------------------------------
   read pars
 ------------------------------------------------- */
  trim("init.giis","init-trim.giis");
  openIn("init-trim.giis");
  scanInt(" data-npatients ", &dta.npa);
  scanInt(" data-nstudies ",&dta.J);
  scanInt(" model-p1 ", &mdl.p1);
  dta.j = ivector(0,dta.npa);
  scanIntArray(" data-study ", dta.j, dta.npa);
  /* this is new compared to meta.c !! */


  scanInt(" mcmc-niter ",&g.niter);
  scanInt(" mcmc-npi ",&g.npi);
  scanInt(" mcmc-nbatch ",&g.nbatch);
  scanInt(" mcmc-ndiscard ",&g.ndiscard);
  scanInt(" mcmc-verbose ",&g.verbose);
  /* r.v. seeds */
  scanLong(" mcmc-seed1", &g.is1);
  scanLong(" mcmc-seed2", &g.is2);
  setseed(g.is1, g.is2);

  scanInt(" mcmc-eps ",&g.sampleeps);


  scanDouble(" model-eps ",&mdl.eps);
  scanDouble(" model-ae ",&mdl.ae);
  scanDouble(" model-be ",&mdl.be);

  scanDouble(" model-pe0 ",&mdl.pe0);
  scanDouble(" model-pe1 ",&mdl.pe1);

  closeIn();

/* -------------------------------------------------
   allocate memory
 ------------------------------------------------- */
  dta.z = dmatrix(0,dta.npa,0,mdl.p1);
  dta.n = ivector(0,dta.J);
  i_zero(dta.n,dta.J);

/* -------------------------------------------------
   read data
 ------------------------------------------------- */

  openIn("data-z.giis");
  scanDoubleMatrix(" data-z ",dta.z,dta.npa,mdl.p1);
  closeIn();
  /* initialize hash table for first obs on patient i */
  for(i=0; i<dta.npa; i++){
    s = dta.j[i];
    dta.n[s] += 1;
  }

  /* check ni's */
  for(npa=0,j=1;j<=dta.J;j++)
    npa+=dta.n[j];
  if (npa != dta.npa){ /* should never happen ... */
    printf("\n *** Error: Sum n[j] != npa.\n");
    exit(1);
  }
    

}

void	finish()
{
  mdp_finish();
}


/* =========================================================== */
/* MDP
/* 
/*
/* =========================================================== */

/* ============================================================= 
   Multivariate MDP model
   for use in gdm simulation
 ============================================================= 
   Model:
   z[i] | pi[i] is N(y[i]; mu[i], S)
                   pi[i] = (mu[i])
   pi[i]        is DP(G0, alpha)
   G0(pi)       is N(mu; mean, B)
   S            is W(1/S; q, 1/R*1/q), i.e. E(S)=R

Notation:
   Some of the pi[i]'s can be identical. The list of distinct
   pi[i]'s is denoted by pi*[0]..pi*[n_class].

   Details as in Escobar & West */


/* ***********************************************************
   all MDP pars
/* *********************************************************** */

FILE *mfile, *Vfile, *Lfile, *Bfile, *Sfile, 
  *parFile, *meanFile;

/* -----------------------------------------------------------
   Gibbs iterations
  ----------------------------------------------------------- */

int 
   m_prior, B_prior, S_prior, alpha_prior,
   p,
   n_predupdate;

int 
   n_class0, *member0, *r0; /* initial config */


/* *************************************************
/* initialization
 ************************************************* */

/* ***********************************************************
   init mdp
/* *********************************************************** */
void mdp_init()
{
  int 
    i,j;
  double 
    alpha;

/* -----------------------------------------------------------
   read in sample size and dim
/* ----------------------------------------------------------- */
  trim("init.mdp","init-trim.mdp");
  openIn("init-trim.mdp");
  mdp.N = dta.npa;
  mdp.J = dta.J;
  scanInt(" data-p ", &p);
  mdp.p = p;

/* -----------------------------------------------------------
   allocate memory for pars & data
/* -----------------------------------------------------------  */
  mdp.Y = dmatrix(0,mdp.N-1,0,p-1);
  mdp.Vt = darray_3(0,mdp.N-1);
  mdp.Vt_cd_inv = darray_3(0,mdp.N-1);
  mdp.mu = darray_2(0,mdp.N);
  mdp.mj = darray_2(0,mdp.N);
  mdp.m = dvector(0,p-1);

  mdp.S = dmatrix(0,p-1,0,p-1);
  mdp.St = dmatrix(0,p-1,0,p-1);
  mdp.St_cd_inv = dmatrix(0,p-1,0,p-1);
  mdp.S_cd_inv = dmatrix(0,p-1,0,p-1);
  mdp.S_inv = dmatrix(0,p-1,0,p-1);
  mdp.R = dmatrix(0,p-1,0,p-1);
  mdp.B = dmatrix(0,p-1,0,p-1);
  mdp.B_inv = dmatrix(0,p-1,0,p-1);
  mdp.B_cd_inv = dmatrix(0,p-1,0,p-1);
  mdp.C = dmatrix(0,p-1,0,p-1);

  mdp.aa = dvector(0,p-1);
  mdp.A_inv = dmatrix(0,p-1,0,p-1);

  mdp.alpha = dvector(0,mdp.J);
  mdp.eta = dvector(0,mdp.J);

  /* note: count, first & upd allow for one extra class =
     "new class" */
  mdp.count = ivector(0,mdp.N);
  mdp.first = ivector(0,mdp.N);
  mdp.member = ivector(0,mdp.N-1);
  mdp.updated = ivector(0,mdp.N-1);
  mdp.prev = ivector(0,mdp.N-1);
  mdp.next = ivector(0,mdp.N);
  mdp.new = ivector(0,mdp.N-1);
  mdp.r = ivector(0,mdp.N-1);
  mdp.j = ivector(0,mdp.N-1);
  mdp.K = ivector(0,mdp.J);
  mdp.n = ivector(0,mdp.J);
  mdp.qv = dvector(0,mdp.N+1);

  i_zero(mdp.count,mdp.N+1);
  i_zero(mdp.n,mdp.J);
  i_zero(mdp.K,mdp.J);

  member0 = ivector(0,mdp.N-1);
  r0 = ivector(0,mdp.N-1);

/* -----------------------------------------------------------
   read in hyperpars and initial values
/* ----------------------------------------------------------- */

  scanInt(" model-m_prior", &m_prior);
  scanInt(" model-B_prior", &B_prior);
  scanInt(" model-S_prior", &S_prior);
  scanInt(" model-alpha_prior", &alpha_prior);

  scanInt(" mcmc-npredupdate ",&n_predupdate);
  /* scanInt(" par-s ", &mdp.s); not used... */

  /* read and invert S */
  scanDoubleMatrix(" par-S ",mdp.S,p,p);
  /* rA(0.1,mdp.S,mdp.S,p,p); */
  d_invert(p,mdp.S,mdp.S_inv);

  scanInt(" par-q ",&mdp.q);
  scanDoubleMatrix(" par-R ",mdp.R, p,p);
  rA(0.1,mdp.R,mdp.R,p,p);

  /* read and invert B */
  scanDoubleMatrix(" par-B ",mdp.B,p,p);
  d_invert(p,mdp.B,mdp.B_inv);
  d_chol_decomp(p,mdp.B_inv,mdp.B_cd_inv);

  /* initialize mdp.St */
  rA_plus_sB(1.0,mdp.S,1.0,mdp.B,mdp.St,p);
  invcd(p,mdp.St,mdp.St_cd_inv);

  scanInt(" par-c ",&mdp.c);
  scanDoubleMatrix(" par-C ",mdp.C, p,p);
  scanDoubleArray(" par-m ",mdp.m, p);
  scanDoubleArray(" par-aa ", mdp.aa, p);
  scanDoubleMatrix(" par-Ainv ", mdp.A_inv, p,p);
  scanDouble(" par-alpha ",&alpha);
  x_assgn_r(mdp.alpha,alpha,mdp.J+1);
  scanDouble(" par-aalpha ", &mdp.a0);
  scanDouble(" par-balpha ", &mdp.b0);

/* -----------------------------------------------------------
   initial config
/* ----------------------------------------------------------- */
  scanInt(" init-k0 ", &n_class0);
  if (n_class0 == 1)		/* one common cluster */
    for(i=0;i<mdp.N;i++){
      member0[i] = 0;
      r0[i] = 0;
    }
  else if (n_class0 == mdp.N)	/* N idiosync clusters */
    for(i=0;i<mdp.N;i++){
      member0[i] = i;
      r0[i] = 1;
    }
  else if (n_class0 == mdp.J){	/* J idiosync clusters */
    for(i=0;i<mdp.N;i++){
      member0[i] = dta.j[i]-1;
      r0[i] = 1;
    }
  }
  else if (n_class0 == mdp.J+1){ /* J idiosync clusters, 1 common cl */
    for(i=0;i<mdp.N;i++){
      if (i%2 == 0){ /* common cluster */
	member0[i] = 0;
	r0[i] = 0;
      } else { /* idiosync cluster */
	member0[i] = dta.j[i];
	r0[i] = 1;
      }
    }
  }
  else{
    scanIntArray(" init-member0 ", member0, mdp.N);
    scanIntArray(" init-r0 ", r0, mdp.N);
  }
  closeIn();

/* -----------------------------------------------------------
   read in data
/* ----------------------------------------------------------- */

  /* consistency check */
  if ((mdl.p1 != mdp.p) | (dta.npa != mdp.N) ){
    printf("\n*** inconsistent mdl.p1/mdp.p or dta.npa/mdp.N!!\n");
    exit(1);
  }
    
  R_assgn_T(mdp.Y,dta.z,dta.npa,mdl.p1);
  mdp_initconfig(); /* after initializing the data ! */

  /* open files */
  mfile = openOut("mj.mdp");
  Vfile = openOut("Vt.mdp");
  Lfile = openOut("L.mdp");
  parFile = openOut("par.mdp");
  meanFile = openOut("mean.mdp");
  Sfile = openOut("S.mdp");
  Bfile = openOut("B.mdp");
}

/* ***********************************************************
   initial configuration
/* *********************************************************** */
void mdp_initconfig()
{
  int i, k;

  /* Note: values for aa and tau are required here! */
  mdp.n_class = 0;
  for (i = 0; i < mdp.N; i++){
    k = member0[i];
    if (mdp.count[k]==0) /* first obs in cluster */
      mdp_newclass(i,r0[i]);
    else
      mdp_addclass(i,k);
  }
  for (k=0;k<n_class0;k++){
    if (mdp.count[k]==0){
      printf("\n *** Error: class %d has no members.\n",k);
      exit(1);
    }
    mdp_setupmu(k,1);
  }
  mdp_checkclass();
}

/* ************************************************* 
   mdp_finish
 ************************************************* */ 
void mdp_finish()
{
  fclose(mfile);
  fclose(Vfile);
  fclose(Lfile);
  fclose(parFile);
  fclose(meanFile);
  fclose(Sfile);
  fclose(Bfile);
}

/* ***********************************************************
   Gibbs Sampler
/* *********************************************************** */
void mdp_gibbs()
{
  int
    start,k;

  /* update latent data */
  R_assgn_T(mdp.Y,dta.z,mdp.N,mdp.p);

  if (S_prior) mdp_sampleS(); 
  mdp_sampleconfig();
  for(k=0; k<mdp.n_class; k++){
    mdp_setupmu(k,1);
  } 
  if (alpha_prior) mdp_samplealpha(); 
  if (B_prior) mdp_sampleB();
  if (m_prior) mdp_samplem();
  if (g.iter > g.ndiscard){
    if (g.iter % g.npi == 0){
      mdp_writepars();
    }
    if (g.iter % n_predupdate == 0) /* update predictive */
      mdp_writepi();
  }

  /* UPDATE latent data */
  R_assgn_T(dta.z,mdp.Y,mdp.N,p);
}

/* ============================================================= 
   Conditional draws
/* ============================================================= */ 
	  
/* ***********************************************************
   sample config
/* ***********************************************************
   does the Gibbs over configurations
   Pr(y[i] in class j | ...) proportional to
   for j=0,..n_class-1:
      count[j]*MV-N(y[i]; mut[j], Vt[j])
      mut[j] and Vt[j] computed in make_mut_Vt
   for j=n_class (i.e. new class): same expression with 
                                   count[n_class] = alpha
*/
void mdp_sampleconfig()
{
  int 
    ji, jl, l, k, k_new, i;
  double 
    *alpha,eps,kd,u,py;

  kd = 1.0*mdp.n_class;
  alpha = mdp.alpha;
  eps = mdl.eps;

  for (i = 0; i < mdp.N; i++){
    x_zero(mdp.qv,mdp.n_class+2);
    mdp_takeout(i);

    /* compute prob pi[i] = pi*[k]       */
    for(l=0;l<mdp.n_class;l++){
      jl = mdp.j[l];
      if (jl == dta.j[i]){ /* idiosync cluster of same  Fj */
	if (mdp.updated[l]==0)
	  mdp_setupmu(l,0);
	mdp.qv[l] = mdp.count[l] * eps/(alpha[jl]+mdp.n[jl])*
	  mvn_pdf(mdp.Y[i],mdp.mj[l],mdp.Vt_cd_inv[l],p);
      } else if (jl == 0){   /* other common cluster */
	if (mdp.updated[l]==0)
	  mdp_setupmu(l,0);
	mdp.qv[l] = mdp.count[l] * (1-eps)/(alpha[0]+mdp.n[0])*
	  mvn_pdf(mdp.Y[i],mdp.mj[l],mdp.Vt_cd_inv[l],p);
      }
    }

    /* Prob of new class */
    /* qv[n_class = alpha * MV-N(y[i],m,V+tau*B)*/
    py = mvn_pdf(mdp.Y[i],mdp.m, mdp.St_cd_inv,p);
    /* new idiosync cluster */
    ji = dta.j[i];
    mdp.qv[mdp.n_class] = alpha[ji]*eps/(alpha[ji]+mdp.n[ji])*py;
    /* new common cluster */
    mdp.qv[mdp.n_class+1] = alpha[0]*(1-eps)/(alpha[0]+mdp.n[0])*py;

    /* draw the index of the new class for ind */
    multinomial(1,mdp.n_class+2,mdp.qv,&k_new);     
    if (k_new == mdp.n_class){ /* create a new idiosync class */
      mdp_newclass(i,1);
    } else if (k_new == mdp.n_class+1){ /* create a new common class */
      mdp_newclass(i,0);
    } else{ 
      /* add to existing class */
      mdp_addclass(i,k_new);
    }
  }
  mdp_checkclass();

}

/* -----------------------------------------------------------
   check class
/* ----------------------------------------------------------- */
/* checks for bogous classes */

void mdp_checkclass(){
  int i,n,k;

  for(k=0;k<mdp.n_class;k++){
    for(n=0,i=mdp.first[k]; i!=-1; n++, i = mdp.next[i]){
      if (n>mdp.count[k])
	printf("\n*** Error: class %d is bogous (too big).\n", k);
      if ( (mdp.j[k] > 0) & (dta.j[i] != mdp.j[k]) )
	printf("\n*** Error: class %d is bogous (inconsistent j).\n",
	       k);
      if ( (mdp.j[k] == 0) & (mdp.r[i] != 0) )
	printf("\n*** Error: class %d is bogous (inconsistent j/r).\n",
	       k);
    }
    if (n != mdp.count[k])
	printf("\n*** Error: class %d is bogous (too small).\n", k);
  }
}

/* -----------------------------------------------------------
   mdp_takeout
/* ----------------------------------------------------------- */
void mdp_takeout(int ind)
{
  int k_old,i,j_old,k;
  

  /* take y[ind] out of it's current class */
  k_old = mdp.member[ind];
  j_old = mdp.j[k_old];
  mdp.updated[k_old] = 0;

  mdp.count[k_old] -= 1;
  mdp.n[j_old] -= 1;
  mdp.member[ind] = -1;
  if (mdp.prev[ind] == -1)
    mdp.first[k_old] = mdp.next[ind];
  else mdp.next[mdp.prev[ind]] = mdp.next[ind];
  mdp.next[ind] = -1;
  mdp.prev[ind] = -1;
  mdp.r[ind] = -1;

  /* if class k=k_old is empty */
  if (mdp.count[k_old] == 0){
    k = mdp.n_class-1;

    /* remove empty class form list */
    free_dvector(mdp.mu[k_old],0,p-1);
    free_dvector(mdp.mj[k_old],0,p-1);
    free_dmatrix(mdp.Vt_cd_inv[k_old],0,p-1,0,p-1);
    free_dmatrix(mdp.Vt[k_old],0,p-1,0,p-1);

    /* relabel class k=n_class-1 as k_old */
    mdp.count[k_old] = mdp.count[k];
    mdp.first[k_old] = mdp.first[k];
    mdp.j[k_old] = mdp.j[k];
    mdp.count[k] = 0;
    mdp.first[k] = -1;
    mdp.j[k] = -1;
    mdp.new[k_old] = mdp.new[k];

    mdp.mu[k_old] = mdp.mu[k];
    mdp.mj[k_old] = mdp.mj[k];
    mdp.Vt_cd_inv[k_old] = mdp.Vt_cd_inv[k];
    mdp.Vt[k_old] = mdp.Vt[k];
    mdp.updated[k_old] = mdp.updated[k];

    /* update member-references to the relabeled class */
    for(i=mdp.first[k_old];i!=-1;i=mdp.next[i])
      mdp.member[i] = k_old;

    /* decrement n-class by one */
    mdp.n_class -= 1;
    mdp.K[j_old] -= 1;
    
    mdp_checkclass();
  }
}
/* -----------------------------------------------------------
   mdp_addclass
/* -----------------------------------------------------------
   adds obs[ind] to class[k]
*/
void mdp_addclass(int ind, int k)
{
  int j;

  if (mdp.first[k]== -1)
    {
      error("add-class", "tried to add to empty class...", 1);
    }
  j = mdp.j[k];
  if ( (j > 0) & (j != dta.j[ind]) ){
    printf(
    "\n *** Error: tried to add obs %d (j=%d) to class %d (j=%d).\n",
    ind,dta.j[ind],k,j);
    exit(1);
  }
  mdp.n[j] += 1;
  mdp.count[k] += 1;
  mdp.member[ind] = k;
  mdp.next[ind] = mdp.first[k];
  mdp.first[k] = ind;
  mdp.prev[mdp.next[ind]] = ind;
  mdp.prev[ind] = -1;
  mdp.r[ind] = (j==0) ? 0 : 1;
  
  mdp.updated[k] = 0;
  mdp_checkclass();
}

/* ***********************************************************
   new-class
/* *********************************************************** 
  open a new class with y[ind] as first member 
*/
void mdp_newclass(int ind, int r)
{ int k, j;

  k = mdp.n_class;
  j = (r==0) ? 0 : dta.j[ind];

  /* allocate memory for the new class */
  mdp.mu[k] = dvector(0,p-1);
  mdp.mj[k] = dvector(0,p-1);
  mdp.Vt_cd_inv[k] = dmatrix(0,p-1,0,p-1);
  mdp.Vt[k] = dmatrix(0,p-1,0,p-1);

  /* update count and member */
  mdp.count[k] = 1;
  mdp.first[k] = ind;
  mdp.new[k] = g.iter;
  mdp.member[ind] = k;
  mdp.next[ind] = -1;
  mdp.prev[ind] = -1;
  mdp.r[ind] = r;
  mdp.j[k] = j;

  mdp.n_class += 1;
  mdp.K[j] += 1;
  mdp.n[j] += 1;
  mdp.updated[k] = 0;
  mdp_checkclass();
}

/* -----------------------------------------------------------
   mdp_ybar
/* ----------------------------------------------------------- */
void mdp_ybar(int k, double *ybar)
{
  int i;
  double nk;

  nk = 1.0*mdp.count[k];
  x_zero(ybar,p); /* initialize ybar=0 */
  for(i=mdp.first[k];i!=-1;i=mdp.next[i]){
    x_plus_y(ybar,mdp.Y[i],ybar,p);
  }
  x_div_r(ybar,nk,ybar,p);
}
/* -----------------------------------------------------------
   mdp_Sk
/* ----------------------------------------------------------- */
void mdp_Sk(int k, double **S)
{
  int i;
  double nk, *z;

  /* allocate memory */
  z = dvector(0,p-1);

  nk = 1.0*mdp.count[k];
  A_zero(S,p);
  for(i=mdp.first[k];i!=-1;i=mdp.next[i]){
    x_min_y(mdp.Y[i],mdp.mu[k],z,p);
    A_plus_rxxt(S,1.0,z,S,p);
  }

  free_dvector(z,0,p-1);
}


/* *************************************************
   mdp_setupmu
 ************************************************* */
void mdp_setupmu(int k, int draw)
{
  double **V_inv, *ybar, nk;
  double **Vj, **Vj_inv;


  /* allocate memory for the auxilary arrays */
  ybar = dvector(0,p-1);
  Vj_inv = dmatrix(0,p-1,0,p-1);
  Vj = dmatrix(0,p-1,0,p-1);

  if (mdp.count[k] == 0){
    printf("\n*** tried to update mu for empty class (k=%d)\n",k);
    exit(1);
  }
  /* get pars for conditional posterior p(mu|...) */
  mdp_ybar(k, ybar);
  nk = 1.0*mdp.count[k];
  
  /* compute posterior mean cov matrix for mu[k] */
  nn_bayes(1.0,mdp.B_inv,mdp.m, 1.0/nk, mdp.S_inv, ybar,
	   Vj, Vj_inv, mdp.mj[k], p);


  /* if draw, then sample mu[k] */
  if (draw==1)
    mvnS_rand(p,mdp.mj[k],Vj,mdp.mu[k]);

  /* Vt = (Vj + V[k]) */
  rA_plus_sB(1.0,Vj,1.0,mdp.S,mdp.Vt[k],p); 
  invcd(p,mdp.Vt[k],mdp.Vt_cd_inv[k]);

  mdp.updated[k] = 1;

  /* free memory */
  free_dvector(ybar, 0,p-1);
  free_dmatrix(Vj_inv,0,p-1,0,p-1);
  free_dmatrix(Vj,0,p-1,0,p-1);
}

/* ***********************************************************
   mdp_sampleS
/* *********************************************************** */
/* S from W(1/S; q+n, (qR + sum S[k])^-1),
   where S[k] = sum (y[i]-mu[k])(y[i]-mu[k])'
*/
void mdp_sampleS()
{ double **Sk, **S1, **S1_inv_cd,
    **S_cd;

  int j;

  /* allocate memory */
  S_cd = dmatrix(0,p-1,0,p-1);
  Sk = dmatrix(0,p-1,0,p-1);
  S1 = dmatrix(0,p-1,0,p-1);
  S1_inv_cd = dmatrix(0,p-1,0,p-1);
  
  /* Sk = sum (y[i]-mu[j])(y[i]-mu[j]) over member[i]=j */
  A_zero(S1,p);
  for(j=0;j<mdp.n_class;j++){
    mdp_Sk(j,Sk);
    A_plus_B(S1,Sk,S1,p,p);
  }

  /* S1 = s*S+Sk */
  rA_plus_sB(mdp.q*1.0, mdp.R, 1.0, S1, S1, p);
  cdinv(p,S1,S1_inv_cd);
  wishart(mdp.q+mdp.N, p, S1_inv_cd, mdp.S_inv);
  
  /* compute V, V_cd, V_inv_cd etc. */
  d_chol_decomp(p,mdp.S_inv, mdp.S_cd_inv);
  d_inv_triang(p,mdp.S_cd_inv,S_cd);
  ABt(S_cd,p,p,S_cd,p,p,mdp.S);

  /* compute Vt_cd_inv */
  rA_plus_sB(1.0,mdp.S,1.0,mdp.B,mdp.St,p);
  invcd(p,mdp.St,mdp.St_cd_inv);

  /* allocate memory for the auxilary arrays */
  free_dmatrix(S_cd,0,p-1,0,p-1);
  free_dmatrix(S1,0,p-1,0,p-1);
  free_dmatrix(Sk,0,p-1,0,p-1);
  free_dmatrix(S1_inv_cd,0,p-1,0,p-1);
}


/* ***********************************************************
   sample m
/* *********************************************************** 
/* m from p(m|..) = N(m; m(mubar), T)
   m(ybar[k]) = T*(B-inv/tau*k*mubar + A-inv*a)
   T-inv = B-inv/tau*k + A-inv
*/
void mdp_samplem()
{
  double *mubar, K;
  int i;

  K = mdp.n_class; /* need double */

  /* allocate memory */
  mubar = dvector(0,p-1);
  
  /* make mubar = 1/K sum mu[k] */
  x_zero(mubar,p); /* initialize mubar=0 */
  for(i=0;i<mdp.n_class;i++){
    x_plus_y(mubar,mdp.mu[i],mubar,p);
  }
  x_div_r(mubar,K,mubar,p);

  /* compute and sample posterior on m */
  nn_bayes_rand(1.0,mdp.A_inv, mdp.aa,
		1.0/K, mdp.B_inv, mubar,
		mdp.m, p);

  /* free memory */
  free_dvector(mubar,0,p-1);
}

/* ***********************************************************
   draw B
/* *********************************************************** */
/* dummy */
/* sample B from p(B|..) =
   W(1/B; k+c, [cC + 1/tau*sum (mu[k]-m)(mu[k]-m)']^-1)
*/
void mdp_sampleB()
{ double **Sk, **S1, **S1_inv_cd,
    **B_cd, *z;
  int k;

  /* allocate memory */
  B_cd = dmatrix(0,p-1,0,p-1);
  Sk = dmatrix(0,p-1,0,p-1);
  S1 = dmatrix(0,p-1,0,p-1);
  S1_inv_cd = dmatrix(0,p-1,0,p-1);
  z = dvector(0,p-1);
  
  /* Sk = sum (mu[k]-m)(mu[k]-m)' */
  A_zero(Sk,p);
  for(k=0;k<mdp.n_class;k++){
    x_min_y(mdp.mu[k],mdp.m,z,p);
    A_plus_rxxt(Sk,1.0,z,Sk,p);
  }

  /* S1 = s*S+Sk */
  rA_plus_sB(mdp.c*1.0, mdp.C, 1.0, Sk, S1, p);
  cdinv(p,S1,S1_inv_cd);
  wishart(mdp.c+mdp.n_class, p, S1_inv_cd, mdp.B_inv);
  
  /* compute B, B_cd, B_inv_cd etc. */
  d_chol_decomp(p,mdp.B_inv, mdp.B_cd_inv);
  d_inv_triang(p,mdp.B_cd_inv,B_cd);
  ABt(B_cd,p,p,B_cd,p,p,mdp.B);

  /* error check */
  if (isnan(mdp.B[0][0])){
    error("mdp_sampleB", "B is NaN.", 1);
  }

  /* allocate memory for the auxilary arrays */
  free_dmatrix(B_cd,0,p-1,0,p-1);
  free_dmatrix(S1,0,p-1,0,p-1);
  free_dmatrix(Sk,0,p-1,0,p-1);
  free_dmatrix(S1_inv_cd,0,p-1,0,p-1);
  free_dvector(z,0,p-1);
}



/* ***********************************************************
   sample alpha
/* *********************************************************** 
   draw alpha and eta from:
   alpha from     pi*Gamma(a0+n_class, b0-log(eta)) +
              (1-pi)*Gamma(a0+n_class-1, b0-log(eta)),
         where pi = (a+n-class-1)/(a0+n_class-1 + n*b0 - n*log(eta)).
   eta   from Beta(alpha+1,n).
*/
void mdp_samplealpha()
{
     double a,  b, pi;       /* parameters in Gamma call */
     double u;
     int j;

     for (j=0;j<=mdp.J;j++){
       if (mdp.n[j]==0){ /* no obs alloc to measure j */
	 a = mdp.a0;
	 b = mdp.b0;
       } else {
	 /* draw eta */
	 mdp.eta[j] = betadev(mdp.alpha[j]+1.0, mdp.n[j]*1.0);

	 /* compute pars for gammas for alpha */
	 b = mdp.b0 - log(mdp.eta[j]);
	 pi = (mdp.a0 + mdp.K[j] - 1)/
	   (mdp.a0 + mdp.K[j] -1 + mdp.n[j]*(mdp.b0 - log(mdp.eta[j])));
	 duniform(1,&u);
 	 a = (u<pi) ? mdp.a0+mdp.K[j] : mdp.a0+mdp.K[j]-1;
	 /* draw alpha with prob pi from Gamma(a1,b) else Gamma(a2,b) */
       }
       mdp.alpha[j] = gamdev_ab(a,b);
     }
}

/* *************************************************
   sample_eps
 ************************************************* */
void sample_eps()
{ /* samples eps from conditional posterior 
     p(eps|..) = 0 with prob c*pi0*I(n0=N)
               = 1 with prob c*pi1*I(n1=N)
	                                            B(ae*,be*)
               ~ Be(ae*,be*), with pr c*(1-pi0-pi1)*----------
	                                            B(ae, be )
     where ae*=ae+n1, be*=be+n0, n1=Sum r[ji], n0=N-n1
     */
  int 
    N,n1,n0,i,j,k;
  double 
    q[3],u,ae1,be1,ae,be;
  N = mdp.N;
  
  /* compute n0, n1 */
  for(n0=n1=0,j=1;j<dta.J;j++)
    for(i=0;i<dta.n[j];i++)
      if (mdp.r[i]==1)
	n1++;
  n0 = N-n1;	
  
  /* compute q0,q1,q2 */
  ae = mdl.ae;
  be = mdl.be;
  mdl.ae1 = ae1 = ae+n1;     
  mdl.be1 = be1 = be+n0;
  if ( (n0==N) | (n1==N) ){
    q[0] = (n0==N) ? mdl.pe0 : 0.0;
    q[1] = (n1==N) ? mdl.pe1 : 0.0;
    q[2] = (1-mdl.pe0-mdl.pe1)*
      exp( lgamma(ae1)+lgamma(be1)-lgamma(ae1+be1)-   /* B(ae1,be1) */
	  (lgamma(ae )+lgamma(be )-lgamma(ae +be ))); /* B(ae, be ) */
    multinomial(1,3,q,&k);
    if (k==0)
      mdl.eps = 0;
    else if (k==1)
      mdl.eps = 1;
    else 
      mdl.eps = betadev(ae1,be1);
  } else
    mdl.eps = betadev(ae1,be1);
}



/* *************************************************
   write pars
 ************************************************* */
void mdp_writepi()
{
  int i,j,k,jj;
  double
    w,eps,*alpha;

  eps = mdl.eps;
  alpha = mdp.alpha;

  if (line % 100 == 0){
    fclose(mfile);
    fclose(Vfile);
    fclose(Lfile);	
    fclose(Sfile);
    fclose(Bfile);
    mfile = openAppend("mj.mdp");
    Vfile = openAppend("Vt.mdp");
    Lfile = openAppend("L.mdp");
    Sfile = openAppend("S.mdp");
    Bfile = openAppend("B.mdp");
  }

  for (jj=0;jj<=mdp.J+1;jj++){
    if (jj==0) w = alpha[0]*(1-eps)/(alpha[0]+mdp.n[0]);
    else if (jj==mdp.J+1) w = eps;
    else w = alpha[jj]*eps/(alpha[jj]+mdp.n[jj]);
    /* write m */
    fprintf(mfile,"%d %d %5.3f ",g.iter, jj, w);
    for(i=0;i<mdp.p;i++)
      fprintf(mfile, " %4.2f ",mdp.m[i]);
    fprintf(mfile,"\n");
    /* write V */
    fprintf(Vfile,"%d %d %4.3f ", g.iter,jj,w);
    for(i=0;i<mdp.p;i++){
      for(k=0;k<mdp.p;k++)
	fprintf(Vfile," %7.5e ",mdp.St[i][k]);
      fprintf(Vfile,"\n\t");
    }
    fprintf(Vfile,"\n");
  }
  for(j=0;j<mdp.n_class;j++){
    jj = mdp.j[j];
    w = (jj==0) ?
      mdp.count[j] * (1-eps)/(alpha[0]+mdp.n[0]) :
      mdp.count[j] * eps/(alpha[jj]+mdp.n[jj]);
    /* write m */
    fprintf(mfile,"%d %d %5.3f ", g.iter, jj, w);
    for(i=0;i<mdp.p;i++)
      fprintf(mfile, " %4.2f ",mdp.mj[j][i]);
    fprintf(mfile,"\n");
    /* write V */
    fprintf(Vfile,"%d %d %4.3f ", g.iter,jj,w);
    for(i=0;i<mdp.p;i++){
      for(k=0;k<mdp.p;k++)
	fprintf(Vfile," %7.5e ",mdp.Vt[j][i][k]);
      fprintf(Vfile,"\n\t");
    }
    fprintf(Vfile,"\n");
  }

/*  Vfile = openAppend("L.mdp");
  fprintf(mfile,"%d %7.5e ",g.iter,mdp.alpha);
  for(i=0;i<mdp.p;i++){
    for(k=0;k<mdp.p;k++)
      fprintf(Vfile," %7.5e ",mdp.St_cd_inv[i][k]);
    fprintf(Vfile,"\n\t");
  }
  fprintf(Vfile,"\n");
  for(j=0;j<mdp.n_class;j++){
    fprintf(mfile,"%d %d ",g.iter,mdp.count[j]);
    for(i=0;i<mdp.p;i++){
      for(k=0;k<mdp.p;k++)
	fprintf(Vfile," %7.5e ",mdp.Vt_cd_inv[j][i][k]);
      fprintf(Vfile,"\n\t");
    }
    fprintf(Vfile,"\n");
  }
  closeOut();
*/

  fprintf(Sfile,"%4d \n", g.iter);
  for(i=0;i<p;i++){
    for(j=0;j<p;j++)
      fprintf(Sfile," %7.5e",mdp.S[i][j]);
    fprintf(Sfile,"\n");
  }

  if (B_prior){
    fprintf(Bfile,"%4d \t", g.iter);
    for(i=0;i<p;i++){
      for(j=0;j<p;j++)
	fprintf(Bfile," %7.5e",mdp.B[i][j]);
      fprintf(Bfile,"\n");
    }
  }
}

/* *************************************************
/*  mdp_writepars 
 ************************************************* */
void mdp_writepars()
{
  int i,j;
  
  line += 1;
  if (line % 100 == 0){
    fclose(parFile);
    parFile = openAppend("par.mdp");
    fclose(meanFile);
    meanFile = openAppend("mean.mdp");
  }
  fprintf(parFile,"%4d %4.3f %4.3f %4.3f ", 
	  g.iter, mdl.eps, mdl.ae1, mdl.be1);
  for(j=0;j<=mdp.J;j++)
    fprintf(parFile," %6.3f ", mdp.alpha[j]);
  fprintf(parFile," %4d %5.3f\n", mdp.n_class,mdl.sig);

  fprintf(meanFile,"%4d \t", g.iter);
  for(i=0;i<p;i++)
    fprintf(meanFile,"%4.2f ",mdp.m[i]);
  fprintf(meanFile,"\n");
}

