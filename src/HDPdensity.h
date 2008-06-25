/* =========================================================== */
/*							       */
/*  meta.h                 				       */
/*							       */
/* =========================================================== */

/* *************************************************
   Model structure 
 ************************************************* */
struct Model{
  int
    p1,            /* dim of random effect plus cov vector 
		      = n columns in data-z */
    nb,            /* n pars in logistic regr/loglog  */
    pltpa;         /* dummy for debugging purposes only */
  int
    model;         /* model indicator
		       0: logistic, 1: complementary log log */
  double 
    sig,           /* error variance */
    as,bs,	   /* IG hyperprior on sig2 */
    eps,	   /* fraction of idiosyncratic measure */
    pe0,pe1,ae,be,   /* hyperprior on eps */
    ae1,be1;	   /* beta par's for complete conditional of
		      eps - used for "rao-blackwellization" */
};

/* *************************************************
   Data structure
 ************************************************* */
struct Data {
  int 
    npa,	     /* total # patients */
    J,		     /* # studies */
    N,		     /* total # obs's */
    *n,		     /* n pat's for each study, j=1..J */
    d,               /* dimension of data vector at time t */
    p1;              /* dim of random effects per patient */
  double **y,        /* data matrix */
    **z;             /* imputed data */
  int 
    *j,		     /* study, j=1..J */
    *first,          /* hash table for first obs on patient i */
    *ni;             /* n obs per patient */
};

/* *************************************************
   Markov chain Monte Carlo pars
 ************************************************* */
struct MCMC {
  int 
    iter,             /* number of the current iteration */
    niter,            /* number iterations */
    nbatch,           /* batch siyes for each step of indep chain */
    npi,              /* size of batches for ergodic avges */
    ndiscard,         /* start accumulating MC sample after ncollect its */
    verbose,          /* be verbose about comments */
    seed,             /* seed for RV generator */
    sampleeps,	      /* indicator for resampling epsilon */
    t0;               /* par for each reparametrization */
  double
    c;                /* scale factor for Metropolis */
  long
    is1,is2;	      /* r.v. seeds */
};


/* ***********************************************************
   all MDP pars
/* *********************************************************** */

  struct MDP_PARS
{
  /* -----------------------------------------------------------
     Data & Parameters */
  /* ----------------------------------------------------------- */
  double 
     **Y;                           /* data = patient random eff's
				              and covariates */
  int 
     *n,                            /* n obs allocated to each
				       random measure, j=0..J */
     N,				    /* total # obs in the mdp */
     J,				    /* number of different random
				       measures (including meas F0)*/
     p;                             /* dimension of each obs vector */
  double  
     *qv,			    /* aux var: resampling prob's */
     **mu,                          /* mu[j], location of cluster j */
     **mj, ***Vt_cd_inv, ***Vt;
  double 
     **S, **S_cd_inv,**St,
     **S_inv, **St_cd_inv,          /* prior variance, chol decomp and 
				      inverse of chol decomp */
     **R, **R_inv;                  /* Matrix par of hyperrior on S */
  int 
     q;                             /* df of Wishhart hyperprior on S */
  double *eta,                      /* latent par for sampling of alpha */
     *alpha,                        /* toal mass par of DP */
     a0,b0;                         /* prior pars for alpha */
  double  
     *m,                            /* mean of base measure */
     *aa, **A_inv;                  /* hyperpars on mean */
  double 
     **B, **B_inv, **B_cd_inv,      /* cov matrix of base measure */
     **C;
  int c;
  
  /* -----------------------------------------------------------
     configuration */
  /* ----------------------------------------------------------- 
     The current configuration, i.e. number of distinct pi[i]'s 
     and classes of identical pi[i]'s is described by count, member
     and n_class.
     Notation: pi*[0],... pi*[n_class] is the list of distinct pi[i]'s.
     Of course n_class <= n_obs.
     n_class:                 number of distinct pi[i]'s.
     count[j]                 number of pi[i]'s equal to pi*[j]
     member[i]:               = j if pi[i] = pi*[j] 
     new[j]                   = r if class j was formed during
     iteration r (for debugging purposes) 
     first[j]                 = index of first obs in class j
     next[i]:                 = index of next obs in same class as obs i
     prev[i]:                 = index of previous obs in same class
     (prev[first] and next[last] = -1) 
     */ 
  int 
     n_class,			/* total # clusters */
     *K,			/* K[j] = # clusters with measure j */
     *count, *first,             
     *member, *next, *prev,
     *new,*updated,
     *j,			    /* j[i]=j iff cluster i is
				       associated with measure j */
     *r;			    /* r[i]=0 iff patient i is
				       allocated to common measure */
};

/* Indep chain */
void par_rand(double *z, double *z1, double *z2, double *z3,
	      double *t1, double *t2, double *b0, double *b1,
	      double *b2);
double et(double *y, double *z, struct Model *mdl, double day);
double etat(double *y, double *z, struct Model *mdl, double day);
double et2(double *y, double *z, struct Model *mdl);
double ft(double *y, double *z, struct Model *mdl, double day);


void init();
void finish();
void write_th();
void write_current_th();
void sample_z1();
void sample_z2();
void sample_z3();
void sample_sig();
void setup_b(int patid, double *mu, double **V_inv, double **L,
	     int *proper);
void sample_b(int nbatch, int *n_jump);
void sample_t12(int nbatch, int *n_jump);
void sample_eps();

void llikl_init(struct Data *dta);
double llikl(struct Model *mdl, struct Data *dta); 
double llikl_it(double *y, double *z, struct Model *mdl);
double llikl_i(int patid, struct Model *mdl, struct Data *dta, double *z);
double llikl_moments(double *y, double *z, struct Model *mdl, double *m,
	     double *s);

void write_th();
void write_current_th();
/* MDP */
void mdp_gibbs();
void mdp_init();
void mdp_initconfig();
void mdp_finish();
void mdp_sampleconfig();
void mdp_checkclass();
void mdp_takeout(int ind);
void mdp_addclass(int ind, int k);
void mdp_newclass(int ind, int r);
void mdp_ybar(int k, double *ybar);
void mdp_Sk(int k, double **S);
void mdp_setupmu(int k, int draw);
void mdp_sampleS();
void mdp_sampleB();
void mdp_samplealpha();
void mdp_samplem();


void mdp_writepars();
void mdp_writepi();
