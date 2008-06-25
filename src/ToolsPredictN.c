/* =========================================================== 
  							       
    predict.c						       
  							       
   =========================================================== */
/* compute predictive for new observations, including predicitives for
   random effects  */

/* shortened version from smooth/meta/predict.c
   without marginal for r.e.'s z1...b1 and nadir, tlo and fnadir */

#include "math.h"
#include "stdio.h"
#include <stdlib.h>
#include "ToolsInterface.h"
#include "ToolsVector.h"
#include "ToolsCMatrix.h"
#include "ToolsRand.h"
#include "ToolsNrutil.h"
#include "HDPdensity.h"


void init_predict();
void finish_predict();


/* =========================================================== 

		 main 

  =========================================================== */

/* dummy fct declarations */
char *nm(char *nm0);

/* *************************************************
   global vars
 ************************************************* */

int 
  n,      /* n MC sample points */
  p,      /* dim of par vector */
  px,	  /* dim of cov vector */
  nx;     /* n grid points on time axis */


double 


  *wij, /* weights for individual simulations */

  /* arrays to accumulate mean process and s.d. for
     predictive profiles */
  *xgrid, /* grid of days over which to evaluate profile */ 

  /* arrays to read in moments of clusters and 
     to compute cond moments */
  *mj,**Vt,**S22,**S22_inv,**S11,*m1,*z1 /* z_b */, *z0 /* z_c */,
  
  /* cov and random effects vector for new patient */
  *z,
  *mu, **V;
 
int
  /* idx of random effects which are conditioned upon etc. */
  idx[10],idx_c[10], nidx, nidx_c,

  *ok,			/* =1 if a term is used, =0 if not */
  iter, ij0,
  jpred,		/* index of subpopulation to be evaluated */
  rpred;  		/* if rpred=0, then mix with common measure
		   if rpred=1, don't mix with common measure */

char 
   name0[80], name1[80]; 
   /* aux to construct filenames */

FILE
  *xfile,
  *mjfile,*Vtfile,
  *zfile,
  *ctrlfile,
  *xgridfile;

long is1=981963, is2=123456789;
struct Model mdl;
struct Data dta;

/********************************************
         main
 ********************************************/
void predictN(int *jpredp, int *rpredp)
// main(int argc, char *argv[])
{
  int
    i, /* index of simulation */
    k,r, iter_last, tlo, h, 
    ii, /* index of repeated simulation */
    jj, jc, I=10,nk,
    j /* index of subpopulation */;
  double 
    wj, sumw, m, s, lp, wj2, pmarg,
    niter, u, maxw, eps /* mixing par for mixture of measures */;
  FILE
    *tmpfile;

  jpred=*jpredp; /* R passes pointers to integers! */
  rpred=*rpredp;
  /* for debugging as main program.... 
  if (argc != 3){
    printf("\n Usage 'predict j r'.\n");
    exit(1);
  }
  sscanf(argv[1]," %d ",&jpred);   
  sscanf(argv[2]," %d ",&rpred);    */


  init_predict();

  if ( (rpred != 0) & (rpred != 1) ){
    printf("\n *** Error: rpred=0 or 1 only.\n");
    exit(1);
  }
  if ( (jpred <0) | (jpred > dta.J+1) ){
    printf("\n *** Error: 0 <= jpred <= J+1 only.\n");
    exit(1);
  }
  /*  loop over all data points */
  for(r=0;r<dta.npa;r++){
    rewind(mjfile);
    rewind(Vtfile);
    /* fill in covariates & z1 */
    y_assgn_x(z,dta.z[r],p);
    printf("%d: ",r);

    nk =0;


    /* 1. compute weights */
    for(i=0,ij0=0, maxw = -999999.9,iter_last=0,niter=0;
	(i<n) & (fscanf(mjfile,"%d %d %lf ", &iter,&j,&wj) != EOF);
	i++){
      fscanDoubleArray(mjfile,mj,p);
      fscanf(Vtfile,"%d %d %lf ", &iter, &j, &wj2);
      fscanDoubleMatrix(Vtfile,Vt,p,p);
      if ( ((j>0) & (j!=jpred)) /* not the desired idiosync mesure */
	   | ((j==0) & (jpred != 0) & (rpred != 0)) ){
	/* no mixing with common measure, not jpred=0 */
	ok[i] = 0;
	wij[i] = -1.0;
	continue;
      }
      ok [i] = 1;
      wj = log(wj); /* keep on log scale until standardization */
      /* 1.1. standardize the terms at the current iteration
	 (ij0 records the index of the first term at this
	 iteration) */
      if ( (i>0) & (iter != iter_last) ){
	for( ii=ij0, sumw=0; ii<i; ii++ ){
	  if (ok[ii]==0){
	    continue;
	  }
	  wij[ii] = exp(wij[ii]-maxw);	/* change to absolute scale */
	  sumw += wij[ii];		
	}
	for( ii=ij0; ii<i; ii++ )		/* standardize */
	  if (ok[ii]==1)
	    wij[ii] /= sumw;
	ij0 = i;
	maxw =  -999999.9;
	niter += 1.0;
      }
      iter_last = iter;

      /* 1.2. compute marginal weight */
      dsubmatrix(Vt,p,p,S11,nidx_c,nidx_c,idx_c,idx_c);
      /* for(ii=0;ii<3;ii++)
	 printf("%6.3f ",sqrt(S11[ii][ii]));
	printf("\n"); */
      dsubvector(mj,p,m1,nidx_c,idx_c);
      dsubvector(z,p,z0,nidx_c,idx_c);
      pmarg = mvnS_logpdf(z0,m1,S11,nidx_c);
      wj += pmarg; /* still on log basis */
      if (wj>maxw) maxw = wj;
      wij[i] = wj; 
    } /* end for i */
    if (i<n) n=i; /* total # records on mj file */
    printf(" n=%d \n",n);
    /* change to absolute scale for last batch.. */
    for(ii=ij0, sumw=0; ii<i; ii++ ){
      if (ok[ii]==1){
	wij[ii] = exp(wij[ii]-maxw);	
	sumw += wij[ii];		
      }
    }
    for(ii=ij0; ii<i; ii++ )		/* standardize */
      if (ok[ii]==1)
	wij[ii] /= sumw;
    niter += 1.0;

    /* 2. post predictive inference for patient i=n+r */
    for (ii=0;ii<I;ii++){ /* repeat I times to get more simuulations */
      rewind(mjfile);
      rewind(Vtfile);
      
      /* compute densities etc */
      for(i=0,iter_last=0;i<n;i++){
	fscanf(mjfile,"%d %d %lf ",&iter, &j, &wj);
	iter_last = iter;
	
	fscanDoubleArray(mjfile,mj,p);
	fscanf(Vtfile,"%d %d %lf ", &iter, &j, &wj2);
	fscanDoubleMatrix(Vtfile,Vt,p,p);
	if (ok[i]==0)
	  continue;

	/* compute marginal weight */
	wj = wij[i]; 
	if (wj==0)
	  continue;

	/* compute cond moments for z2,z3,t,b */
	dsubmatrix(Vt,p,p,S22,nidx_c,nidx_c,idx_c,idx_c);
	d_invert(nidx_c,S22,S22_inv);
	cond_mvnmoments(p,nidx,z,mj,Vt,S22_inv,idx,idx_c,V,mu);
	mvnS_rand(nidx,mu,V,z1);
	for(jj=0;jj<nidx;jj++)  z[idx[jj]] = z1[jj];
	  
	u = runif();
	if (u<wj){
	  /* update predicted profile */
 	  fprintf(zfile,"%d \t",r); /* pat number */
	  for(jj=0;jj<p;jj++)
	    fprintf(zfile,"%5.3f ",z[jj]);
	  fprintf(zfile,"\n");
	}

	/* 3. Marginal posterior predictive dist for rand effects */
      } /* end for i */
      /* linefeed on nadir file etc */
    }/* for ii */
    fprintf(xfile,"%5.3f ",r);
    for(jc=0;jc<px;jc++)
      fprintf(xfile,"%5.3f ",dta.z[r][idx_c[jc]]);
    fprintf(xfile,"\n");
  } /* end for r*/
  

  finish_predict();
}

/* =========================================================== 
    
  		 initialization & finish
  
   =========================================================== */
/********************************************
           init
  *******************************************/

void  init_predict()
{
  int
    i,j,jj,jb,jc,incl;
  double
    from,to;
/* -------------------------------------------------
   allocate memory
 ------------------------------------------------- */
  setseed(is1,is2);

/* -------------------------------------------------
   read initial values & MC points
 ------------------------------------------------- */
  trim("init.p","init-trim.p");
  ctrlfile = openIn("init-trim.p");
  fscanInt(ctrlfile," n ",&n);	 /* n MCMC simulations to use */
  fscanInt(ctrlfile," p ",&p);	 /* dim of rand eff vector (?) */
  /* jpred and rpred are now command line arg's ... 
     fscanInt(ctrlfile," jpred ",&jpred); 
     jpred indicates the random meas
     for which pred inf is desired 
     fscanInt(ctrlfile," rpred ",&rpred); 
     if rpred=0, mix with common measure
     if rpred=1, don't mix with common m. */
  wij = dvector(0,n);
  ok = ivector(0,n);

  
  /* dta.z will be used to store the 
     covariates of the predicted cases */
  scanInt(" npa ",&dta.npa);
  scanInt(" px ",&px);
  scanIntArray(" idx_c ",idx_c, px);

/* -------------------------------------------------
   allocate memory
 ------------------------------------------------- */
  dta.z =    dmatrix(0,dta.npa,0,p);

/* -------------------------------------------------
   read *fake* data
 ------------------------------------------------- */
  scanDoubleMatrix(" data-z ",dta.z, dta.npa, p);
  closeIn();

/* -------------------------------------------------
   alloc mem
 ------------------------------------------------- */
  /* indices for cond moments */
  for(jj=0,jb=0;jj<p;jj++){
    incl=1;
    for(jc=0;jc<px;jc++)
      if (idx_c[jc]==jj) incl=0;
    if (incl==1){
      idx[jb]=jj;
      jb++;
    }
  }//jj
  nidx_c = px; /*z1,ctx and gm (and wr) */
  nidx = p-nidx_c;

  /* allocate mem */
  z = dvector(0,p);
  mj = dvector(0,p);
  Vt = dmatrix(0,p,0,p);
  mu = dvector(0,p);
  V = dmatrix(0,p,0,p);

  S22 = dmatrix(0,nidx_c,0,nidx_c);
  S22_inv = dmatrix(0,nidx_c,0,nidx_c);
  S11 = dmatrix(0,nidx,0,nidx);
  m1 = dvector(0,nidx);
  z1 = dvector(0,nidx);
  z0 = dvector(0,nidx_c);

/* -------------------------------------------------
   open files 
 ------------------------------------------------- */

  ctrlfile = openIn("init.p");

  xfile = openOut(nm("x"));
  zfile = openOut(nm("z"));


  /* loop over all MC points */
  mjfile = openIn("mj.mdp");
  Vtfile = openIn("Vt.mdp");

}

char *nm(char *nm0){
  sprintf(name1,"%s-%d%d.p",nm0,jpred,rpred);
  return(name1);
}


/* *************************************************
   init dta and mdl
 ************************************************* */


/* *************************************************
   finish
 ************************************************* */
void finish_predict()
{
  fclose(  xfile );
  fclose(  zfile );
  fclose( mjfile );
  fclose( Vtfile );
}
