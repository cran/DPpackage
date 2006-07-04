#include <R.h>
#include <Rmath.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////
//CDF
///////////////////////////////////////////////////////////////////////////////////////////////////////

double F77_SUB(cdfcauchy)(double *x, double *location, double *scale, int *lower_tail, int *give_log)
{
	return pcauchy(*x, *location, *scale, *lower_tail, *give_log);
}

double F77_SUB(cdfnorm)(double *x, double *mu, double *sigma, int *lower_tail, int *give_log)
{
	return pnorm(*x, *mu, *sigma, *lower_tail, *give_log);
}

double F77_SUB(cdflogis)(double *x, double *location, double *scale, int *lower_tail, int *give_log)
{
	return plogis(*x, *location, *scale, *lower_tail, *give_log);
}

double F77_SUB(cdflnorm)(double *x, double *logmean, double *logsd, int *lower_tail, int *give_log)
{
	return plnorm(*x, *logmean, *logsd, *lower_tail, *give_log);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//Quantile
///////////////////////////////////////////////////////////////////////////////////////////////////////


double F77_SUB(invcdflogis)(double *p, double *location, double *scale, int *lower_tail, int *log_p)
{
	return qlogis(*p, *location, *scale, *lower_tail, *log_p);
}


double F77_SUB(invcdfcauchy)(double *p, double *location, double *scale, int *lower_tail, int *log_p)
{
	return qcauchy(*p, *location, *scale, *lower_tail, *log_p);
}


double F77_SUB(invcdfnorm)(double *p, double *mu, double *sigma, int *lower_tail, int *log_p)
{
	return qnorm(*p, *mu, *sigma, *lower_tail, *log_p);
}


int F77_SUB(invcdfbinom)(double *p, int *size, double *prob, int *lower_tail, int *log_p)
{
	return qbinom(*p, *size, *prob, *lower_tail, *log_p);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//Density
///////////////////////////////////////////////////////////////////////////////////////////////////////

double F77_SUB(dnrm)(double *x, double *mu, double *sigma, int *give_log)
{
	return dnorm(*x, *mu, *sigma, *give_log);
}


double F77_SUB(dlnrm)(double *x, double *logmean, double *logsd, int *give_log)
{
	return dlnorm(*x, *logmean, *logsd, *give_log);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//Math functions
///////////////////////////////////////////////////////////////////////////////////////////////////////

double F77_SUB(dgamlog)(double *x)
{
	return lgammafn(*x);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//Print on screen
///////////////////////////////////////////////////////////////////////////////////////////////////////

int F77_SUB(sprint)(int *iter, int *tot_iter, double *sec)
{
	Rprintf("\nMCMC scan %i of %i (CPU time: %.3f s)", *iter, *tot_iter, *sec);
	return *iter;
}

