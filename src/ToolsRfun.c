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


double F77_SUB(cdfbetas)(double *x, double *a, double *b, int *lower_tail, int *give_log)
{
	return pbeta(*x, *a, *b, *lower_tail, *give_log);
}


double F77_SUB(cdfchisq)(double *x, double *df, int *lower_tail, int *log_p)
{
	return pchisq(*x, *df, *lower_tail, *log_p);
}

double F77_SUB(cdfpoiss)(double *x, double *lambda, int *lower_tail, int *log_p)
{
	return ppois(*x, *lambda, *lower_tail, *log_p);
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


double F77_SUB(invcdfbetas)(double *p, double *a, double *b, int *lower_tail, int *log_p)
{
	return qbeta(*p, *a, *b, *lower_tail, *log_p);
}


double F77_SUB(invcdfchisq)(double *p, double *df, int *lower_tail, int *log_p)
{
	return qchisq(*p, *df, *lower_tail, *log_p);
}

double F77_SUB(invcdfpoiss)(double *p, double *lambda, int *lower_tail, int *log_p)
{
	return qpois(*p, *lambda, *lower_tail, *log_p);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//Density
///////////////////////////////////////////////////////////////////////////////////////////////////////


double F77_SUB(dbin)(double *x, double *n, double *p, int *give_log)
{
	return dbinom(*x, *n, *p, *give_log);
}


double F77_SUB(dnrm)(double *x, double *mu, double *sigma, int *give_log)
{
	return dnorm(*x, *mu, *sigma, *give_log);
}


double F77_SUB(dlnrm)(double *x, double *logmean, double *logsd, int *give_log)
{
	return dlnorm(*x, *logmean, *logsd, *give_log);
}


double F77_SUB(dlogit)(double *x, double *location, double *scale, int *give_log)
{
	return dlogis(*x, *location, *scale, *give_log);
}

double F77_SUB(dcauch)(double *x, double *location, double *scale, int *give_log)
{
	return dcauchy(*x, *location, *scale, *give_log);
}


double F77_SUB(dbet)(double *x, double *a, double *b, int *give_log)
{
	return dbeta(*x, *a, *b, *give_log);
}


double F77_SUB(dgamm)(double *x, double *shape, double *scale, int *give_log)
{
	return dgamma(*x, *shape, *scale, *give_log);
}

double F77_SUB(dpoiss)(double *x, double *lambda, int *give_log)
{
	return dpois(*x, *lambda, *give_log);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//Math functions
///////////////////////////////////////////////////////////////////////////////////////////////////////

double F77_SUB(dgamlog)(double *x)
{
	return lgammafn(*x);
}

double F77_SUB(lbetaf)(double *a, double *b)
{
	return lbeta(*a, *b);
}

double F77_SUB(trigamm)(double *x)
{
	return trigamma(*x);
}

double F77_SUB(powerxy)(double *x, double *y)
{
    return R_pow(*x, *y);
}

double F77_SUB(powerxi)(double *x, int *i)
{
    return R_pow_di(*x, *i);
}

double F77_SUB(log1px)(double *x)
{
    return log1p(*x);
}

double F77_SUB(log1pxmx)(double *x)
{
    return log1pmx(*x);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//Print on screen
///////////////////////////////////////////////////////////////////////////////////////////////////////

int F77_SUB(sprint)(int *iter, int *tot_iter, double *sec)
{
	Rprintf("\nMCMC scan %i of %i (CPU time: %.3f s)", *iter, *tot_iter, *sec);
	return *iter;
}

