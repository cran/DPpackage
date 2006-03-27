#include <R.h>
#include <Rmath.h>

//CDF

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

//Density


//Math functions

double F77_SUB(dgamlog)(double *x)
{
	return lgammafn(*x);
}


//Print on screen

int F77_SUB(sprint)(int *iter, int *tot_iter)
{
	Rprintf("\nMCMC scan %i of %i", *iter, *tot_iter);
	return *iter;
}



