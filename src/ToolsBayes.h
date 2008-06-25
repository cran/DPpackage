void nn_bayes(double, double **, double *, 
	      double, double **, double *,
	      double **, double **, double *, int);
void nn_bayes2(double **Spr_inv, double *mpr, 
	       double **Slik_inv, double *y,
	       double **Spo, double **Spo_inv, double *mpo, 
	       int p);
void nn_bayes_rand(double, double **, double *, 
		   double, double **, double *,
		   double *, int);
void nn_bayes1_rand(double, double, double , 
		   double, double , double ,
		   double *);
void nn_bayes1(double r1, double Spr_inv, double mpr, 
		  double r2, double Slik_inv, double y,
		  double *Spo, double *mpo);

