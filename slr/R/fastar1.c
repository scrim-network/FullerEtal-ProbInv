/*   Copyright 2010 Nathan M. Urban                                           */
/*   build with:  system("rm fastar1.o; R CMD SHLIB fastar1.c")               */

/*   Simulates AR(1) red noise                                                */

/*   n = length of time series to generate                                    */
/*   rhoin = (pointer to) lag-1 autocorrelation coefficient                   */
/*   sigmain = (pointer to) innovation standard deviation                     */
/*   w = time series array of iid standard normal white noise (of length n)   */
/*   r = time series array of AR1 red noise (of length n)                     */

#include <math.h>

void ar1sim(int* nin, double* rhoin, double* sigmain, double* w, double* r)
{
  double rho, sigma, sigma_proc;
  int n, i;
  
  n = nin[0]; rho = rhoin[0]; sigma = sigmain[0];

  sigma_proc = sigma/sqrt(1.0-rho*rho); /* process variance */
  r[0] = sigma_proc*w[0]; /* draw initial value from stationary distribution */

  for(i=1; i<n; i++)
    r[i] = rho*r[i-1] + sigma*w[i];
}
