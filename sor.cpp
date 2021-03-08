// SOR Gauss-Seidel method for 1D Poisson eqn
// $ ./sor -n 100 -NIT 500 -omega 1.94
// $ ./sor -n 1000 -NIT 5000 -omega 1.994
// $ ./sor -n 10000 -NIT 50000 -omega 1.9994
/*
      matrix-free SOR method for (non-symmetric) Gauss-Seidel iteration

            u_xx = f(x) = -1, 0 < x < 1, u(0) = u(1) = 0

            u[0] = 0, u[n+1] = 0   [so allocate n+2 values]

            1 <= omega < 2 is (empirical) relaxation factor

            standard gauss-seidel if omega = 1
*/

#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include <math.h>

int main(int argc, char *argv[]) 
{
  Timer t;
  long n = read_option<long>("-n", argc, argv);
  long NIT = read_option<long>("-NIT", argc, argv, "1");
  float omega = read_option<float>("-omega", argc, argv, "2");
  long np  = n+1;
  long npp = n+2;

  printf("n = %ld, NITmax = %ld, \n",n, NIT);
  printf("over-relaxation factor = %12.5e \n", omega);

  long NITmax = NIT;   // do at most NITmax iterations
  long nitTot = NIT;   // actual number of iterations when rmax < rtol
  int  num_print = 15; // number of lines of output (iterations)
  long mod_print = NITmax/num_print;
  int  num_plot = 150; // number of data points to plot (iterations)
  long mod_plot = NITmax/num_plot;


  double emax = 1.0e-12;  // maximum error in solution [ans = x (1 - x) /2]
  double rmax = 1.0e-12;  // maximum residual |u_xx + f|
  double rtol = 1.0e-14;  // terminate when max residual is less than rtol

  double h = 1.0/np;   // mesh width for unit interval [0 < x < 1]
  double hsq = h*h; 
  printf("h = %12.5e, h^2 = %12.5e \n", h, hsq);

  t.tic();
  double* u    = (double*) malloc(npp * sizeof(double)); // unknown at current time step
  double* f    = (double*) malloc(npp * sizeof(double)); // right hand side
  double* ans  = (double*) malloc(npp * sizeof(double)); // exact solution

  printf("time to malloc = %12.5e s\n", t.toc());

  // initialize ...

  t.tic();
  for (long i = 0; i <= np; i++) 
  {
      double x = i*h;
      u[i]     = 0.0;
      f[i]     = 1.0;
      ans[i]   = 0.5*x*(1.0 - x);
  }
  printf("time to initialize = %12.5e s\n", t.toc());

  t.tic();

  printf("          nit         rmax         emax\n");

  // double omega = 1.2;   // over-relaxation factor > 1

  FILE *fptr;

  fptr = fopen("resid.txt","w");

  fprintf(fptr,"          nit         rmax         emax\n");

  t.tic();

  //  ************* main iteration loop ***********

  for (long nit = 1; nit <= NITmax; nit++) 
  { 
     double emax = 0.0;
     double rmax = 0.0;
     for (long i = 1; i <= n; i++) 
     {
         double Ru = (u[i+1] - 2*u[i] + u[i-1])/hsq + f[i];

         // GS is: solve Ru = 0 for new u[i] = 0.5*(u[i+1] + u[i-1] + hsq * f[i])
         // use latest values in Ru as soon as they are updated (not like Jacobi)

         double utemp = u[i] + 0.5*hsq * Ru;      // provisional Gauss-Seidel update (pre-SOR)
         u[i] = omega*utemp + (1.0 - omega)*u[i]; // overwrite old by new value ("lexicographic")
         rmax = fmax(hsq * fabs(Ru),rmax);
         emax = fmax(fabs(u[i] - ans[i]),emax);
     }

     if (nit%mod_plot == 0) 
     {
        fprintf(fptr,"%13ld %12.5e %12.5e \n", nit, rmax, emax); 
     }

     if (nit == 1 || nit%mod_print == 0)
     {
         printf("%13ld %12.5e %12.5e \n", nit, rmax, emax);
     }
     if (rmax <+ rtol)
     { 
        nitTot = nit;
        printf("%13ld %12.5e %12.5e \n", nit, rmax, emax);
        break;
     }
  }

  //  ************* end of main iteration loop ***********

  fclose(fptr);

  double time = t.toc();
  double time_perstep = time/nitTot;
  double time_perpt = time_perstep/np;
  printf("\n cpu_time (s) = %12.5e\n time/step    = %12.5e\n time/step/pt = %12.5e\n", 
          time, time_perstep, time_perpt);


  printf(" \n            i              x              u          error \n");

  for (long i = 0; i < np; i = i + n/10)
  {
     double x = i*h;
     double err = fabs(u[i] - ans[i]);
     printf("%13ld   %12.5e   %12.5e   %12.5e \n",i, x, u[i], err);
  }
  long i = np;
  double x = i*h;
  double err = abs(u[i] - ans[i]);
  printf("%13ld   %12.5e   %12.5e   %12.5e \n",i, x, u[i], err);
  
  t.tic();

  free(u);
  free(f);
  free(ans);

  printf("time to free = %12.5e s\n", t.toc());

  return 0;
}
