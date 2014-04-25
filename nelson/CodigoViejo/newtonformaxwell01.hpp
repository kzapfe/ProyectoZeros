//Newtons Method for root searching
// Free ware robado de internet, modificado por karel

using namespace arma;

double newton(double x_0, double tol, int max_iters, 
          int* iters_p, int* converged_p) {
   double x = x_0;
   double x_prev;
   int    iter = 0;

   do {
      iter++;
      x_prev = x;
      x = x_prev - f(x_prev)/f_prime(x_prev);
   } while (fabs(x - x_prev) > tol && iter < max_iters);

   if (fabs(x - x_prev) <= tol)
      *converged_p = 1;
   else
      *converged_p = 0;
   *iters_p = iter;
   
   return x;
}  /* newton algorithm */


double f(double x) {
   return x*x-2;
}  /* f */

double f_prime(double x) {
   return 2*x; //the derivative
}  /* f_prime */

