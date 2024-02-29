#include <stdio.h>
#include <math.h>

extern "C"
{
  double integrand(int n, double *params){
    double eps = params[0];
    double x   = params[1];
    double delta = params[2];
    double tau = params[3];
    double res =
      (1 - (delta / (2*delta - 1)) *
       pow((x - eps), (-(1 - delta) / delta)) +
       ((1 - delta) / (2*delta - 1)) * pow((x - eps), -1)) *
      exp(-eps*eps/(2*tau*tau))/(tau*sqrt(2*M_PI));
    return res;
  }
}
