#include <cstddef>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <math.h>
#include <stdio.h>

extern "C" {

  struct my_f_params {
    double x;
    double delta;
    double tau;
  };

  double integrand(double eps, void * p) {
    struct my_f_params *params = (struct my_f_params *)p;
    double x = (params->x);
    double delta = (params->delta);
    double tau = (params->tau);

    double res =
      (1 - (delta / (2 * delta - 1)) * pow((x - eps), (-(1 - delta) / delta)) +
       ((1 - delta) / (2 * delta - 1)) * pow((x - eps), -1)) *
      exp(-eps * eps / (2 * tau * tau)) / (tau * sqrt(2 * M_PI));
    return res;
  }

  double f_gsl(double x, double delta, double tau) {
    double result = 0.0;
    double error = 0.0;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
    gsl_set_error_handler_off();
    struct my_f_params params = {x, delta, tau};
    int err = 0;

    gsl_function F;
    F.function = &integrand;
    F.params = &params;

    double sig_pts[2];
    if (x-1 <= -100) {
      sig_pts[0] = x - 3;
      sig_pts[1] = x - 2;
    } else if (x-1 <= 100) {
      sig_pts[0] = -100;
      sig_pts[1] = x - 2;
    } else {
      sig_pts[0] = -10;
      sig_pts[1] = 100;
    }
    double pts[4] = {-1000000.0, sig_pts[0], sig_pts[1], x-1.0};

    err = gsl_integration_qagp(&F, pts, 4,
                               1.49e-8, 1.49e-9, 50,
                               w, &result,&error);
    // Error handling: decrease epsabs if diverging
    if (err == GSL_EDIVERGE) {
      printf("Quadrature didn't converge.\n");
      exit(EXIT_FAILURE);
    }

    gsl_integration_workspace_free(w);
    return result;
  }
}
