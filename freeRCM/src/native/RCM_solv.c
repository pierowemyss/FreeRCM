//
// Python-callable C script to draw a residue curve
// Author: Piero Wemyss
// Date: January 16, 2025
//
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern void nrtl_(double *x, double *Tcel, const double *nrtlA,
                  const double *nrtlB, const double *nrtlC, const int *Ncomps,
                  double *gamma);
extern void srk_(double *x, double *T, const double *P, const double *TcCel,
                 const double *Pc, const double *omega, const int *Ncomps,
                 double *phi);

// Passing parameter struct for VLE and Antoine functions
typedef struct {
  double *x0;
  const double P;
  const double *antProps;
  const double (*nrtlA)[];
  const double (*nrtlB)[];
  const double (*nrtlC)[];
  const double *TcCel;
  const double *Pc;
  const double *omega;
  const int Ncomps;
  const int antMethod;
  const int actMethod;
  const double dxi;
  const int n_it;
  const int maxiter;
  const double ftol;
  const double xtol;
} params_t;

// Return structure for RCM curves
typedef struct {
  double(*x);
  double(*y);
  double *T;
} curves_t;

// For flipping arrays before iterating in other direction
void flipNreverse(int n_it, const int Ncomps, double x_arr1[n_it][Ncomps],
                  double *x_arr2) {

  for (int i = 0; i < n_it; i++) {
    for (int j = 0; j < Ncomps; j++) {
      x_arr2[i * Ncomps + j] = x_arr1[n_it - 1 - i][j];
    }
  }
}

void flipT(int n_it, const int Ncomps, double T_arr1[n_it], double *T_arr2) {
  for (int i = 0; i < n_it; i++) {
    T_arr2[i] = T_arr1[n_it - 1 - i];
  }
}

// Definition of Antoine function
void antoineCalc(const double *T, void *params, double *Psat) {

  params_t *p = (params_t *)params;

  const int Ncomps = p->Ncomps;
  const int antMethod = p->antMethod;

  if (antMethod == 1) {
    double logP[Ncomps];
    const double *antProps = p->antProps;
    for (int i = 0; i < Ncomps; i++) {
      logP[i] = antProps[i + Ncomps * 0] -
                antProps[i + Ncomps * 1] / (antProps[i + Ncomps * 2] + *T);
      Psat[i] = pow(10, logP[i]); // need to fix indexing here
    }
  }

  if (antMethod == 2) {
    const double Tk = *T + (double)273.15;
    double lnP[Ncomps];
    const double *antProps = p->antProps;
    for (int i = 0; i < Ncomps; i++) {
      lnP[i] = (antProps[i * 7 + 0] +
                antProps[i * 7 + 1] / (antProps[i * 7 + 2] + Tk) +
                antProps[i * 7 + 3] * Tk + antProps[i * 7 + 4] * log(Tk) +
                antProps[i * 7 + 5] * pow(Tk, antProps[i * 7 + 6]));
      Psat[i] = exp(lnP[i]);
    }
  }
}

// Wrapper for using Antoine to calculate pressure residual for saturated
// temperature estimate
int antPressResid(const gsl_vector *T_vec, void *params, gsl_vector *f_vec) {
  params_t *p = (params_t *)params;

  double T = gsl_vector_get(T_vec, 0);

  const int Ncomps = p->Ncomps;

  double Psat[Ncomps];

  antoineCalc(&T, params, Psat);

  double Psat_sum = 0.0;
  for (int i = 0; i < Ncomps; i++) {
    Psat_sum += Psat[i];
  }

  double residual = Psat_sum - p->P;

  gsl_vector_set(f_vec, 0, residual);

  return GSL_SUCCESS;
}

// Wrapper for solving antPressResid
int antsolve(double T_initial, params_t *params, gsl_vector *T_sol) {
  const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc(T, 1);

  gsl_multiroot_function f = {&antPressResid, 1, params};

  gsl_vector *T_vec = gsl_vector_alloc(1);
  gsl_vector_set(T_vec, 0, T_initial);

  gsl_multiroot_fsolver_set(s, &f, T_vec);

  int status;
  size_t iter = 0;

  do {
    iter++;
    status = gsl_multiroot_fsolver_iterate(s);
    if (status)
      break;

    status = gsl_multiroot_test_residual(s->f, params->ftol);
  } while (status == GSL_CONTINUE && iter < params->maxiter);

  if (status == GSL_SUCCESS) {
    gsl_vector_set(T_sol, 0, gsl_vector_get(s->x, 0));
  } else {
    printf("antsolver did not converge.\n");
  }

  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(T_vec);
  return status;
}

// Define the residual function for VLE
int VLEfunc(const gsl_vector *Y, void *params, gsl_vector *f) {
  params_t *p = (params_t *)params;

  double *x = p->x0;
  const double P = p->P;
  const int Ncomps = p->Ncomps;
  const double(*nrtlA)[Ncomps] = p->nrtlA;
  const double(*nrtlB)[Ncomps] = p->nrtlB;
  const double(*nrtlC)[Ncomps] = p->nrtlC;
  const double *TcCel = p->TcCel;
  const double *Pc = p->Pc;
  const double *omega = p->omega;
  const int actMethod = p->actMethod;

  double y[Ncomps];
  for (int i = 0; i < Ncomps; i++) {
    y[i] = gsl_vector_get(Y, i);
  }
  double T = gsl_vector_get(Y, Ncomps);

  double Psat[Ncomps], gamma[Ncomps], phi[Ncomps], raoult[Ncomps];
  double sumy = 0.0;

  antoineCalc(&T, params, Psat);
  if (actMethod >= 2) {
    nrtl_(x, &T, &nrtlA[0][0], &nrtlB[0][0], &nrtlC[0][0], &Ncomps, gamma);
  } else {
    for (int i = 0; i < Ncomps; i++)
      gamma[i] = 1.0;
  }
  if (actMethod >= 3) {
    srk_(y, &T, &P, &TcCel[0], &Pc[0], &omega[0], &Ncomps, phi);
  } else {
    for (int i = 0; i < Ncomps; i++)
      phi[i] = 1.0;
  }

  for (int i = 0; i < Ncomps; i++) {
    raoult[i] = x[i] * gamma[i] * Psat[i] - y[i] * phi[i] * P;
    gsl_vector_set(f, i, raoult[i]);
    sumy += fabs(y[i]);
  }

  gsl_vector_set(f, Ncomps, sumy - (double)1.00000000);

  return GSL_SUCCESS;
}

// Wrapper for VLE solver
int VLEsolve(gsl_vector *Y_init, void *params, gsl_vector *Y_sol) {
  params_t *p = (params_t *)params;
  const int Ncomps = p->Ncomps;

  const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc(T, Ncomps + 1);

  gsl_multiroot_function f = {&VLEfunc, Ncomps + 1, params};

  gsl_multiroot_fsolver_set(s, &f, Y_init);

  int status;
  size_t iter = 0;

  do {
    iter++;
    status = gsl_multiroot_fsolver_iterate(s);
    if (status)
      break;

    status = gsl_multiroot_test_residual(s->f, p->ftol);
  } while (status == GSL_CONTINUE && iter < p->maxiter);

  if (status == GSL_SUCCESS) {
    for (int i = 0; i < Ncomps + 1; i++) {
      gsl_vector_set(Y_sol, i, gsl_vector_get(s->x, i));
    }
  } else {
    printf("VLEsolver did not converge.\n");
  }

  gsl_multiroot_fsolver_free(s);
  return status;
}

// Parent solver
curves_t *RCM(void *params) {
  params_t *p = (params_t *)params;

  int status = 1;
  const int Ncomps = p->Ncomps;
  double P = p->P;
  gsl_vector *T0 = gsl_vector_alloc(1);
  double dxi = p->dxi;
  int n_it = p->n_it;
  double *x0 = p->x0;

  gsl_vector *Y_init = gsl_vector_alloc(Ncomps + 1);
  gsl_vector *Y_sol = gsl_vector_alloc(Ncomps + 1);
  double x_arr1[n_it][Ncomps];
  double y_arr1[n_it][Ncomps];
  double T_arr1[n_it];

  while (status != GSL_SUCCESS) {
    status = antsolve(50 + (3 * P), p, T0);
    if (status != 0) {
      int inc = 20;
      while (status != GSL_SUCCESS || inc < 80) {
        status = antsolve(inc + (3 * P), p, T0);
        inc += 0.5;
      }
      if (inc >= 80) {
        printf("Initial Antoine guess failed\n");
      }
    }
  }

  for (int i = 0; i < Ncomps; i++) {
    x_arr1[0][i] = x0[i];
    gsl_vector_set(Y_init, i, x0[i]);
  }
  gsl_vector_set(Y_init, Ncomps, gsl_vector_get(T0, 0));

  for (int i = 0; i < n_it; i++) {
    VLEsolve(Y_init, params, Y_sol);
    for (int j = 0; j < Ncomps; j++) {
      y_arr1[i][j] = gsl_vector_get(Y_sol, j);
    }
    T_arr1[i] = gsl_vector_get(Y_sol, Ncomps);
    gsl_vector_memcpy(Y_init, Y_sol);
    for (int j = 0; j < Ncomps; j++) {
      x_arr1[i + 1][j] = x_arr1[i][j] + dxi * (x_arr1[i][j] - y_arr1[i][j]);
    }
    p->x0 = x_arr1[i + 1];
  }

  curves_t *c = (curves_t *)malloc(sizeof(curves_t));

  c->x = (double *)malloc(2 * n_it * Ncomps * sizeof(double));
  c->y = (double *)malloc(2 * n_it * Ncomps * sizeof(double));
  c->T = (double *)malloc(2 * n_it * sizeof(double));
  flipNreverse(n_it, Ncomps, x_arr1, c->x);
  flipNreverse(n_it, Ncomps, y_arr1, c->y);
  flipT(n_it, Ncomps, T_arr1, c->T);

  for (int i = 0; i < Ncomps; i++) {
    gsl_vector_set(Y_init, i, x0[i]);
  }
  gsl_vector_set(Y_init, Ncomps, gsl_vector_get(T0, 0));

  p->x0 = x0;

  for (int i = n_it - 1; i < 2 * n_it; i++) {
    VLEsolve(Y_init, params, Y_sol);
    for (int j = 0; j < Ncomps; j++) {
      c->y[i * Ncomps + j] = gsl_vector_get(Y_sol, j);
    }
    c->T[i] = gsl_vector_get(Y_sol, Ncomps);
    gsl_vector_memcpy(Y_init, Y_sol);
    for (int j = 0; j < Ncomps; j++) {
      c->x[(i + 1) * Ncomps + j] =
          c->x[i * Ncomps + j] -
          dxi * (c->x[i * Ncomps + j] - c->y[i * Ncomps + j]);
    }
    p->x0 = &c->x[(i + 1) * Ncomps];
  }

  return c;
}

void freeCurveMem(curves_t *c) {
  free(c->x);
  free(c->y);
  free(c->T);
  free(c);
}
