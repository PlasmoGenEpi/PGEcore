#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>

void McCOIL_categorical(int *max, int *iterations, int *n0, int *k0,
                        double *sampleS2, int *M0, double *P0, double *error1,
                        double *error2, char **file_index, char **path,
                        int *err_method0);

void McCOIL_prop(int *max, int *iterations, int *n0, int *k0, double *A1,
                 double *A2, int *M0, double *P0, double *A, double *B,
                 double *c0, char **file_index, char **path, int *err_method0);

static const R_CMethodDef CEntries[] = {
    {"McCOIL_categorical", (DL_FUNC)&McCOIL_categorical, 12},
    {"McCOIL_prop", (DL_FUNC)&McCOIL_prop, 14},
    {NULL, NULL, 0}};

void R_init_PGEcore(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
