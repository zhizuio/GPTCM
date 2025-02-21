#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _GPTCM_arms_gibbs_beta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GPTCM_arms_gibbs_xi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GPTCM_ars(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GPTCM_ars_debug(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GPTCM_ars_gibbs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GPTCM_sampleGamma(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_GPTCM_arms_gibbs_beta", (DL_FUNC) &_GPTCM_arms_gibbs_beta, 21},
    {"_GPTCM_arms_gibbs_xi",   (DL_FUNC) &_GPTCM_arms_gibbs_xi,   16},
    {"_GPTCM_ars",             (DL_FUNC) &_GPTCM_ars,             12},
    {"_GPTCM_ars_debug",       (DL_FUNC) &_GPTCM_ars_debug,       13},
    {"_GPTCM_ars_gibbs",       (DL_FUNC) &_GPTCM_ars_gibbs,       11},
    {"_GPTCM_sampleGamma",     (DL_FUNC) &_GPTCM_sampleGamma,      4},
    {NULL, NULL, 0}
};

void R_init_GPTCM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}