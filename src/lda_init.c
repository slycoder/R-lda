#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP collapsedGibbsSampler(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cvb0(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP nubbi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rtm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"collapsedGibbsSampler", (DL_FUNC) &collapsedGibbsSampler, 18},
    {"cvb0",                  (DL_FUNC) &cvb0,                   7},
    {"nubbi",                 (DL_FUNC) &nubbi,                  9},
    {"rtm",                   (DL_FUNC) &rtm,                   10},
    {NULL, NULL, 0}
};

void R_init_lda(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
