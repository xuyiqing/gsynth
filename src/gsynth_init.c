#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP gsynth_beta_iter(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gsynth_inter_fe(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gsynth_panel_beta(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gsynth_panel_est(SEXP, SEXP, SEXP);
extern SEXP gsynth_panel_factor(SEXP, SEXP);
extern SEXP gsynth_XXinv(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"gsynth_beta_iter",    (DL_FUNC) &gsynth_beta_iter,    6},
    {"gsynth_inter_fe",     (DL_FUNC) &gsynth_inter_fe,     6},
    {"gsynth_panel_beta",   (DL_FUNC) &gsynth_panel_beta,   5},
    {"gsynth_panel_est",    (DL_FUNC) &gsynth_panel_est,    3},
    {"gsynth_panel_factor", (DL_FUNC) &gsynth_panel_factor, 2},
    {"gsynth_XXinv",        (DL_FUNC) &gsynth_XXinv,        1},
    {NULL, NULL, 0}
};

void R_init_gsynth(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
	R_RegisterCCallable("gsynth", "gsynth_XXinv", (DL_FUNC) gsynth_XXinv);
    R_RegisterCCallable("gsynth", "gsynth_panel_est", (DL_FUNC) gsynth_panel_est);
    R_RegisterCCallable("gsynth", "gsynth_panel_beta", (DL_FUNC) gsynth_panel_beta);
    R_RegisterCCallable("gsynth", "gsynth_panel_factor", (DL_FUNC) gsynth_panel_factor);
    R_RegisterCCallable("gsynth", "gsynth_beta_iter", (DL_FUNC) gsynth_beta_iter);
    R_RegisterCCallable("gsynth", "gsynth_inter_fe", (DL_FUNC) gsynth_inter_fe);
}
