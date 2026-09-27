.onAttach <-
function (libname, pkgname){
   # echo output to screen
   packageStartupMessage("## gsynth v1.5.0 (requires fect >= 2.4.7).")
   packageStartupMessage("## *gsynth* is a wrapper of the *fect* package; see ?gsynth.")
   packageStartupMessage("## NOTE: estimator='ife' and estimator='mc' are soft-deprecated in v1.5.0;")
   packageStartupMessage("##   for IFE-EM or matrix completion, use fect::fect() directly.")
   packageStartupMessage("##   See https://yiqingxu.org/packages/gsynth/02-ife-mc.html.")
}
