.onAttach <-
function (libname, pkgname){
   # echo output to screen
   packageStartupMessage("## v1.4.0: estimator='ife' now uses the EM algorithm; default is 'gsynth'. See ?gsynth.")
   packageStartupMessage("## Since v.1.3.0, *gsynth* is a wrapper of the *fect* package.")
}
