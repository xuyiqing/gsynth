# gsynth

## Generalized Synthetic Control Method
---

**Authors:** Yiqing Xu [<yiqingxu@ucsd.edu>]; Licheng Liu [<liulch.16@sem.tsinghua.edu.cn>] 

**How to Uses:** [Examples](http://yiqingxu.org/software/gsynth/gsynth_examples.html)

**Reference:**  Yiqing Xu. 2017. "Generalized Synthetic Control Method: Causal Inference  with Interactive Fixed Effects Models." Political Analysis, Vol. 25, Iss. 1, January 2017, pp. 57-76. Available at: <http://dx.doi.org/10.1017/pan.2016.2>

**Note:**

Rcpp, RcppArmadillo and MacOS "-lgfortran" and "-lquadmath" error, see: http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/

Installation failture related to OpenMP on MacOS, see:
http://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/

To fix these problems, consider installing: 
gfortran 6.1 from https://gcc.gnu.org/wiki/GFortranBinaries#MacOS
clang4 R Binaries from https://github.com/coatless/r-macos-clang
