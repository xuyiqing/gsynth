# gsynth 1.1.7

Using normal approximation instead of the percentile method to obtain confidence intervals based on bootstrapped standard errors. 

# gsynth 1.1.4

1. Import function *felm* from **lfe** to fit two-way fixed effects model as the 
starting value for estimation of interactive fixed effects model with unbalanced 
panel data.
2. Add cluster bootstrap option for uncertainty estimates.
3. Add jackknife uncertainty estimates.
3. Add a new function *cumuEff* for calculation of sub-group and cumulative 
treatment effects.
 
# gsynth 1.0.9

1. Function `panelView()` is removed from **gsynth** and becomes an independent package [**panelview**](https://yiqingxu.org/packages/panelview/). 
2. Implement the matrix completion method.
2. Fix bugs.
3. Change the color scheme.

# gsynth 1.0.8

1. Add a function `panelView()` to visualize raw data and data structure before estimation.
2. Fix bugs.

# gsynth 1.0.7

1. Add "implied weights" of control units for each treated unit to the output of the main function (`wgt.implied`).
2. Add a plot to visualize missing data and treatment status (`type = "missing"`).
3. Accommodate unbalanced panels.
