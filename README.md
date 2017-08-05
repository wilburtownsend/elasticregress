# elasticreg
Stata implementation of the Friedman, Hastie and Tibshirani (2010, JStatSoft) coordinate descent algorithm for elastic net regression and its famous special cases: LASSO and ridge regression.

```

Title

    elasticreg -- Elastic net regression

    lassoreg -- LASSO regression

    ridgereg -- Ridge regression


Syntax

        elasticreg depvar [indepvars] [if] [in] [weight] , alpha(#) [options]

        lassoreg depvar [indepvars] [if] [in] [weight] [, options]

        ridgereg depvar [indepvars] [if] [in] [weight] [, options]

    options               Description
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Main
      alpha               weight placed on the LASSO (L1) norm, one minus weight placed on the ridge (L2) norm;
      lambda              penalty placed on larger coefficients — by default found by cross-validation.
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Options which only matter when lambda is found through cross-validation
      numlambda           number of lambda tested when lambda is found by cross-validation;
      lambdamin           lambda selected is that which minimises cross-validation MSE (default);
      lambda1se           lambda selected is largest within a standard error of that selected under lambdamin;
      numfolds            number of folds used when cross-validating lambda — default is 10;
      epsilon             ratio of the smallest lambda tested to the largest — default is 0.001, the user is prompted when this constraint is binding.
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Technical options which you can hopefully ignore
      tol                 tolerance on the Euclidean norm of beta.
    ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    by is allowed; see [D] by.
    aweights are allowed; see weight.


Description

    elasticreg calculates an elastic net-regularized regression: an estimator of a linear model in which larger parameters are discouraged.  This estimator nests the LASSO and the ridge regression, which can be estimated by setting
    alpha equal to 1 and 0 respectively.

    lassoreg estimates the LASSO; it is a convenience command equivalent to elasticreg with the option alpha(1).

    ridgereg estimates a ridge regression; it is a convenience command equivalent to elasticreg with the option alpha(0).


Remarks

    elasticreg implements the coordinate descent algorithm in Friedman, Hastie and Tibshirani (2008) to produce very efficient estimates of the LASSO, of the ridge regression and of a generalization of the two. The algorithm used is
    the 'covariance updates' algorithm presented in Section 2 of that paper, except that it does not exploit sparsity in the covariates.

    When lambda is not supplied, it is calculated using K-fold cross-validation. In each cross-validation sample beta(lambda) is calculated using pathwise coordinate descent — starting at some large lambda, and moving down a
    sequence of lambda using the previous solution as the starting vector. Generally the starting lambda is sufficiently large to guarantee that beta=0, but when alpha is small no such lambda exists — you will be warned if the upper
    limit on lambda seems binding. The sequence of lambda are equidistant in log space and thus the smallest lambda must be positive. The user will also be warned if this lower limit seems to be binding (it often is when alpha is
    small and thus the maximal lambda is very large).

    To submit bugs or request additional features, please post an issue on the project Github.


Example


. * Set the random number seed for consistent cross-validation.
. set seed 1
. * Generate some simple data.
. clear
. set obs 100000
number of observations (_N) was 0, now 100,000
. gen     x1 = rnormal(0,2) > 1. gen x2 = rnormal(4,1)
. gen x3 = rnormal(5,3)
. gen x4 = rnormal(0,1)
. gen u  = rnormal(0,4)
. gen y  = 0.2*x1 + 0.5*x2  + 1*x4 + u + 6

. * Calculate a LASSO model.
. lassoreg y x*
------------------------------------------------------------------------------
           y |      Coef.
-------------+----------------------------------------------------------------
          x1 |   .0930302
          x2 |   .4437604
          x3 |          0
          x4 |   .9417223
       _cons |   6.256568
------------------------------------------------------------------------------

. * Calculate a ridge-regression model.
. ridgereg y x*
Warning: the smallest λ tested was the MSE-minimising λ.
Consider re-running estimation with a smaller epsilon.
------------------------------------------------------------------------------
           y |      Coef.
-------------+----------------------------------------------------------------
          x1 |   .1026057
          x2 |   .2477366
          x3 |   .0024702
          x4 |   .4987638
       _cons |   7.022845
------------------------------------------------------------------------------

. * Calculate a ridge-regression model, testing smaller lambda.
. ridgereg y x*, epsilon(0.00001)
------------------------------------------------------------------------------
           y |      Coef.
-------------+----------------------------------------------------------------
          x1 |   .1880368
          x2 |   .4506527
          x3 |   .0048642
          x4 |   .9032322
       _cons |   6.175198
------------------------------------------------------------------------------

. * Calculate OLS — equivalent to lasso or ridge with lambda=0.
. lassoreg y x*, lambda(0)
------------------------------------------------------------------------------
           y |      Coef.
-------------+----------------------------------------------------------------
          x1 |    .207751
          x2 |   .4970736
          x3 |    .005458
          x4 |   .9952633
       _cons |   5.980998
------------------------------------------------------------------------------


Stored Results

    elasticreg, lassoreg and ridgereg store the following in e():

    Scalars        
      e(N)                number of observations
      e(lambda)           lambda provided or selected by cross-validation
      e(r2)               R-squared
      e(alpha)            alpha provided

    Macros         
      e(cmd)              the command used: elasticreg, lassoreg or ridgereg
      e(varlist_nonzero)  list of covariates with non-zero coefficients
      e(properties)       b
      e(depvar)           name of dependent variable

    Matrices       
      e(b)                coefficient vector

    Functions      
      e(sample)           marks estimation sample


References

    Jerome Friedman, Trevor Hastie and Rob Tibshirani. (2008). Regularization Paths for Generalized Linear Models via Coordinate Descent.  Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010.


```
