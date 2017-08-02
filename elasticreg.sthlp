{smcl}
{* *! version 1.2.1  07mar2013}{...}
{findalias asfradohelp}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "elasticreg" "help elasticreg"}{...}
{viewerjumpto "Syntax"         "elasticreg##syntax"}{...}
{viewerjumpto "Description"    "elasticreg##description"}{...}
{viewerjumpto "Remarks"        "elasticreg##remarks"}{...}
{viewerjumpto "Examples"       "elasticreg##examples"}{...}
{viewerjumpto "Stored results" "elasticreg##results"}{...}
{viewerjumpto "References"     "elasticreg##references"}{...}

{title:Title}

{phang} {bf:elasticreg} {hline 2} Elastic net regression

{phang} {bf:lassoreg}   {hline 2} LASSO regression	

{phang} {bf:ridgereg}   {hline 2} Ridge regression


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}{opt elasticreg} {depvar} [{indepvars}] {ifin} {weight} {cmd:,} alpha(#) [{it:options}]

{p 8 17 2}{opt lassoreg}   {depvar} [{indepvars}] {ifin} {weight} [{cmd:,} {it:options}]

{p 8 17 2}{opt ridgereg}   {depvar} [{indepvars}] {ifin} {weight} [{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt alpha}}weight placed on the LASSO (L1) norm, one minus weight placed on the ridge (L2) norm;{p_end}
{synopt:{opt lambda}}penalty placed on larger coefficients — by default found by cross-validation.{p_end}
{synoptline}
{syntab:Options which only matter when lambda is found through cross-validation}
{synopt:{opt numlambda}}number of lambda tested when lambda is found by cross-validation;{p_end}
{synopt:{opt lambdamin}}lambda selected is that which minimises cross-validation MSE (default); {p_end}
{synopt:{opt lambda1se}}lambda selected is largest within a standard error of that selected under lambdamin; {p_end}
{synopt:{opt numfolds}}number of folds used when cross-validating lambda — default is 10; {p_end}
{synopt:{opt epsilon}}ratio of the smallest lambda tested to the largest — default is 0.001, the user is prompted when this constraint is binding. {p_end}
{synoptline}
{syntab:Technical options which you can hopefully ignore}
{synopt:{opt tol}}tolerance on the Euclidean norm of beta. {p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{cmd:by} is allowed; see {manhelp by D}.{p_end}
{p 4 6 2}
{cmd:aweight}s are allowed; see {help weight}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:elasticreg} calculates an elastic net-regularized regression: an estimator 
				of a linear model in which larger parameters are discouraged. 
				This estimator nests the LASSO and the ridge regression, which
				can be estimated by setting alpha equal to 1 and 0 respectively.
	
{pstd}			
{cmd:lassoreg}  estimates the LASSO; it is a convenience command equivalent to 
				elasticreg with the option alpha(1).
		
{pstd}		
{cmd:ridgereg}  estimates a ridge regression; it is a convenience command equivalent to 
				elasticreg with the option alpha(0).


{marker remarks}{...}
{title:Remarks}

{pstd}
{cmd:elasticreg} implements the coordinate descent algorithm in Friedman, Hastie
and Tibshirani (2008) to produce very efficient estimates of the LASSO, of the ridge 
regression and of a generalization of the two. The algorithm used is the 
'covariance updates' algorithm presented in Section 2 of that paper, except that
it does not exploit sparsity in the covariates.

{pstd}
When lambda is not supplied, it is calculated using K-fold cross-validation. In
each cross-validation sample beta(lambda) is calculated using pathwise coordinate
descent — starting at some large lambda, and moving down a sequence of lambda
using the previous solution as the starting vector. Generally the starting lambda
is sufficiently large to guarantee that beta=0, but when alpha is small no such lambda
exists — you will be warned if the upper limit on lambda seems binding. The 
sequence of lambda are equidistant in log space and thus the smallest
lambda must be positive. The user will also be warned if this lower limit seems
to be binding (it often is when alpha is small and thus the maximal lambda
is very large).

{pstd}
To submit bugs or request additional features, please post an issue on the
{browse "http://github.com/wilbur-t/elasticreg/issues/":project Github}.


{marker examples}{...}
{title:Example}


{com}. * Set the random number seed for consistent cross-validation.
. set seed 1
{com}. * Generate some simple data.
. clear
{com}. set obs 100000
{txt}{p}
number of observations (_N)  was 0,
now 100,000
{p_end}
{com}. gen     x1 = rnormal(0,2) > 1{com}. gen x2 = rnormal(4,1)
{com}. gen x3 = rnormal(5,3)
{com}. gen x4 = rnormal(0,1)
{com}. gen u  = rnormal(0,4)
{com}. gen y  = 0.2*x1 + 0.5*x2  + 1*x4 + u + 6

{com}. * Calculate a LASSO model.
. lassoreg y x*
{res}{txt}{hline 13}{c TT}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}
{col 1}           y{col 14}{c |}      Coef.
{hline 13}{c +}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}
{space 10}x1 {c |}{col 14}{res}{space 2} .0930302
{txt}{space 10}x2 {c |}{col 14}{res}{space 2} .4437604
{txt}{space 10}x3 {c |}{col 14}{res}{space 2}        0
{txt}{space 10}x4 {c |}{col 14}{res}{space 2} .9417223
{txt}{space 7}_cons {c |}{col 14}{res}{space 2} 6.256568
{txt}{hline 13}{c BT}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}

{com}. * Calculate a ridge-regression model.
. ridgereg y x*
{res}Warning: the smallest λ tested was the MSE-minimising λ.
Consider re-running estimation with a smaller epsilon.
{txt}{hline 13}{c TT}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}
{col 1}           y{col 14}{c |}      Coef.
{hline 13}{c +}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}
{space 10}x1 {c |}{col 14}{res}{space 2} .1026057
{txt}{space 10}x2 {c |}{col 14}{res}{space 2} .2477366
{txt}{space 10}x3 {c |}{col 14}{res}{space 2} .0024702
{txt}{space 10}x4 {c |}{col 14}{res}{space 2} .4987638
{txt}{space 7}_cons {c |}{col 14}{res}{space 2} 7.022845
{txt}{hline 13}{c BT}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}

{com}. * Calculate a ridge-regression model, testing smaller lambda.
. ridgereg y x*, epsilon(0.00001)
{res}{txt}{hline 13}{c TT}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}
{col 1}           y{col 14}{c |}      Coef.
{hline 13}{c +}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}
{space 10}x1 {c |}{col 14}{res}{space 2} .1880368
{txt}{space 10}x2 {c |}{col 14}{res}{space 2} .4506527
{txt}{space 10}x3 {c |}{col 14}{res}{space 2} .0048642
{txt}{space 10}x4 {c |}{col 14}{res}{space 2} .9032322
{txt}{space 7}_cons {c |}{col 14}{res}{space 2} 6.175198
{txt}{hline 13}{c BT}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}

{com}. * Calculate OLS — equivalent to lasso or ridge with lambda=0.
. lassoreg y x*, lambda(0)
{res}{txt}{hline 13}{c TT}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}
{col 1}           y{col 14}{c |}      Coef.
{hline 13}{c +}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}
{space 10}x1 {c |}{col 14}{res}{space 2}  .207751
{txt}{space 10}x2 {c |}{col 14}{res}{space 2} .4970736
{txt}{space 10}x3 {c |}{col 14}{res}{space 2}  .005458
{txt}{space 10}x4 {c |}{col 14}{res}{space 2} .9952633
{txt}{space 7}_cons {c |}{col 14}{res}{space 2} 5.980998
{txt}{hline 13}{c BT}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}


{marker results}{...}
{title:Stored Results}

{pstd}
{cmd:elasticreg}, {cmd:lassoreg} and {cmd:ridgereg} store the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(lambda)}}lambda provided or selected by cross-validation{p_end}
{synopt:{cmd:e(r2)}}R-squared{p_end}
{synopt:{cmd:e(alpha)}}alpha provided{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}the command used: {cmd:elasticreg}, {cmd:lassoreg} or {cmd:ridgereg}{p_end}
{synopt:{cmd:e(varlist_nonzero)}}list of covariates with non-zero coefficients{p_end}
{synopt:{cmd:e(properties)}}{cmd:b}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{pstd}
Jerome Friedman, Trevor Hastie and Rob Tibshirani. (2008). Regularization Paths for Generalized Linear Models via Coordinate Descent.
Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010.
