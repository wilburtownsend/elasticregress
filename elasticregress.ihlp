{smcl}
{* *! version 1.3  16apr2018}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "elasticregress" "help elasticregress"}{...}
{viewerjumpto "Syntax"         "elasticregress##syntax"}{...}
{viewerjumpto "Description"    "elasticregress##description"}{...}
{viewerjumpto "Remarks"        "elasticregress##remarks"}{...}
{viewerjumpto "Examples"       "elasticregress##examples"}{...}
{viewerjumpto "Stored results" "elasticregress##results"}{...}
{viewerjumpto "Author"         "elasticregress##author"}{...}
{viewerjumpto "References"     "elasticregress##references"}{...}

{title:Title}

{phang} {bf:elasticregress} {hline 2} Elastic net regression

{phang} {bf:lassoregress}   {hline 2} LASSO regression  

{phang} {bf:ridgeregress}   {hline 2} Ridge regression


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}{opt elasticregress} {depvar} [{indepvars}] {ifin} {weight} [{cmd:,} alpha(#) {it:options}]

{p 8 17 2}{opt lassoregress}   {depvar} [{indepvars}] {ifin} {weight} [{cmd:,} {it:options}]

{p 8 17 2}{opt ridgeregress}   {depvar} [{indepvars}] {ifin} {weight} [{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt alpha}}weight placed on the LASSO (L1) norm, one minus weight placed on the ridge (L2) norm — by default found by cross-validation;{p_end}
{synopt:{opt lambda}}penalty placed on larger coefficients — by default found by cross-validation; {p_end}
{synopt:{opt numfolds}}number of folds used when cross-validating lambda or alpha — default is 10. {p_end}
{synoptline}
{syntab:Options which only matter when alpha is found through cross-validation}
{synopt:{opt numalpha}}number of alpha tested when alpha is found by cross-validation.{p_end}
{synoptline}
{syntab:Options which only matter when lambda is found through cross-validation}
{synopt:{opt numlambda}}number of lambda tested when lambda is found by cross-validation;{p_end}
{synopt:{opt lambdamin}}lambda selected is that which minimises cross-validation MSE (default); {p_end}
{synopt:{opt lambda1se}}lambda selected is largest within a standard error of that selected under lambdamin; {p_end}
{synopt:{opt epsilon}}ratio of the smallest lambda tested to the largest — default is 0.001, the user is prompted when this constraint is binding. {p_end}
{synoptline}
{syntab:Technical options which you can hopefully ignore}
{synopt:{opt tol}}tolerance on the Euclidean norm of beta; {p_end}
{synopt:{opt collinear}}retain collinear independent variables. {p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{cmd:by} is allowed; see {manhelp by D}.{p_end}
{p 4 6 2}
{cmd:aweight}s are allowed; see {help weight}.


{marker description}{...}
{title:Description}

{pstd}
{cmd:elasticregress} calculates an elastic net-regularized regression: an estimator 
                of a linear model in which larger parameters are discouraged. 
                This estimator nests the LASSO and the ridge regression, which
                can be estimated by setting alpha equal to 1 and 0 respectively.
    
{pstd}          
{cmd:lassoregress}  estimates the LASSO; it is a convenience command equivalent to 
                elasticregress with the option alpha(1).
        
{pstd}      
{cmd:ridgeregress}  estimates a ridge regression; it is a convenience command equivalent to 
                elasticregress with the option alpha(0).


{marker remarks}{...}
{title:Remarks}

{pstd}
{cmd:elasticregress} implements the coordinate descent algorithm in Friedman, Hastie
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
elasticregress extends the Friedman, Hastie and Tibshirani (2008) by allowing
the elastic-net mixing parameter alpha to be also found using cross-validation.
Consistent cross-validation samples are estimated for a series of alpha, and the
alpha selected is that which minimises CV-MSE(alpha, lambda_min(alpha)), where
CV-MSE() is the mean error in the cross validation samples and lambda_min(alpha)
is the lambda which minimises cross-validation mean squared error given alpha.

{pstd}
To submit bugs or request additional features, please post an issue on the
{browse "http://github.com/wilbur-t/elasticregress/issues/":project Github}.


{marker examples}{...}
{title:Example}

{pstd}
Load our favourite data:
{p_end}{pstd}
{stata sysuse auto, clear}

{pstd}
Set the seed, for consistent cross-validation:
{p_end}{pstd}
{stata set seed 1}

{pstd}
Calculate a LASSO model: {p_end}{pstd}
{stata lassoregress mpg weight foreign} {p_end}

{txt}LASSO regression{col 40}Number of observations{col 67}= {res}        74
{col 40}{txt}R-squared{col 67}= {res}    0.6466
{col 40}{txt}alpha{col 67}= {res}    1.0000
{col 40}{txt}lambda{col 67}= {res}    0.4034
{col 40}{txt}Cross-validation MSE{col 67}= {res}   12.0339
{col 40}{txt}Number of folds{col 67}= {res}        10
{col 40}{txt}Number of lambda tested{col 67}= {res}       100
{txt}{hline 13}{c TT}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}
{col 1}         mpg{col 14}{c |}      Coef.
{hline 13}{c +}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}
{space 6}weight {c |}{col 14}{res}{space 2}-.0054861
{txt}{space 5}foreign {c |}{col 14}{res}{space 2}        0
{txt}{space 7}_cons {c |}{col 14}{res}{space 2} 37.86231
{txt}{hline 13}{c BT}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}


{pstd}
Calculate a ridge-regression model: {p_end}{pstd}
{stata ridgeregress mpg weight foreign}  {p_end}

{res}Warning: the smallest λ tested was the MSE-minimising λ.
Consider re-running estimation with a smaller epsilon.

{txt}Ridge regression{col 40}Number of observations{col 67}= {res}        74
{col 40}{txt}R-squared{col 67}= {res}    0.2343
{col 40}{txt}alpha{col 67}= {res}    0.0000
{col 40}{txt}lambda{col 67}= {res}    4.6383
{col 40}{txt}Cross-validation MSE{col 67}= {res}   27.1411
{col 40}{txt}Number of folds{col 67}= {res}        10
{col 40}{txt}Number of lambda tested{col 67}= {res}       100
{txt}{hline 13}{c TT}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}
{col 1}         mpg{col 14}{c |}      Coef.
{hline 13}{c +}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}
{space 6}weight {c |}{col 14}{res}{space 2}-.0010224
{txt}{space 5}foreign {c |}{col 14}{res}{space 2} .6956361
{txt}{space 7}_cons {c |}{col 14}{res}{space 2} 24.17757
{txt}{hline 13}{c BT}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}


{pstd}
Calculate OLS — equivalent to lasso or ridge with lambda=0: {p_end}{pstd}
{stata lassoregress mpg weight foreign, lambda(0)}  {p_end}

{txt}LASSO regression{col 40}Number of observations{col 67}= {res}        74
{col 40}{txt}R-squared{col 67}= {res}    0.6627
{col 40}{txt}alpha{col 67}= {res}    1.0000
{col 40}{txt}lambda{col 67}= {res}    0.0000
{txt}{hline 13}{c TT}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}
{col 1}         mpg{col 14}{c |}      Coef.
{hline 13}{c +}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}
{space 6}weight {c |}{col 14}{res}{space 2}-.0065875
{txt}{space 5}foreign {c |}{col 14}{res}{space 2}-1.649645
{txt}{space 7}_cons {c |}{col 14}{res}{space 2} 41.67843
{txt}{hline 13}{c BT}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}

{pstd}
Calculate an elastic-net regression: {p_end}{pstd}
{stata elasticregress mpg weight foreign}  {p_end}

{res}Warning: the smallest λ tested was the MSE-minimising λ.
Consider re-running estimation with a smaller epsilon.

{txt}Elastic-net regression{col 40}Number of observations{col 67}= {res}        74
{col 40}{txt}R-squared{col 67}= {res}    0.6622
{col 40}{txt}alpha{col 67}= {res}    0.2000
{col 40}{txt}lambda{col 67}= {res}    0.0232
{col 40}{txt}Cross-validation MSE{col 67}= {res}   12.8399
{col 40}{txt}Number of folds{col 67}= {res}        10
{col 40}{txt}Number of alpha tested{col 67}= {res}         6
{col 40}{txt}Number of lambda tested{col 67}= {res}       100
{txt}{hline 13}{c TT}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}
{col 1}         mpg{col 14}{c |}      Coef.
{hline 13}{c +}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}
{space 6}weight {c |}{col 14}{res}{space 2}-.0063764
{txt}{space 5}foreign {c |}{col 14}{res}{space 2}-1.402116
{txt}{space 7}_cons {c |}{col 14}{res}{space 2} 40.96739
{txt}{hline 13}{c BT}{hline 11}{hline 11}{hline 9}{hline 8}{hline 13}{hline 12}


{marker results}{...}
{title:Stored Results}

{pstd}
{cmd:elasticregress}, {cmd:lassoregress} and {cmd:ridgeregress} store the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(alpha)}}alpha provided or selected by cross-validation{p_end}
{synopt:{cmd:e(lambda)}}lambda provided or selected by cross-validation{p_end}
{synopt:{cmd:e(r2)}}R-squared{p_end}
{synopt:{cmd:e(cvmse_minimal)}}minimal cross-validation mean squared error{p_end}
{synopt:{cmd:e(cvmse_actual)}}cross-validation mean squared error for selected lambda (only differs from cvmse_minimal if option lambda1se is used){p_end}
{synopt:{cmd:e(numfolds)}}number of folds used in cross-validation{p_end}
{synopt:{cmd:e(numalpha)}}number of alpha tested in cross-validation{p_end}
{synopt:{cmd:e(numlambda)}}number of lambda tested in cross-validation{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}the command used: {cmd:elasticregress}, {cmd:lassoregress} or {cmd:ridgeregress}{p_end}
{synopt:{cmd:e(varlist_nonzero)}}list of covariates with non-zero coefficients{p_end}
{synopt:{cmd:e(properties)}}{cmd:b}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}mark estimation sample{p_end}
{p2colreset}{...}


{marker author}{...}
{title:Author}

{pstd}Wilbur Townsend{p_end}
{pstd}wilbur.townsend@gmail.com{p_end}


{marker references}{...}
{title:References}

{pstd}
Jerome Friedman, Trevor Hastie and Rob Tibshirani. (2008). Regularization Paths for Generalized Linear Models via Coordinate Descent.
Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010.
