

/* To do before releaasing to the OpLab:
	+ should we demean y in cross-validation?? and standardise x??
	+ export more scalars.
	+ check remaining xx's.
	+ comment lots.
	+ use cross() more.
	+ check on glmnet data.
	+ get code-review.
*/

/* To do before public release:
	+ ensure can handle factor variables.
	+ add capacity to handle logistic regressions.
	+ add option to sort results by absolute coefficient size.
*/


*! version 0.1
program define elasticreg, eclass byable(recall)
	version 15

syntax varlist(min=3 numeric) [if] [in] [aweight], alpha(real) [ ///
	lambda(real -1) numlambda(integer 100) lambda1se lambdamin   ///
	numfolds(integer 10) epsilon(real 0.001) tol(real 0.001) ] 

/*
  EXPLAIN SYNTAX XX
  if lambda is -1 (the default) it is found by cross-validation
*/
marksample touse

local depvar : word 1 of `varlist'
local indvars = substr("`varlist'",length("`depvar'")+1,length("`varlist'"))


* If a weight is not provided, set all weights equal.
if "`weight'" == "" {
	tempvar weight
	generate `weight' = 1
}

* Count observations.
summarize `touse', meanonly
local N = _N*`r(mean)'

* Assert that alpha is in [0,1].
if (`alpha' < 0) | (`alpha' > 1) {
	display as error `"alpha must be between 0 and 1."'
	exit
}

* Assert that at most one of lambdamin and lambda1se is specified, and set the
* heuristic on the basis of the option provided.
if ("`lambdamin'" == "lambdamin") & ("`lambda1se'" == "lambda1se") {
	display as error `"only one of lambdamax and lambda1se can be specified"'
	exit
}
else if ("`lambdamin'" == "lambdamin") local heuristic min
else if ("`lambda1se'" == "lambda1se") local heuristic 1se
else                                   local heuristic min

* Assert that lambda is non-negative.
if (`lambda' < 0 & `lambda' != -1) {
	display as error `"lambda cannot be negative."'
	exit
}

* Assert that numlambda is > 1 when lambda is not provided.
if (`numlambda' < 2 & `lambda' == -1) {
	display as error `"numlambda must be greater than 1."'
	exit
}

* Assert that epsilon is > 0. (Note that we cannot have a lambda series ending 
* in 0 as we require the lambda to be log-equidistant. When lambda is provided,
* lambda = 0 is acceptable.)
if (`epsilon' <= 0) {
	display as error `"epsilon must be positive."'
	exit
}

* Assert that tol is >= 0.
if (`tol' < 0) {
	display as error `"tol cannot be negative."'
	exit
}

* Assert that numfolds is >=2 and <= N.
if !inrange(`numfolds', 2, `N') {
	display as error `"The number of folds must be between 2 and N."'
	exit
}

* Format the data: 
* ... setting beta_0 equal to the dependent variable's mean and making the
*     dependent variable mean zero,
tempvar depvar_demeaned 
summarize `depvar' if `touse', meanonly
local ymean = r(mean)
generate `depvar_demeaned' = `depvar' - r(mean)
* ... making the weights sum to 1,
tempvar weight_sum1
summarize `weight' if `touse', meanonly
generate `weight_sum1' = `weight'/(`N' * r(mean))
* ... and standardising the x's (storing their standard deviations and means).
tempname mean_x
tempname sd_x
local K = wordcount("`indvars'")
matrix `sd_x'   = J(`K', 1, .)
matrix `mean_x' = J(`K', 1, .)
forvalues j = 1/`K' {
	local xname : word `j' of `indvars'
	quietly summarize `xname' if `touse'
	matrix `mean_x'[`j',1] = `r(mean)'
	matrix `sd_x'[`j',  1] = `r(sd)'
	tempname `xname'_std
	generate ``xname'_std' = (`xname' - `r(mean)')/`r(sd)'
	local indvars_std `indvars_std' ``xname'_std'
}

* Estimate the regression within Mata, storing outputs in temporary matrices.
tempname beta_handle lambda_handle beta0
mata: notEstimation("`depvar_demeaned'", "`indvars_std'", "`weight_sum1'",      ///
					`alpha', `numfolds', `numlambda', `lambda',	"`heuristic'",  ///
					`epsilon', `tol',					                        ///
					"`beta_handle'", "`lambda_handle'")
* Replace the estimated beta with one corresponding to the unstandardised
* variables and note the list of the non-zero covariates.
forvalues j = 1/`K' {
	matrix `beta_handle'[`j',1] = `beta_handle'[`j',1] / `sd_x'[`j',1]
	if `beta_handle'[`j',1] != 0 {
		local xname : word 1 of `indvars'
		local varlist_nonzero `varlist_nonzero' `xname'
	}
}
* Calculate the intercept as mean(Y) - beta' * mean(X) xx check this is standard.
matrix `beta0' = `ymean' - `beta_handle''*`mean_x'
matrix `beta_handle' = (`beta_handle' \ `beta0')
* Set beta rownames.
matrix rownames `beta_handle' = `indvars' _cons

* Return the covariates of interest.
matrix `beta_handle' = `beta_handle''
ereturn post `beta_handle' , depname(`depvar') obs(`N') esample(`touse')
ereturn scalar lambda = `lambda_handle'
ereturn scalar alpha  = `alpha'
ereturn local varlist_nonzero `varlist_nonzero'
ereturn display
	
end




version 15
mata:

// This function load our data into Mata and calls our estimation subroutine. It
// then stores result matrices available to Stata.
void notEstimation(
	string scalar y_var, string scalar x_varlist, string scalar weight_var,
	real scalar alpha, real scalar numfolds, real scalar numlambda, 
	real scalar lambda, string scalar heuristic,
	real scalar epsilon,  real scalar tol,
	string scalar beta_handle, string scalar lambda_handle)
{
//

	// Import data into Mata.
	real matrix x
	real colvector y
	real colvector weight
	st_view(y,      ., y_var)     
	st_view(x,      ., x_varlist) 
	st_view(weight, ., weight_var) 
	// Calculate the full sample weighted covariance between each independent 
	// variable and the dependent variable. (This is used both when calculating
	// the series of lambda and, after cross-validation, when estimating the
	// final beta, so for efficiency we calculate it only once).
	cov_xy = x' * (weight :* y)
	// Select the series of lambda for which we will estimate beta. If lambda
	// is provided, this is trivial. If not, the series depends on the data.
	if (lambda==-1) lambda_vec = findLambda(cov_xy, alpha, numlambda, epsilon, tol)
	else            lambda_vec = (lambda)
	// If lambda is not provided, select the MSE-minimising lambda using 
	// cross-validation.
	if (lambda==-1) lambda = crossValidateLambda(numfolds, heuristic, 
		x, y, weight, alpha, tol, lambda_vec)
	// Estimate the beta on the full data, given the lambda selected.
	beta = findAllBeta(x, y, weight, alpha, tol, lambda, cov_xy)	
	// Store lambda, beta.)
	st_matrix(beta_handle, beta) 
	st_numscalar(lambda_handle, lambda) 
	
}

// This function calculates the series of lambda for which our model will be
// estimated when no lambda has been provided. It calculates the largest lambda
// such that all beta are guaranteed zero and then creates a decreasing sequence
// of lambda which are equidistant in log space.
real colvector findLambda(real colvector cov_xy, 
						  real scalar alpha, real scalar numlambda,
						  real scalar epsilon, real scalar tol)
{
	lambda_max      = max(cov_xy)/max((alpha, 0.001))   //xx note this is a variation on the paper
	lambda_min      = epsilon * lambda_max
	loglambda_delta = (log(lambda_max) - log(lambda_min))/(numlambda-1)
	lambda      = lambda_max :/ exp((0..(numlambda-1))*loglambda_delta)'
	return(lambda)
}

// This function calculates the optimal lambda using cross-validation.
real scalar crossValidateLambda(
	real scalar numfolds, string scalar heuristic,
	real matrix x, real colvector y, real colvector weight,
	real scalar alpha, real scalar tol, real colvector lambda_vec)
{
	N = length(y)
	numlambda = length(lambda_vec)
	// Divide the data into numfolds equally-sized (+- one observation) cross-
	// validation subsamples.
	Npersample = floor(N/numfolds)
	cv = J(N, 4, .)
	cv[,1] = 1::N
	cv[,2] = runiform(N, 1)
	_sort(cv, 2)
	cv[,3] = 1::N
	cv[,4] = ceil(cv[,3]:/Npersample)
	cv[,4] = cv[,4] + (cv[,4] :> numfolds) :*  (cv[,3] :- numfolds*Npersample :- cv[,4])
	_sort(cv, 1)
	cv = cv[,4]
	// Now for each subsample...
	MSE = J(numfolds, numlambda, .)
	for (s=1; s<=numfolds; s++) {
		// estimate beta for all lambda, excluding that subsample,
		beta = findAllBeta(
				select(x, cv:!=s), select(y, cv:!=s), select(weight, cv:!=s), 
				alpha, tol, lambda_vec, .)
		// extrapolate beta to that subsample for all lambda and store within
		// an N_s x numlambda matrix,
		r = select(y, cv:==s) :- select(x, cv:==s)*beta
		// and calculate the weighted MSE for each lambda.
		for (l=1; l<=numlambda; l++) MSE[s,l] = 
										r[,l]'*(r[,l]:*select(weight, cv:==s))
	}
	// Then collapse MSE(k,lambda) to mean MSE(lambda) and se MSE(lambda).
	MSEmean = mean(MSE)
	// xx cf p27 of Athie/Imbens NBER slides. should we use the sample s.d.? I think so?
	MSEse   = sqrt(mean((MSE:-MSEmean):^2):/numlambda)
	// Then apply the heuristic by which we select lambda -- the traditional
	// choice is the MSE-minimising lambda but a more severe choice is the
	// maximal lambda within a standard error of the MSE-minimising lambda.
	minMSE = selectindex(MSEmean :== min(MSEmean))
	if (heuristic == "min")       lambda = lambda_vec[minMSE]
	else if (heuristic == "1se")  lambda = lambda_vec[max((MSEmean :< (MSEmean+MSEse)[minMSE]) :* (1..numlambda))]
	else _error("Heuristic must be 'min' or '1se'.")
	// Warn the user if the MSE-minimising lambda is the smallest lambda tested.
	if (lambda_vec[minMSE] == lambda_vec[length(lambda_vec)]) { 
		display("Warning: the smallest λ tested was the MSE-minimising λ.")
		display("Consider re-running estimation with a smaller epsilon.")
	}
	return(lambda)
}

// This function finds beta for a series of lambda.
real matrix findAllBeta(
	real matrix x, real colvector y, real colvector weight,
	real scalar alpha, real scalar tol, real colvector lambda,
	real colvector cov_xy)
{
	K         = cols(x)
	numlambda = length(lambda)
	// Calculate wx2: the K vector with kth element equal to the inner
	// product of weight and x_k^2. 
	wx2    = (x:^2)' * weight
	// If cov_xy (the K vector with kth element equal to the 'triple inner
	// product' of weight, x_k and y) has not been provided, calculate that too.
	if (cov_xy == .) cov_xy = x'      * (weight :* y)
	// Instantiate a matrix storing the (weighted) covariances of the
	// x's -- but leave this empty, as we only fill it as necessary.
	cov_x = J(K, K, 0)
	// Store the series of beta in a K x numlamda matrix, with Kth column equal
	// to the beta found from the Kth lambda.
	beta = J(K, numlambda, 0)
	// If estimating for a series of lambda, beta should remain at zero when
	// calculated given lambda_max = lambda[1]. (This doesn't hold for the 
	// limiting case alpha = 0). Each subsequent column of beta is calculated
	// starting from the previous.
	if (numlambda > 1) {
		if (alpha > 0.001) assert(findBeta(x, weight, J(K,1,0), lambda[1],
												alpha, tol, cov_x, cov_xy, wx2) 
										== J(K,1,0))
		for (l=2; l<=numlambda; l++) {
					beta[,l] = findBeta(x, weight, beta[,l-1], lambda[l],
												alpha, tol, cov_x, cov_xy, wx2)
		}
	}
	// Otherwise we estimate beta(lambda) starting from beta = 0. Note that this
	// isn't necessarily efficient -- sometimes hot starts are quicker.
	else beta = findBeta(x, weight, J(K,1,0), lambda[1],
												alpha, tol, cov_x, cov_xy, wx2)
	
	// Return the the beta matrix.
	return(beta)
}

// This function continues refining beta until the set of variables is stable
// and some tolerance is reached.
real colvector findBeta(
	real matrix x, real colvector weight, real colvector beta_start,
	real scalar lambda, real scalar alpha, real scalar tol,
	real matrix cov_x, real colvector cov_xy, real colvector wx2)
{
	beta          = beta_start
	beta_previous = .
	elements      = .
	// We loop,
	while (1) {	
		// in each loop, starting with the full set of variables,
		elements_previous = elements
		elements          = 1::length(beta)
		// improving beta but not re-checking invariant elements of beta, 
		// until beta stabilises.
		while (norm(beta:-beta_previous) > tol) {
			beta_previous = beta
			covUpdate(x, weight, beta, elements, lambda, alpha, cov_x, cov_xy, wx2)
		}
		// and we only stop looping when the set of variables has stabilised.
		if (elements == elements_previous) break
	}
 	// This function modifies the cov_x matrix and returns beta.
	return(beta)	
}

// This function updates certain elements of the beta vector (those indexed by 
// the elements vector) using the covariance update formula. It modifies the
// beta, elements and cov_x matrices and returns null.
void covUpdate(
	real matrix x, real colvector weight,
	real colvector beta, real colvector elements,
	real scalar lambda, real scalar alpha,
	real matrix cov_x, real colvector cov_xy, real colvector wx2)
{	// Loop over the selected elements of beta,
	K = length(elements)
	for (k=1; k<=K; k++) {
		// finding the element of beta,
		j = elements[k]
		// determining if the element changes,
		beta_old_j = beta[j]
		z = cov_xy[j] - cov_x[j,]*beta + beta_old_j
		beta_new_j = softThreshold(z, lambda*alpha) / (wx2[j] + lambda*(1-alpha))
		// if so, updating the beta vector and adding the corresponding
		// covariate to the cov_x matrix,
		 if (beta_old_j != beta_new_j) {
			beta[j] = beta_new_j
			if (beta_old_j == 0) updateCovX(x, weight, j, cov_x)
		 }
		 // and if not, removing the element from the elements vector.
		 else elements[j] = 0
	}
	elements = select(elements, elements)
}

// This is the soft-thresholding operator.
real scalar softThreshold(real scalar z, real scalar gamma)
{
	if      ((z > 0) & (gamma < z))  return(z-gamma)
	else if ((z < 0) & (gamma < -z)) return(z+gamma)
	else                             return(0)
}


// This function adds to the covariance matrix of covariates when beta[j] has 
// been made non-zero.
void updateCovX(
	real matrix x, real colvector w, real scalar j, real matrix cov_x)
{	
	this_cov = (w:*x[,j])'*x
	cov_x[j,] = this_cov
	cov_x[,j] = this_cov'
	assert(issymmetric(cov_x))
}

// xx delete this later -- it's useful for debugging.
void matlist(
    real matrix X,
    | string scalar fmt
    )
{
    real scalar     i, j, wd, rw, cw
    string scalar   sfmt

    if (fmt=="") fmt = "%g"
    wd = strlen(sprintf(fmt,-1/3))

    if (length(X)==0) return

    rw = trunc(log10(rows(X))) + 1
    cw = trunc(log10(cols(X))) + 1
    wd = max((cw,wd)) + 2
    sfmt = "%"+strofreal(wd)+"s"

    printf("{txt}"+(2+rw+1+1)*" ")
    for (j=1;j<=cols(X);j++) {
        printf(sfmt+" ", sprintf("%g", j))
    }
    printf("  \n")
    printf((2+rw+1)*" " + "{c TLC}{hline " +
        strofreal((wd+1)*cols(X)+1) + "}{c TRC}\n")
    for (i=1;i<=rows(X);i++) {
        printf("{txt}  %"+strofreal(rw)+"s {c |}{res}", sprintf("%g", i))
        for (j=1;j<=cols(X);j++) {
            printf(sfmt+" ",sprintf(fmt, X[i,j]))
        }
        printf(" {txt}{c |}\n")
    }
    printf((2+rw+1)*" " + "{c BLC}{hline " +
        strofreal((wd+1)*cols(X)+1) + "}{c BRC}\n")
}


end
