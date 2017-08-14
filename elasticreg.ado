*! version 0.2
program define elasticreg, eclass byable(recall)
	version 14

syntax varlist(min=2 numeric fv) [if] [in] [aweight], [             ///
	alpha(real -1)  numalpha(integer 6)                             ///
	lambda(real -1) numlambda(integer 100) lambdamin lambda1se      ///
	numfolds(integer 10) epsilon(real 0.001)                        ///
	tol(real 0.001) collinear] 

/*
  alpha is the weight placed on the L1 (LASSO) constraint and (1-alpha) is the 
	weight placed on the L2 (ridgereg) constraint. if alpha is -1 (the
	default) it is found by cross-validation.
  numalpha is the number of alpha for which beta is calculated when alpha is
	being found via cross-validation.
  lambda is the penalty for placed on larger coefficients. if lambda is -1 (the
	default) it is found by cross-validation.
  numlambda is the number of lambda for which beta is calculated when lambda is
	being found via cross-validation.
  If lambdamin, the lambda selected by cross-validation is that which minimises
	MSE in the cross-validation samples. This is the default.
  If lambda1se, the lambda selected by cross-validation is the largest within
	a standard error of that selected under lambdamin.
  numfolds is the number of folds used when cross-validating lambda. Its default
	is 10.
  epsilon is ratio of the smallest lambda tested to the largest. Its default is
	0.001. The user will receive a warning if epsilon seems to be binding.
  tol is the tolerance used when optimising beta.
  If collinear, the program is estimated including collinear variables.
*/
marksample touse

* Seperate the dependant variable and the independant variable, and assert that
* the dependant variable isn't a factor variable.
local depvar : word 1 of `varlist'
local indvars = substr("`varlist'",length("`depvar'")+1,length("`varlist'"))
_fv_check_depvar `depvar'

* Expand factor variables if they are specified, and remove collinear variables
* unless the user has decided not to.
if "`collinear'" == "collinear"  fvexpand `indvars' if `touse'
else                             _rmcoll  `indvars' if `touse', expand
local indvars_uncoded  `r(varlist)'
* Create temporary variables for collinear and factor variables.
fvrevar `indvars_uncoded'
local indvars_coded    `r(varlist)'

* If a weight is not provided, set all weights equal.
if "`weight'" == "" {
	tempvar weight
	quietly generate `weight' = 1
}

* Count observations.
summarize `touse', meanonly
local N = _N*`r(mean)'

* Assert that alpha is in [0,1].
if ((`alpha' < 0) | (`alpha' > 1)) & (`alpha' != -1) {
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

* Assert that numalpha is > 1 when alpha is not provided.
if (`numalpha' < 2 & `alpha' == -1) {
	display as error `"numalpha must be greater than 1."'
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
quietly generate `depvar_demeaned' = `depvar' - r(mean)
* ... making the weights sum to 1,
tempvar weight_sum1
summarize `weight' if `touse', meanonly
quietly generate `weight_sum1' = `weight'/(`N' * r(mean))
* ... and standardising the x's (storing their standard deviations and means).
tempname mean_x
tempname sd_x
local K = wordcount("`indvars_coded'")
matrix `sd_x'   = J(`K', 1, .)
matrix `mean_x' = J(`K', 1, .)
forvalues j = 1/`K' {
	local xname : word `j' of `indvars_coded'
	quietly summarize `xname' if `touse'
	matrix `mean_x'[`j',1] = `r(mean)'
	local sd_unless0 = cond(`r(sd)'==0, 1, `r(sd)')
	matrix `sd_x'[`j',  1] = `sd_unless0'
	tempname x_std
	quietly generate `x_std' = (`xname' - `r(mean)')/`sd_unless0'
	local indvars_std `indvars_std' `x_std'
}

* Estimate the regression within Mata, storing outputs in temporary matrices.
tempname beta_handle lambda_handle alpha_handle r2_handle ///
	     cvmse_minimal_handle cvmse_actual_handle beta0
mata: notEstimation(                                                        ///
	"`depvar_demeaned'", "`indvars_std'", "`weight_sum1'", "`touse'",       ///
	`numfolds', `alpha', `numalpha', `lambda', `numlambda', "`heuristic'",  ///
	`epsilon', `tol',					                                    ///
	"`beta_handle'", "`lambda_handle'", "`alpha_handle'",                   ///
	"`r2_handle'", "`cvmse_minimal_handle'", "`cvmse_actual_handle'")
* Replace the estimated beta with one corresponding to the unstandardised
* variables and note the list of the non-zero covariates.
forvalues j = 1/`K' {
	matrix `beta_handle'[`j',1] = `beta_handle'[`j',1] / `sd_x'[`j',1]
	if `beta_handle'[`j',1] != 0 {
		local xname : word `j' of `indvars_uncoded'
		local varlist_nonzero `varlist_nonzero' `xname'
	}
}
* Calculate the intercept as mean(Y) - beta' * mean(X)
matrix `beta0' = `ymean' - `beta_handle''*`mean_x'
matrix `beta_handle' = (`beta_handle' \ `beta0')
* Set beta rownames.
matrix rownames `beta_handle' = `indvars_uncoded' _cons

* Return the covariates of interest.
matrix `beta_handle' = `beta_handle''
ereturn post `beta_handle' , depname(`depvar') obs(`N') esample(`touse')
ereturn scalar alpha  = `alpha_handle'
ereturn scalar lambda = `lambda_handle'
ereturn scalar r2     = `r2_handle'
ereturn scalar cvmse_minimal  = `cvmse_minimal_handle'
ereturn scalar cvmse_actual   = `cvmse_actual_handle'
if (`lambda' != -1 & `alpha' != -1) ereturn scalar numfolds       = .
else                                ereturn scalar numfolds       = `numfolds'
if `alpha' != -1  ereturn scalar numalpha      = .
else              ereturn scalar numalpha      = `numalpha'
if `lambda' != -1 ereturn scalar numlambda      = .
else              ereturn scalar numlambda      = `numlambda'
ereturn local varlist_nonzero `varlist_nonzero'
ereturn local cmd "elasticreg" 

* Display the results.
if      `alpha' == 1 local title LASSO regression
else if `alpha' == 0 local title Ridge regression
else                 local title Elastic-net regression
display _newline as text "`title'" _continue
local textpreface _col(40) as text 
local intpreface  _col(67) "= " as res %10.0fc 
local realpreface _col(67) "= " as res %10.4f 
display `textpreface' "Number of observations" `intpreface'  e(N)
display `textpreface' "R-squared"              `realpreface' e(r2)
display `textpreface' "alpha"                  `realpreface' e(alpha)
display `textpreface' "lambda"                 `realpreface' e(lambda)
if (`lambda' == -1 | `alpha' == -1 )  {
  display `textpreface' "Cross-validation MSE"    `realpreface' e(cvmse_actual)
  display `textpreface' "Number of folds"         `intpreface'  e(numfolds)
}
if (`alpha'  == -1) display `textpreface' "Number of alpha tested"  `intpreface'  e(numalpha)
if (`lambda' == -1) display `textpreface' "Number of lambda tested" `intpreface'  e(numlambda)
_coef_table

end




version 14
mata:

// This function load our data into Mata and calls our estimation subroutine. It
// then stores result matrices available to Stata.
void notEstimation(
	string scalar y_var, string scalar x_varlist,
	string scalar weight_var, string scalar touse_var,
	real scalar numfolds,
	real scalar alpha, real scalar numalpha,
	real scalar lambda, real scalar numlambda, string scalar heuristic,
	real scalar epsilon,  real scalar tol,
	string scalar beta_handle, string scalar lambda_handle, string scalar alpha_handle, string scalar r2_handle,
	string scalar cvmse_minimal_handle, string scalar cvmse_actual_handle	
	)
{
	// Import data into Mata.
	real matrix x
	real colvector y
	real colvector weight
	st_view(y,      ., y_var,      touse_var)     
	st_view(x,      ., x_varlist,  touse_var) 
	st_view(weight, ., weight_var, touse_var) 	
	// Calculate the full sample weighted covariance between each independent 
	// variable and the dependent variable. (This is used both when calculating
	// the series of lambda and, after cross-validation, when estimating the
	// final beta, so for efficiency we calculate it only once).
	cov_xy = cross(x, weight, y)
	// If alpha and lambda are both provided, we select them and move on. If not
	// we loop first over a vector of possible alphas (which will be a singleton
	// when alpha is provided but lambda is not) and, for each alpha, calculate
	// the optimal lambda.
	if (lambda!=-1 & alpha!=-1) {
		CVMSE = (.\.)
		lambda_found = lambda
		alpha_found  = alpha
	}
	else {
		// If alpha is not provided, we loop over a series of alpha. Otherwise
		// this loop is trivial. Across each interation we store CVMSE[1,2] and
		// lambda_found.
		if (alpha == -1) alpha_vec = rangen(0, 1, numalpha)
		else             alpha_vec = (alpha)
		lambda_found_vec = J(length(alpha_vec), 1, .)
		CVMSE_vec        = J(length(alpha_vec), 2, .)
		// To ensure cross-validation samples are consistent we reset the seed
		// in each loop. To ensure that we get varying estimates when the seed
		// has not been set by the user, the value we seed to is random.
		seed = rdiscrete(1,1,J(10000,1,0.0001))
		for (alphaindex=1; alphaindex<=length(alpha_vec); alphaindex++) {
			thisalpha = alpha_vec[alphaindex]
			rseed(seed)
			// Select the series of lambda for which we will estimate beta. If
			// lambda is provided, this is trivial. If not, the series depends
			// on the data.
			if (lambda==-1) lambda_vec = findLambda(cov_xy, thisalpha, numlambda, epsilon, tol)
			else            lambda_vec = (lambda)
			// We now run a function which finds the MSE-minimising lambda using 
			// cross-validation. (CVMSE is a vector which crossValidateLambda 
			// alters to include the minimal cross-validation mean error and the
			// cross-validation mean error corresponding to the selected lambda.)
			// If lambda is provided but alpha isn't, this function is still
			// useful because it provides the CVMSE. In practice it'd be pretty
			// weird to provide lambda but not alpha.
			thisCVMSE = (.\.)
			lambda_found_vec[alphaindex] = crossValidateLambda(
				  numfolds, heuristic, x, y, weight, thisalpha, tol, lambda_vec,
				  thisCVMSE)
			CVMSE_vec[alphaindex, ] = thisCVMSE'
		}
		// We select alpha on the basis of CVMSE[1, ], regardless of heuristic.
		minMSEindex  = selectindex(CVMSE_vec[,1] :== min(CVMSE_vec[,1]))
		alpha_found  = alpha_vec[       minMSEindex]
		lambda_found = lambda_found_vec[minMSEindex]
		CVMSE        = CVMSE_vec[       minMSEindex, ]'
	}
	// Estimate the beta on the full data, given the lambda selected.
	beta = findAllBeta(x, y, weight, alpha_found, tol, lambda_found, cov_xy)	
	// Calculate the weighted r2.
	r2  = 1 - norm(y - x*beta)^2/norm(y)^2
	// Store lambda, beta, the minimal cross-validation MSE and the selected
	// cross-validation MSE.
	st_matrix(beta_handle, beta) 
	st_numscalar(lambda_handle, lambda_found) 
	st_numscalar(alpha_handle, alpha_found) 
	st_numscalar(r2_handle, r2) 
	st_numscalar(cvmse_minimal_handle, CVMSE[1]) 
	st_numscalar(cvmse_actual_handle,  CVMSE[2]) 
	
}


// This function calculates the series of lambda for which our model will be
// estimated when no lambda has been provided. It calculates the largest lambda
// such that all beta are guaranteed zero and then creates a decreasing sequence
// of lambda which are equidistant in log space. (When alpha is 0, there is no
// largest lambda and so we calculate the one corresponding to alpha = 0.001,
// we will later warn the user if this seems to low.)
real colvector findLambda(real colvector cov_xy, 
						  real scalar alpha, real scalar numlambda,
						  real scalar epsilon, real scalar tol)
{
	lambda_max      = max(cov_xy)/max((alpha, 0.001)) 
	lambda_min      = epsilon * lambda_max
	// We allow for the trivial case when the empirical correlation between all
	// of the x's and y is 0 by setting lambda equal to zero.
	if (lambda_max == 0) lambda = 0
	else {
		loglambda_delta = (log(lambda_max) - log(lambda_min))/(numlambda-1)
		lambda      = lambda_max :/ exp((0..(numlambda-1))*loglambda_delta)'
	}
	return(lambda)
}

// This function calculates the optimal lambda using cross-validation.
real scalar crossValidateLambda(
	real scalar numfolds, string scalar heuristic,
	real matrix x, real colvector y, real colvector weight,
	real scalar alpha, real scalar tol, real colvector lambda_vec,
	real colvector CVMSE)
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
		selectedweights = weightNorm(select(weight, cv:!=s))
		beta = findAllBeta(
				standardise(select(x, cv:!=s), selectedweights),
				demean(select(y, cv:!=s), selectedweights), 
				selectedweights, 
				alpha, tol, lambda_vec, .)
		// extrapolate beta to that subsample for all lambda and store within
		// an N_s x numlambda matrix,
		unselectedweights = weightNorm(select(weight, cv:==s))
		r = demean(select(y, cv:==s), unselectedweights) :- 
							standardise(select(x, cv:==s), unselectedweights)*beta
		// and calculate the weighted MSE for each lambda.
		for (l=1; l<=numlambda; l++) MSE[s,l] = 
									 cross(r[,l], unselectedweights, r[,l])
	}
	// Then collapse MSE(k,lambda) to mean MSE(lambda) and se MSE(lambda).
	MSEmean = mean(MSE)
	MSEse   = sqrt(mean((MSE:-MSEmean):^2):/(numlambda-1))
	// Then apply the heuristic by which we select lambda -- the traditional
	// choice is the MSE-minimising lambda but a more severe choice is the
	// maximal lambda within a standard error of the MSE-minimising lambda.
	minMSE = selectindex(MSEmean :== min(MSEmean))
	if (heuristic == "min")       lambda = lambda_vec[minMSE]
	else if (heuristic == "1se")  lambda = lambda_vec[
							selectindex(MSEmean :< (MSEmean+MSEse)[minMSE])[1]]
	else _error("Heuristic must be 'min' or '1se'.")
	// Warn the user if the MSE-minimising lambda is the smallest lambda tested.
	if (lambda_vec[minMSE] == lambda_vec[length(lambda_vec)]) { 
		display("Warning: the smallest λ tested was the MSE-minimising λ.")
		display("Consider re-running estimation with a smaller epsilon.")
	}
	// Warn the user if the MSE-minimising lambda is the largest lambda tested
	// if the maximal lambda was truncated.
	if (lambda_vec[minMSE] == lambda_vec[1] & alpha < 0.001) { 
		display("Warning: the largest λ tested was the MSE-minimising λ.")
		display("elasticreg includes a very large λ for ridge regression but")
		display("here it appears that this λ constraint is binding.")
	}
	// Return the lambda selected and (implicitly) a vector composed of the
	// minimal cross-validation mean error and the cross-validation mean error
	// corresponding to the selected lambda -- these only differ if heuristic
	// == "1se".
	CVMSE = (MSEmean[minMSE] \ MSEmean[selectindex(lambda :== lambda_vec)])
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
	if (cov_xy == .) cov_xy = cross(x, weight, y)
	// Instantiate a matrix storing the (weighted) covariances of the
	// x's -- but leave this empty, as we only fill it as necessary.
	cov_x = J(K, K, 0)
	// Store the series of beta in a K x numlamda matrix, with Kth column equal
	// to the beta found from the Kth lambda.
	beta = J(K, numlambda, 0)
	// If estimating for a series of lambda, beta should remain at zero when
	// calculated given lambda_max = lambda[1]. (This doesn't hold for the 
	// limiting case alpha = 0). Each subsequent column of beta is calculated
	// starting from the previous. (Note that beta(lambda_max) can be non-zero
	// in cross-validation samples because lambda_max takes into account the
	// total sample cov_xy, not the subsample cov_xy. However given that 
	// beta(lambda_max) will be zero in the full sample, for the sake of cross-
	// validation it makes sense to leave to at zero in that case too.)
	if (numlambda > 1) {
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
		if (wx2[j] == 0) beta_new_j = 0
		else beta_new_j = softThreshold(z, lambda*alpha) / (wx2[j] + lambda*(1-alpha))
		if (beta_new_j == .) ///
			_error("Beta became missing. Please report this bug to " +
									  "https://github.com/wilbur-t/elasticreg.")
		// if so, updating the beta vector and adding the corresponding
		// covariate to the cov_x matrix,
		 if (beta_old_j != beta_new_j) {
			beta[j] = beta_new_j
			if (beta_old_j == 0) updateCovX(x, weight, j, cov_x)
		 }
		 // and if not, removing the element from the elements vector.
		 else elements[k] = 0
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

// This produces a de-meaned vector.
real colvector demean(real colvector vec, real colvector w)
{
	return(vec :- vec'*weightNorm(w))
}

// This standardises a matrix of variables (constant variables are left unstandardised).
real matrix standardise(real matrix x, real colvector w)
{	
	xmean = x'*weightNorm(w)
	xsd   = J(cols(x),1,.)
	for (k=1;k<=cols(x);k++) xsd[k,1] = sqrt(((x[,k] :- xmean[k]):^2)'*weightNorm(w))
	if (sum((xsd:==0)) > 0) ///
			   xsd[selectindex(xsd:==0)] = J(length(selectindex(xsd:==0)), 1, 1)
	return((x :- xmean'):/xsd')
}

// This renormalises weights so that they sum to 1.
real colvector weightNorm(real colvector w)
{	
	return(w:/sum(w))
}

// This function adds to the covariance matrix of covariates when beta[j] has 
// been made non-zero.
void updateCovX(
	real matrix x, real colvector w, real scalar j, real matrix cov_x)
{	
	this_cov = cross(x[,j], w, x)
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
