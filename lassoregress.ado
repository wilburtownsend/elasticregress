*! version 1.2

* This command calculates the LASSO estimator of a linear regression. It is a
* wrapper for elasticregress.
program define lassoregress, eclass byable(onecall)
	version 13

syntax varlist(min=2 numeric fv) [if] [in] [aweight], [             ///
	lambda(real -1) numlambda(integer 100) lambda1se lambdamin      ///
	numfolds(integer 10) epsilon(real 0.001) tol(real -1) collinear ] 	
	
if _by() local byprefix by `_byvars': 

`byprefix' elasticregress `varlist' `if' `in' [`weight'`exp'], alpha(1)      ///
	lambda(`lambda') numlambda(`numlambda') `lambda1se' `lambdamin'          ///
	numfolds(`numfolds') epsilon(`epsilon') tol(`tol')  `collinear'
	
ereturn local cmd "lassoregress" 
	
end
