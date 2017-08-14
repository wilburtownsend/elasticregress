*! version 0.1

* This command calculates the LASSO estimator of a linear regression. It is a
* wrapper for elasticreg.
program define lassoreg, eclass byable(onecall)
	version 14

syntax varlist(min=2 numeric fv) [if] [in] [aweight], [             ///
	lambda(real -1) numlambda(integer 100) lambda1se lambdamin      ///
	numfolds(integer 10) epsilon(real 0.001) tol(real 0.001) collinear ] 	
	
if _by() local byprefix by `_byvars': 

`byprefix' elasticreg `varlist' `if' `in' [`aweight'], alpha(1)      ///
	lambda(`lambda') numlambda(`numlambda') `lambda1se' `lambdamin'  ///
	numfolds(`numfolds') epsilon(`epsilon') tol(`tol')  `collinear'
	
ereturn local cmd "lassoreg" 
	
end
