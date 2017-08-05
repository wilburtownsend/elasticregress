

set seed 1
clear
discard

/*
use "C:\Users\wilb\GitHub\elasticreg\rscfp2013.dta" 
lassoreg income hhsex age educ married kids lf lifecl race 
*/

set obs 100000
gen     x1 = rnormal(0,2) > 1
gen x2 = rnormal(4,1)
gen x3 = rnormal(5,3)
gen x4 = rnormal(0,1)
gen u  = rnormal(0,4)
gen y  = 0.2*x1 + 0.5*x2  + 1*x4 + u + 6
* Calculate a LASSO model.
lassoreg y x*
* Calculate a ridge-regression model.
ridgereg y x*
* Calculate a ridge-regression model, testing smaller lambda.
ridgereg y x*, epsilon(0.00001)
* Calculate OLS  equivalent to lasso or ridge with lambda=0.
lassoreg y x*, lambda(0)


* Allow for missing data.
replace y  = . if rnormal(0,1) > 2
replace x1 = . if rnormal(0,2) > 3
* Calculate OLS  equivalent to lasso or ridge with lambda=0.
lassoreg y x*, lambda(0)
