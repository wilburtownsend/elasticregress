# elasticregress
lasso, ridge regression and elastic net regression

To install elasticregress, type	`ssc install elasticregress` into the Stata terminal.
 
The elasticregress package is a Stata implementation of the Friedman,
Hastie and Tibshirani (2010, JStatSoft) coordinate descent algorithm
for elastic net regression and its famous special cases: lasso and
ridge regression. These regression estimates tend to have better 
out-of-sample fit than those from ordinary least squares. When tuning 
parameters are not provided by the user, elasticregress can find them
with K-fold cross-validation.

Keywords:
lasso, ridge regression, elastic net, regularized regression, machine learning, linear regression

Requires: Stata version 13.0

Distribution date: 20170809

Author: Wilbur Townsend, Stanford University

Support: https://github.com/wilbur-t/elasticregress
