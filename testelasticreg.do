clear
set obs 1000000
gen     x1 = rnormal(0,2) > 1
replace x1 = 2 if rnormal(0,2) > 1
gen x2 = rnormal(4,1)
gen x3 = rnormal(5,3)
gen u  = rnormal(0,4)
gen y  = 0.2*(x1==1) + 0.4*(x1==2) + 0.5*x2 +u+6

reg y i.x1 x2 x3
elasticreg y x1 x2 x3, alpha(1)
