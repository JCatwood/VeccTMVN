# VeccTMVN 1.3.1

* Change the optimization method used to optimize psi w.r.t. x, beta to 'BFGS' from 'L-BFGS-B'. 'BFGS' appears more suitable for convex optimization.
* Add the optimization control of above optimization to use the absolute tolerance of 1e-6.
* Increase the max iteration number of the above optimization to 1500.
* Remove the second stage optimization.
* Increase the limit for TMVN/TMVT sampling dimension from 1000 to 2000.
