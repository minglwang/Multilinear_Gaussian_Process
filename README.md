# Multilinear_Gaussian_Process

The source code for the proposed multilinear Gaussian Process is implemented with Matlab. 

The code is the realization of Experiment 1 in the paper. However, one can easily adapt it to Experiment 2 and the Fed-batch experiments.

The source code includes 4 files:

1. "main.m" corresponds to our proposed algorithm.

2. "gibbs estimate.m" is the implementation multilinear GP.

3. cv_lambda.m" is the L1-regularization with penalty lambda selected by the cross validation.

4. ec_lambda.m" is the L1-regularization with an empirical choice of lambda.

The codes are well documented. 

The variables and functions used in the codes agree with those defined in the paper. 
