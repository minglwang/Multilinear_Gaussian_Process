# Multilinear_Gaussian_Process


> The [poster](https://github.com/minglwang/Multilinear_Gaussian_Process/blob/master/poster.pdf) for the "European Research Network System Identification (ERNSI) 2018" is available in the repository which provide clears demonstrations of the proposed method.
<!-- 
# Table of Content
- [Description](#description)
- [Multilinear Gaussian Process](multi-linear-gp)
- [L1 Regularization](#l1-regularization)
- [How to use](#how-to-use)-->

<!--# Description-->

<!--# Kinetic model-->

<!--# Multilinear Gaussian Process-->

# How to use

The source code for the proposed algorithm is implemented with Matlab. The code is the realization of Experiment 1. However, one can easily adapt it to Experiment 2 and the Fed-batch experiments. It includes 4 files: 
  - ['main.m'](#) corresponds to the proposed algorithm.
  - ['gibbs_estimate.m'](#) is the implementation of multilinear GP which includes 3 nested functions: 1) expectation maximization,  2) Gibbs sampling, 3) grid search optimization.
  - ['cv_lambda.m'](#) is the L1-regularization with lambda selected by the cross validation in Algorithm 4.
  - ['ec_lambda.m'](#) is the L1-regularization with an empirical choice of lambda.

The codes are well noted. The variables and functions used in the codes agree with those defined in the paper. 

