# Multilinear_Gaussian_Process
===
![project_figure](https://user-images.githubusercontent.com/45757826/57381582-82d73900-71ab-11e9-8cfd-2a2e1f31ea47.png)
> The poster for the "European Research Network System Identification (ERNSI) 2018" is available in the repository which provide clear demonstrations of the proposed method.
# Table of Content
- Gibbs sampling
- Gaussian process
- cv
- "ec_lambda"

# How to use

The source code for the proposed multilinear Gaussian Process is implemented with Matlab. 

The code is the realization of Experiment 1 in the paper. However, one can easily adapt it to Experiment 2 and the Fed-batch experiments.

The source code includes 4 files:

1. "main.m" corresponds to our proposed algorithm.

2. "gibbs estimate.m" is the implementation multilinear GP.

3. cv_lambda.m" is the L1-regularization with penalty lambda selected by the cross validation.

4. ec_lambda.m" is the L1-regularization with an empirical choice of lambda.

The Matlab codes are well documented.


The variables and functions used in the codes agree with those defined in the paper. 
