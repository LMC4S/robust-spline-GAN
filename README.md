# R code for robust logit fGAN and hinge GAN
[Tractable and Near-Optimal Adversarial Algorithms for Robust Estimation in Contaminated Gaussian Models](https://arxiv.org/abs/2112.12919)

------
* main.R:
Controls input arguments and output RData/log files. Call main.R with Rscript front end for Linux or MacOS or with R CMD BATCH for Windows. For example, for Linux or MacOS try
```{shell}
Rscript --vanilla <path_to_"main.R">  0.1 5 1000 30 1 JS MOSEK .5 0.025 A TRUE <path_to_result_folder> 0 
```
The argument names are contamination proportion *eps*, dimension *p*, sample size *n*, iteration times *iter*, *repeat_times*, *Loss*,  *optim*, *lr*, *pen*, *Q*, *l1*, *wd*, and *seed*. For details see the comments in main.R. Argument p > 20 is not recommanded if using default *ECOS* solver.

* utilities.R:
Utility functions.

* train.R:
Compiles the model settings  and generates experiment data correspondingly. Organizes output. 

* fGAN.R:
Core functions of robust logit fGAN and hinge GAN, including the spline basis generating functions, the generator update function, discriminator update function, and the fGAN itself.

