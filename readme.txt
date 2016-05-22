Matlab code for M. Zhou, O. H. M. Padilla and J. G. Scott, "Priors for Random Count Matrices Derived from a Family of Negative Binomial Processes," arXiv:1404.3331, to appear in the Journal of the American Statistical Association (Theory and Methods), 2015.

(1)Demo code for naive Bayes classifiers using the negative binomial process, the gamma-negative binomial process, and the beta-negative binomial process:

Demo_NBP_samples.m

(2) Demo code for the 20-newsgroups dataset with the default by-date training/testing partition:

Demo_Laplace_20newsTrainTest.m
Demo_NBP_samples_20newsTrainTest.m

(3) Infer the parameters for a random count matrix:

NBP_Train.m
GNBP_Train.m
BNBP_Train.m

(5) Calculate the predictive likelihood of a count vector under a random count matrix:

predict_NBP_Par.m
predict_GNBP_Par.m
predict_BNBP_Par.m

(6) Calculate the log likelihoods of various distributions:
poisspdf_log.m
nbinpdf_log.m
Logrithmic_pdf_log.m
gammaNB_pdf_log.m
sumStirling_log.m
betaNB_pdf_log.m
digamma_pdf_log_unnormalized.m

(7) Calculate Stirling numbers of the first kind in the log space:
LogFmatrix.m

(8) Generate a logbeta random variable:
logBeta_rnd
psi_complex

(9) C mex function to sample Chinese restaurant table random variables:
CRT_sum_mex_matrix.c
cokus.c
cokus.h

(10) Example datasets:
20news-bydate data set: http://qwone.com/~jason/20Newsgroups/20news-bydate-matlab.tgz
        
Top 30 categories in TDT2: http://www.cad.zju.edu.cn/home/dengcai/Data/TDT2/TDT2.mat
    
CNAE-9: https://archive.ics.uci.edu/ml/machine-learning-databases/00233/CNAE-9.data


Copyright (c), 2015, Mingyuan Zhou
http://mingyuanzhou.github.io/

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.