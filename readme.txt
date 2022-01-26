APPENDIX

The pdf file AppendixSupp.pdf contains the Appendix to the article to be published as online supplementary material. This contains all proofs and further numerical results. 

DATA

The excel datafile microbiome_data contains counts of 18 bacterial phyla from 92 individuals who were infected with helminths in 2008 (two observations were deleted from the file because they were deemed to be outliers in our analysis). For further details see Martin et al. (2018) which is a cited reference in the article.


CODE

The code calculates the score matching estimates of the parameters for the PPI model. Standard error estimates of the parameters are also calculated for the PPI model with beta fixed in advance. The code is in R and we include all relevant functions. 

1. functions.R
This file contains all the R functions. It needs to be run first because all the other files below use these functions.

2. figure_1_and_2.R
This file produces Figures 1 and 2 in the article.

3. h_calculations.R
This file calculates the percentiles of the empirical distribution of the h weight functions under the PPI model in Section 2.3 with beta=(-0.8,-0.8,-0.5). The values of the cap a_c that were used to produce Tables 1 and 2 in the article can be extracted from the output of this file. With minor adjustments to the code (i.e simulating under other models) the percentiles in brackets in Table 11 in Appendix A.11 (SM) can also be calculated. 

4. modelA.R
This file simulates a single sample of size n=100 from the PPI model in Section 2.3 with beta=(-0.8,-0.8,-0.5). This file uses the weight function h given by (12) in the article and the cap a_c is set to 20%. First the score matching estimator is calculated for the full model with only beta_p fixed at -0.5. This estimate is denoted by EST in Tables 1 and 2. Then the score matching estimator is calculated assuming beta is fixed in advance. This estimate is denoted by FIX in Tables 1 and 2. If this code is repeated 1000 times with different seeds, then the simulation results for n=100, beta=(-0.8,-0.8,-0.5) and a_c=20% in Tables 1 and 2 (and Tables 7 and 8 in Appendix A.10 (SM)) are easily generated. With minor adjustments to the parameters, sample sizes and a_c choices, all the results in Tables 1 and 2 (and Tables 7 and 8 in Appendix A.10 (SM)) can be generated.

5. model1.R
This file simulates a single sample of size n=92 from either Model 1 or Model 13; see Appendix A.11 (SM) in the article. It then calculates Score1ac (plus standard errors), Score2ac (plus standard errors), Score2 and Score1. See Appendix A.11 (SM) in the article for the definitions of these estimators. If this code is repeated 1000 times with different seeds, then the simulation results in Table 12 are easily generated. The Table 13 results can also be obtained with minor modification to the parameters and dimensions. Note: ScoreMult for p=3 is given at the end of the code (it is currently commented out).
  

6. model4.R
This file simulates a single sample of size n=1000 from Model 4; see Appendix A.11 (SM) in the article. It then calculates Score1ac (plus standard errors), Score2ac (plus standard errors), Score2 and Score1. See Appendix A.11 (SM) in the article for the definitions of these estimators. If the code is repeated 1000 times with different seeds, then the simulation results for Model 4 in Table 14 are easily generated. With minor modification to the parameters, Models 5 and Model 6 results in Table 14 can also be obtained. 


7. model7.R
This file simulates a single sample of size n=92 from Model 7; see Appendix A.11 (SM) in the article. It then calculates Score1, Score2, Score1ac, Score2ac, MLdir, MLDirMult and MomDir. See Appendix A.11 (SM) in the article for the definitions of these estimators. If the code is repeated 1000 times with different seeds, then the simulation results for model 7 in Table 15 are easily generated. With minor modification to the parameters and dimensions, Models 16 and Models 8-12 results in Tables 15 and 16 can also be obtained. 

8. real_data_analysis.R
This file computes all the numbers in Section 7 of the article. This includes the parameter estimates and standard errors for the PPI model given in Tables 3 and 4, all the goodness of fit test output and Figures 3 and 4. 






