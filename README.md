# Nonlinear-fitness-with-TRM
We show how transient responses to pulse disturbances accumulate to determine a press disturbances. For stage-structured models, this cumulative change is given by a new Transient Response Matrix (TRM) that can be readily computed. TRM also yields the second derivatives of population growth rate with respect to matrix elements

Much work in ecology has developed methods to help predict how natural populations respond to disturbances, but analyses of pulse (one-time) and press (persistent) disturbances have been largely disconnected. Here we show how transient responses to pulse disturbances accumulate to determine the long-term response to press disturbances. For stage-structured models, this cumulative change is given by a new Transient Response Matrix (TRM) that can be readily computed. Strikingly, the TRM also yields the second derivatives of population growth rate with respect to matrix elements. Thus there is an intimate but unexpected relationship between nonlinear selection pressures on demographic rates, and the transient dynamics of populations. This relationship yields a strong correlation between nonlinear selection and generation time across 439 unique plant and animal species (2690 population models). We also show that the TRM is directly related to Cohen's cumulative distance measure for populations converging to stability. Finally, we indicate how our method is generalizable to stochastic disturbances, and to equilibria in nonlinear models.

The code for running the analysis is attached and the setup is as follows:

Data file for Figure 2 and 4:

1) mpm1.csv - This is matrix used for Figure 2 and Figure 4.

Data files for Figure 5:

2) final_tree_corrected_newversion v2.tre
3) full animal and plant data v2.csv
4) match_list_uncorrected_newversion v2.csv

These files contain phylogenetically corrected tree for the data from Compadre and Comadre databases. We use these to then regress the traits generation time against dominant eigenvalue of the Transient Response Matrix (J0).

Please load all files in the same working directory to run the analysis in the Rscript titled "code for pulse and press paper figures.R"
