# SVAR Identification From Higher Moments

This repo contains the MATLAB code for Montiel Olea, Plagborg-Møller, Qian (2022): "SVAR Identification From Higher Moments: Has the Simultaneous Causality Problem Been Solved?" ([paper](https://scholar.princeton.edu/mikkelpm/svar_higher_moments)).

**Development note**: Tested with MATLAB R2020b on macOS Monterey 12.1. The MATLAB code requires the Optimization, Global Optimization, Parallel Computing, and Datafeed toolboxes. The R code requires the steadyICA package.


## Contents

* [empirical](https://github.com/eric-qian/higher_moments/tree/main/empirical): Code for the empirical application. 
	* ``makeData.m``: Fetches and saves data from FRED.
	* ``mainEmpirical.m``: Runs the bootstrap procedure and produces output.
	* ``test_DN.R``: Run the test featured in Davis and Ng (2021).
* [simulation](https://github.com/eric-qian/higher_moments/tree/main/simulation): Code for the simulation study. 	
	* ``higher_moments_simul.m``: Defines the class and methods for simulating the data.
	*  ``mainSimul.m``: Main file for running the simulations.
	*  ``plotSimul.m``: Produces figures using the output of ``mainSimul.m``.
* [functions](https://github.com/eric-qian/higher_moments/tree/main/functions): Routines used in the empirical application and simulation study.
	* ``var_test_indep.m``: Function that executes the bootstrap test. 
	* Also contains other estimation routines and subroutines.

## Replication
1.  **Empirical application**: To replicate the results of the empirical application (as featured in Table 1 of the Appendix), run the following scripts:
	1. Run ``empirical/makeData.m`` to fetch and save the data from FRED.
	2. Next, run ``empirical/mainEmpirical.m`` to produce our test's 5% critical value, 10% critical value, and test statistic.
	3. Run ``empirical/test_DN.R`` to produce the *p*-values for the procedure featured in Davis and Ng (2021).
2. **Simulation study**: After completing Step 1, run the following scripts to replicate Figure 1:
	1. Run ``simulation/mainSimul.m`` to produce the results of the simulation study. The DGP is taken directly from the results of the empirical application. Here, we use MATLAB's Parallel Computing toolbox. 
	2. Equipped with the output of the previous step, ``simulation/plotSimul.m`` generates Figure 1.



