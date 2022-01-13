# Identification from Higher Moments

This repo contains the MATLAB code for Montiel Olea, Plagborg-Møller, Qian (2022): "SVAR Identification From Higher Moments: Has the Simultaneous Causality Problem Been Solved?" 

**Development note**: Tested with MATLAB R2020b on macOS Monterey 12.1. 

## Contents

* [empirical](https://github.com/eric-qian/higher_moments/tree/main/empirical): Code for the empirical application. 
	* ``makeData.m``: Fetches and saves data from FRED.
	* ``mainEmpirical.m``: Runs the bootstrap procedure and produces output.
* [simulation](https://github.com/eric-qian/higher_moments/tree/main/simulation): Code for the simulation study. 	
	* ``higher_moments_simul.m``: Defines the class and methods for simulating the data.
	*  ``mainSimul.m``: Main file for running the simulations.
	*  ``plotSimul.m``: Produces figures using the output of ``mainSimul.m``.
* [functions](https://github.com/eric-qian/higher_moments/tree/main/functions): Routines and subroutines used in the empirical application and simulation study.