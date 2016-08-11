"Neglecting cliff instability in Antarctic ice sheet models can reduce melting during warming periods" Codes

Reference:
Ruckert, KL, Shaffer, G, Pollard, D, Forest, FE, and Keller, K. Neglecting cliff instability in Antarctic ice sheet   models can reduce melting during warming periods, (in prep.).

This program is distributed in the hope that it will be useful,
but WITH NO WARRANTY (NEITHER EXPLICITNOR IMPLICIT). We are not liable
for the behavior of these codes in your own application. You are free
to share this code so long as the authors(s) and version history remain
intact. 

Kelsey Ruckert, klr324@psu.edu
Yawen Guan, yig5031@psu.edu

 ======================================================================
| For the impatient:(NOTE:This datasets take time to load)             |
|  1. Open terminal, open seperate screen, and navigate to directory   |
|  2. Type: Rscript LHS_plots.R                                        |
|  3. Type: Rscript MCMC_plots.R                                       |
 ======================================================================

Required software:
    R (3.2.1)
    matlab

Required libraries:
    coda
    lhs
    mcmc (matlab script included)
    R.matlab
    RColorBrewer
    ash
    fields
    pscl
    DEoptim
    compiler

Please note that the folder directory MUST be in the same format as when downloaded
otherwise the scripts will not locate the files/scripts needed to run.
 =============================================================================
| How to run MCMC analysis:
|(NOTE:MCMC analysis takes a day to run)
|
|  1. Open terminal, open seperate screen, and navigate to directory
|  2. Type: Rscript Initial_Parameters.R (optional, can use previously created "OptimizedInitialParameters.txt)
|  3. Type: cd DAIS_matlab
|  4. Type: matlab -nodisplay -nosplash < DAIScali_hetero_model_iid_mcmc.m      
|  5. Type: cd..
|  6. Type: Rscript DAIS_convergence_plotsR.R (to check for convergence)
|  7. Type: Rscript DAIScali_hetero_model_iid_mcmcmat.R
|  8. Type: Rscript MCMC_plots.R
|
| How to run LHS analysis:
|(NOTE:LHS analysis takes ~2hrs to run)
|
|  1. Open terminal, open seperate screen, and navigate to directory
|  2. Type: Rscript DAIS_precali_LHS.R
|  3. Type: Rscript LHS_plots.R
 =============================================================================


 =============================================================================
| Short description of main scripts:
|  The main scripts were used in the analysis of the paper.
|
|  1. DAIS_datal.R & DAIS_data.m: Input forcing data, standards, and other important info.
|  2. DAIS_IceFlux_model.R & DAIS_IceFlux_model.m: Antarctic ice sheet model.
|  3. DAIS_precali_LHS.R: Latin hypercube sampling of the AIS model.
|  4. DAIScali_hetero_model_iid_mcmc.m & DAIScali_hetero_model_iid_mcmcRversion.R: Heteroskedastic 
|     IID MCMC fitting of the AIS model (matlab version used in analysis; R version will run for months).      
|  5. DAIS_convergence_plots.R: Tests MCMC convergence based on the Potential Scale Reduction 
|     factor, Trace plots, and Heidelberger and Welch's convergence diagnostic.
|  6. DAIScali_hetero_model_iid_mcmcmat.R: hindcast and projects AIS melt using MCMC data from matlab.
|  7. DAISobs_likelihood_AR.R: Estimates the likelihood functions for the MCMC methods (if using R version).
|  8. log_lik_calibration_copy.m: Estimates the log-likelihood functions for the MCMC methods.
|  9. log_pri_copy.m: Estimates the prior for the MCMC methods.
|  10. LHS_plots.R & MCMC_plots.R: Generates figures shown in paper.
 =============================================================================


Credits:
    AIS Model: 
    Shaffer, G.: Formulation, calibration and validation of the DAIS model (version 1), a simple Antarctic ice sheet
    model sensitive to variations of sea level and ocean subsurface temperature, Geoscientific Model Development, 7(4),
    1803–1818, 2014.


