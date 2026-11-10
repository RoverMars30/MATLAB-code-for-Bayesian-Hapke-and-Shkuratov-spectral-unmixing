# Bayesian Spectral Unmixing With Hapke and Shkuratov Models

This repository contains the MATLAB source code supporting the study:

**Bayesian Spectral Unmixing of Al-Bearing Phyllosilicates and Sulfates in Martian Analog Mixtures With Radiative Transfer Models**

The code implements Hapke and Shkuratov radiative transfer calculations, derives the imaginary refractive index \(k(\lambda)\) from endmember reflectance spectra, and performs Bayesian retrievals of endmember abundances and effective grain sizes from laboratory mixture spectra.


## Repository Contents

### Main inversion program

- `MCMC_RTM_20260827_v6_BatchExcel.m`  
  Batch inversion of mixture spectra using the Hapke and Shkuratov radiative transfer models. The program retrieves endmember abundances and effective grain sizes and saves the posterior samples, maximum a posteriori solutions, modeled spectra, residuals, and spectral root-mean-square errors.

### Endmember optical-constant retrieval

- `run_independent_k_from_excel.m`  
  Reads the grain-size-specific endmember spectra and derives separate \(k(\lambda)\) spectra for each measured size fraction.

- `load_endmember_grain_excel.m`  
  Imports and organizes endmember spectra and their assigned grain-size ranges from the Excel dataset.

- `fit_endmember_k_by_fraction.m`  
  Derives \(k(\lambda)\) independently for each endmember size fraction and calculates the mean of the selected estimates.

- `invert_endmember_k_size.m`  
  Jointly optimizes \(k(\lambda)\) and effective particle diameter within specified bounds.

- `compute_endmember_k_bounded.m`  
  Provides the endmember \(k(\lambda)\) inputs required by the mixture inversion workflow.

### Radiative transfer calculations

- `Hapke_forward_model.m`
- `Hapke_nk_to_ssa.m`
- `Hapke_ssa2reff.m`
- `Hapke_likelihood.m`
- `Shkuratov_forward_model.m`
- `Shkuratov_Calc_rb_rf.m`
- `sk_likelihood.m`

These functions calculate endmember and mixture spectra and evaluate the likelihoods used in the Bayesian inversions.

## Software Requirements

The code was developed in MATLAB and requires:

- MATLAB R2022a or later
- Optimization Toolbox
- The CATMIP/MCMC routines used by Lapôtre et al. (2017)

## Third-Party MCMC Dependency

The directory `MCMC_unmixing_routines` is not distributed in this repository because it contains third-party CATMIP routines originating from the code associated with:

Lapôtre, M. G. A., Ehlmann, B. L., & Minson, S. E. (2017). A probabilistic approach to remote compositional analysis of planetary surfaces. *Journal of Geophysical Research: Planets, 122*, 983–1009.  
https://doi.org/10.1002/2016JE005248

Consequently, the Bayesian batch-inversion program cannot be executed without a compatible implementation of these routines. Users must obtain permission and access to the original routines separately or supply compatible MCMC functions.

## Files Not Included

The following locally generated or third-party files are intentionally excluded:

- `MCMC_unmixing_routines/`: third-party CATMIP routines.
- `Independent_k_Results/`: generated \(k(\lambda)\) spectra, diagnostic figures, and intermediate calibration files.
- `dependency_manifest.csv`: local dependency inventory.

The files in `Independent_k_Results` are analysis outputs rather than source code. They can be regenerated from `Endmember_Various_Grain_Size.xlsx` using `run_independent_k_from_excel.m`.

## General Workflow

1. Download the spectral dataset from Zenodo.
2. Use `run_independent_k_from_excel.m` to derive the endmember \(k(\lambda)\) spectra.
3. Place a compatible implementation of the required CATMIP routines in a directory named `MCMC_unmixing_routines`.
4. Update the local data paths and inversion settings in `MCMC_RTM_20260827_v6_BatchExcel.m`.
5. Run `MCMC_RTM_20260827_v6_BatchExcel.m` and select the directory containing the mixture Excel files.

The principal output for each input workbook is a MATLAB file containing retrieved abundances, effective grain sizes, posterior samples, modeled spectra, residuals, and fit statistics.

## Current Release Status

This repository is an archival source-code snapshot prepared to document the computational methods used in the associated manuscript. The third-party MCMC dependency and generated intermediate products are not included in the present release. A complete executable release may be provided after publication, subject to the redistribution conditions of the third-party routines.


