# kurtosis-fieldtrip

This repository contains FieldTrip/MATLAB scripts for performing kurtosis beamforming on MEG data.

The workflow is designed for source localisation of interictal epileptiform activity in epilepsy MEG recordings. It includes MRI–MEG co-registration, MEG preprocessing, construction of the head model and source model, LCMV beamforming, estimation of source-level kurtosis, and visualisation of kurtosis peaks.

The workflow follows the general approach described for kurtosis beamforming in epilepsy MEG, where kurtosis peaks in source space are used to guide candidate source localisation and downstream visual inspection of virtual electrode activity.

Written and maintained by Caroline Witton and [Kirandeep Kaur](https://github.com/kkauraiims).

## Workflow

### 01_ft_coregister_meg_mri.m

Co-registers Elekta Neuromag MEG data with structural MRI and prepares the geometrical and covariance inputs required for source analysis.

Main outputs:
- preprocessed MEG data covariance matrix
- head model
- source model
- figures for checking MRI, headshape, sensor and source-model alignment

### 02_ft_kurtosis.m

Performs LCMV beamformer source analysis, estimates source-level kurtosis, identifies kurtosis peaks, and visualises outputs in FieldTrip, MRIcro and AnyWave.

Main outputs:
- leadfield matrix
- LCMV source analysis output
- source-level kurtosis estimates
- kurtosis peak visualisations
- exported files for further visual inspection

## Reference

Hall MBH, Nissen IA, van Straaten ECW, Furlong PL, Witton C, Foley E, Seri S, Hillebrand A. An evaluation of kurtosis beamforming in magnetoencephalography to localize the epileptogenic zone in drug resistant epilepsy patients. *Clinical Neurophysiology*. 2018;129(6):1221–1229. doi:10.1016/j.clinph.2017.12.040.
