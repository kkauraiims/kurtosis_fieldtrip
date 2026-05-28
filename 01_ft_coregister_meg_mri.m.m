function ft_coregister_meg_mri (patient_id)
% This function co-registers Elekta neuromag MEG data with the subject's structural MRI 
% the code is executed in FieldTrip
% Input: 
%  The function only requires the patient/subject id as input
%  it is assumed that ft_defaults is already initiated and completed
%  it is assumed that both the MEG (.fif) and MRI (.nii) are in the same patient/subject directory
%  the name of the patient/subject directory is the same as the patient/subject ID
% Output: 
%  (The function automatically generates an output directory from the patient/subject_ID)
%  i. Computes the covariance matrix of the MEG data for creating the downstream leadfield matrix
%  ii. Headmodel 
%  iii. Sourcemodel 
%  iv. Figure depicting the alignment of all geometrical data
% Notes: 
%   1. currently coded for '.nii' MRI files only 
%   2. requires the following manual input: 
%     Do you want to change the anatomical labels for the axes [Y, n]? y
%     What is the anatomical label for the positive X-axis [r, l, a, p, s, i]? r
%     What is the anatomical label for the positive Y-axis [r, l, a, p, s, i]? a
%     What is the anatomical label for the positive Z-axis [r, l, a, p, s, i]? s
%     Is the origin of the coordinate system at the a(nterior commissure), i(nterauricular), n(ot a landmark)? i
%   3. It is important to check the end alignment of the geometric data, misaligned geometry is likely to generate spurious sources
%   4. The pre-processing for kurtosis differs from similar analytical procedures (e.g. computing Hurst exponent) with respect to bp filter.
%   The lp for kurtosis in this function was kept high as 150Hz. 
% Authors: CW, KK June 2022

% specify data and output directories 
patient_dir = strcat ('/path/to/', patient_id); 
patient_output_dir = strcat ('/path/to/', patient_id);
mkdir(patient_output_dir)
cd (patient_dir)

% specify and read mri file from patient folder
PFMRI =dir ('*.nii');
mri_file = PFMRI.name; 
mri_orig = ft_read_mri(mri_file);

% specify and read MEG file 
PFMEG= dir ('*.fif'); 
meg_file = PFMEG.name; 


% automatically detect headshape from .fif file
headshape = ft_read_headshape(meg_file);
% convert dimensions of headshape for further analysis
headshape = ft_convert_units(headshape, 'mm');

% check axis of coordinate system
ft_determine_coordsys(mri_orig, 'interactive', 'no') % x-axis should be right
ft_plot_headshape(headshape);

cfg = [];
cfg.method = 'headshape';
cfg.headshape.interactive = 'yes';
cfg.headshape.icp = 'yes';
cfg.headshape.headshape = headshape;
cfg.coordsys = 'neuromag';
cfg.spmversion = 'spm12';
mri_realigned = ft_volumerealign(cfg, mri_orig);
mri_realigned.coordsys = 'neuromag';

% Do you want to change the anatomical labels for the axes [Y, n]? y
% What is the anatomical label for the positive X-axis [r, l, a, p, s, i]? r
% What is the anatomical label for the positive Y-axis [r, l, a, p, s, i]? a
% What is the anatomical label for the positive Z-axis [r, l, a, p, s, i]? s
% Is the origin of the coordinate system at the a(nterior commissure), i(nterauricular), n(ot a landmark)? i

ft_determine_coordsys(mri_realigned, 'interactive', 'no')
ft_plot_headshape(headshape);
cd (patient_output_dir)
save ('headshape');

%% pre-process MEG 

dataset = strcat (PFMEG(1).folder, '/', PFMEG(1).name);

cfg = [];
cfg.dataset = dataset;
cfg.hpfilter = 'yes';
cfg.hpfreq = 1; % can be changed to 10
cfg.lpfilter = 'yes';
cfg.lpfreq = 150;
cfg.channel = {'megmag', 'meggrad'};
cfg.coilaccuracy = 0;
data = ft_preprocessing(cfg);

% Resample data file to save memory space
cfg = [];
cfg.resamplefs = 500;
data_resampled = ft_resampledata(cfg, data);

%% compute data covariance window 
cfg = [];
cfg.channel = 'MEG';
cfg.covariance = 'yes';
cov_matrix = ft_timelockanalysis(cfg, data_resampled);

cd (patient_output_dir);
save ('cov_matrix');

%% create headmodel 
cfg = [];
cfg.tissue = 'brain';
cfg.spmversion = 'spm12';
seg = ft_volumesegment(cfg, mri_realigned);

save ('seg'); 

cfg = [];
cfg.tissue = 'brain';
cfg.spmversion = 'spm12';
brain_mesh = ft_prepare_mesh(cfg, seg);

cfg = [];
cfg.method = 'singleshell'; % previously was singleshell 
cfg.grad = data_resampled.grad; 
headmodel = ft_prepare_headmodel(cfg, brain_mesh);

save ('headmodel')
clear seg

%% construction of source model 
cfg = [];
cfg.resolution = 10; % for clinical purpose should be 5mm or less
cfg.unit = 'mm';
cfg.headmodel = headmodel;
cfg.grad = data_resampled.grad;
sourcemodel = ft_prepare_sourcemodel(cfg);

save ('sourcemodel') 

%% plot all geometrical data to check their alignment
figure
ft_plot_axes([], 'unit', 'mm', 'coordsys', 'neuromag');
ft_plot_headmodel(headmodel, 'unit', 'mm'); % this is the brain shaped head model volume
ft_plot_sens(data_resampled.grad, 'unit', 'mm', 'coilsize', 10); % the sensor locations
ft_plot_mesh(sourcemodel.pos, 'unit', 'mm'); % the source model is a cubic grid of points
ft_plot_ortho(mri_realigned.anatomy, 'transform', mri_realigned.transform, 'style', 'intersect');
alpha 0.5 % make the anatomical MRI slices a bit transparent

savefig (gcf, 'alignment.fig', 'compact');
close all
