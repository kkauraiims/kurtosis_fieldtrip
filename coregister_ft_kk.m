function coregister_ft_kk (patient_id)

% A custom function to co-register Elekta MEG data with individual MRI 
% in fieldtrip 
% requires manual inputs 
% currently coded for '.nii' MRI files only 

% KK, June 2022

%% coregister MEG and MRI 

% specify data and output directories 
patient_dir = strcat ('/Users/neelbazro/Desktop/he_db/input','/', patient_id); 
patient_output_dir = strcat ('/Users/neelbazro/Desktop/he_db/output', '/', patient_id);

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

















