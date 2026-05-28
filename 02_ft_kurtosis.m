function ft_kurtosis (patient_id) 
%% This function computes source-level kurtosis from MEG data data using 
% LCMV beamformer in FieldTrip
% This function is implemented in continuation with 01_ft_coregister_meg_mri.m
% Input: 
%       (i) The function only requires the patient/subject id as input
%       (ii) it is assumed that the output directory of 01_ft_coregister_meg_mri.m is the 
%       patient/input directory for this function i.e. the input directory contains the 
%       headmodel, cov_matrix and the sourcemodel 
%       (iii) It is assumed that ft_defaults is already intiated and completed
% Output: 
%       (i) LCMV beamformer source analysis 
%       (ii) extraction and visualization of kurtosis peaks in ft, MRIcro and AnyWave
% Notes: 
%        no manual inputs required
%        this function uses Brainnetome template atlas (available in ft) to plot kurtosis peaks
% Authors: CW, KK; July 2022

% Specify the patient/input directory
patient_dir = strcat ('/Path/to/input_dir/', patient_id);
patient_output_dir = strcat ('/Path/to/output_dir/',patient_id);
cd (patient_dir)

load sourcemodel.mat 

% compute leadfield 
cfg = [];
cfg.channel = 'MEG';
cfg.headmodel = headmodel;
cfg.sourcemodel = sourcemodel;
cfg.normalize = 'yes'; % normalisation avoids power bias towards centre of head
cfg.reducerank = 2;
cfg.resolution = 2; 
cfg.unit      = 'cm';
leadfield = ft_prepare_leadfield(cfg, cov_matrix);

save ('leadfield'); 

% plot svd of covariance matrix 

[u,s,v] = svd(cov_matrix.cov);

figure;
semilogy(diag(s),'o-');
savefig (gcf, 'svd_covariance.fig', 'compact');  
close all

% compute the LCMV beamformer 
cfg = [];
cfg.method = 'lcmv';
cfg.sourcemodel = leadfield;
cfg.headmodel = headmodel;
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori = 'yes'; % project on axis of max variance using SVD
cfg.lcmv.lambda = '5%';
cfg.lcmv.kappa = 69;
cfg.lcmv.projectmom = 'yes'; % project dipole time series in direction of maximal power 
cfg.lcmv.kurtosis = 'yes';
cfg.lcmv.keepmom = 'yes';
source = ft_sourceanalysis(cfg, cov_matrix);

save ('source'); 
       
%% interpolate kurtosis into the mri 
cfg = [];
cfg.parameter = 'kurtosis';
source_interp = ft_sourceinterpolate(cfg, source, mri_realigned);

%% load brainnetome atlas
atlas_brainnetome = ft_read_atlas ('/path/to/fieldtrip-20220104/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii');

% plot kurtosis output in 'ortho'
cfg = [];
cfg.funparameter = 'kurtosis';
cfg.method = 'ortho'; % orthogonal slices with crosshairs at peak (default anyway if not specified)
cfg.atlas = atlas_brainnetome; 
ft_sourceplot(cfg, source_interp);

% plot kurtosis output in 'slices'
cfg = [];
cfg.funparameter = 'kurtosis';
cfg.method = 'slice'; % plot slices
ft_sourceplot(cfg, source_interp);

% find regions of max kurtosis
array = reshape(source.avg.kurtosis, source.dim);
array(isnan(array)) = 0;
ispeak = imregionalmax(array); % findpeaksn is an alternative that does not require the image toolbox
peakindex = find(ispeak(:));
[peakval, i] = sort(source.avg.kurtosis(peakindex), 'descend'); % sort on the basis of kurtosis value
peakindex = peakindex(i);

npeaks = 5; %we are interested in the top 5 peaks
disp(source.pos(peakindex(1:npeaks),:));% output positions
poi = (source.pos(peakindex(1:npeaks),:));

% plot peaks 
for i = 1:npeaks
  cfg = [];
  cfg.funparameter = 'kurtosis';
  cfg.location = source.pos(peakindex(i),:)*1000; % convert from m to mm
  cfg.atlas = atlas_brainnetome; 
  ft_sourceplot(cfg, source_interp);
  savefig (gcf, strcat ('kurtosis','_atlas','_',num2str (i), '.fig'), 'compact');
end 
close all

% visualize the kurtosis in MRIcro 
cfg = [];
cfg.filename = strcat (patient_id, '.nii');
cfg.parameter = 'anatomy';
cfg.format = 'nifti';
ft_volumewrite(cfg, source_interp);

cfg = [];
cfg.filename = strcat(patient_id, '_kurtosis.nii');
cfg.parameter = 'kurtosis';
cfg.format = 'nifti';
cfg.datatype = 'float'; % integer datatypes will be scaled to the maximum, floating point datatypes not
ft_volumewrite(cfg, source_interp);

%% visualize results in Anywave 
load source.mat
dat = ft_fetch_data(data_resampled);
hdr = ft_fetch_header(data_resampled);

%
npeaks = 5;
for i = 1:npeaks
  dat(end+1,:) = source.avg.mom{peakindex(i),:}; % see comment below about scaling
  hdr.label{end+1}= ['S' num2str(i)];
  hdr.chantype{end+1} = 'Source';
  hdr.chanunit{end+1} = 'T' ; % see note below about scaling
end
hdr.nChans = hdr.nChans+npeaks;
ft_write_data(patient_id, dat, 'header', hdr, 'dataformat', 'ades');

%
fid = fopen(fullfile(patient_output_dir, [patient_id, '.mrk']), 'w');
fprintf(fid,'%s\r\n','// AnyWave Marker File ');
k = 1;
for i = 1:npeaks
  dat = source.avg.mom{peakindex(i),:};
  sd = std(dat);
  tr = zeros(size(dat));
  tr(dat>6*sd)=1;
  [tmp, peaksample] = findpeaks(tr, 'MinPeakDistance', 300); % peaks have to be separated by 300 sample points to be treated as separate
  for j = 1:length(peaksample)
    fprintf(fid, 'S%d_%02d\t', i, j); % marker name
    fprintf(fid, '%d\t', dat(peaksample(j))); % marker value
    fprintf(fid, '%d\t',source.time(peaksample(j)) ); % marker time
    fprintf(fid, '%d\t', 0); % marker duration
    fprintf(fid, 'S%d\r\n', i); % marker channel
    k = k + 1;
  end
end
fclose(fid);
