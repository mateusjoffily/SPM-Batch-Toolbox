%=======================================================================
% Define Pre-processing, Stats & Util params 
%--------------------------------------------------------------------------

% initialize param - do not change this line
param = struct('preproc', [], 'model', [], 'util', []);

%=======================================================================
% General variables
%-----------------------------------------------------------------------
% Volumes prefix and extension
%--------------------------------------------------------------------------
param.ext     = '.img';        % files extension ('.img' or '.nii')
param.sprefix = 's';           % 3d volume prefix: usualy 's' in SPM
param.fprefix = 'f';           % functional volume prefix: usualy 'f' in SPM
param.mprefix = '';            % mean image prefix: any letter that 
                               % comes before the word 'mean' in the
                               % file name

% Volumes acquisition parameters
%--------------------------------------------------------------------------
param.nslices     = 28;        % Number of slices
param.TR          = 2.65;      % TR in seconds
param.TA          = [];        % Leave it empty, if TA=TR-(TR/nslices)
param.slice_order = 'interleaved(siemens)';     % 'ascending',
                                                % 'descending',
                                                % 'interleaved(middle-top)',
                                                % 'interleaved(bottom-up)',
                                                % 'interleaved(top-down)'
                                                % 'interleaved(siemens)'
param.refslice = [];        % Leave this field empty if you want the program 
                            % to automaticaly select the slice at TR/2 as
                            % the reference slice, taking into account the
                            % slice order. This information will be used by
                            % the slice timing pre-processing and model
                            % specification routines.
                                                            

%=======================================================================
% Pre-processing specifics
%-----------------------------------------------------------------------
param.preproc.cleandir = false;      % boolean: true or false
                                     % CAUTION! It will delete every file 
                                     % that doesn't match 'f*' and '*.mat'
                                     % in the 'data.preproc' path.
                                                 
                                                  
% Realign
%--------------------------------------------------------------------------
% realign - reslice
param.preproc.realign.which = [0 1];               % [2 1] = Write All Images + Mean Image
                                                   % [0 1] = Write Mean Image only
                                                   
% Realign - unwarp                                                   
%--------------------------------------------------------------------------
param.preproc.realignunwarp.which = [2 1];         % [2 1] = Write All Images + Mean Image
                                                   % [2 0] = Write All Images
param.preproc.realignunwarp.pmscan = false;        % true  - use precalculated phase map (vdm*)
                                                   % false - no phase correction.
                                                   
                                                  
% Coregister
%--------------------------------------------------------------------------
% coregister - estimate&reslice and estimate
param.preproc.coregister.source_ref = '3d_mean';   % 'mean_3d': mean image (source) is matched to 3d (reference)
                                                   % '3d_mean': 3d (source) is matched to mean image (reference)
% coregister - reslice
param.preproc.coregister.ref = '3d';               % '3d'   : images are resliced to match 3d image
                                                   % 'mean' : images are resliced to match mean image
                                                   % 'other': images are resliced to match user specified image
param.preproc.coregister.other = '';               % 'other image', required if ref is 'other'
param.preproc.coregister.source = 'func';          % '3d'   : reslice 3d image
                                                   % 'func' : reslice functional images (All Images + Mean Image)

% coregister - estimate&reslice and reslice
param.preproc.coregister.interp = 2;               % interpolation method:
                                                   % 0 = Nearest Neighbour
                                                   % 1 = Trilinear
                                                   % 2..7 = Nth Degree B-Spline

% Segment
%--------------------------------------------------------------------------
param.preproc.segment.data='mean';               % '3d' or 'mean': image modality to be segmented
param.preproc.segment.output.gm  = [0 0 0];      % [Modulated Normalised Native]
param.preproc.segment.output.wm  = [0 0 0];      % [Modulated Normalised Native]
param.preproc.segment.output.csf = [0 0 0];      % [Modulated Normalised Native]
param.preproc.segment.output.biascor = 0;        % save bias corrected image: 0-no / 1-yes
param.preproc.segment.output.cleanup = 0;        % clean up any patitions: 0-no / 1-yes
param.preproc.segment.mask = {''};               % masking images

% Normalise
%--------------------------------------------------------------------------
% normalise - estimate & write
param.preproc.normalise.estwrite.source = '3d';      % '3d' or 'mean': source image to be matched to templates
param.preproc.normalise.estwrite.templates = {'/usr/local/lib/spm/spm5/templates/T1.nii,1'}; %template(s) to match source image
% normalise - write only
param.preproc.normalise.write.params = 'mean_seg_sn';       % use normalisation parameters estimated from: 
                                                        % '3d_sn'       - 3d image with normalise routine
                                                        % 'mean_sn'     - functional image with normalise routine
                                                        % '3d_seg_sn'   - 3d image with segment routine
                                                        % 'mean_seg_sn' - functional image with segment routine
% normalise - all (estimate&write and write)
param.preproc.normalise.vox_3d   = [1 1 1];        % Voxel size (mm) to write - 3d images
param.preproc.normalise.vox_func = [3 3 3];        % Voxel size (mm) to write - Functional images

% Smooth
%--------------------------------------------------------------------------
param.preproc.smooth.fwhm = [8.0 8.0 8.0];         % FWHM

% Pre-processing tasks (enter steps in the desired order):
%--------------------------------------------------------------------------
%   'slice_timing', 
%   'realign_estimate', 'realign_estimate&reslice', 'realign_unwarp', 
%   'coregister_estimate', 'coregister_reslice', 'coregister_estimate&reslice', 
%   'segment', 
%   'normalise_estimate&write', 'normalise_write', 
%   'smooth'
%
% N.B. Leave this field empty if you don't want to run pre-processing
% (i.e. param.preproc.order={};)
param.preproc.order={'slice_timing' 'realign_estimate&reslice' 'coregister_estimate' 'segment' 'normalise_write' 'smooth'};  


%=======================================================================
% Stats specifics
%-----------------------------------------------------------------------
param.model.cleandir = false;        % boolean: true or false
                                     % CAUTION! It will delete every file 
                                     % in the 'data.model' path (SPM.mat).
                                     
param.model.regress_estimated_movement=false; % include movement regressors
                                              % into the model? true/false

% Specification
%--------------------------------------------------------------------------                      
% model specification
param.model.specification.units = 'secs';                 % units for design: 'secs' or 'scans'
param.model.specification.T = param.nslices;              % Microtime resolution: usualy equal to nslices
param.model.specification.T0 = [];                        % Microtime onset: leave it empty for automaticaly setting at TR/2 
param.model.specification.hpf = 128;                      % High-pass filter (s)
param.model.specification.bases.type = 'hrf';             % 'hrf, 'fourier', 'hanning', 'gamma', 'fir'
param.model.specification.bases.hrf.derivatives = [0 0];      % [Time Dispersion], only 'hrf'
param.model.specification.bases.others.length = [];           % only 'fourier', 'hanning', 'gamma', 'fir'
param.model.specification.bases.others.order = [];            % only 'fourier', 'hanning', 'gamma', 'fir'
param.model.specification.global = 'None';                    % 'None' or 'Scaling'
param.model.specification.mask = {''};                        % explicit mask images

% % factorial design specification
% % !! NOT WORKING YET !! See batch_builder.m for more information.
% param.model.factorial_specification.des = 1;            % 1 = "One-sample t-test"
%                                                         % 2 = "Two-sample t-test"
%                                                         % 3 = "Paired t-test"
%                                                         % 4 = "Multiple regression"
%                                                         % 5 = "Full factorial"
%                                                         % 6 = "Flexible factorial"
% param.model.factorial_specification.mask.tm = -1;       % Threshold masking: -1 = 'None, ]1:Inf] = 'Absolute' or [0:1] = 'Relative'
% param.model.factorial_specification.mask.im = 1;        % Implicit mask: 0 = no / 1 = yes
% param.model.factorial_specification.mask.em = {''};     % Explicit mask images
% param.model.factorial_specification.globalc = -2;                % Global Calculation: -2 = 'Omit', -1 = 'Mean', or vector of global values = 'User'
% param.model.factorial_specification.globalm.gmsca = 0;           % Overall Grand Mean Scaling: 0 = 'No', or grand mean scaled value = 'Yes'
% param.model.factorial_specification.globalm.normalisation = 1;   % Normalisation: 1 = 'None', 2 = 'Proportional' or 3 = 'ANCOVA'

% Contrast
%--------------------------------------------------------------------------
param.model.contrast.delete  = 0;        % delete existing contrasts? 0 = no / 1 = yes


% Stats tasks (enter steps in the desired order):
%--------------------------------------------------------------------------
%        'specification'
%        'estimation_classical' 
%        'estimation_bayesian_1st_level' 
%        'contrast'
% N.B. Leave this field empty if you don't want to run statistics
% (i.e. param.model.order={};)
param.model.order={'specification' 'estimation_classical' 'contrast'};
        

%==========================================================================
% Util specifics
%--------------------------------------------------------------------------
% NOTE: At the the moment, only SPM5 'image calculator' Util task 
%       is emplemented.

param.util.images = '^f.*\.img$';

% Image Calculator                                                  
%--------------------------------------------------------------------------
param.util.imcalc.expression     = 'var(X)';
param.util.imcalc.options.dmtx   = 1;
param.util.imcalc.options.mask   = 0;
param.util.imcalc.options.interp = 1;
param.util.imcalc.options.dtype  = 4;

% Util tasks (enter steps in the desired order):
%--------------------------------------------------------------------------
%   'imcalc'
%
% N.B. Leave this field empty if you don'twant to run util
% (i.e. param.util.order={};)
param.util.order={'imcalc'};  

%-END======================================================================