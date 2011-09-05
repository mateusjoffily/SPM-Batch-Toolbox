function jobs = batch_builder(data, param, bfile, dispOK)
%
% BATCH_BUILDER Builds SPM5 compatible jobs.
%
%   BATCH_BUILDER (DATA, PARAM, BFILE, DISPOK) load user defined input
%     structures DATA and PARAM, create SPM5 'jobs' structure  
%     and, eventualy, launch spm_jobman('run', jobs).
%
%     DATA   - (struct) see 'batch_data_example.m'.
%     PARAM  - (struct) see 'batch_param_example.m'.
%     BFILE  - (string) 'jobs' structure is saved in BFILE (optional).
%              If BFILE is not passed as input argument or is empty, the 
%              program will automaticaly create an unique file name like 
%              'bb_yyyy-mm-dd-time'.
%     DISPOK - (boolean) If matlab program is running in '-nodisplay' 
%              mode (see matlab(UNIX) documentation), 'dispOk' should  be 
%              set to false. So, dialog windows will not, eventualy, block 
%              the program execution. This is always true for remote nodes 
%              submission.
%
%-----------------------------------------------------------------------
% Mateus Joffily - CNC/ISC - Jan/06 - v0.04
%

% Copyright (C) 2006, 2010 Mateus Joffily, mateusjoffily@gmail.com.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% Load spm defaults
spm('defaults','fmri');

if nargin < 2
    disp('Usage: jobs = BATCH_BUILDER(data, param, bfile, dispOK)');
    return;
end

if nargin < 3 ||  isempty(bfile)
    % Create unique filename
    bfile = ['bb_' datestr(now,30)];  
end

if nargin < 4
    % Set dispOK true
    dispOK = true;   
end

%-----------------------------------------------------------------------
% Reset virtual filenames
%-----------------------------------------------------------------------
spm_select('clearvfiles')

%-----------------------------------------------------------------------
% Initialize variables
%-----------------------------------------------------------------------
jobs = {};       % jobs struct
njc =  0;        % jobs counter

%-----------------------------------------------------------------------
% Backup files prefix
%-----------------------------------------------------------------------
bak.sprefix = param.sprefix;         % 3d file prefix
bak.fprefix = param.fprefix;         % functional files prefix
bak.mprefix = param.mprefix;         % mean file prefix

%-----------------------------------------------------------------------
% Pre-processing
%-----------------------------------------------------------------------
% Confirm file removal
if dispOK && isfield(param, 'preproc') && ...
   ~isempty(param.preproc.order) && param.preproc.cleandir
    ButtonName=questdlg('Old pre-processed images will be removed!', ...
                        'Warning', 'Continue', 'Cancel', 'Cancel');
    if strcmp(ButtonName, 'Cancel')
        return;
    end
end

% Loop over subjects/sessions
if isfield(param, 'preproc') && ~isempty(param.preproc.order)
    for Ns=1:length(data.preproc)
        
        if isempty(data.preproc(Ns).fdir) 
            continue;
        end
        
        disp(sprintf('%d - pre-processing: %s', Ns, ...
            data.preproc(Ns).fdir{1}));
        for k = 2:length(data.preproc(Ns).fdir)
            disp(sprintf('                    %s', ...
                data.preproc(Ns).fdir{k}));
        end
        
        % Clean directories, if requested
        %---------------------------------------------------------------
        if param.preproc.cleandir
            % Delete old pre processed functional data: every filename that
            % doesn't begin with 'f'
            disp('cleanning functional directory');
            for fdir = data.preproc(Ns).fdir
                % clean image files
                [ffiles fdirs] = select_files(char(fdir), '^[^f]', false);
                if ~isempty(fdirs)
                    ffiles = strcat(fdirs, ffiles);
                    for i = 1:size(ffiles,1)
                        delete(deblank(ffiles(i,:)));
                    end
                end
                % clean mat files
                delete([char(fdir) filesep '*.mat']);
            end
        end
        
        % Reset files prefix
        %---------------------------------------------------------------
        param.sprefix = bak.sprefix;         % 3d file prefix
        param.fprefix = bak.fprefix;         % functional files prefix
        param.mprefix = bak.mprefix;         % mean file prefix

        % Loop over preprocessing steps
        %---------------------------------------------------------------
        for i = 1:length(param.preproc.order)
            njc = njc + 1;
            switch param.preproc.order{i}
                case 'slice_timing'
                    [temporal param.fprefix] = ...
                        slice_timing(data.preproc(Ns), param);
                    jobs{njc}.temporal{1} = temporal;
                case 'realign_estimate'
                    spatial = ...
                        realign_estimate(data.preproc(Ns), param);
                    jobs{njc}.spatial{1} = spatial;
                case 'realign_estimate&reslice'
                    [spatial param.fprefix] = ...
                        realign_estimate_reslice(data.preproc(Ns), param);
                    jobs{njc}.spatial{1} = spatial;
                case 'realign_unwarp'
                    [spatial param.fprefix] = ...
                        realign_unwarp(data.preproc(Ns), param);
                    jobs{njc}.spatial{a} = spatial;
                case 'coregister_estimate'
                    spatial = ...
                        coregister_estimate(data.preproc(Ns), param);
                    jobs{njc}.spatial{1} = spatial;
                case 'coregister_reslice'
                    [spatial param.fprefix param.sprefix param.mprefix] = ...
                        coregister_reslice(data.preproc(Ns), param);
                    jobs{njc}.spatial{1} = spatial;
                case 'coregister_estimate&reslice'
                    [spatial param.fprefix param.sprefix param.mprefix] = ...
                        coregister_estimate_reslice(data.preproc(Ns), param);
                    jobs{njc}.spatial{1} = spatial;
                case 'segment'
                    [spatial param.sprefix param.mprefix] = ...
                        segment(data.preproc(Ns), param);
                    jobs{njc}.spatial{1} = spatial;
                case 'normalise_write'
                    [spatial param.fprefix param.sprefix param.mprefix] =  ...
                        normalise_write(data.preproc(Ns), param);
                    jobs{njc}.spatial{1} = spatial;
                case 'normalise_estimate&write'
                    [spatial param.fprefix param.sprefix param.mprefix] =  ...
                        normalise_estimate_write(data.preproc(Ns), param);
                    jobs{njc}.spatial{1} = spatial;
                case 'smooth'
                    [spatial param.fprefix param.mprefix] = ...
                        smooth(data.preproc(Ns), param);
                    jobs{njc}.spatial{1} = spatial;
                otherwise
                    disp(['Unknown pre-processing procedure: ' param.preproc.order{i}]);
                    njc = njc - 1;
                    break
            end
        end
    end
end

%-----------------------------------------------------------------------
% Statistical model
%-----------------------------------------------------------------------
% Confirm file removal
if dispOK && isfield(param, 'model') && ...
   ~isempty(param.model.order) && param.model.cleandir
    ButtonName=questdlg('Old model files will be removed!', ...
                        'Warning', 'Continue', 'Cancel', 'Cancel');
    if strcmp(ButtonName, 'Cancel')
        return;
    end
end

% Loop over subjects/sessions
if isfield(param, 'model') && ~isempty(param.model.order)
    for Ns = 1:length(data.model)
        
        if isempty(data.model(Ns).fdir) || ...
           isempty(data.model(Ns).log) || ...
           isempty(data.model(Ns).mdir) 
            continue;
        end
        
        disp(sprintf('%d - model: %s', Ns, char(data.model(Ns).mdir)));
        
        % Force full directory path
        mpath  = spm_select('CPath', data.model(Ns).mdir);
        % Add file separator at the end of directory path
        mpath  = deblank(fullfile(mpath, ' '));
        
        % Check if Model directory already exists, otherwise create it
        %----------------------------------------------------------------
        if ~exist(mpath, 'dir')
            disp('making model directory');
            [dir_root, dir_model] = fileparts(mpath(1:end-1));
            mkdir(dir_root, dir_model);
        end
        
        % Clean directories, if requested
        %-----------------------------------------------------------------
        if param.model.cleandir
            disp('cleanning model directory');
            % Delete old model data
            delete([mpath '*.*']);
        end
        
        % Loop over model steps
        %-----------------------------------------------------------------
        for i = 1:length(param.model.order)
            njc = njc + 1;
            switch param.model.order{i}
                case 'specification'
                    stats = specification(data.model(Ns), param);
                    jobs{njc}.stats{1} = stats;
                case 'specification_factorial'
                    stats = specification_factorial(data.model(Ns), param);
                    jobs{njc}.stats{1} = stats;
                case 'estimation_classical'
                    stats = estimation_classical(data.model(Ns));
                    jobs{njc}.stats{1} = stats;
                case 'estimation_bayesian_1st_level'
                    stats = estimation_bayesian_1st_level(data.model(Ns), param);
                    jobs{njc}.stats{1} = stats;
                case 'contrast'
                    stats = contrast(data.model(Ns), param);
                    jobs{njc}.stats{1} = stats;
                otherwise
                    disp(['Unknown model procedure: ' param.model.order{i}]);
                    njc = njc - 1;
                    break
            end
        end
    end
end

%--------------------------------------------------------------------------
% Util
%--------------------------------------------------------------------------
% Loop over subjects/sessions
if isfield(param, 'util') && isfield(data, 'util') && ...
   ~isempty(param.util.order)
    for Ns = 1:length(data.util)
        
        if isempty(data.util(Ns).idir) 
            continue;
        end
        
        disp(sprintf('%d - util: %s', Ns, ...
            data.util(Ns).idir{1}));
        for k = 2:length(data.util(Ns).idir)
            disp(sprintf('                    %s', ...
                data.util(Ns).idir{k}));
        end
        
        % Loop over preprocessing steps
        %------------------------------------------------------------------
        for i = 1:length(param.util.order)
            njc = njc + 1;
            switch param.util.order{i}
                case 'imcalc'
                    util = imcalc(data.util(Ns), param);
                    jobs{njc}.util{1} = util;
                otherwise
                    disp(['Unknown util procedure: ' param.util.order{i}]);
                    njc = njc - 1;
                    break
            end
        end
    end
end

%--------------------------------------------------------------------------
% Clear virtual filenames
%--------------------------------------------------------------------------
spm_select('clearvfiles')

%--------------------------------------------------------------------------
% Save spm job
%--------------------------------------------------------------------------
save(bfile, 'jobs');                      % save jobs
disp(['jobs saved in ' bfile '.']);

disp('done.');
%--------------------------------------------------------------------------
% End of main function
%=========================================================================

end

%==========================================================================
% Temporal Pre-Processing: slice_timing
%
function [temporal, fprefix] = slice_timing(data, param)

disp('temporal pre-processing: slice_timing');

% Rename some variables
%--------------------------------------------------------------------------
fpath       = data.fdir;
fprefix     = param.fprefix;
nslices     = param.nslices;
TR          = param.TR;
if isfield(param, 'TA') && ~isempty(param.TA)
    TA      = param.TA;
else
    TA      = TR-(TR/nslices);
end
slice_order = param.slice_order;
fext        = param.ext;
refslice    = param.refslice;

% Initialize variables
%--------------------------------------------------------------------------
temporal = [];    % temporal structure
vfiles   = {};    % virtual files list

for i=1:length(fpath)
    [ffiles fdirs] = select_files(fpath{i}, ['^' fprefix '.*\' fext '$']);    
    temporal.st.scans{i} = cellstr(strcat(fdirs, ffiles, ',1'));
    
    % update virtual files list
    vfiles = cellstr(strvcat(char(vfiles), strcat(fdirs, 'a', ffiles)));
end

temporal.st.nslices = nslices;
temporal.st.tr      = TR;
temporal.st.ta      = TA;
temporal.st.so      = calc_slice_order(nslices, slice_order);
if isempty(refslice)
    temporal.st.refslice = temporal.st.so(fix(nslices/2));
else
    temporal.st.refslice = refslice;
end
    
% update functional files prefix
fprefix = sprintf('a%s', fprefix);

% add virtual files
%--------------------------------------------------------------------------
spm_select('addvfiles', vfiles);

end

%==========================================================================
% Spatial Pre-Processing: realign_estimate (estimate only)
%
% Note: In the coregistration step, the sessions are first realigned to
%       each other, by aligning the first scan from each session to the 
%       first scan of the first session.  Then the images within each session 
%       are aligned to the first image of the session. 
function spatial = realign_estimate(data, param)

disp('spatial pre-processing: realign_estimate (estimate only)');

global defaults;

% Rename some variables
%--------------------------------------------------------------------------
fpath   = data.fdir;
fprefix = param.fprefix;
fext    = param.ext;

% Initialize variables
%--------------------------------------------------------------------------
spatial = [];     % spatial structure
vfiles  = {};     % virtual files list

for i=1:length(fpath)
    [ffiles fdirs] = select_files(fpath{i}, ['^' fprefix '.*\' fext '$']);    
    spatial.realign{1}.estimate.data{i} = cellstr(strcat(fdirs, ffiles, ',1'));

    % realigment parameters
    vfiles = cellstr(strvcat(char(vfiles), [fdirs(1,:) 'rp_' ffiles(1,1:end-4) '.txt']));
end

spatial.realign{1}.estimate.eoptions.quality = defaults.realign.estimate.quality;
spatial.realign{1}.estimate.eoptions.sep     = defaults.realign.estimate.sep;
spatial.realign{1}.estimate.eoptions.fwhm    = defaults.realign.estimate.fwhm;
spatial.realign{1}.estimate.eoptions.rtm     = defaults.realign.estimate.rtm;
spatial.realign{1}.estimate.eoptions.interp  = defaults.realign.estimate.interp;
spatial.realign{1}.estimate.eoptions.wrap    = defaults.realign.estimate.wrap;
if ~iscell(defaults.realign.estimate.weight)
    spatial.realign{1}.estimate.eoptions.weight = {};
else
    spatial.realign{1}.estimate.eoptions.weight = defaults.realign.estimate.weight;
end

% add virtual files
%--------------------------------------------------------------------------
spm_select('addvfiles', vfiles);

end

%==========================================================================
% Spatial Pre-Processing: realign_estimate_reslice (estimate & reslice)
%
% Note: In the coregistration step, the sessions are first realigned to
%       each other, by aligning the first scan from each session to the 
%       first scan of the first session.  Then the images within each session 
%       are aligned to the first image of the session. 
function [spatial, fprefix] = realign_estimate_reslice(data, param)

disp('spatial pre-processing: realign_estimate_reslice (estimate & reslice)');

global defaults;

% Rename some variables
%--------------------------------------------------------------------------
fpath   = data.fdir;
fprefix = param.fprefix;
fext    = param.ext;
rwhich  = param.preproc.realign.which;    % [2 1] = Write All Images + Mean Image
                                           % [0 1] = Write Mean Image only

% Initialize variables
%--------------------------------------------------------------------------
spatial = [];    % spatial structure
vfiles  = {};    % virtual files list

for i=1:length(fpath)
    [ffiles fdirs] = select_files(fpath{i}, ['^' fprefix '.*\' fext '$']);    
    spatial.realign{1}.estwrite.data{i} = cellstr(strcat(fdirs, ffiles, ',1'));
    
    % update virtual files list
    if rwhich(1)~=0
        vfiles = cellstr(strvcat(char(vfiles), strcat(fdirs, 'r', ffiles)));
    end
    if i==1 && rwhich(2)~=0
        vfiles = cellstr(strvcat(char(vfiles), [fdirs(1,:) 'mean' ffiles(1,:)]));   
    end
    
    % realigment parameters
    vfiles = cellstr(strvcat(char(vfiles), [fdirs(1,:) 'rp_' ffiles(1,1:end-4) '.txt']));
end

if rwhich(1)~=0
    fprefix=sprintf('r%s', fprefix);    
end

spatial.realign{1}.estwrite.eoptions.quality = defaults.realign.estimate.quality;
spatial.realign{1}.estwrite.eoptions.sep     = defaults.realign.estimate.sep;
spatial.realign{1}.estwrite.eoptions.fwhm    = defaults.realign.estimate.fwhm;
spatial.realign{1}.estwrite.eoptions.rtm     = defaults.realign.estimate.rtm;
spatial.realign{1}.estwrite.eoptions.interp  = defaults.realign.estimate.interp;
spatial.realign{1}.estwrite.eoptions.wrap    = defaults.realign.estimate.wrap;
if ~iscell(defaults.realign.estimate.weight)
    spatial.realign{1}.estwrite.eoptions.weight = {};
else
    spatial.realign{1}.estwrite.eoptions.weight = defaults.realign.estimate.weight;
end

spatial.realign{1}.estwrite.roptions.which  = rwhich;
spatial.realign{1}.estwrite.roptions.interp = defaults.realign.write.interp;
spatial.realign{1}.estwrite.roptions.wrap   = defaults.realign.write.wrap;
spatial.realign{1}.estwrite.roptions.mask   = defaults.realign.write.mask;

% add virtual files
%--------------------------------------------------------------------------
spm_select('addvfiles', vfiles);

end

%==========================================================================
% Spatial Pre-Processing: realign_unwarp
%
% Note: In the coregistration step, the sessions are first realigned to
%       each other, by aligning the first scan from each session to the 
%       first scan of the first session.  Then the images within each session 
%       are aligned to the first image of the session. 
function [spatial, fprefix] = realign_unwarp(data, param)

disp('spatial pre-processing: realign_unwarp');

global defaults;

% Rename some variables
%--------------------------------------------------------------------------
fpath   = data.fdir;
fprefix = param.fprefix;
fext    = param.ext;
rwhich  = param.preproc.realignunwarp.which;    % [2 1] = Write All Images + Mean Image
                                                % [2 0] = Write All Images
pmscan  = param.preproc.realignunwarp.pmscan;   % true  - use precalculated phase map (vdm*)
                                                % false - no phase correction.


% Initialize variables
%--------------------------------------------------------------------------
spatial = [];    % spatial structure
vfiles  = {};    % virtual files list

for i=1:length(fpath)
    [ffiles fdirs] = select_files(fpath{i}, ['^' fprefix '.*\' fext '$']);    
    spatial.realignunwarp.data(i).scans = cellstr(strcat(fdirs, ffiles, ',1'));
    if pmscan
        [pmfile pmdir] = select_files(fpath{i}, ['^vdm5_.*\' fext '$']);
        if ~isempty(pmdir)
            spatial.realignunwarp.data(i).pmscan = ...
                           cellstr(strcat(pmdir(1,:), pmfile(1,:), ',1'));
        else
            spatial.realignunwarp.data(i).pmscan = {};
        end
    else
        spatial.realignunwarp.data(i).pmscan = {};
    end
    
    % update virtual files list
    if rwhich(1)~=0
        vfiles = cellstr(strvcat(char(vfiles), strcat(fdirs, 'u', ffiles)));
    end
    if i==1 && rwhich(2)~=0
        vfiles = cellstr(strvcat(char(vfiles), [fdirs(1,:) 'meanu' ffiles(1,:)]));    
    end
    
    % realigment parameters
    vfiles = cellstr(strvcat(char(vfiles), [fdirs(1,:) 'rp_' ffiles(1,1:end-4) '.txt']));
end

if rwhich(1)~=0
    fprefix=sprintf('u%s', fprefix);
end

spatial.realignunwarp.eoptions.quality = defaults.realign.estimate.quality;
spatial.realignunwarp.eoptions.sep     = defaults.realign.estimate.sep;
spatial.realignunwarp.eoptions.fwhm    = defaults.realign.estimate.fwhm;
spatial.realignunwarp.eoptions.rtm     = defaults.realign.estimate.rtm;
spatial.realignunwarp.eoptions.einterp = defaults.realign.estimate.interp;
spatial.realignunwarp.eoptions.ewrap   = defaults.realign.estimate.wrap;
if ~iscell(defaults.realign.estimate.weight)
    spatial.realignunwarp.eoptions.weight = {};
else
    spatial.realignunwarp.eoptions.weight = defaults.realign.estimate.weight;
end

spatial.realignunwarp.uweoptions.basfcn   = defaults.unwarp.estimate.basfcn;
spatial.realignunwarp.uweoptions.regorder = defaults.unwarp.estimate.regorder;
spatial.realignunwarp.uweoptions.lambda   = defaults.unwarp.estimate.regwgt;
spatial.realignunwarp.uweoptions.jm       = defaults.unwarp.estimate.jm;
spatial.realignunwarp.uweoptions.fot      = [4 5];
spatial.realignunwarp.uweoptions.sot      = [];
spatial.realignunwarp.uweoptions.uwfwhm   = defaults.unwarp.estimate.fwhm;
spatial.realignunwarp.uweoptions.rem      = defaults.unwarp.estimate.rem;
spatial.realignunwarp.uweoptions.noi      = defaults.unwarp.estimate.noi;
spatial.realignunwarp.uweoptions.expround = defaults.unwarp.estimate.expround;

spatial.realignunwarp.uwroptions.uwwhich = rwhich;
spatial.realignunwarp.uwroptions.rinterp = defaults.realign.write.interp;
spatial.realignunwarp.uwroptions.wrap    = defaults.realign.write.wrap;
spatial.realignunwarp.uwroptions.mask    = defaults.realign.write.mask;

% add virtual files
spm_select('addvfiles', vfiles);

end

%==========================================================================
% Spatial Pre-Processing: coregister_estimate (estimate only)
%
function spatial = coregister_estimate(data, param)

disp('spatial pre-processing: coregister_estimate (estimate only)');

global defaults;

% Rename some variables
%--------------------------------------------------------------------------
fpath      = data.fdir;
spath      = char(data.sdir);
fprefix    = param.fprefix;
sprefix    = param.sprefix;
mprefix    = param.mprefix;
fext       = param.ext;
source_ref = param.preproc.coregister.source_ref;

% Initialize variables
%--------------------------------------------------------------------------
spatial = [];    % spatial structure

[mfile mdir] = select_files(fpath{1}, ['^' mprefix 'mean.*\' fext '$']);
[sfile sdir] = select_files(spath, ['^' sprefix '.*\' fext '$']);

if strcmp(source_ref, '3d_mean')
    spatial.coreg{1}.estimate.ref = {[mdir(1,:) mfile(1,:) ',1']};
    spatial.coreg{1}.estimate.source = {[sdir(1,:) sfile(1,:) ',1']};
    spatial.coreg{1}.estimate.other{1} = '';
    
elseif strcmp(source_ref, 'mean_3d')
    spatial.coreg{1}.estimate.ref = {[sdir(1,:) sfile(1,:) ',1']};
    spatial.coreg{1}.estimate.source = {[mdir(1,:) mfile(1,:) ',1']};
        
    Nimg=0;
    for i=1:length(fpath)
        [ffiles fdirs] = select_files(fpath{i}, ['^' fprefix '.*\' fext '$']);    
        for j=1:size(ffiles,1)
            Nimg = Nimg + 1;
            spatial.coreg{1}.estimate.other{Nimg,1}= [fdirs(j,:) ffiles(j,:) ',1'];
        end
    end
else
    disp('coregister_estimate: unknown source reference matching operation');
end

spatial.coreg{1}.estimate.eoptions.cost_fun = defaults.coreg.estimate.cost_fun;
spatial.coreg{1}.estimate.eoptions.sep      = defaults.coreg.estimate.sep;
spatial.coreg{1}.estimate.eoptions.tol      = defaults.coreg.estimate.tol;
spatial.coreg{1}.estimate.eoptions.fwhm     = defaults.coreg.estimate.fwhm;

end

%==========================================================================
% Spatial Pre-Processing: coregister_estimate_reslice (estimate & reslice)
%
function [spatial, fprefix, sprefix, mprefix] = coregister_estimate_reslice(data, param)

disp('spatial pre-processing: coregister_estimate_reslice (estimate & reslice)');

global defaults;

% Rename some variables
%--------------------------------------------------------------------------
fpath      = data.fdir;
spath      = char(data.sdir);
fprefix    = param.fprefix;
sprefix    = param.sprefix;
mprefix    = param.mprefix;
fext       = param.ext;
source_ref = param.preproc.coregister.source_ref;
cinterp    = param.preproc.coregister.interp;

% Initialize variables
%--------------------------------------------------------------------------
spatial = [];   % spatial structure
vfiles  = {};   % virtual files list

[mfile mdir] = select_files(fpath{1}, ['^' mprefix 'mean.*\' fext '$']);
[sfile sdir] = select_files(spath, ['^' sprefix '.*\' fext '$']);

if strcmp(source_ref, '3d_mean')
    spatial.coreg{1}.estwrite.ref = {[mdir(1,:) mfile(1,:) ',1']};
    spatial.coreg{1}.estwrite.source = {[sdir(1,:) sfile(1,:) ',1']};
    spatial.coreg{1}.estwrite.other{1} = '';
    
    sprefix=sprintf('r%s', sprefix);
    
    % update virtual files list
    vfiles = cellstr(strvcat(char(vfiles), [sdir(1,:) 'r' sfile(1,:)]));
    
elseif strcmp(source_ref, 'mean_3d')
    spatial.coreg{1}.estwrite.ref = {[sdir(1,:) sfile(1,:) ',1']};
    spatial.coreg{1}.estwrite.source = {[mdir(1,:) mfile(1,:) ',1']};
    
    % update mean file prefix
    mprefix=sprintf('r%s', mprefix);
    % update virtual files list
    vfiles = cellstr(strvcat(char(vfiles), [mdir(1,:) 'r' mfile(1,:)]));
    
    Nimg=0;
    for i=1:length(fpath)
        [ffiles fdirs] = select_files(fpath{i}, ['^' fprefix '.*\' fext '$']);    
        for j=1:size(ffiles,1)
            Nimg = Nimg + 1;
            spatial.coreg{1}.estwrite.other{Nimg,1}= [fdirs(j,:) ffiles(j,:) ',1'];
        end
        % update virtual files list
        vfiles = cellstr(strvcat(char(vfiles), strcat(fdirs, 'r', ffiles)));
    end
    % update functional files prefix
    fprefix=sprintf('r%s', fprefix);
    
else
    disp('coregister_estimate_reslice: unknown source reference matching operation');
end

spatial.coreg{1}.estwrite.eoptions.cost_fun = defaults.coreg.estimate.cost_fun;
spatial.coreg{1}.estwrite.eoptions.sep      = defaults.coreg.estimate.sep;
spatial.coreg{1}.estwrite.eoptions.tol      = defaults.coreg.estimate.tol;
spatial.coreg{1}.estwrite.eoptions.fwhm     = defaults.coreg.estimate.fwhm;

spatial.coreg{1}.estwrite.roptions.interp   = cinterp;
spatial.coreg{1}.estwrite.roptions.wrap     = defaults.coreg.write.wrap;
spatial.coreg{1}.estwrite.roptions.mask     = defaults.coreg.write.mask;

% add virtual files
%--------------------------------------------------------------------------
spm_select('addvfiles', vfiles);

end

%==========================================================================
% Spatial Pre-Processing: coregister_reslice (reslice only)
%
function [spatial, fprefix, sprefix, mprefix] = coregister_reslice(data, param)

disp('spatial pre-processing: coregister_reslice (reslice only)');

global defaults;

% Rename some variables
%--------------------------------------------------------------------------
fpath      = data.fdir;
spath      = char(data.sdir);
fprefix    = param.fprefix;
sprefix    = param.sprefix;
mprefix    = param.mprefix;
fext       = param.ext;
ref        = param.preproc.coregister.ref;
other      = param.preproc.coregister.other;
source     = param.preproc.coregister.source;
cinterp    = param.preproc.coregister.interp;

% Initialize variables
%--------------------------------------------------------------------------
spatial = [];   % spatial structure
vfiles  = {};   % virtual files list

[mfile mdir] = select_files(fpath{1}, ['^' mprefix 'mean.*\' fext '$']);
[sfile sdir] = select_files(spath, ['^' sprefix '.*\' fext '$']);

switch ref
    case '3d'
        spatial.coreg{1}.write.ref = {[sdir(1,:) sfile(1,:) ',1']};        
    case 'mean'
        spatial.coreg{1}.write.ref = {[mdir(1,:) mfile(1,:) ',1']};
    case 'other'
        other = spm_select('CPath',other);
        spatial.coreg{1}.write.ref = {[other ',1']};
end

switch source
    case '3d'
        spatial.coreg{1}.write.source = {[sdir(1,:) sfile(1,:) ',1']};
        
        sprefix = sprintf('r%s', sprefix);
    
        % update virtual files list
        vfiles = cellstr(strvcat(char(vfiles), [sdir(1,:) 'r' sfile(1,:)]));
    case 'func'
        if ~isempty(mfile)
            % If mean image exists, reslice it with other functional images
            spatial.coreg{1}.write.source = {[mdir(1,:) mfile(1,:) ',1']};

            % update mean file prefix
            mprefix=sprintf('r%s', mprefix);
            % update virtual files list
            vfiles = cellstr(strvcat(char(vfiles), [mdir(1,:) 'r' mfile(1,:)]));
            
            % Increment images counter
            Nimg = 1;
        else
            % Set images counter to zero
            Nimg = 0;
        end
        
        % loop over functional images
        for i=1:length(fpath)
            [ffiles fdirs] = select_files(fpath{i}, ['^' fprefix '.*\' fext '$']);
            for j=1:size(ffiles,1)
                Nimg = Nimg + 1;
                spatial.coreg{1}.write.source{Nimg,1}= [fdirs(j,:) ffiles(j,:) ',1'];
            end
            % update virtual files list
            vfiles = cellstr(strvcat(char(vfiles), strcat(fdirs, 'r', ffiles)));
        end
        % update functional files prefix
        fprefix=sprintf('r%s', fprefix);
end

spatial.coreg{1}.write.roptions.interp   = cinterp;
spatial.coreg{1}.write.roptions.wrap     = defaults.coreg.write.wrap;
spatial.coreg{1}.write.roptions.mask     = defaults.coreg.write.mask;

% add virtual files
%--------------------------------------------------------------------------
spm_select('addvfiles', vfiles);

end

%==========================================================================
% Spatial Pre-Processing: segment
%
function [spatial, sprefix, mprefix] = segment(data, param)

disp('spatial pre-processing: segment');

global defaults;

% Rename some variables
%--------------------------------------------------------------------------
fpath    = data.fdir;
spath    = char(data.sdir);
fprefix  = param.fprefix;
sprefix  = param.sprefix;
mprefix  = param.mprefix;
fext     = param.ext;
seg_data = param.preproc.segment.data;
if isfield(param.preproc.segment, 'mask')
    mask = param.preproc.segment.mask;
else
    mask = {''};
end

% Initialize variables
%--------------------------------------------------------------------------
spatial = [];   % spatial structure
vfiles  = {};   % virtual files list

% input images
switch(seg_data)
    case '3d'
        [sfile sdir] = select_files(spath, ['^' sprefix '.*\' fext '$']);
        spatial.preproc.data = {[sdir(1,:) sfile(1,:) ',1']};
        
        % update virtual files list
        vfiles = cellstr(strvcat(char(vfiles), [sdir(1,:) sfile(1,1:end-4) '_seg_sn.mat']));            
        vfiles = cellstr(strvcat(char(vfiles), [sdir(1,:) sfile(1,1:end-4) '_seg_inv_sn.mat']));            
    case 'mean'
        [mfile mdir] = select_files(fpath{1}, ['^' mprefix 'mean.*\' fext '$']);
        spatial.preproc.data  = {[mdir(1,:) mfile(1,:) ',1']};
        
        % update virtual files list
        vfiles = cellstr(strvcat(char(vfiles), [mdir(1,:) mfile(1,1:end-4) '_seg_sn.mat']));            
        vfiles = cellstr(strvcat(char(vfiles), [mdir(1,:) mfile(1,1:end-4) '_seg_inv_sn.mat']));   
end

% output images
spatial.preproc.output.GM      = param.preproc.segment.output.gm;
spatial.preproc.output.WM      = param.preproc.segment.output.wm;
spatial.preproc.output.CSF     = param.preproc.segment.output.csf;
spatial.preproc.output.biascor = param.preproc.segment.output.biascor;
spatial.preproc.output.cleanup = param.preproc.segment.output.cleanup;

% segment options
spatial.preproc.opts.tpm      = cellstr(defaults.preproc.tpm);
spatial.preproc.opts.ngaus    = defaults.preproc.ngaus;
spatial.preproc.opts.regtype  = defaults.preproc.regtype;
spatial.preproc.opts.warpreg  = defaults.preproc.warpreg;
spatial.preproc.opts.warpco   = defaults.preproc.warpco;
spatial.preproc.opts.biasreg  = defaults.preproc.biasreg;
spatial.preproc.opts.biasfwhm = defaults.preproc.biasfwhm;
spatial.preproc.opts.samp     = defaults.preproc.samp;
for i=1:length(mask)
    if ~isempty(mask{i})
        spatial.preproc.opts.msk{i} = sprintf('%s,1', mask{i});
    else
        spatial.preproc.opts.msk{i} = '';
    end
end

if param.preproc.segment.output.biascor == 1
    switch(seg_data)
        case '3d'
            % update 3d file prefix
            sprefix=sprintf('m%s', sprefix);
            % update virtual files list
            vfiles = cellstr(strvcat(char(vfiles), [sdir(1,:) 'm' sfile(1,:)]));
        case 'mean'
            % update virtual files list
            vfiles = cellstr(strvcat(char(vfiles), [mdir(1,:) 'm' mfile(1,:)]));            
            % update mean file prefix
            mprefix=sprintf('m%s', mprefix);
    end
end

% add virtual files
spm_select('addvfiles', vfiles);
    
end

%==========================================================================
% Spatial Pre-Processing: normalise_write (write only)
%
function [spatial, fprefix, sprefix, mprefix] = normalise_write(data, param)

disp('spatial pre-processing: normalise_write (write only)');

global defaults;

% Rename some variables
%--------------------------------------------------------------------------
fpath   = data.fdir;
spath   = char(data.sdir);
fprefix = param.fprefix;
sprefix = param.sprefix;
mprefix = param.mprefix;
fext    = param.ext;
if isfield(param.preproc.normalise, 'vox_3d')
    img3d_vox = param.preproc.normalise.vox_3d;
else
    % use default value
    img3d_vox = defaults.normalise.write.vox;
end
if isfield(param.preproc.normalise, 'vox_func')
    imgfunc_vox = param.preproc.normalise.vox_func;
else
    % use default value
    imgfunc_vox = defaults.normalise.write.vox;
end
norm_params = param.preproc.normalise.write.params;

% Initialize variables
%--------------------------------------------------------------------------
spatial = [];    % spatial structure
vfiles  = {};    % virtual files list

%Normalise: Write
%Note: we don't normalise 3D and functional images together because they might be written 
%      with different voxel sizes

% Select normalisation parameters file
%--------------------------------------------------------------------------
switch norm_params
    case '3d_sn'
        [sfile sdir] = select_files(spath, ...
            ['^.*[^_][^s][^e][^g]_sn\.mat$']);
    case '3d_seg_sn'
        [sfile sdir] = select_files(spath, ...
            ['^.*_seg_sn\.mat$']);
    case 'mean_sn'
        [sfile sdir] = select_files(fpath{1}, ...
            ['^.*mean.*[^_][^s][^e][^g]_sn\.mat$']);
    case 'mean_seg_sn'
        [sfile sdir] = select_files(fpath{1}, ...
            ['^.*mean.*_seg_sn\.mat$']);
end
norm_params_file = [sdir(1,:) sfile(1,:)];

% Normalise 3D image
%--------------------------------------------------------------------------
[sfile sdir] = select_files(spath, ['^' sprefix '.*\' fext '$']);
% update 3d file prefix
sprefix=sprintf('w%s', sprefix);
% update virtual files list
vfiles = cellstr(strvcat(char(vfiles), [sdir(1,:) 'w' sfile(1,:)]));

spatial.normalise{1}.write.subj.matname       = {norm_params_file};
spatial.normalise{1}.write.subj.resample{1,1} = [sdir(1,:) sfile(1,:) ',1'];
spatial.normalise{1}.write.roptions.preserve  = defaults.normalise.write.preserve;
spatial.normalise{1}.write.roptions.bb        = defaults.normalise.write.bb;
spatial.normalise{1}.write.roptions.interp    = defaults.normalise.write.interp;
spatial.normalise{1}.write.roptions.vox       = img3d_vox;
spatial.normalise{1}.write.roptions.wrap      = defaults.normalise.write.wrap;

% Normalise functional & mean images
%--------------------------------------------------------------------------
spatial.normalise{2}.write.subj.matname{1} = norm_params_file;

% check for mean image
[mfile mdir] = select_files(fpath{1}, ['^' mprefix 'mean.*\' fext '$']);
if ~isempty(mfile)
    spatial.normalise{2}.write.subj.resample{1,1} = [mdir(1,:) mfile(1,:) ',1'];
    % update mean files prefix
    mprefix=sprintf('w%s', mprefix);
    % update virtual files list
    vfiles = cellstr(strvcat(char(vfiles), [mdir(1,:) 'w' mfile(1,:)]));
    Nimg = 1;
else
    Nimg = 0;
end

% check for functional images
for i=1:length(fpath)
    [ffiles fdirs] = select_files(fpath{i}, ['^' fprefix '.*\' fext '$']);    
    for j=1:size(ffiles,1)
        Nimg = Nimg + 1;
        spatial.normalise{2}.write.subj.resample{Nimg,1} = [fdirs(j,:) ffiles(j,:) ',1'];
    end
    % update virtual files list
    vfiles = cellstr(strvcat(char(vfiles), strcat(fdirs, 'w', ffiles)));
end
% update functional files prefix
fprefix=sprintf('w%s', fprefix);

spatial.normalise{2}.write.roptions.preserve = defaults.normalise.write.preserve;
spatial.normalise{2}.write.roptions.bb       = defaults.normalise.write.bb;
spatial.normalise{2}.write.roptions.interp   = defaults.normalise.write.interp;
spatial.normalise{2}.write.roptions.vox      = imgfunc_vox;
spatial.normalise{2}.write.roptions.wrap     = defaults.normalise.write.wrap;

% add virtual files
%--------------------------------------------------------------------------
spm_select('addvfiles', vfiles);

end

%==========================================================================
% Spatial Pre-Processing: normalise_estimate_write (estimate & write)
%
% Note: First calculate normalization parameters from source image.
%       Then apply transformation into functional images.
function [spatial, fprefix, sprefix, mprefix] = normalise_estimate_write(data, param)

disp('spatial pre-processing: normalise_estimate_write (estimate & write)');

global defaults;

% Rename some variables
%--------------------------------------------------------------------------
fpath          = data.fdir;
spath          = data.sdir;
fprefix        = param.fprefix;
sprefix        = param.sprefix;
mprefix        = param.mprefix;
fext           = param.ext;
source_img     = param.preproc.normalise.estwrite.source;
template_files = param.preproc.normalise.estwrite.templates;
if isfield(param.preproc.normalise, 'vox_3d')
    img3d_vox = param.preproc.normalise.vox_3d;
else
    % use default value
    img3d_vox = defaults.normalise.write.vox;
end
if isfield(param.preproc.normalise, 'vox_func')
    imgfunc_vox = param.preproc.normalise.vox_func;
else
    % use default value
    imgfunc_vox = defaults.normalise.write.vox;
end

% Initialize variables
%--------------------------------------------------------------------------
spatial = [];    % spatial structure
vfiles  = {};    % virtual files list

[sfile sdir] = select_files(spath, ['^' sprefix '.*\' fext '$']);
[mfile mdir] = select_files(fpath{1}, ['^' mprefix 'mean.*\' fext '$'], false);

%Normalise: Estimate & Write
%Note: we don't resample 3D and functional images together becayse they might be written 
%      with different voxel sizes
if strcmp(source_img, '3d')
    spatial.normalise{1}.estwrite.subj.source = {[sdir(1,:) sfile(1,:) ',1']};
    spatial.normalise{1}.estwrite.subj.wtsrc = {};
    spatial.normalise{1}.estwrite.roptions.vox = img3d_vox;
    
    % resample 3d
    spatial.normalise{1}.estwrite.subj.resample = {[sdir(1,:) sfile(1,:) ',1']};
    
    % update virtual files
    vfiles = cellstr(strvcat(char(vfiles), [sdir(1,:) sfile(1,1:end-4) '_sn.mat']));            
    vfiles = cellstr(strvcat(char(vfiles), [sdir(1,:) sfile(1,1:end-4) '_inv_sn.mat']));            
elseif strcmp(source_img, 'mean')
    spatial.normalise{1}.estwrite.subj.source = {[mdir(1,:) mfile(1,:) ',1']};
    spatial.normalise{1}.estwrite.subj.wtsrc = {};
    spatial.normalise{1}.estwrite.roptions.vox = imgfunc_vox;
    
    % resample mean image
    spatial.normalise{1}.estwrite.subj.resample{1,1} = [mdir(1,:) mfile(1,:) ',1'];
    
    Nimg = 1;
    % resample functional images
    for i=1:length(fpath)
        [ffiles fdirs] = select_files(fpath{i}, ['^' fprefix '.*\' fext '$']);    
        for j=1:size(ffiles,1)
            Nimg = Nimg + 1;
            spatial.normalise{1}.estwrite.subj.resample{Nimg,1} = [fdirs(j,:) ffiles(j,:) ',1'];
        end
        % update virtual files list
        vfiles = cellstr(strvcat(char(vfiles), strcat(fdirs, 'w', ffiles)));
    end
    % update functional files prefix
    fprefix=sprintf('w%s', fprefix);
    
    % update virtual files
    vfiles = cellstr(strvcat(char(vfiles), [mdir(1,:) mfile(1,1:end-4) '_sn.mat']));            
    vfiles = cellstr(strvcat(char(vfiles), [mdir(1,:) mfile(1,1:end-4) '_inv_sn.mat']));            
else
    % not implemented yet
    return;
end

spatial.normalise{1}.estwrite.eoptions.template = template_files;

spatial.normalise{1}.estwrite.eoptions.weight{1} = defaults.normalise.estimate.weight;
spatial.normalise{1}.estwrite.eoptions.smosrc    = defaults.normalise.estimate.smosrc;
spatial.normalise{1}.estwrite.eoptions.smoref    = defaults.normalise.estimate.smoref;
spatial.normalise{1}.estwrite.eoptions.regtype   = defaults.normalise.estimate.regtype;
spatial.normalise{1}.estwrite.eoptions.cutoff    = defaults.normalise.estimate.cutoff;
spatial.normalise{1}.estwrite.eoptions.nits      = defaults.normalise.estimate.nits;
spatial.normalise{1}.estwrite.eoptions.reg       = defaults.normalise.estimate.reg;

spatial.normalise{1}.estwrite.roptions.preserve  = defaults.normalise.write.preserve;
spatial.normalise{1}.estwrite.roptions.bb        = defaults.normalise.write.bb;
spatial.normalise{1}.estwrite.roptions.interp    = defaults.normalise.write.interp;
spatial.normalise{1}.estwrite.roptions.wrap      = defaults.normalise.write.wrap;

% Resample other images using previoulsy estimated parameters: Write only
if strcmp(source_img, '3d')
    spatial.normalise{2}.write.subj.matname = {[sdir(1,:) sfile(1,1:end-4) '_sn.mat']};
    spatial.normalise{2}.write.roptions.vox = imgfunc_vox;
    
    % If mean image exists...
    if ~isempty(mfile)
        % resample mean image
        spatial.normalise{2}.write.subj.resample{1,1} = [mdir(1,:) mfile(1,:) ',1'];    
        Nimg = 1;
    else
        % Otherwise, don't do anything
        Nimg = 0;
    end
    
    % resample functional images
    for i=1:length(fpath)
        [ffiles fdirs] = select_files(fpath{i}, ['^' fprefix '.*\' fext '$']);    
        for j=1:size(ffiles,1)
            Nimg = Nimg + 1;
            spatial.normalise{2}.write.subj.resample{Nimg,1} = [fdirs(j,:) ffiles(j,:) ',1'];
        end
        % update virtual files list
        vfiles = cellstr(strvcat(char(vfiles), strcat(fdirs, 'w', ffiles)));
    end
    
elseif strcmp(source_img, 'mean')
    spatial.normalise{2}.write.subj.matname = {[mdir(1,:) mfile(1,1:end-4) '_sn.mat']};
    spatial.normalise{2}.write.roptions.vox = img3d_vox;    
   
    % resample 3d image
    spatial.normalise{2}.write.subj.resample = {[sdir(1,:) sfile(1,:) ',1']};

else
    % not implemented yet
    return;
end

spatial.normalise{2}.write.roptions.preserve = defaults.normalise.write.preserve;
spatial.normalise{2}.write.roptions.bb       = defaults.normalise.write.bb;
spatial.normalise{2}.write.roptions.interp   = defaults.normalise.write.interp;
spatial.normalise{2}.write.roptions.wrap     = defaults.normalise.write.wrap;

% update functional files prefix
fprefix=sprintf('w%s', fprefix);
% update 3d file prefix
sprefix=sprintf('w%s', sprefix);
% update virtual files list
vfiles = cellstr(strvcat(char(vfiles), [sdir(1,:) 'w' sfile(1,:)]));

% If mean image exists...
if ~isempty(mfile)
    % update mean file prefix
    mprefix=sprintf('w%s', mprefix);
    % Add mean files to the virtual list
    vfiles = cellstr(strvcat(char(vfiles), [mdir(1,:) 'w' mfile(1,:)]));
end

% add virtual files
%--------------------------------------------------------------------------
spm_select('addvfiles', vfiles);

end

%==========================================================================
% Spatial Pre-Processing: smooth
%
function [spatial, fprefix, mprefix] = smooth(data, param)

disp('spatial pre-processing: smooth');

% Rename some variables
%--------------------------------------------------------------------------
fpath   = data.fdir;
fprefix = param.fprefix;
mprefix = param.mprefix;
fext    = param.ext;
fwhm    = param.preproc.smooth.fwhm;

% Initialize variables
%--------------------------------------------------------------------------
spatial = [];    % spatial structure
vfiles  = {};    % virtual files list

[mfile mdir] = select_files(fpath{1}, ['^' mprefix 'mean.*\' fext '$'], false);

% If mean image exists...
if ~isempty(mfile)
    % smooth mean image
    spatial.smooth.data{1,1} = [mdir(1,:) mfile(1,:) ',1'];
    Nimg = 1;
    
    % update mean file prefix
    mprefix=sprintf('s%s', mprefix);
    % update virtual files list
    vfiles = cellstr(strvcat(char(vfiles), [mdir(1,:) 's' mfile(1,:)]));
else
    % Otherwise, don't do anything
    Nimg = 0;
end

% smooth functional images
for i=1:length(fpath)
    [ffiles fdirs] = select_files(fpath{i}, ['^' fprefix '.*\' fext '$']);    
    for j=1:size(ffiles,1)
        Nimg = Nimg + 1;
        spatial.smooth.data{Nimg,1} = [fdirs(j,:) ffiles(j,:) ',1'];
    end
    % update virtual files list
    vfiles = cellstr(strvcat(char(vfiles), strcat(fdirs, 's', ffiles)));
end

spatial.smooth.fwhm = fwhm;

% update functional files prefix
fprefix=sprintf('s%s', fprefix);

% add virtual files
%--------------------------------------------------------------------------
spm_select('addvfiles', vfiles);

end

%==========================================================================
% Model: specification
%
function stats = specification(data, param)

disp('model: specification');

% Rename some variables
%--------------------------------------------------------------------------
mpath       = spm_select('CPath',data.mdir);
fpath       = data.fdir;
log         = data.log;
fprefix     = param.fprefix;
TR          = param.TR;
nslices     = param.nslices;
slice_order = param.slice_order;
refslice    = param.refslice;
fext        = param.ext;

if isfield(param, 'TA') && ~isempty(param.TA)
    TA      = param.TA;
else
    TA      = TR-(TR/nslices);
end

if isempty(refslice)
    rslice = fix(nslices/2);
else
    so     = calc_slice_order(nslices, slice_order);
    rslice = find(so == refslice);
end

fmri_t = param.model.specification.T;

if isfield(param.model.specification, 'T0') && ...
                                    ~isempty(param.model.specification.T0)
    fmri_t0 = param.model.specification.T0;
else
    t = TA * (rslice-1) / (nslices-1);     % refslice acquisition time
    fmri_t0 = round(1 + (fmri_t * t / TR));
    % If TA = TR-(TR/nslices):
    % => fmri_t0 = round(1+(fmri_t*(rslice-1)/nslices));
end

units  = param.model.specification.units;
if ~ (strcmp(units, 'scans') || strcmp(units, 'secs') )
     disp(['error (model specification): unknown units for design: ' units]);
end

regress_estimated_movement = param.model.regress_estimated_movement;
if isfield(param.model.specification, 'mask')
    mask = param.model.specification.mask;
else
    mask = {''};
end

% 'hrf, 'fourier', 'hanning', 'gamma', 'fir'
bases.type = param.model.specification.bases.type;  
% [Time[0/1] Dispersion[0/1]], only 'hrf'
bases.hrf.derivatives = param.model.specification.bases.hrf.derivatives; 
% only 'fourier', 'hanning', 'gamma', 'fir'
bases.others.length = param.model.specification.bases.others.length; 
% only 'fourier', 'hanning', 'gamma', 'fir'
bases.others.order = param.model.specification.bases.others.order; 

hpf         = param.model.specification.hpf;
global_norm = param.model.specification.global;

% Initialize variables
%--------------------------------------------------------------------------
stats = [];    % stats structure

stats.fmri_spec.timing.units   = units;
stats.fmri_spec.timing.RT      = TR;
stats.fmri_spec.timing.fmri_t  = fmri_t;
stats.fmri_spec.timing.fmri_t0 = fmri_t0;


for i=1:length(fpath)
    [ffiles fdirs] = select_files(fpath{i}, ['^' fprefix '.*\' fext '$']);
    stats.fmri_spec.sess(i).scans = cellstr(strcat(fdirs, ffiles, ',1'));
    
    % Check if time modulation or parametric modulation is defined in logfile
    logparam=load(log{i});
    if isfield(logparam, 'tmod') || isfield(logparam, 'modulation')
        for j=1:numel(logparam.names)
            stats.fmri_spec.sess(i).cond(j).name     = logparam.names{j};
            stats.fmri_spec.sess(i).cond(j).onset    = logparam.onsets{j};
            stats.fmri_spec.sess(i).cond(j).duration = logparam.durations{j};
            if isfield(logparam, 'tmod')
                stats.fmri_spec.sess(i).cond(j).tmod = logparam.tmod(j);
            else
                stats.fmri_spec.sess(i).cond(j).tmod = 0;
            end
            stats.fmri_spec.sess(i).cond(j).mod = struct('name',{},'param',{},'poly',{});
            if isfield(logparam, 'modulation') && ~isempty(logparam.modulation(j).names)
                stats.fmri_spec.sess(i).cond(j).mod(1).name  = ...
                    logparam.modulation(j).names;
                stats.fmri_spec.sess(i).cond(j).mod(1).param = ...
                    reshape(logparam.modulation(j).param, numel(logparam.modulation(j).param), 1);
                stats.fmri_spec.sess(i).cond(j).mod(1).poly  = ...
                    logparam.modulation(j).poly;
            end
            stats.fmri_spec.sess(i).multi={''};
        end
    else   % Otherwise, use multiple coditions structure
        stats.fmri_spec.sess(i).cond = struct([]);
        stats.fmri_spec.sess(i).multi{1} = spm_select('CPath',log{i});
    end
    
    if regress_estimated_movement
        d = dir([fpath{i} 'rp_*.txt']);
        if length(d) == 1
            mov_regress_name = {'x' 'y' 'z' 'pitch' 'roll' 'yaw'};
            % [x y z pitch roll yaw] = textread([fpath{i} d(1).name], '%f %f %f %f %f %f');
            for c=1:6
                stats.fmri_spec.sess(i).regress(c).name = mov_regress_name{c};
                stats.fmri_spec.sess(i).regress(c).val  = eval(mov_regress_name{c});
            end
        else
            stats.fmri_spec.sess(i).regress = struct([]);
            disp(['Could not find rp_*.txt file, or more than one file was found, for session ' fpath{i}]);
        end
    else
        stats.fmri_spec.sess(i).regress = struct([]);
    end
    
    stats.fmri_spec.sess(i).multi_reg{1} = '';
    stats.fmri_spec.sess(i).hpf = hpf;
end

if isfield(data, 'fact') && ~isempty(data.fact)
    % Load factors from .mat file
    load(data.fact);
    stats.fmri_spec.fact = fact;
else
    stats.fmri_spec.fact = struct([]);
end

switch bases.type
    case 'hrf'
        stats.fmri_spec.bases.hrf.derivs         = bases.hrf.derivatives;
    case 'fourier'
        stats.fmri_spec.bases.fourier.length     = bases.others.length;
        stats.fmri_spec.bases.fourier.order      = bases.others.order;
    case 'hanning'
        stats.fmri_spec.bases.fourier_han.length = bases.others.length;
        stats.fmri_spec.bases.fourier_han.order  = bases.others.order;
    case 'gamma'
        stats.fmri_spec.bases.gamma.length       = bases.others.length;
        stats.fmri_spec.bases.gamma.order        = bases.others.order;
    case 'fir'
        stats.fmri_spec.bases.fir.length         = bases.others.length;
        stats.fmri_spec.bases.fir.order          = bases.others.order;       
end
    
stats.fmri_spec.volt = 1;
stats.fmri_spec.dir{1} = mpath;
stats.fmri_spec.global = global_norm;
for i=1:length(mask)
    if ~isempty(mask{i})
        stats.fmri_spec.mask{i} = sprintf('%s,1', mask{i});
    else
        stats.fmri_spec.mask{i} = '';
    end
end
stats.fmri_spec.cvi = 'AR(1)';

end

%==========================================================================
% Model: factorial design specification
%        Not operational yet, still needs to implement "des" (design) field
function stats = specification_factorial(data, param)

disp('model: factorial design specification');

% Rename some variables
%--------------------------------------------------------------------------
mpath       = spm_select('CPath',data.mdir);
fpath       = data.fdir;
fext        = param.ext;

if isfield(param.model.factorial_specification.mask, 'em')
    mask = param.model.factorial_specification.mask.em;
else
    mask = {''};
end

% Initialize variables
%--------------------------------------------------------------------------
stats = [];    % stats structure

% Covariates and nuissance variables
if isfield(data, 'cov') && ...
   ~isempty(data.cov)
    % Load covariates from .mat file
    load(data.cov);
    stats.factorial_design.cov = cov;
    return
end

% Threshold masking: 0='None, ]1:Inf]='Absolute' or [0:1]='Relative'
if param.model.factorial_specification.mask.tm == -1
    stats.factorial_design.masking.tm.tm_none = [];
elseif param.model.factorial_specification.mask.tm >= 0 && ...
       param.model.factorial_specification.mask.tm <= 1
    stats.factorial_design.globalc.tm.tmr.rthresh = ...
        param.model.factorial_specification.mask.tm;
else
    stats.factorial_design.globalc.tm.tma = ...
        param.model.factorial_specification.mask.tm;
end

% Implicit mask: 0=no / 1=yes
stats.factorial_design.masking.im = ...
    param.model.factorial_specification.mask.im;

% Explicit mask images file
for i=1:length(mask)
    if ~isempty(mask{i})
        stats.factorial_design.masking.em{i} = sprintf('%s,1', mask{i});
    else
        stats.factorial_design.masking.em{i} = '';
    end
end

% Global Calculation: -2='Omit', -1='Mean' or vector of global values='User'
if param.model.factorial_specification.globalc == -2
    stats.factorial_design.globalc.g_omit = [];
elseif param.model.factorial_specification.globalc == -1
    stats.factorial_design.globalc.g_mean = [];
else
    stats.factorial_design.globalc.g_user.global_uval = ...
        param.model.factorial_specification.globalc;
end

% Overall Grand Mean Scaling: 0='No' or grand mean scaled value='Yes'
if param.model.factorial_specification.globalm.gmsca <= 0
    stats.factorial_design.globalm.gmsca.gmsca_no = [];
else
    stats.factorial_design.globalm.gmsca.gmsca_yes.gmscv = ...
        param.model.factorial_specification.globalm.gmsca;
end

% Normalisation: 'None', 'Proportional' or 'ANCOVA'
stats.factorial_design.globalm.glonorm = ...
    param.model.factorial_specification.globalm.glonorm;

% SPM.mat dir
stats.factorial_design.dir{1} = mpath;

end

%==========================================================================
% Model: estimation (classical)
%
function stats = estimation_classical(data)

disp('model: estimation_classical');

% Rename some variables
%--------------------------------------------------------------------------
mpath = spm_select('CPath',data.mdir);

% Initialize variables
%--------------------------------------------------------------------------
stats = [];    % stats structure

stats.fmri_est.spmmat{1} = [mpath filesep 'SPM.mat'];
stats.fmri_est.method.Classical= 1;

end

%==========================================================================
% Model: estimation (bayesian 1st level)
%
function stats = estimation_bayesian_1st_level(data, param)

disp('model: estimation_bayesian_1st_level');

% Rename some variables
%--------------------------------------------------------------------------
mpath = spm_select('CPath',data.mdir);
log   = data.log;
regress_estimated_movement = param.model.regress_estimated_movement;

% Initialize variables
%--------------------------------------------------------------------------
stats = [];    % stats structure

stats.fmri_est.spmmat{1} = [mpath filesep 'SPM.mat'];
stats.fmri_est.method.Bayesian.space.Volume=1;
stats.fmri_est.method.Bayesian.signal='GMRF';
stats.fmri_est.method.Bayesian.ARP=3;
stats.fmri_est.method.Bayesian.noise.GMRF=1;
stats.fmri_est.method.Bayesian.anova.first='No';
stats.fmri_est.method.Bayesian.anova.second='Yes';

addmov=[];
if regress_estimated_movement
    % Add 6 additional columns after each session
    addmov=[0 0 0 0 0 0];
end

% Define automaticaly the contrasts
%-They will be automaticaly defined by ANOVA Second Level
names = {};      % define 'name', but it will be loaded from log file
load([log{1}]);
n=0;
Ncond=length(names);
for i=1:Ncond
    for j=1:Ncond
        if i~=j
            n=n+1;
            stats.fmri_est.method.Bayesian.gcon(n).name = sprintf('%s-%s', names{i}, names{j});
            v=[zeros(1,Ncond) addmov zeros(1,Ncond)];
            v([i j])=[1 -1];
            stats.fmri_est.method.Bayesian.gcon(n).convec = repmat(v,1,length(log))';
        else
            n=n+1;
            stats.fmri_est.method.Bayesian.gcon(n).name = sprintf('%s', names{i});
            v=[zeros(1,Ncond) addmov zeros(1,Ncond)];
            v(i)=1;
            stats.fmri_est.method.Bayesian.gcon(n).convec = repmat(v,1,length(log))';                                              
        end
    end
end

end

%==========================================================================
% Model: contrast
%
function stats = contrast(data, param)

disp('model: contrast');

% Rename some variables
mpath = spm_select('CPath',data.mdir);
log   = data.log;
if isfield(param.model, 'contrast') && ...
   isfield(param.model.contrast, 'delete')
    condel = param.model.contrast.delete;
else
    condel = 0;
end
if param.model.regress_estimated_movement
% If movement parameters are modelled, six additional
% regreesors have been included in the design matrix.
    regmov = 6;
else
    regmov = 0;
end

% Initialize variables
%--------------------------------------------------------------------------
stats = [];           % stats structure
Nlog  = numel(log);   % Number of logfiles

stats.con.spmmat{1} = [mpath filesep 'SPM.mat'];
stats.con.delete    = condel;

if isfield(data, 'con') && ...
   ~isempty(data.con)
    % Load contrasts from .mat file
    load(data.con);
    stats.con.consess = consess;
    return
end

% Get conditions for each session
sessrep = 'both';  % 'repl', 'sess' or 'both'
for k = 1:Nlog   % Loop over logfiles (i.e. sessions)
    load([log{k}]);
    
    [sess(k).names{1:numel(names)}] = deal(names{:});
    
    if k > 1
        % Check if conditions are identical across sessions
        if ~all(strcmp(sess(k-1).names, sess(k).names))
            sessrep = 'none';
        end
    end
end    

% Set contrasts
nc  = 0;    % contrast counter
ccs = 0;    % conditions cumulative counter
for k=1:numel(sess)  % Loop over sessions
    
    nm = numel(sess(k).names);   % number of conditions
    
    % F-contrast
    nc = nc + 1;                 % increment contrast counter
    convec = [zeros(nm,ccs) eye(nm)];
    if ~strcmp(sessrep, 'none')
        name = 'main effect';
    else
        name = sprintf('main effect - Session %d', k);
    end
    stats.con.consess{nc}.fcon = struct('name'    , name, ...
                                        'convec'  , {{convec}}, ...
                                        'sessrep' , sessrep);
    
    % T-contrast
    for i = 1:nm
        nc = nc + 1;                 % increment contrast counter
        convec = zeros(1,ccs+nm);
        convec(ccs+i) = 1;
        if ~strcmp(sessrep, 'none')
            name = sess(k).names{i};
        else
            name = sprintf('%s - Session %d',sess(k).names{i}, k);
        end
        stats.con.consess{nc}.tcon = struct('name'    , name, ...
                                            'convec'  , convec, ...
                                            'sessrep' , sessrep);
        for ii = 1:nm
            if i == ii, continue, end
            nc = nc + 1;                 % increment contrast counter
            convec = zeros(1,ccs+nm);
            convec([ccs+i ccs+ii]) = [1 -1];
            if ~strcmp(sessrep, 'none')
                name = sprintf('%s - %s', ...
                                  sess(k).names{i}, sess(k).names{ii});
            else
                name = sprintf('%s - %s - Session %d', ...
                                  sess(k).names{i}, sess(k).names{ii}, k);
            end
            stats.con.consess{nc}.tcon = struct('name'    , name, ...
                                                'convec'  , convec, ...
                                                'sessrep' , sessrep);
        end
    end

    % Increment conditions cumulative counter
    ccs = ccs + nm + regmov;
    
    % If conditions are identical across sessions, contrasts will be 
    % replicated over sessions. So, we stop after first iteration.
    if ~strcmp(sessrep, 'none')
        break;
    end
end

end

%==========================================================================
% Util: imcalc
%
function util = imcalc(data, param)

disp('util: imcalc');

% Rename some variables
%--------------------------------------------------------------------------
regimg     = param.util.images;
ipath      = data.idir;
expression = param.util.imcalc.expression;
dmtx       = param.util.imcalc.options.dmtx;
mask       = param.util.imcalc.options.mask;
uinterp    = param.util.imcalc.options.interp;
dtype      = param.util.imcalc.options.dtype;

% Initialize variables
%--------------------------------------------------------------------------
util = [];       % util structure
vfiles  = {};    % virtual files list

% get images
Nimg = 0;
for i=1:length(ipath)
    [ifiles idirs] = select_files(ipath{i}, regimg);    
    for j=1:size(ifiles,1)
        Nimg = Nimg + 1;
        util.imcalc.input{Nimg,1} = [idirs(j,:) ifiles(j,:) ',1'];
    end
end

% Output filename
[pout fout eout] = fileparts(util.imcalc.input{1});
util.imcalc.output          = [pout filesep 'imcalc_' fout eout(1:end-2)];

util.imcalc.expression      = expression;
util.imcalc.options.dmtx    = dmtx;
util.imcalc.options.mask    = mask;
util.imcalc.options.uinterp = uinterp;
util.imcalc.options.dtype   = dtype;

% update virtual files list
vfiles = cellstr(strvcat(char(vfiles), util.imcalc.output));

% add virtual files
%--------------------------------------------------------------------------
spm_select('addvfiles', vfiles);

end

%==========================================================================
% Automaticaly select files
%
function [myfiles, mydirs] = select_files(mydir, myfilt, warnON)

if nargin < 3
    warnON = true;
end

% Get 'mydir' full directory path
mydir   = spm_select('CPath',mydir);
% Add file separator at the end of directory path
mydir   = deblank(fullfile(mydir, ' '));
% Get 'myfiles'
myfiles = spm_select('List',mydir,myfilt);

if isempty(myfiles), 
    mydirs = [];
    if warnON
        disp(['warning: Files (''' myfilt ''') not found in ' mydir]);
    end
else
    % mydirs has the same length as myfiles
    mydirs  = repmat(mydir,size(myfiles,1),1);  
end

end

%==========================================================================
% Calculate slice order vector
%
function so = calc_slice_order(nslices, slice_order)

switch slice_order
    case 'ascending'
        so = 1:1:nslices;
    case 'descending'
        so = nslices:-1:1;
    case  'interleaved(middle-top)'
        for k = 1:nslices
            nord = round( (nslices-k)/2 + ( rem( (nslices-k),2 ) * (nslices-1) / 2 ) ) + 1;
            so(k) = nord;
        end
    case 'interleaved(bottom-up)'
        so = [1:2:nslices 2:2:nslices];
    case 'interleaved(top-down)'
        so = [nslices:-2:1, nslices-1:-2:1];        
    case  'interleaved(siemens)'
        if rem(nslices,2) == 1
            % Odd number of slices
            so = [1:2:nslices 2:2:nslices];
        else
            % Even number of slices
            so = [2:2:nslices 1:2:nslices];
        end
    otherwise
        disp('error (slice_timing): unknown slice order');
end

end