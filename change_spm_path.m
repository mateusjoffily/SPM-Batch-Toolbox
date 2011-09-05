function change_spm_path(SPMmat, old_fdir, new_fdir)
%
% CHANGE_SPM_PATH Change SPM.mat images path.
%
%   SPMmat   - SPM.mat full path name
%   old_fdir - Functional images old directories (cell array)
%   new_fdir - Functional images new directories (cell array). Each old
%              directory in old_fdir must match a new directory in new_fdir.
%
%--------------------------------------------------------------------------
% Mateus Joffily - CNC/CNRS - July/2007

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

% Initialize SPM
SPM = [];

% Get SPM.mat full path name
SPMmat   = spm_select('CPath', SPMmat);

% Append .mat to SPMmat filename, if it doesn't exist
SPMmat = [spm_str_manip(SPMmat, 's') '.mat'];

% Get current SPM.mat directory
SPMnew = spm_str_manip(SPMmat, 'H');

% Load SPM.mat
if exist(SPMmat, 'file')
    load(SPMmat, 'SPM');
else
    disp(['SPM.mat file not found: ' SPMmat]);
    return
end

% Change statistic images path
%--------------------------------------------------------------------------
if isfield(SPM, 'swd')
    % Get old SPM path
    SPMold  = SPM.swd;

    % Set new SPM path
    SPM.swd = SPMnew;

    if isfield(SPM, 'Vbeta')
        new_fnames = strrep({SPM.Vbeta.fname}, SPMold, SPMnew);
        for i = 1:numel(new_fnames)
            SPM.Vbeta(i).fname = new_fnames{i};
        end
    end

    if isfield(SPM, 'VResMS')
        new_fnames = strrep({SPM.VResMS.fname}, SPMold, SPMnew);
        SPM.VResMS.fname = new_fnames{1};
    end

    if isfield(SPM, 'VM')
        new_fnames = strrep({SPM.VM.fname}, SPMold, SPMnew);
        SPM.VM.fname = new_fnames{1};
    end

    if isfield(SPM, 'xCon')
        for i = 1:numel(SPM.xCon)
            new_fnames = strrep({SPM.xCon(i).Vcon.fname}, SPMold, SPMnew);
            SPM.xCon(i).Vcon.fname = new_fnames{1};
            new_fnames = strrep({SPM.xCon(i).Vspm.fname}, SPMold, SPMnew);
            SPM.xCon(i).Vspm.fname = new_fnames{1};
        end
    end
end

% Change functional images path
%--------------------------------------------------------------------------
% current funtional images' path
fdir  = spm_str_manip(SPM.xY.P, 'H');
% current functional image's name
fname = spm_str_manip(SPM.xY.P, 't');

for i = 1:size(fdir,1)
    idx = find(strcmp(old_fdir, fdir(i,:)));
    new_fname = fullfile(new_fdir{idx}, fname(i,:));
    SPM.xY.P(i,1:length(new_fname)) = new_fname;
    SPM.xY.VY(i).fname = new_fname;
end
SPM.xY.P = deblank(SPM.xY.P);

% Save SPM.mat
save(SPMmat, 'SPM');

end

