function batch_remote(data, param, jobID, runmode)
%
% BATCH_REMOTE Remote job submission.
%
%   data    - see 'batch_data_example.m' for DATA structure definition
%   param   - see 'batch_param_example.m' for PARAM structure definition
%   jobID   - job ID
%   runmode - 'remote_submit'  - submit job remotely
%             'remote_build'   - build remote jobs struct locally
%                                (do not launch spm_jobman.m)
%
%-----------------------------------------------------------------------
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

if nargin < 4
    runmode = 'remote_submit';
end

% Force full pathnames in data structure
data = batch_data_fullpath(data);

% Create new diretory and launch job from there
jobdir = fullfile(pwd, sprintf('job%04d', jobID));
[status,msg, msgid] = mkdir(jobdir);
if status == 0 
    error(msgid, msg); 
end

% Keep track of current directory
cwd = pwd;

% Change to new directory
cd(jobdir);

% Create node specific matfile
% This file will contain all matlab variables required for node start a 
% SPM processing in a specific node. It includes 'jobs', 'param' and 
% 'data' structures
fmat = fullfile(pwd, sprintf('job%04d.mat', jobID));
save(fmat, 'data', 'param');

% Create node specific shell script
% This script will contain all the system variables required to start
% a job in a specific node.
fbsh = fullfile(pwd, sprintf('job%04d.bsh', jobID));
fjob = sprintf('job%04d_bb.mat', jobID);  
batch_remote_shell(fbsh, fjob, fmat);

% Launch job
jobname = sprintf('job%04d', jobID);
switch runmode
    case 'remote_submit'
        disp(sprintf('Submitting job%04d...\n',jobID));
        
        cmdstr = sprintf('qsub -N %s %s', jobname, fbsh);
        [status,msg] = unix(cmdstr);
        if status ~= 0,
            error([mfilename ': ' msg]);
        end

    case 'remote_build'
        disp(sprintf('Making job%04d...\n',jobID));
        
        % Create node diretcory
        % Set environment variable TMPBATCH
        setenv('TMPBATCH', jobdir);
        
        % Launch batch_remote_node
        batch_remote_node(fmat, fjob, 'emulate');
        
        disp(sprintf('\njob%04d done.\n\n',jobID));
end

% Go back to current working directory
cd(cwd);

end

%==========================================================================
function fulldata = batch_data_fullpath(data)
%
% BATCH_DATA_FULLPATH Force full pathnames in data structure.
%
%--------------------------------------------------------------------------
% Mateus Joffily - CNC/CNRS - July/2007

% Loop over preprocessing data
for n = 1:numel(data.preproc)
    % Functional files
    for d = 1:numel(data.preproc(n).fdir)
        fulldata.preproc(n).fdir{d} = ...
                            spm_select('CPath',data.preproc(n).fdir{d});
    end
    % Structural file
    fulldata.preproc(n).sdir = spm_select('CPath',data.preproc(n).sdir);
end

% Loop over statistics data
for n = 1:numel(data.model)
    % Functional files
    for d = 1:numel(data.model(n).fdir)
    	fulldata.model(n).fdir{d} = ...
                   spm_select('CPath',data.model(n).fdir{d});
    end
    % Log files
    for d = 1:numel(data.model(n).log)
    	fulldata.model(n).log{d} = ...
                   spm_select('CPath',data.model(n).log{d});
    end
    % SPM.mat file
    fulldata.model(n).mdir = spm_select('CPath',data.model(n).mdir);
end

end