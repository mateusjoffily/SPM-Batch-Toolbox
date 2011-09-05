function batch_launcher(data, param, runmode)
%
% BATCH_LAUNCHER Launch SPM5 compatible jobs.
%
%   data    - See 'data_example.m' for DATA structure definition
%   param   - See 'param_example.m' for PARAM structure definition
%   runmode - 'local_submit', 'local_build', 
%             'remote_submit' or 'remote_build'
%
%-----------------------------------------------------------------------
% Mateus Joffily - CNC/ISC - July/06
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

% Initialize jobs
jobs = {};

% Launch spm jobs
%--------------------------------------------------------------------------
if strcmp(runmode, 'local_submit') || ...
   strcmp(runmode, 'local_build')
    % Submit jobs locally
    %----------------------------------------------------------------------
    
    % Create jobs struct (see spm_jobman.m)
    for n = 1:numel(data)
        if numel(param) > 1
            np = n;
        else
            np = 1;
        end
        if iscell(param)
            pp = param{np};
        else
            pp = param(np);
        end
        if iscell(data)
            dd = data{n};
        else
            dd = data(n);
        end
        jobs{n} = batch_builder(dd, pp);
    end
    
    % Launch jobs
    if strcmp(runmode, 'local_submit')
        spm('defaults','fmri');          % spm defaults
        spm_jobman('run', [jobs{:}]);    % launch jobs
    end
    
elseif strcmp(runmode, 'remote_submit') || ...
       strcmp(runmode, 'remote_build')
    % Submit jobs remotely
    %----------------------------------------------------------------------
    
    % Parses jobs to each node
    for n = 1:numel(data)
        if numel(param) > 1
            np = n;
        else
            np = 1;
        end
        if iscell(param)
            pp = param{np};
        else
            pp = param(np);
        end
        if iscell(data)
            dd = data{n};
        else
            dd = data(n);
        end
        % Get job ID available between 1-9999
        for jobID = 1:9999
            cmdstr = sprintf('qjob job%04d', jobID);
            [status,msg] = unix(cmdstr);
            if status ~= 0,
                error([mfilename ': ' msg]);
                break;
            end
            if strfind(msg, 'No job found')
                break;
            end
        end
        batch_remote(dd, pp, jobID, runmode);
    end
    
    if strcmp(runmode, 'remote_submit')
        % Display job status
        pause(1);
        unix('qjob');   
    end
    
    disp('Done.');
    disp('Type "!qjob" for job monitoring.');
else
    % Unknown option
    %----------------------------------------------------------------------
    
    disp([mfilename ': unknown <runmode>: ' runmode]);
end

end