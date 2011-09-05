function status = batch_remote_node(fmat, fjob, runmode)
%
% BATCH_REMOTE_NODE Matlab routine for remote nodes.
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

% Get TMPBATCH pathname
rootdir = getenv('TMPBATCH');
if isempty(deblank(rootdir))
    error('System''s environment variable $TMPBATCH not found.');
end

% Load 'param' and 'data' structures
load(fmat);

% Create mirror of local directory tree in remote node
[rdata pdir] = create_remote(rootdir, data);

% Copy local data files to remote node
status = copy_data(pdir, [1 2], [1 1 0 1 1], param);
if status == 0
    error('Error copying files from local to remote directory'); 
end

% Create jobs structure, using remote node directory tree
jobs = batch_builder(rdata, param, fjob, false);

% If run mode is 'emulation', do not launch SPM job and return...
if strcmp(runmode, 'emulate')
    disp([mfilename ': done.']);
    return
end

% Load spm defaults
spm('defaults','fmri');

% Launch spm job
spm_jobman('run', jobs);

% Copy back files from remote node to local disk
status = copy_data(pdir, [2 1], [1 1 1 0 0]);
if status == 0
    disp('Error copying files from remote to local directory');
end

% Change SPM.mat paths to match local directories
for n = 1:numel(rdata.model)
    fSPM     = fullfile(data.model(n).mdir, 'SPM.mat');
    old_fdir = rdata.model(n).fdir;    % old functional images path
    new_fdir = data.model(n).fdir;     % new fucntional images path
    change_spm_path(fSPM, old_fdir, new_fdir);
end
   
% Done
disp([mfilename ': done.']);

end

%==========================================================================
function [status tmpdir] = create_tmpdir(rootdir)

% create new temporary directory name
tmpdir = fullfile(rootdir, tempname);

% create new directory
[status,msg, msgid] = mkdir(tmpdir);
if status == 0 
    error(msgid, msg); 
end

end

%==========================================================================
function [rdata pdir] = create_remote(rootdir, data)

rdata  = [];
pdir   = struct('f', [], 's', [], 'm', [], 'l', [], 'c', []);
rn     = [0 0 0 0 0];

% Preprocessing
%--------------------------------------------------------------------------
for n = 1:numel(data.preproc)
    % Create remote directories for functional images
    %----------------------------------------------------------------------
    for d = 1:numel(data.preproc(n).fdir)
        rn(1) = rn(1) + 1;
        pdir.f{1,rn(1)} = data.preproc(n).fdir{d};
        [status pdir.f{2,rn(1)}] = create_tmpdir(rootdir);
        rdata.preproc(n).fdir{d} = pdir.f{2,rn(1)};
        if status == 0
            error('error creating directory or copying files to node');
        end
    end
    % Create remote directory for strcutural images
    %----------------------------------------------------------------------
    rn(2) = rn(2) + 1;
    pdir.s{1,rn(2)} = data.preproc(n).sdir;
    [status pdir.s{2,rn(2)}] = create_tmpdir(rootdir);
    rdata.preproc(n).sdir = pdir.s{2,rn(2)};
    if status == 0
        error('error creating directory or copying files to node');
    end
end

% Statistics
%--------------------------------------------------------------------------
for n = 1:numel(data.model)
    % Create remote directories for functional images
    %----------------------------------------------------------------------
    for d = 1:numel(data.model(n).fdir)
        % Check if directory has already been created for preprocessing
        rdata.model(n).fdir{d} = '';
        for f = 1:size(pdir.f,2)
            if strcmp(data.model(n).fdir{d}, pdir.f{1,f})
                rdata.model(n).fdir{d} = pdir.f{2,f};
                break
            end
        end
        
        if isempty(rdata.model(n).fdir{d})
            % If directory doesn't exist, create it
            rn(1) = rn(1) + 1;
            pdir.f{1,rn(1)} = data.model(n).fdir{d};
            [status pdir.f{2,rn(1)}] = create_tmpdir(rootdir);
            rdata.model(n).fdir{d} = pdir.f{2,rn(1)};
            if status == 0
                error('error creating directory or copying files to node');
            end
        end
    end
    % Create remote directory for statistics (SPM.mat)
    %----------------------------------------------------------------------
    rn(3) = rn(3) + 1;
    pdir.m{1,rn(3)} = data.model(n).mdir;
    [status pdir.m{2,rn(3)}] = create_tmpdir(rootdir);
    rdata.model(n).mdir = pdir.m{2,rn(3)};
    if status == 0
        error('error creating directory in node');
    end
    % Create remote directories for logfiles
    %----------------------------------------------------------------------
    for d = 1:numel(data.model(n).log)
        % Check if directory has already been created
        logdir  = spm_str_manip(data.model(n).log{d}, 'H');
        logfile = spm_str_manip(data.model(n).log{d}, 't');
        rdata.model(n).log{d} = '';
        for f = 1:size(pdir.l,2)
            if strcmp(logdir, pdir.l{1,f})
                rdata.model(n).log{d} = fullfile(pdir.l{2,f}, logfile);
                pdir.l{3,f} = {pdir.l{3,f}{:} logfile};
                break
            end
        end
        
        if isempty(rdata.model(n).log{d})
            % If directory doesn't exist, create it
            rn(4) = rn(4) + 1;
            pdir.l{1,rn(4)} = logdir;
            pdir.l{3,rn(4)} = {logfile};
            [status pdir.l{2,rn(4)}] = create_tmpdir(rootdir);
            rdata.model(n).log{d} = fullfile(pdir.l{2,rn(4)}, logfile);
            if status == 0
                error('error creating directory or copying files to node');
            end
        end
    end 
    % Create remote directory for contrasts files
    %----------------------------------------------------------------------
    % Check if 'con' field exists
    if ~isfield(data.model(n), 'con')
        % If 'con' field doesn''t exist, create it empty. Contrasts will be
        % automaticaly created based on the conditions defined in the logfile. 
        rdata.model(n).con = '';
        continue;
    end
    
    condir  = spm_str_manip(data.model(n).con, 'H');
    confile = spm_str_manip(data.model(n).con, 't');
    rdata.model(n).con = '';
    
    % If confile is empty, we don't need to do anything
    if isempty(confile)
        continue
    end
    
    % Check if directory has already been created
    for f = 1:size(pdir.c,2)
        if strcmp(condir, pdir.c{1,f})
            rdata.model(n).con = fullfile(pdir.c{2,f}, confile);
            pdir.c{3,f} = {pdir.c{3,f}{:} confile};
            break
        end
    end

    if isempty(rdata.model(n).con)
        % If directory doesn't exist, create it
        rn(5) = rn(5) + 1;
        pdir.c{1,rn(5)} = condir;
        pdir.c{3,rn(5)} = {confile};
        [status pdir.c{2,rn(5)}] = create_tmpdir(rootdir);
        rdata.model(n).con = fullfile(pdir.c{2,rn(5)}, confile);
        if status == 0
            error('error creating directory or copying files to node');
        end
    end
end

end
%==========================================================================
function status = copy_data(pdir, direc, filetype, param)

if nargin < 4 
    param = struct('sprefix', '', 'fprefix', '', 'mprefix', '');
end

% Initialize status variable
status = 1;

% Copy functional files
if filetype(1)
    for i = 1:size(pdir.f,2)
        src = fullfile(pdir.f{direc(1),i}, [param.fprefix '*']);
        dst = pdir.f{direc(2),i};
        [status, msg] = copyfile(src, dst);
        if status == 0
        	disp(sprintf('Error using ==> %s>copy_data\n%s', mfilename, msg));
        end
        if ~isempty(param.fprefix) && ~isempty(param.mprefix)
            src = fullfile(pdir.f{direc(1),i}, [param.mprefix '*']);
            dst = pdir.f{direc(2),i};
            [status, msg] = copyfile(src, dst);
            if status == 0
                disp(sprintf('Error using ==> %s>copy_data\n%s', mfilename, msg));
            end
        end
    end
end

% Copy structural files
if filetype(2)
    for i = 1:size(pdir.s,2)
        src = fullfile(pdir.s{direc(1),i}, [param.sprefix '*']);
        dst = pdir.s{direc(2),i};
        [status, msg] = copyfile(src, dst);
        if status == 0
        	disp(sprintf('Error using ==> %s>copy_data\n%s', mfilename, msg));
        end
        if ~isempty(param.fprefix) && ~isempty(param.mprefix)
            src = fullfile(pdir.s{direc(1),i}, [param.mprefix '*']);
            dst = pdir.s{direc(2),i};
            [status, msg] = copyfile(src, dst);
            if status == 0
                disp(sprintf('Error using ==> %s>copy_data\n%s', mfilename, msg));
            end
        end
    end
end

% Copy model files
if filetype(3)
    for i = 1:size(pdir.m,2)
        % If source directory doesn't exist, don't do anything
        if ~exist(pdir.m{direc(1),i}, 'dir'), break; end
        % If destination directory doesn't exist, create it
        % This may be necessary if destination is in local disk.
        if ~exist(pdir.m{direc(2),i}, 'dir')
            [status,msg] = mkdir(pdir.m{direc(2),i});
            if status == 0
                disp(sprintf('Error using ==> %s>copy_data\n%s', mfilename, msg));
            end
        end
        src = fullfile(pdir.m{direc(1),i}, '*');
        dst = pdir.m{direc(2),i};
        [status, msg] = copyfile(src, dst);
        if status == 0
        	disp(sprintf('Warning! ==> %s>copy_data\n%s', mfilename, msg));
        end
    end
end

% Copy log files
if filetype(4)
    for i = 1:size(pdir.l,2)
        for j = 1:numel(pdir.l{3,i})
            src = fullfile(pdir.l{direc(1),i}, pdir.l{3,i}{j});
            dst = fullfile(pdir.l{direc(2),i}, pdir.l{3,i}{j});
            [status, msg] = copyfile(src, dst);
            if status == 0
                disp(sprintf('Error using ==> %s>copy_data\n%s', mfilename, msg));
            end
        end
    end
end

% Copy contrasts files
if filetype(5)
    for i = 1:size(pdir.c,2)
        for j = 1:numel(pdir.c{3,i})
            src = fullfile(pdir.c{direc(1),i}, pdir.c{3,i}{j});
            dst = fullfile(pdir.c{direc(2),i}, pdir.c{3,i}{j});
            [status, msg] = copyfile(src, dst);
            if status == 0
                disp(sprintf('Error using ==> %s>copy_data\n%s', mfilename, msg));
            end
        end
    end
end

end

