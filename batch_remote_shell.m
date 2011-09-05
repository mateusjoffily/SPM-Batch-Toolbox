function batch_remote_shell(fbsh, fjob, fmat)
%
% BATCH_REMOTE_SHELL Create shell script for remote job submission.
%
%-----------------------------------------------------------------------
% Adpated from codes written by Anne Cheylus and Sylvain Maurin 
% - Mateus Joffily - CNC/CNRS - July/2007

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

% Create shell script
[fid,msg] = fopen(fbsh,'w');
if fid == -1,
    error(['Error openning ' fbsh ' : ' msg]);
end;

% Shell type
fprintf(fid,'#!/bin/bash\n');
% BQS options
fprintf(fid,'#PBS -q T\n');
fprintf(fid,'#PBS -l platform=LINUX # Plateforme d''execution\n');
fprintf(fid,'#PBS -l model=Xeon\n');
fprintf(fid,'#PBS -l T=320000 # Nombre d''unites normalisees (consommation cpu)\n');
fprintf(fid,'#PBS -l M=2048MB # Memoire en MB\n');
fprintf(fid,'#PBS -l scratch=5GB # Taille du scratch en GB\n');
fprintf(fid,'#PBS -l spool=25MB # Taille du spool en MB\n');
fprintf(fid,'#PBS -l u_sps_isc\n');
fprintf(fid,'#PBS -eo # Redirection du flux stderr vers stdout\n');
fprintf(fid,'#PBS -mb # Mail begin\n');
fprintf(fid,'#PBS -me # Mail end\n');

% Virtual framebuffer X server
fprintf(fid,'echo localhost >$HOME/Xauth.localhost\n');
fprintf(fid,'( Xvfb :6666 -screen 0 800x600x24 -auth $HOME/Xauth.localhost >/dev/null 2>&1 ) & \n');
fprintf(fid,'DISPLAY=localhost:6666\n');
fprintf(fid,'export DISPLAY\n');

% Matlab command
mpath   = spm_str_manip(mfilename('fullpath'), 'H');    % Add batch_toolbox path
spmpath = spm_str_manip(which('spm'), 'H');             % Add spm path
fprintf(fid,'matlab -r "addpath %s %s; batch_remote_node(''%s'',''%s'', ''run''); quit;"\n', mpath, spmpath, fmat, fjob);

% Copy node's 'jobs' structure (created by batch_builder.m)
fprintf(fid,'cp %s %s\n', fjob, fullfile(pwd, '.'));

% Copy any *.ps created by spm
fprintf(fid,'cp spm_*.ps %s\n', fjob, fullfile(pwd, '.'));

% Clean directory
fprintf(fid,'rm -R $TMPBATCH/*\n');

fclose(fid);

% Set shell script executable
[s,msg,msgid] = fileattrib(fbsh,'+x','a');
if s == 0 
    error(msgid, msg);
end