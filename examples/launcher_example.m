%==========================================================================
% Set job's 'data', 'param' and 'runmode' variables to be parsed by
% batch_launcher.m
%--------------------------------------------------------------------------

runmode = 'local_submit';    % 'local_submit', 'local_build'
                             % 'remote_submit', 'remote_build'

% Load same param for all subjects/sessions
param_example;

% - Subject 1
data_example_s1;
dataS{1} = data;

% - Subject 2
data_example_s2;
dataS{2} = data;

% - Subject 3
data_example_s3;
dataS{3} = data;

% - All subjects 
%   If runmode is 'remote_*', this line should not be included, otherwise
%   the data from all subjetcs defined in 'data_example_all' will be submitted to 
%   the same node. The correct way of entering data from several 
%   subjects for remote job submission is separating them into different
%   data.m files.
data_example_all;
dataS{4} = data;

%--------------------------------------------------------------------------
% Launch spm jobs
batch_launcher(dataS, param, runmode);
