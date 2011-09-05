%=======================================================================
% Define data for batch jobs.
%

% initialize data - do not change this line
data = struct('preproc', [], 'model', [], 'util', []);

%=======================================================================
% Pre-processing data

% Session/Subjet 1
%-----------------------------------------------------------------------
% Directories where functional images are saved
data.preproc(1).fdir{1}='./s1/run1';
data.preproc(1).fdir{2}='./s1/run2';
data.preproc(1).fdir{3}='./s1/run3';
% Directory where 3D image is saved
data.preproc(1).sdir='./s1/3d';

% Session/Subjet 2
%-----------------------------------------------------------------------
% Directories where functional images are saved
data.preproc(2).fdir{1}='./s2/run1';
data.preproc(2).fdir{2}='./s2/run2';
data.preproc(2).fdir{3}='./s2/run3';
% Directory where 3D image is saved
data.preproc(2).sdir='./s2/3d';

% Session/Subjet 3
%-----------------------------------------------------------------------
% Directories where functional images are saved
data.preproc(3).fdir{1}='./s3/run1';
data.preproc(3).fdir{2}='./s3/run2';
data.preproc(3).fdir{3}='./s3/run3';
% Directory where 3D image is saved
data.preproc(3).sdir='./s3/3d';


%=======================================================================
% Model data

% Fixed effect models:
%-----------------------------------------------------------------------

% Example 1: Create a design matrix for each subject/session
%-----------------------------------------------------------------------
% Functional images directories
data.model(1).fdir{1}='./s1/run1';
data.model(1).fdir{2}='./s1/run2';
data.model(1).fdir{3}='./s1/run3';
data.model(2).fdir{1}='./s2/run1';
data.model(2).fdir{2}='./s2/run2';
data.model(2).fdir{3}='./s2/run3';
data.model(3).fdir{1}='./s3/run1';
data.model(3).fdir{2}='./s3/run2';
data.model(3).fdir{3}='./s3/run3';
% Logfiles: .mat file containing 'names', 'onsets' and 'durations' cells
data.model(1).log{1}='./s1/log/run1.mat';
data.model(1).log{2}='./s1/log/run2.mat';
data.model(1).log{3}='./s1/log/run3.mat';
data.model(2).log{1}='./s2/log/run1.mat';
data.model(2).log{2}='./s2/log/run2.mat';
data.model(2).log{3}='./s2/log/run3.mat';
data.model(3).log{1}='./s3/log/run1.mat';
data.model(3).log{2}='./s3/log/run2.mat';
data.model(3).log{3}='./s3/log/run3.mat';
% SPM.mat directories
data.model(1).mdir='./s1/model2';
data.model(2).mdir='./s2/model2';
data.model(3).mdir='./s3/model2';
% Contrast files: .mat files containing user defined contrasts.
% Leave it empty if you want the contrasts to be automaticaly created
% based on the conditions defined in the logfile.                                         
% The contrast files must contain a cell array of 'tcon' and/or 'fcon' 
% structures named 'consess'. The 'tcon' and 'fcon' structures 
% must contain 'name', 'convec' and 'sessrep' subfields.
% See SPM5 'Contrast Manager' documentation for more information
data.model(1).con='';
data.model(2).con='';
data.model(3).con='';

% Example 2: Create an unique design matrix with all subjects/sessions
%-----------------------------------------------------------------------
% Functional images directories
data.model(4).fdir{1}='./s1/run1';
data.model(4).fdir{2}='./s1/run2';
data.model(4).fdir{3}='./s1/run3';
data.model(4).fdir{4}='./s2/run1';
data.model(4).fdir{5}='./s2/run2';
data.model(4).fdir{6}='./s2/run3';
data.model(4).fdir{7}='./s3/run1';
data.model(4).fdir{8}='./s3/run2';
data.model(4).fdir{9}='./s3/run3';
% Logfiles (matfile containing 'names', 'onsets' and 'durations' cells)
data.model(4).log{1}='./s1/log/run1.mat';
data.model(4).log{2}='./s1/log/run2.mat';
data.model(4).log{3}='./s1/log/run3.mat';
data.model(4).log{4}='./s2/log/run1.mat';
data.model(4).log{5}='./s2/log/run2.mat';
data.model(4).log{6}='./s2/log/run3.mat';
data.model(4).log{7}='./s3/log/run1.mat';
data.model(4).log{8}='./s3/log/run2.mat';
data.model(4).log{9}='./s3/log/run3.mat';
% Contrast files
data.model(4).con='';
% SPM.mat directories
data.model(4).mdir='./model_all';

%=======================================================================
