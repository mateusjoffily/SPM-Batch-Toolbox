%=======================================================================
% Define data for batch jobs.
%

% initialize data - do not change this line
data = struct('preproc', [], 'model', [], 'util', []);

%=======================================================================
% Pre-processing data
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
% - In this example, we have one single subject/session (number of 
%   elements in data.preproc) and three runs (number of elements in 
%   data.preproc.fdir). Only one structural (3d) image can be specified 
%   for each subject/session.

% Session/Subjet 2
%-----------------------------------------------------------------------
% Directories where functional images are saved:
data.preproc(1).fdir{1}='./s2/run1';
data.preproc(1).fdir{2}='./s2/run2';
data.preproc(1).fdir{3}='./s2/run3';
% Directory where 3D image is saved
data.preproc(1).sdir='./s2/3d';


%=======================================================================
% Model data
%-----------------------------------------------------------------------
% In total, four different models will be created in this set of 
% examples (number of elements in data.model):
 
%-----------------------------------------------------------------------
% - In this example, a single design matrix will be created to model  
%   the three runs.

% Functional images directories
data.model(1).fdir{1}='./s2/run1';
data.model(1).fdir{2}='./s2/run2';
data.model(1).fdir{3}='./s2/run3';
% Logfiles: .mat files containing 'names', 'onsets' and 'durations' cells
% See SPM5 documentation for 'Multiple Conditions'
data.model(1).log{1}='./s2/log/run1.mat';
data.model(1).log{2}='./s2/log/run2.mat';
data.model(1).log{3}='./s2/log/run3.mat';
% SPM.mat directory
data.model(1).mdir='./s2/model1a';
% Contrast files: .mat files containing user defined contrasts.
% Leave it empty if you want the contrasts to be automaticaly created
% based on the conditions defined in the logfile.                                         
% The contrast files must contain a cell array of 'tcon' and/or 'fcon' 
% structures named 'consess'. The 'tcon' and 'fcon' structures 
% must contain 'name', 'convec' and 'sessrep' subfields.
% See SPM5 'Contrast Manager' documentation for more information
data.model(1).con='';

%-----------------------------------------------------------------------
% - In this example, a design matrix will be created to model  
%   each run separately.

% Run 1
%-----------------------------------------------------------------------
% Functional images directories
data.model(2).fdir{1}='./s2/run1';
% Logfiles: .mat files containing 'names', 'onsets' and 'durations' cells
% See SPM5 documentation for 'Multiple Conditions'
data.model(2).log{1}='./s2/log/run1.mat';
% SPM.mat directory
data.model(2).mdir='./s2/model1b';
% Contrast files: .mat files containing user defined contrasts.
% See SPM5 'Contrast Manager' documentation for more information
data.model(2).con='';

% Run 2
%-----------------------------------------------------------------------
% Functional images directories
data.model(3).fdir{1}='./s2/run2';
% Logfiles: .mat files containing 'names', 'onsets' and 'durations' cells
% See SPM5 documentation for 'Multiple Conditions'
data.model(3).log{1}='./s2/log/run2.mat';
% SPM.mat directory
data.model(3).mdir='./s2/model1c';
% Contrast files: .mat files containing user defined contrasts.
% See SPM5 'Contrast Manager' documentation for more information
data.model(3).con='';

% Run 3
%-----------------------------------------------------------------------
% Functional images directories
data.model(4).fdir{1}='./s2/run3';
% Logfiles: .mat files containing 'names', 'onsets' and 'durations' cells
% See SPM5 documentation for 'Multiple Conditions'
data.model(4).log{1}='./s2/log/run3.mat';
% SPM.mat directory
data.model(4).mdir='./s2/model1d';
% Contrast files: .mat files containing user defined contrasts.
% See SPM5 'Contrast Manager' documentation for more information
data.model(4).con='';


%==========================================================================
% Util data
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
% - Define the set of images to which SPM5 util task defined in 
%   script para_example.m will be applied.(see para_example.m)

% Session/Subjet 1
%-----------------------------------------------------------------------
% Directories where images are saved
data.util(1).idir{1}='./s2/run1';
data.util(2).idir{1}='./s2/run2';
data.util(3).idir{1}='./s2/run3';
