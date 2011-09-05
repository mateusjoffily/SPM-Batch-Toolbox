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

% Session/Subjet 1
%-----------------------------------------------------------------------
% Directories where functional images are saved:
data.preproc(1).fdir{1}='./s1/run1';
data.preproc(1).fdir{2}='./s1/run2';
data.preproc(1).fdir{3}='./s1/run3';
% Directory where 3D image is saved
data.preproc(1).sdir='./s1/3d';


%=======================================================================
% Model data
%-----------------------------------------------------------------------
% In total, four different models will be created in this set of 
% examples (number of elements in data.model):
 
%-----------------------------------------------------------------------
% - In this example, a single design matrix will be created to model  
%   the three runs.

% Functional images directories
data.model(1).fdir{1}='./s1/run1';
data.model(1).fdir{2}='./s1/run2';
data.model(1).fdir{3}='./s1/run3';
% Logfiles: .mat files containing 'names', 'onsets' and 'durations' cells
% See SPM5 documentation for 'Multiple Conditions'
data.model(1).log{1}='./s1/log/run1.mat';
data.model(1).log{2}='./s1/log/run2.mat';
data.model(1).log{3}='./s1/log/run3.mat';
% SPM.mat directory
data.model(1).mdir='./s1/model1a';
% Contrast files: .mat files containing user defined contrasts.
% The contrast files must contain a cell array of 'tcon' and/or 'fcon' 
% structures named 'consess'. The 'tcon' and 'fcon' structures 
% must contain 'name', 'convec' and 'sessrep' subfields.
% See SPM5 'Contrast Manager' documentation for more information
data.model(1).con='';
% Factors files: .mat files containing user defined factors for factorial 
% design.
% The factors files must contain an array of 'fact' structures. 
% The 'fact' structure must contain 'name' and 'levels' subfields.
% See SPM5 'fMRI model specification' documentation for more information
data.model(1).fact='';
% Covariates file: .mat files containing user defined covariates and
% nuissance variables (only used by factorial design specification).
% !! NOT WORKING YET !! See batch_builder.m for more information.
% Leave it empty if you don't want to specify covariates and nuissance 
% variables.                                         
% The covariates file must contain a cell array of covariates
% structure named 'cov'. The covariates structure must contain 
% 'c', 'cname', 'iCFI' and 'iCC' subfields:
% 'c'     is a vector of covariate values
% 'cname' is the covariate's name
% 'iCFI'  creates an additional regressor that is the interaction 
%         between the covariate and a chosen experimental factor. It
%         can be 0='None', 1='With Factor 1', 2='With Factor 2' or
%         3='With Factor 3'
% 'iCC'   it can be 0='Overall mean', 1='Factor 1 mean', 
%         2='Factor 2 mean', 3='Factor 3 mean', 4='No centering',
%         5='User specified value', 6='As implied by ANCOVA' or 7='GM'
% See SPM5 'Factorial Design Specification' documentation for more information
%data.model(1).cov='';   

%-----------------------------------------------------------------------
% - In this example, a design matrix will be created to model  
%   each run separately.

% Run 1
%-----------------------------------------------------------------------
% Functional images directories
data.model(2).fdir{1}='./s1/run1';
% Logfiles: .mat files containing 'names', 'onsets' and 'durations' cells
% See SPM5 documentation for 'Multiple Conditions'
data.model(2).log{1}='./s1/log/run1.mat';
% SPM.mat directory
data.model(2).mdir='./s1/model1b';
% Contrast files: .mat files containing user defined contrasts.
% See SPM5 'Contrast Manager' documentation for more information
data.model(2).con='';
% Covariates files: .mat files containing user defined covariates and
% nuissance variables (only used by factorial design specification).
% See SPM5 'Factorial Design Specification' documentation for more information
data.model(2).cov='';

% Run 2
%-----------------------------------------------------------------------
% Functional images directories
data.model(3).fdir{1}='./s1/run2';
% Logfiles: .mat files containing 'names', 'onsets' and 'durations' cells
% See SPM5 documentation for 'Multiple Conditions'
data.model(3).log{1}='./s1/log/run2.mat';
% SPM.mat directory
data.model(3).mdir='./s1/model1c';
% Contrast files: .mat files containing user defined contrasts.
% See SPM5 'Contrast Manager' documentation for more information
data.model(3).con='';
% Covariates files: .mat files containing user defined covariates and
% nuissance variables (only used by factorial design specification).
% See SPM5 'Factorial Design Specification' documentation for more information
data.model(3).cov='';

% Run 3
%-----------------------------------------------------------------------
% Functional images directories
data.model(4).fdir{1}='./s1/run3';
% Logfiles: .mat files containing 'names', 'onsets' and 'durations' cells
% See SPM5 documentation for 'Multiple Conditions'
data.model(4).log{1}='./s1/log/run3.mat';
% SPM.mat directory
data.model(4).mdir='./s1/model1d';
% Contrast files: .mat files containing user defined contrasts.
% See SPM5 'Contrast Manager' documentation for more information
data.model(4).con='';
% Covariates files: .mat files containing user defined covariates and
% nuissance variables (only used by factorial design specification).
% See SPM5 'Factorial Design Specification' documentation for more information
data.model(4).cov='';


%==========================================================================
% Util data
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
% - Define the set of images to which SPM5 util task defined in 
%   script para_example.m will be applied.(see para_example.m)

% Session/Subjet 1
%-----------------------------------------------------------------------
% Directories where images are saved
data.util(1).idir{1}='./s1/run1';
data.util(2).idir{1}='./s1/run2';
data.util(3).idir{1}='./s1/run3';
