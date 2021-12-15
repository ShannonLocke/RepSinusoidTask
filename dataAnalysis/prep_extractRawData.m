function [] = prep_extractRawData(all_sID)
% This script will run the prepare the raw data for the Repeated-Sinusoids
% task (Locke, Goettker, Gegenfurtner, & Mamassian, 2021).
%
% INPUTS:
% 1) all_sID - vector of subject IDs to be extracted
%
% OUTPUTS (per participant):
% 1) main-task raw data file for further processing in Matlab
% 2) training & main-task raw eye data file for OSF (CSV format)
%
% Data: https://osf.io/j9asn/
% Code: https://github.com/ShannonLocke/RepSinusoidTask
%
% Created by SML Dec 2021
% License: CC-By Attribution 4.0 International

%% Preamble:

% Inputs:
if ~isvector(all_sID)
    error('Error: the input all_sID must be a vector.')
end
nSs = length(all_sID); % number of participants

% Directories:
expName = 'SPC';
dataFromPath = '../data/';
dataToPath_matFiles = 'output_data/raw/';
dataToPath_osfFiles = 'output_data/forOSF/';

% Add paths to necessary files:
try
    addpath(genpath('/Users/shannonlocke/Documents/MATLAB/edf-converter-master/'))
catch
    warning('EDF converter not found in assumed place. Add to your Matlab path manually, otherwise the code will not run.')
end
% addpath(genpath([cd '../experimentToolbox'])) % <== WHY THIS?

%% Extract Raw Data:

for nn = 1:nSs % EACH participant
    
    %% Get stimulus and confidence data:
    
    % Load data:
    sID = all_sID(nn); % participants' ID
    dataPath = [dataFromPath, 's', num2str(sID), '/'];
    dataFilePattern = ['/results_' expName '_*.mat'];
    dataOnFile = dir([dataPath, dataFilePattern]);
    fname = [dataPath, '/', dataOnFile.name];
    load(fname, 'expData')
    
    % Extract key info:
    % pixPerDeg = expData.calibration.pixPerDeg; % pixels per degree for conversion
    nTrainTrials = length(expData.expDesign.designMat{1}); % number of training trials
    nTrials = length(expData.expDesign.designMat{end}); % number of test trials
    idxR = strcmp(expData.expDesign.modePhases,'training'); % training sessions
    idxT = find(~cellfun(@isempty,expData.res.resp)); % non-empty test sessions
    nSessions = length(idxT); % number of sessions
    
    % Experiment data:
    trial_training = repmat((1:nTrainTrials)',[nSessions,1]);
    trial = repmat((1:nTrials)',[nSessions,1]);
    N_training = length(trial_training);
    N = length(trial);
    session_training = repelem((1:nSessions)',nTrainTrials);
    session = repelem((1:nSessions)',nTrials);
    
    % Stimulus data:
    getData = cell2mat(expData.expDesign.designMat(idxR)');
    signedTrajectoryID_training = getData(:,2) .* getData(:,1); % (only for OSF)
    getData = cell2mat(expData.expDesign.designMat(idxT)');
    signedTrajectoryID = getData(:,2) .* getData(:,1);
    
    % Prepare structure for .mat file:
    eyeData.sID = all_sID(nn);
    eyeData.session = session;
    eyeData.trial = trial;
    eyeData.trajectory = signedTrajectoryID;
    
    %% Get eye data:
    
    % Setup info:
    sCenter = expData.hardware.sCenter; % screen center coordinates
    pixPerDeg = 80; % testing room
    
    % Preallocate:
    eyeData.eyePosX = cell(N,1);
    eyeData.eyePosY = cell(N,1);
    eyeData.pupilSize = cell(N,1);
    eyePosX_training = {};
    eyePosY_training = {};
    pupilSize_training = {};
    
    for ss = 1:nSessions
        
        % Get training-task data:
        fname = [dataPath, 'R', num2str(ss), '_', num2str(sID), '.edf'];
        [eyePosX, eyePosY, pupilSize] = readEyeData(fname, sCenter, pixPerDeg);
        eyePosX_training = [eyePosX_training; eyePosX];
        eyePosY_training = [eyePosY_training; eyePosY];
        pupilSize_training = [pupilSize_training; pupilSize];
        
        % Get main-task data:
        fname = [dataPath, 'T', num2str(ss), '_', num2str(sID), '.edf'];
        [eyePosX, eyePosY, pupilSize, sf] = readEyeData(fname, sCenter, pixPerDeg);
        idx = (1:nTrials) + (ss-1)*nTrials;
        eyeData.eyePosX(idx) = eyePosX;
        eyeData.eyePosY(idx) = eyePosY;
        eyeData.pupilSize(idx) = pupilSize;
        
    end
    
    %% Export main-task file for further Matlab processing:
    eyeData.sf = sf;
    fname = [dataToPath_matFiles, 'eyeData_s', num2str(sID), '.mat'];
    save(fname, 'eyeData')
    
    %% Export all eye data to csv for OSF:
    
    % Create table:
    for ii = 1:N_training
        
    end
    for ii = 1:N
        
    end
    
    % Save .csv file:
    dataPath = [dataToPath_osfFiles, 's', num2str(sID), '/'];
    % ==> Raw eye data
        
end

end

function [eyePosX, eyePosY, pupilSize, sf] = readEyeData(fname, sCenter, pixPerDeg)

% Convert EDF file to MAT:
Data = Edf2Mat(fname);

%% Find message markers:

% nMessages = length(Data.Events.Messages.time); % number of mesasges recorded

% Start of stimulus presentation indices:
idx_begin = strcmp(Data.Events.Messages.info, 'TRIAL_STIMSTART'); % index of stimulus presentation end events
time_begin = Data.Events.Messages.time(idx_begin); % timing of stimulus presentation end events
time_begin = sort(time_begin,'ascend'); % Make sure stimulus presentation endings are in order

% End of stimulus presentation indices:
idx_stimStop = strcmp(Data.Events.Messages.info, 'TRIAL_STIMSTOP'); % index of stimulus presentation end events
time_stimStop = Data.Events.Messages.time(idx_stimStop); % timing of stimulus presentation end events
time_stimStop = sort(time_stimStop,'ascend'); % Make sure stimulus presentation endings are in order

%% Extract the eye data:

% Find recording eye:
idx_eye = Data.Samples.gx(1,:) ~= Data.MISSING_DATA_VALUE;

% Eye-related data:
N = length(time_begin);
eyePosX = cell([N,1]);
eyePosY = cell([N,1]);
pupilSize = cell([N,1]);
for trial = 1:N % EACH trial
    idx_trial = Data.Samples.time >= time_begin(trial) &  Data.Samples.time <= time_stimStop(trial);
    getX = Data.Samples.gx(idx_trial,idx_eye);
    getY = Data.Samples.gy(idx_trial,idx_eye);
    eyePosX{trial} = (getX - sCenter(1)) / pixPerDeg;
    eyePosY{trial} = (getY - sCenter(2)) / pixPerDeg;
    pupilSize{trial} = Data.Samples.pupilSize(idx_trial);
end
sf = Data.RawEdf.RECORDINGS.sample_rate;  % sampling frequency
sf = double(sf);

end