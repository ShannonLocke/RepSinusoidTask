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
    trial = repmat((1:nTrials)',[nSessions,1]);
    session = repelem((1:nSessions)',nTrials);
    
    % Stimulus and confidence data:
    getData = cell2mat(expData.expDesign.designMat(idxR)');
    trajectoryID_training = getData(:,1); % (only for OSF)
    direction_training = getData(:,2); % (only for OSF)
    getData = cell2mat(expData.expDesign.designMat(idxT)');
    trajectoryID = getData(:,1);
    direction = getData(:,2);
    conf = cell2mat(expData.res.resp(idxT)');
    RT = cell2mat(expData.res.RT(idxT)');
    
    %% Get eye data:
    
    % Pre-allocate results vectors:
    
    for ss = 1:nSessions
        
        % Get training-task data:
        pidx = (ss-1)*2 + 1; % phase index
        respYN = false;
        fname = [dataPath, 'R', num2str(ss), '_', num2str(sID), '.edf'];
        getData_training = readEyeData(fname, nTrainTrials, respYN);
        
        
        % Get main-task data:
        pidx = ss*2; % phase index
        targX = expData.res.targX{pidx}; % target position
        fname = [dataPath, 'T', num2str(ss), '_', num2str(sID), '.edf'];
        respYN = true;
        getData = readEyeData(fname, nTrainTrials, respYN);
        
        
        
    end
    
    %% Export main-task file for further Matlab processing:
    fname = [dataToPath_matFiles, 'rawData_s', num2str(sID), '.mat'];
    load(fname, 'rawData')
    
    %% Export all eye data to csv for OSF:
    dataPath = [dataToPath_osfFiles, 's', num2str(sID), '/'];
    % ==> Raw trial data (goes into 
    % ==> Raw eye data
        
end

end

function [eyeData] = readEyeData(fname, nTrials, respYN)

% Convert EDF file to MAT:
Data = Edf2Mat(fname);

%% Find message markers:
nMessages = length(Data.Events.Messages.time); % number of mesasges recorded

% Start of stimulus presentation indices:
idx_begin = strcmp(Data.Events.Messages.info, 'TRIAL_STIMSTART'); % index of stimulus presentation end events
time_begin = Data.Events.Messages.time(idx_begin); % timing of stimulus presentation end events
time_begin = sort(time_begin,'ascend'); % Make sure stimulus presentation endings are in order

% End of stimulus presentation indices:
idx_stimStop = strcmp(Data.Events.Messages.info, 'TRIAL_STIMSTOP'); % index of stimulus presentation end events
time_stimStop = Data.Events.Messages.time(idx_stimStop); % timing of stimulus presentation end events
time_stimStop = sort(time_stimStop,'ascend'); % Make sure stimulus presentation endings are in order

% Response indices:
idx_resp = strcmp(Data.Events.Messages.info, 'TRIAL_RESP'); % index of trial response events
if ~any(idx_resp)
    respYN = false;
    time_resp = [];
else
    respYN = true;
    time_resp = Data.Events.Messages.time(idx_resp); % timing of trial response events
    time_resp = sort(time_resp,'ascend'); % Make sure responses are in order
end

% End of trial indices:
idx_end = strcmp(Data.Events.Messages.info, 'TRIAL_END'); % index of trial end events
time_end = Data.Events.Messages.time(idx_end); % timing of trial end events
time_end = sort(time_end,'ascend'); % Make sure endings are in order

% Blinks and saccades:
time_sBlink = Data.Events.Eblink.start; % start of blink
time_eBlink = Data.Events.Eblink.end; % end of blink
time_sSacc = Data.Events.Esacc.start; % start of saccade
time_eSacc = Data.Events.Esacc.end; % end of saccade
if length(time_sBlink) > length(time_eBlink)
    warning('Mismatch in blink indices!')
elseif length(time_sBlink) < length(time_eBlink)
    warning('Mismatch in blink indices!')
end

%% Extract the eye data:

% Find recording eye:
idx_eye = Data.Samples.gx(1,:) ~= Data.MISSING_DATA_VALUE;

% Sampling properties:
eyeData.sf = Data.RawEdf.RECORDINGS.sample_rate; eyeData.sf = double(eyeData.sf); % sampling frequency
eyeData.nSamplesPerTrial = (time_end - time_begin) + 1; % How many samples per trial
maxN = max(eyeData.nSamplesPerTrial);

% Trial event timing:
eyeData.timeInTrial = (1/eyeData.sf) * (0:1:(maxN-1));
eyeData.timeOfStimStop = (1/eyeData.sf) * (time_stimStop - time_begin);
eyeData.timeOfResp = (1/eyeData.sf) * (time_resp - time_begin);
eyeData.timeOfTrialEnd = (1/eyeData.sf) * (time_end - time_begin);

% Blinks and saccades:
eyeData.nBlinks = NaN([1,nTrials]);
eyeData.nSaccades = NaN([1,nTrials]);
eyeData.BlinksYN = false([maxN,nTrials]);
eyeData.SaccadeYN = false([maxN,nTrials]);
for trial = 1:nTrials
    time_trial = time_begin(trial):time_end(trial);
    trial_blinks = ismember(time_sBlink,time_trial);
    trial_saccades = ismember(time_sSacc,time_trial);
    eyeData.nBlinks(trial) = sum(trial_blinks);
    eyeData.nSaccades(trial) = sum(trial_saccades);
    if eyeData.nBlinks(trial) > 1 % blinks detected
        aa = find(trial_blinks);
        for ii = 1:eyeData.nBlinks(trial)
            idxS = time_sBlink(aa(ii)) - time_begin(trial) + 1;
            idxE = time_eBlink(aa(ii)) - time_begin(trial) + 1;
            eyeData.BlinkYN(idxS:idxE,trial) = true;
            eyeData.BlinkDur(ii,trial) = (1/eyeData.sf) * (idxE - idxS);
        end
    end
    if eyeData.nSaccades(trial) > 1 % saccades detected
        aa = find(trial_saccades);
        for ii = 1:eyeData.nSaccades(trial)
            idxS = time_sSacc(aa(ii)) - time_begin(trial) + 1;
            idxE = time_eSacc(aa(ii)) - time_begin(trial) + 1;
            eyeData.SaccadeYN(idxS:idxE,trial) = true;
            eyeData.SaccadeDur(ii,trial) = (1/eyeData.sf) * (idxE - idxS);
        end
    end
end

% Eye-related data:
eyePosX = NaN([maxN,nTrials]);
eyePosY = NaN([maxN,nTrials]);
pupilSize = NaN([maxN,nTrials]);
% pupilSizeNorm = NaN([maxN,nTrials]);
for trial = 1:nTrials % EACH trial
    idx_trial = Data.Samples.time >= time_begin(trial) &  Data.Samples.time <= time_end(trial);
    idx_samples = 1:eyeData.nSamplesPerTrial(trial);
    eyePosX(idx_samples,trial) = Data.Samples.gx(idx_trial,idx_eye);
    eyePosY(idx_samples,trial) = Data.Samples.gy(idx_trial,idx_eye);
    getPupil = Data.Samples.pupilSize(idx_trial);
    pupilSize(idx_samples,trial) = getPupil;
    % pupilSizeNorm(idx_samples,trial) = (getPupil - nanmean(getPupil)) / nanstd(getPupil);
end
% eyeData.eyePosX_pix = eyePosX;
% eyeData.eyePosY_pix = eyePosY;
eyeData.eyePosX_deg = (eyePosX - expData.hardware.sCenter(1)) / pixPerDeg;
eyeData.eyePosY_deg = (eyePosY - expData.hardware.sCenter(2)) / pixPerDeg;
eyeData.pupilSize = pupilSize;
% eyeData.pupilSizeNorm = pupilSizeNorm;

% Velocity and Acceleration:
pOrder = 2;
fLength = 51;
dx = denoiseSG(eyeData.eyePosX_deg,eyeData.sf,pOrder,fLength,0);
eyeData.eyeVelX = dx(1:(end-1),:,2);
eyeData.eyeVelX(end,:) = 0;
eyeData.eyeAccX = dx(1:(end-1),:,3);
eyeData.eyeAccX((end-3):end,:) = 0;

%% Add the stimulus and confidence data:

% Basic stimulus trace:
stimDur = 6;
eyeData.sf_stim = 60;
eyeData.conf = expData.res.resp{end};
eyeData.targT = 0:(1/eyeData.sf_stim):(stimDur-(1/eyeData.sf_stim));
eyeData.targX = expData.res.targX{end};

% Interpolate:
eyeData.tSteps_error = 0:(1/eyeData.sf):(stimDur-(1/eyeData.sf));
N = length(eyeData.tSteps_error);
eyeData.targX_int = NaN([N,nTrials]);
for trial = 1:nTrials
    eyeData.targX_int(:,trial) = interp1(eyeData.targT, eyeData.targX(:,trial), eyeData.tSteps_error);
end

% Velocity and Acceleration:
dx = denoiseSG(eyeData.targX_int,eyeData.sf,pOrder,fLength,0);
eyeData.targVelX = dx(1:(end-1),:,2);
eyeData.targVelX(end,:) = 0;
eyeData.targAccX = dx(1:(end-1),:,3);
eyeData.targAccX((end-3):end,:) = 0;

% Calculate error
eyeData.errorX = eyeData.eyePosX_deg(1:N,:) - eyeData.targX_int;
eyeData.errorY = eyeData.eyePosY_deg(1:N,:) - 0;
eyeData.errorEuclid = sqrt(eyeData.errorX.^2 + eyeData.errorY.^2);
eyeData.RMSE_X = sqrt(nanmean(eyeData.errorX.^2));
eyeData.RMSE_Euclid = sqrt(nanmean(eyeData.errorEuclid.^2));

%% Calculate cross-correlations

% Position:
tdat = eyeData.targX_int;
keep_idx = ~isnan(tdat(:,1));
tdat = tdat(keep_idx,:);
N = sum(keep_idx);
edat = eyeData.eyePosX_deg(1:N,:);
truncateBy = N - 1000;
[r,lag] = slidingCrossCorrelationCoefficient(tdat,edat,truncateBy);
eyeData.lag_pos = 1/eyeData.sf * (lag);
eyeData.xcorr_pos = r;

% Velocity:
tdat = eyeData.targVelX;
keep_idx = ~isnan(tdat(:,1));
tdat = tdat(keep_idx,:);
N = sum(keep_idx);
edat = eyeData.eyeVelX(1:N,:);
truncateBy = N - 1000;
[r,lag] = slidingCrossCorrelationCoefficient(tdat,edat,truncateBy);
eyeData.lag_vel = 1/eyeData.sf * (lag);
eyeData.xcorr_vel = r;

end