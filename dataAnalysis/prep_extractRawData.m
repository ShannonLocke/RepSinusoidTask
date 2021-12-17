function [] = prep_extractRawData(all_sID)
% This script will run the prepare the raw eye data for the Repeated-Sinusoids
% task (Locke, Goettker, Gegenfurtner, & Mamassian, 2021).
%
% INPUTS:
% 1) all_sID - vector of subject IDs to be extracted
%
% FILES GENERATED (per participant):
% 1) main-task raw eye data file for further processing in Matlab (*.mat)
% 2) training & main-task raw eye data file for OSF (*.csv)
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
    
    sID = all_sID(nn); % participants' ID
    disp(['Extracting raw data of s', num2str(sID)])
    
    %% Get stimulus and confidence data:
    
    % Load data:
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
    eyeData.nSamples = NaN([N,1]);
    eyeData.eyePosX = cell(N,1);
    eyeData.eyePosY = cell(N,1);
    eyeData.pupilSize = cell(N,1);
    eyeData.saccadeYN = cell(N,1);
    eyePosX_training = {};
    eyePosY_training = {};
    pupilSize_training = {};
    
    for ss = 1:nSessions
        disp(['... Session ', num2str(ss), ' ...'])
        
        % Get training-task data:
        disp('... ... Training data ... ...')
        fname = [dataPath, 'R', num2str(ss), '_', num2str(sID), '.edf'];
        [eyePosX, eyePosY, pupilSize] = readEyeData(fname, sCenter, pixPerDeg);
        eyePosX_training = [eyePosX_training; eyePosX];
        eyePosY_training = [eyePosY_training; eyePosY];
        pupilSize_training = [pupilSize_training; pupilSize];
        
        % Get main-task data:
        disp('... ... Test data ... ...')
        fname = [dataPath, 'T', num2str(ss), '_', num2str(sID), '.edf'];
        [eyePosX, eyePosY, pupilSize, saccadeYN, sf] = readEyeData(fname, sCenter, pixPerDeg);
        idx = (1:nTrials) + (ss-1)*nTrials;
        eyeData.eyePosX(idx) = eyePosX;
        eyeData.eyePosY(idx) = eyePosY;
        eyeData.pupilSize(idx) = pupilSize;
        eyeData.saccadeYN(idx) = saccadeYN;
        eyeData.nSamples(idx) = cellfun(@length,eyePosX);
        
    end
    
    %% Export data:
    
    % Main-task file for further Matlab processing:
    disp('... Exporting mat file ...')
    eyeData.sf = sf;
    fname = [dataToPath_matFiles, 'rawEyeData_s', num2str(sID), '.mat'];
    save(fname, 'eyeData')
    
    % Create table of training eye data for OSF csv:
    disp('... Exporting csv file ...')
    for ii = 1:N_training
        ns = length(eyePosX_training{ii}); % number of samples
        tmpT = table(repmat(session_training(ii), [ns,1]), ... % session
            ones([ns,1]), ... % trainingYN
            repmat(trial_training(ii), [ns,1]), ... % trial
            repmat(signedTrajectoryID_training(ii), [ns,1]), ... % trajectory
            (1:ns)', ... % sample number
            eyePosX_training{ii}, ... % horizontal eye position
            eyePosY_training{ii}, ... % vertical eye position
            pupilSize_training{ii}); % pupil size
        if ii == 1
            T = tmpT;
        else
            T = [T; tmpT];
        end
    end
    
    % Create table of test eye data for OSF csv:
    for ii = 1:N
        ns = length(eyeData.eyePosX{ii}); % number of samples
        tmpT = table(repmat(session(ii), [ns,1]), ... % session
            zeros([ns,1]), ... % trainingYN
            repmat(trial(ii), [ns,1]), ... % trial
            repmat(signedTrajectoryID(ii), [ns,1]), ... % trajectory
            (1:ns)', ... % sample number
            eyeData.eyePosX{ii}, ... % horizontal eye position
            eyeData.eyePosY{ii}, ... % vertical eye position
            eyeData.pupilSize{ii}); % pupil size
        T = [T; tmpT];
    end
    T.Properties.VariableNames = {'session', 'trainingYN', 'trial', ...
        'directionSignedTrajectoryID', 'sampleNumberAt1000HzSF', ...
        'eyePosXdeg', 'eyePosYdeg', 'pupilSizePix'};
    
    % Save .csv file:
    fname = [dataToPath_osfFiles, 'rawEyeData_s', num2str(sID), '.csv'];
    writetable(T,fname);
    
end

end

function [eyePosX, eyePosY, pupilSize, saccadeYN, sf] = readEyeData(fname, sCenter, pixPerDeg)

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

% Save Eyelink-detected saccades:
% <== DO THIS:
% time_sSacc = Data.Events.Esacc.start; % start of saccade
% time_eSacc = Data.Events.Esacc.end; % end of saccade
% elseif length(time_sBlink) < length(time_eBlink)
%     warning('Mismatch in blink indices!')
% end
% eyeData.nBlinks = NaN([1,nTrials]);
% eyeData.nSaccades = NaN([1,nTrials]);
% eyeData.BlinksYN = false([maxN,nTrials]);
% eyeData.SaccadeYN = false([maxN,nTrials]);
% for trial = 1:nTrials
%     time_trial = time_begin(trial):time_end(trial);
%     trial_blinks = ismember(time_sBlink,time_trial);
%     trial_saccades = ismember(time_sSacc,time_trial);
%     eyeData.nBlinks(trial) = sum(trial_blinks);
%     eyeData.nSaccades(trial) = sum(trial_saccades);
%     if eyeData.nBlinks(trial) > 1 % blinks detected
%         aa = find(trial_blinks);
%         for ii = 1:eyeData.nBlinks(trial)
%             idxS = time_sBlink(aa(ii)) - time_begin(trial) + 1;
%             idxE = time_eBlink(aa(ii)) - time_begin(trial) + 1;
%             eyeData.BlinkYN(idxS:idxE,trial) = true;
%             eyeData.BlinkDur(ii,trial) = (1/eyeData.sf) * (idxE - idxS);
%         end
%     end
%     if eyeData.nSaccades(trial) > 1 % saccades detected
%         aa = find(trial_saccades);
%         for ii = 1:eyeData.nSaccades(trial)
%             idxS = time_sSacc(aa(ii)) - time_begin(trial) + 1;
%             idxE = time_eSacc(aa(ii)) - time_begin(trial) + 1;
%             eyeData.SaccadeYN(idxS:idxE,trial) = true;
%             eyeData.SaccadeDur(ii,trial) = (1/eyeData.sf) * (idxE - idxS);
%         end
%     end
% end
saccadeYN = cell([N,1]);

end