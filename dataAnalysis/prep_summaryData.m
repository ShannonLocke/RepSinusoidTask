function [] = prep_summaryData(all_sID, testDataYN, resDir)
% This script will run the prepare the summary data for the Repeated-
% Sinusoids task (Locke, Goettker, Gegenfurtner, & Mamassian, 2021).
%
% INPUTS:
% 1) all_sID - vector of subject IDs to be extracted
%
% OUTPUTS (whole group):
% 1) main-task summary data file for further processing in Matlab
% 2) training and main-task summary data file for OSF (CSV format)
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
if testDataYN
    dataFromPath_conf = '../data/';
else
    dataFromPath_conf = '../data_pilot/';
end
dataFromPath_eye = [resDir 'processed/'];
dataToPath_matFiles = resDir;
dataToPath_osfFiles = [resDir 'forOSF/'];
tableVarNamesForOSF = {'subjectID', 'session', 'trial', ...
    'directionSignedTrajectoryID', 'confidence', 'RT', 'RMSE', ...
    'RMSEsec1', 'RMSEsec2', 'RMSEsec3', 'RMSEsec4', 'RMSEsec5', ...
        'RMSEsec6'};

%% Collate summary data:

for nn = 1:nSs % EACH participant
    
    sID = all_sID(nn); % participants' ID
    disp(['Extracting data of s', num2str(sID)])
    
    % Get stimulus and confidence data:
    dataPath = [dataFromPath_conf, 's', num2str(sID), '/'];
    dataFilePattern = ['/results_SPC_*.mat'];
    dataOnFile = dir([dataPath, dataFilePattern]);
    fname = [dataPath, '/', dataOnFile.name];
    load(fname, 'expData')
    
    % Get eye-tracking data:
    fname = [dataFromPath_eye, 'processedEyeData_s', num2str(sID), '.mat'];
    load(fname, 'eyeDataPro')
    N = length(eyeDataPro.trajectory);
    nTrials = max(eyeDataPro.trial);
    nSessions = max(eyeDataPro.session);
    
    % Preallocate and/or update:
    if nn == 1 % IS first subject
        summaryData.sID = all_sID;
        summaryData.session = NaN([N,nSs]);
        summaryData.trial = NaN([N,nSs]);
        summaryData.trajectory = NaN([N,nSs]);
        summaryData.conf = NaN([N,nSs]);
        summaryData.RT = NaN([N,nSs]);
        summaryData.RMSE = NaN([N,nSs]);
        summaryData.binnedRMSE = NaN([N,6,nSs]);
    end
    tmpN = length(eyeDataPro.session);
    if tmpN > size(summaryData.session,1) % more trials than budgeted for...
        diffN = tmpN - size(summaryData.session,1); % how many to add
        summaryData.session((end+1):(end+diffN),:) = NaN([diffN,nSs]);
        summaryData.trial((end+1):(end+diffN),:) = NaN([diffN,nSs]);
        summaryData.trajectory((end+1):(end+diffN),:) = NaN([diffN,nSs]);
        summaryData.conf((end+1):(end+diffN),:) = NaN([diffN,nSs]);
        summaryData.RT((end+1):(end+diffN),:) = NaN([diffN,nSs]);
        summaryData.RMSE((end+1):(end+diffN),:) = NaN([diffN,nSs]);
        summaryData.binnedRMSE((end+1):(end+diffN),:,:) = NaN([diffN,6,nSs]);
    end
    summaryData.session(1:N,nn) = eyeDataPro.session;
    summaryData.trial(1:N,nn) = eyeDataPro.trial;
    summaryData.trajectory(1:N,nn) = eyeDataPro.trajectory;
    
    % Extract key info for indexing:
    nTrainTrials = length(expData.expDesign.designMat{1}); % number of training trials
    N_training = nTrainTrials * nSessions; % total number of training trials
    idxR = strcmp(expData.expDesign.modePhases,'training'); % training sessions
    idxT = find(~cellfun(@isempty,expData.res.resp)); % non-empty test sessions

    % Confidence data:
    summaryData.conf(1:N,nn) = cell2mat(expData.res.resp(idxT)');
    summaryData.RT(1:N,nn) = cell2mat(expData.res.RT(idxT)');
    
    % RMSE data:
    for ii = 1:N
        getError = eyeDataPro.errorEuclid(:,ii);
        summaryData.RMSE(ii,nn) = sqrt(mean(getError.^2));
        ns = length(getError);
        tbin = ceil(linspace(0,6,ns));
        tbin(1) = 1; % replace leading 0
        for bb = 1:6 % EACH time bin
            idx = tbin == bb;
            summaryData.binnedRMSE(ii,bb,nn) = sqrt(mean(getError(idx).^2));
        end
    end
    
    % Create table of test eye data for OSF csv:
    tmpT = table(sID*ones([N,1]), ... % subjectID
        summaryData.session(1:N,nn), ... % session
        summaryData.trial(1:N,nn), ... % trial
        summaryData.trajectory(1:N,nn), ... % trajectory signed by direction
        summaryData.conf(1:N,nn), ... % confidence
        summaryData.RT(1:N,nn), ... % RT
        summaryData.RMSE(1:N,nn), ... % RMSE <== not computed because unused
        summaryData.binnedRMSE(1:N,1,nn), ... % RMSE in 1 sec bins
        summaryData.binnedRMSE(1:N,2,nn), ... % ...
        summaryData.binnedRMSE(1:N,3,nn), ... % ...
        summaryData.binnedRMSE(1:N,4,nn), ... % ...
        summaryData.binnedRMSE(1:N,5,nn), ... % ...
        summaryData.binnedRMSE(1:N,6,nn)); % RMSE in 1 sec bins
    tmpT.Properties.VariableNames = tableVarNamesForOSF;
    if nn == 1
        T = tmpT;
    else
        T = [T; tmpT];
    end
    
end

%% Perform data-exclusion checks:

% Trials with RMSE >3SD or <-3SD from the mean:
meanRMSE = repmat(mean(summaryData.RMSE),[N,1]);
sdRMSE = repmat(std(summaryData.RMSE),[N,1]);
summaryData.keepTrialYN = true(size(summaryData.RMSE));
summaryData.keepTrialYN(summaryData.RMSE < (meanRMSE - 3*sdRMSE)) = false;
summaryData.keepTrialYN(summaryData.RMSE > (meanRMSE + 3*sdRMSE)) = false;

% No more than 75% of trials sharing the same confidence response:
confBias = max(mean(summaryData.conf==1), mean(summaryData.conf==-1));
summaryData.passChecksYN = confBias <= 0.75;

%% Export main-task file for further Matlab processing:
disp('... Exporting mat file ...')
fname = [dataToPath_matFiles, 'trialSummaryDataSPC.mat'];
save(fname, 'summaryData')

%% Export all eye data to csv for OSF:
disp('... Exporting csv file ...')
fname = [dataToPath_osfFiles, 'trialSummaryData.csv'];
writetable(T,fname);

end