function [] = prep_summaryData(all_sID)
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
dataFromPath_conf = '../data/';
dataFromPath_eye = 'output_data/processed/';
dataToPath_matFiles = 'output_data/';
dataToPath_osfFiles = 'output_data/forOSF/';
tableVarNamesForOSF = {'subjectID', 'session', 'trainingYN', 'trial', ...
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
    summaryData.session(:,nn) = eyeDataPro.session;
    summaryData.trial(:,nn) = eyeDataPro.trial;
    summaryData.trajectory(:,nn) = eyeDataPro.trajectory;
    
    % Extract key info for indexing:
    nTrainTrials = length(expData.expDesign.designMat{1}); % number of training trials
    N_training = nTrainTrials * nSessions; % total number of training trials
    idxR = strcmp(expData.expDesign.modePhases,'training'); % training sessions
    idxT = find(~cellfun(@isempty,expData.res.resp)); % non-empty test sessions

    % Confidence data:
    summaryData.conf = cell2mat(expData.res.resp(idxT)');
    summaryData.RT = cell2mat(expData.res.RT(idxT)');
    
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
    
    % Create table of training  data for OSF csv:
    trajectory_training = cell2mat(expData.expDesign.designMat(idxR)');
    trajectory_training = trajectory_training(:,1) .* trajectory_training(:,2);
    tmpT = table(sID*ones([N_training,1]), ... % subjectID
        repelem((1:nSessions)', nTrainTrials), ... % session
        ones([N_training,1]), ... % trainingYN
        repmat((1:nTrainTrials)', [nSessions,1]), ... % trial
        trajectory_training, ... % trajectory signed by direction
        NaN([N_training,1]), ... % confidence
        NaN([N_training,1]), ... % RT
        NaN([N_training,1]), ... % RMSE <== not computed because unused
        NaN([N_training,1]), ... % RMSE in 1 sec bins
        NaN([N_training,1]), ... % ...
        NaN([N_training,1]), ... % ...
        NaN([N_training,1]), ... % ...
        NaN([N_training,1]), ... % ...
        NaN([N_training,1])); % RMSE in 1 sec bins
    if nn == 1
        T = tmpT;
        T.Properties.VariableNames = tableVarNamesForOSF;
    else
        T = [T; tmpT];
    end
    
    % Create table of test eye data for OSF csv:
    tmpT = table(sID*ones([N,1]), ... % subjectID
        summaryData.session(:,nn), ... % session
        zeros([N,1]), ... % trainingYN
        summaryData.trial(:,nn), ... % trial
        summaryData.trajectory(:,nn), ... % trajectory signed by direction
        summaryData.conf(:,nn), ... % confidence
        summaryData.RT(:,nn), ... % RT
        summaryData.RMSE(:,nn), ... % RMSE <== not computed because unused
        summaryData.binnedRMSE(:,1,nn), ... % RMSE in 1 sec bins
        summaryData.binnedRMSE(:,2,nn), ... % ...
        summaryData.binnedRMSE(:,3,nn), ... % ...
        summaryData.binnedRMSE(:,4,nn), ... % ...
        summaryData.binnedRMSE(:,5,nn), ... % ...
        summaryData.binnedRMSE(:,6,nn)); % RMSE in 1 sec bins
    tmpT.Properties.VariableNames = tableVarNamesForOSF;
    T = [T; tmpT];
    
end

    %% Export main-task file for further Matlab processing:
    disp('... Exporting mat file ...')
    fname = [dataToPath_matFiles, 'trialSummaryDataSPC.mat'];
    save(fname, 'summaryData')
    
    %% Export all eye data to csv for OSF:
    disp('... Exporting csv file ...')
    fname = [dataToPath_osfFiles, 'trialSummaryData.csv'];
    writetable(T,fname);

end