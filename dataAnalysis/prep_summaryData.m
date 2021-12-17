function [] = prep_summaryData(all_sID)
% This script will run the prepare the summary data for the Repeated-
% Sinusoids task (Locke, Goettker, Gegenfurtner, & Mamassian, 2021).
%
% INPUTS:
% 1) all_sID - vector of subject IDs to be extracted
%
% OUTPUTS (per participant):
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
dataToPath_matFiles = 'output_data/summary/';
dataToPath_osfFiles = 'output_data/forOSF/';
tableVarNamesForOSF = {'session', 'trainingYN', 'trial', ...
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
    
    % Create table of training  data for OSF csv:
    trajectory_training = cell2mat(expData.expDesign.designMat(idxR)');
    trajectory_training = trajectory_training(:,1) .* trajectory_training(:,2);
    tmpT = table(repelem((1:nSessions)', nTrainTrials), ... % session
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
    tmpT = table(summaryData.session(:,nn), ... % session
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
    fname = [dataToPath_matFiles, 'rawData_s', num2str(sID), '.mat'];
    load(fname, 'rawData')
    
    %% Export all eye data to csv for OSF:
    disp('... Exporting csv file ...')
    T.Properties.VariableNames = {'session', 'trainingYN', 'trial', ...
        'directionSignedTrajectoryID', 'confidence', 'RT', 'RMSE', ...
        'RMSEsec1', 'RMSEsec2', 'RMSEsec3', 'RMSEsec4', 'RMSEsec5', ...
        'RMSEsec6'};
    dataPath = [dataToPath_osfFiles, 's', num2str(sID), '/'];
    % ==> Raw trial data (goes into 
    % ==> Raw eye data

end

% 
% %% cleanData.m
% %
% % Created by SML Nov 2021
% 
% % expName = 'replicate'; 
% expName = 'sinusoid';
% 
% % Find participants on file:
% dataPath = ['data_', expName, '/'];
% dataFilePattern = 's*';
% dataOnFile = dir([dataPath, dataFilePattern]);
% nSs = size(dataOnFile,1);
% 
% %% Collate data:
% 
% for nn = 1:nSs
%     
%     groupData.sID(nn) = str2double(dataOnFile(nn).name(2:end));
%     
%     % Load data files:
%     dataPath = ['data_', expName, '/' dataOnFile(nn).name, '/'];
%     fName = ['T1_', dataOnFile(nn).name(2:end), '_eyeData.mat'];
%     load([dataPath, fName],'eyeData')
%     
%     % Store data of interest:
% %     if strcmp(expName, 'sinusoid')
% %         groupData.trajID{nn} = expData.expDesign.designMat{2}(:,1);
% %         groupData.trajDir{nn} = expData.expDesign.designMat{2}(:,2);
% %     end
%     groupData.errorX{nn} = eyeData.errorX;
%     groupData.errorEuclid{nn} = eyeData.errorEuclid;
%     groupData.tSteps_error = eyeData.tSteps_error;
%     groupData.RMSE_X{nn} = eyeData.RMSE_X;
%     groupData.RMSE_Euclid{nn} = eyeData.RMSE_Euclid;
%     groupData.conf{nn} = eyeData.conf;
%     
%     % Compute latency shifted RMSE:
%     nTrials = size(eyeData.xcorr_pos,2);
%     nSamples = size(eyeData.targX_int,1);
%     latShiftRMSE = eyeData.RMSE_X;
%     for tt = 1:nTrials
%         get_xcorr = eyeData.xcorr_pos(:,tt);
%         if ~any(isnan(get_xcorr))
%             latency = eyeData.lag_pos(get_xcorr == max(get_xcorr));
%             nShift = -round(latency * eyeData.sf);
%             if nShift > 0
%                 getErr = eyeData.targX_int(:,tt) - eyeData.eyePosX_deg(nShift:(nSamples+nShift-1),tt);
%                 latShiftRMSE(tt) = sqrt(nanmean(getErr.^2));
%             end
%         end
%     end
%     groupData.latShiftRMSE{nn} = latShiftRMSE;
% end
% 
% %% Save data:
% 
% fName = ['data_', expName, '/groupData.mat']; 
% save(fName,'groupData')