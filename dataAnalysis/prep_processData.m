function [] = prep_processData(all_sID)
% This script will run the prepare the processed eye data for the Repeated-
% Sinusoids task (Locke, Goettker, Gegenfurtner, & Mamassian, 2021).
%
% INPUTS:
% 1) all_sID - vector of subject IDs to be extracted
%
% FILES GENERATED (per participant):
% 1) main-task processed eye data file for further processing in Matlab (*.mat)
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
dataFromPath = 'output_data/raw/';
dataToPath = 'output_data/processed/';

%% Process data:

% Load trajectory info:
load('output_data/trajectoryInformation.mat','trajInfo')

for nn = 1:nSs % EACH participant
    
    sID = all_sID(nn); % participants' ID
    disp(['Extracting raw data of s', num2str(sID)])
    
    % Load raw data:
    fname = [dataFromPath, 'rawEyeData_s', num2str(sID), '.mat'];
    load(fname, 'eyeData')
    [nTrials,nSessions] = size(eyeData.eyePosX);
    N = nTrials * nSessions;
    
    % Keep relevant trial info:
    eyeDataPro.sID = eyeData.sID;
    eyeDataPro.session = eyeData.session;
    eyeDataPro.trial = eyeData.trial;
    eyeDataPro.sf = eyeData.sf;
    
    % Preallocate for processed data:
    eyeDataPro.errorEuclid = cell([nTrials,nSessions]);
    eyeDataPro.filtEyeVelocity = cell([nTrials,nSessions]);
    eyeDataPro.noSaccadeEyeVelocity = cell([nTrials,nSessions]);
    eyeDataPro.velocityError = cell([nTrials,nSessions]);
    eyeDataPro.eyeMovClassified = cell([nTrials,nSessions]);
    eyeDataPro.xcorr = cell([nTrials,nSessions]);
    eyeDataPro.lagShiftedError = cell([nTrials,nSessions]);
    eyeDataPro.cleanedPupilSignal = cell([nTrials,nSessions]);
    
    for tt = 1:N % EACH trial
        
        % Compute Euclidean Error signal:
        %     eyeData.errorX = eyeData.eyePosX_deg(1:N,:) - eyeData.targX_int;
        %     eyeData.errorY = eyeData.eyePosY_deg(1:N,:) - 0;
        %     eyeData.errorEuclid = sqrt(eyeData.errorX.^2 + eyeData.errorY.^2);
        %     eyeData.RMSE_X = sqrt(nanmean(eyeData.errorX.^2));
        %     eyeData.RMSE_Euclid = sqrt(nanmean(eyeData.errorEuclid.^2));
        eyeDataPro.errorEuclid{N} = [];
        
        % Compute eye velocity (filtered):
        %     pOrder = 2;
        %     fLength = 51;
        %     dx = denoiseSG(eyeData.eyePosX_deg,eyeData.sf,pOrder,fLength,0);
        %     eyeData.eyeVelX = dx(1:(end-1),:,2);
        %     eyeData.eyeVelX(end,:) = 0;
        %     eyeData.eyeAccX = dx(1:(end-1),:,3);
        %     eyeData.eyeAccX((end-3):end,:) = 0;
        eyeDataPro.filtEyeVelocity{N} = [];
        
        % Compute eye velocity (saccades removed):
        eyeDataPro.noSaccadeEyeVelocity{N} = [];
        
        % Velocity error signal:
        eyeDataPro.velocityError{N} = [];
        
        % Classify samples (blink, fixation, saccade, smooth-pursuit):
        eyeDataPro.eyeMovClassified{N} = [];
        
        % Tracking lag cross-correlation:
        %     % Position:
        %     tdat = eyeData.targX_int;
        %     keep_idx = ~isnan(tdat(:,1));
        %     tdat = tdat(keep_idx,:);
        %     N = sum(keep_idx);
        %     edat = eyeData.eyePosX_deg(1:N,:);
        %     truncateBy = N - 1000;
        %     [r,lag] = slidingCrossCorrelationCoefficient(tdat,edat,truncateBy);
        %     eyeData.lag_pos = 1/eyeData.sf * (lag);
        %     eyeData.xcorr_pos = r;
        %     % Velocity:
        %     tdat = eyeData.targVelX;
        %     keep_idx = ~isnan(tdat(:,1));
        %     tdat = tdat(keep_idx,:);
        %     N = sum(keep_idx);
        %     edat = eyeData.eyeVelX(1:N,:);
        %     truncateBy = N - 1000;
        %     [r,lag] = slidingCrossCorrelationCoefficient(tdat,edat,truncateBy);
        %     eyeData.lag_vel = 1/eyeData.sf * (lag);
        %     eyeData.xcorr_vel = r;
        eyeDataPro.xcorr{N} = [];
        
        % Lag-shifted error signal:
        eyeDataPro.lagShiftedError{N} = [];
        
        % Cleaned pupil signal (z-score?, norm by beginning?):
        eyeDataPro.cleanedPupilSignal{N} = [];
    end
    
    %% Export data for further Matlab processing:
    disp('... Exporting mat file ...')
    fname = [dataToPath, 'processedEyeData_s', num2str(sID), '.mat'];
    save(fname, 'eyeDataPro')
    
end

end