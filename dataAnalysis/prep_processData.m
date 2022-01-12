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

% Visualisation settings:
plotYN = true;
plotEvery = 100; 

% Load trajectory info:
load('output_data/trajectoryInformation.mat','trajInfo')
ns = length(trajInfo.t_1000Hz);

%% Process data:

for nn = 1:nSs % EACH participant
    
    sID = all_sID(nn); % participants' ID
    disp(['Extracting raw data of s', num2str(sID)])
    
    % Load raw data:
    fname = [dataFromPath, 'rawEyeData_s', num2str(sID), '.mat'];
    load(fname, 'eyeData')
    eX = eyeData.eyePosX;
    eY = eyeData.eyePosY;
    [nTrials,nSessions] = size(eyeData.eyePosX);
    N = nTrials * nSessions;
    
    % Keep relevant trial info:
    eyeDataPro.sID = eyeData.sID;
    eyeDataPro.session = eyeData.session;
    eyeDataPro.trial = eyeData.trial;
    eyeDataPro.trajectory = eyeData.trajectory;
    eyeDataPro.nSamples = eyeData.nSamples;
    eyeDataPro.sf = eyeData.sf;
    
    % Preallocate for processed data:
    eyeDataPro.t = trajInfo.t_1000Hz;
    eyeDataPro.eyeX = NaN([ns,nTrials]);
    eyeDataPro.eyeY = NaN([ns,nTrials]);
    eyeDataPro.errorEuclid = NaN([ns,nTrials]); 
    eyeDataPro.filtEyeVelocity = NaN([ns,nTrials]);
    eyeDataPro.noSaccadeEyeVelocity = NaN([ns,nTrials]);
    eyeDataPro.velocityError = NaN([ns,nTrials]);
    eyeDataPro.eyeMovClassified = NaN([ns,nTrials]);
    eyeDataPro.xcorr = NaN([ns,nTrials]);
    eyeDataPro.lagShiftedError = NaN([ns,nTrials]);
    eyeDataPro.cleanedPupilSignal = NaN([ns,nTrials]);
    
    for tt = 1:N % EACH trial
        
        % Find target info:
        tID = abs(eyeData.trajectory(tt)); % trajectory number
        tDir = sign(eyeData.trajectory(tt)); % trajectory direction
        targPosX = tDir * trajInfo.targ_1000Hz(:,tID); % target position
        
        % Interpolate eye data signals for identical sample number:
        nsE = eyeData.nSamples(tt); % number of samples, eye data
        t_eye = linspace(0,6,nsE); % stretched time samples
        idx = ~isnan(eX{tt}); % find non-blink samples
        get_eX = interp1(t_eye(idx), eX{tt}(idx), trajInfo.t_1000Hz);
        get_eY = interp1(t_eye(idx), eY{tt}(idx), trajInfo.t_1000Hz);
        if isnan(get_eX(1)) || isnan(get_eY(1)); disp(['Trial' num2str(tt)]); warning('Missing beginning values.'); end
        if isnan(get_eX(end)) || isnan(get_eY(end)); disp(['Trial' num2str(tt)]); warning('Missing end values.'); end
        eyeDataPro.eyeX(:,tt) = get_eX;
        eyeDataPro.eyeY(:,tt) = get_eY;
        if plotYN && (mod(tt,plotEvery)==0) % visualisation check
            figure; hold on
            plot(t_eye, eX{tt},'-')
            plot(trajInfo.t_1000Hz,get_eX,'o')
        end
        
        % Compute Euclidean Error signal:
        errorX = get_eX - targPosX;
        errorY = get_eY - 0;
        eyeDataPro.errorEuclid(:,tt) = sqrt(errorX.^2 + errorY.^2);
        if plotYN && (mod(tt,plotEvery)==0) % visualisation check
            figure; hold on
            plot(errorX); plot(errorY)
            plot(eyeDataPro.errorEuclid(:,tt))
        end
        
        % Compute eye velocity (filtered):
        % <== DO THIS...
        %     pOrder = 2;
        %     fLength = 51;
        %     dx = denoiseSG(eyeData.eyePosX_deg,eyeData.sf,pOrder,fLength,0);
        %     eyeData.eyeVelX = dx(1:(end-1),:,2);
        %     eyeData.eyeVelX(end,:) = 0;
        %     eyeData.eyeAccX = dx(1:(end-1),:,3);
        %     eyeData.eyeAccX((end-3):end,:) = 0;
        % eyeDataPro.filtEyeVelocity(:,tt) = [];
        
        % Compute eye velocity (saccades removed):
        % <== DO THIS...
        % eyeDataPro.noSaccadeEyeVelocity(:,tt) = [];
        
        % Velocity error signal:
        % <== DO THIS...
        % eyeDataPro.velocityError(:,tt) = [];
        
        % Classify samples (blink, fixation, saccade, smooth-pursuit):
        % <== DO THIS...
        % eyeDataPro.eyeMovClassified(:,tt) = [];
        
        % Tracking lag cross-correlation:
        % <== DO THIS...
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
        % eyeDataPro.xcorr(:,tt) = [];
        
        % Lag-shifted error signal:
        % <== DO THIS...
        % eyeDataPro.lagShiftedError(:,tt) = [];
        
        % Cleaned pupil signal (z-score?, norm by beginning?):
        % <== DO THIS...
        % eyeDataPro.cleanedPupilSignal(:,tt) = [];
    end
    
    %% Export data for further Matlab processing:
    disp('... Exporting mat file ...')
    fname = [dataToPath, 'processedEyeData_s', num2str(sID), '.mat'];
    save(fname, 'eyeDataPro')
    
end

end