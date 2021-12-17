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

%% Process data:

% Load trajectory info:
load('output_data/trajectoryInformation.mat','trajInfo')
nsT = length(trajInfo.t_1000Hz);

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
    eyeDataPro.errorEuclid = cell([nTrials,nSessions]);
    eyeDataPro.filtEyeVelocity = cell([nTrials,nSessions]);
    eyeDataPro.noSaccadeEyeVelocity = cell([nTrials,nSessions]);
    eyeDataPro.velocityError = cell([nTrials,nSessions]);
    eyeDataPro.eyeMovClassified = cell([nTrials,nSessions]);
    eyeDataPro.xcorr = cell([nTrials,nSessions]);
    eyeDataPro.lagShiftedError = cell([nTrials,nSessions]);
    eyeDataPro.cleanedPupilSignal = cell([nTrials,nSessions]);
    
    for tt = 1:N % EACH trial
        
        % Strech trajctory signals to match number of samples:
        tID = abs(eyeData.trajectory(tt)); % trajectory number
        tDir = sign(eyeData.trajectory(tt)); % trajectory direction
        nsE = eyeData.nSamples(tt); % number of samples, eye data
        targPosX = interp1(1:nsT, trajInfo.targ_1000Hz(:,tID), 1:nsE)';
        idx1 = find(diff(isnan(targPosX))==-1); % last beginning NaN
        idx2 = find(diff(isnan(targPosX))==1)+1; % first ending NaN
        targPosX(1:idx1) = trajInfo.targ_1000Hz(1,tID);  % replace initial NaNs with first value
        targPosX(idx2:end) = trajInfo.targ_1000Hz(end,tID); % replace final NaNs with last value
        targPosX = tDir * targPosX; % account for direction
        if plotYN && (mod(tt,plotEvery)==0) % visualisation check
            figure; hold on
            plot(tDir * trajInfo.targ_1000Hz(:,tID),'-')
            plot(linspace(1,nsT,nsE),targPosX,'o')
        end
        
        % Interpolate across blinks:
        blinkYN = isnan(eX{tt}) | isnan(eY{tt});
        if any(blinkYN) % IS blinks, do interpolation
            nBlinks = sum(diff(blinkYN)==1);
            bidx = [find(diff(blinkYN)==-1), find(diff(blinkYN)==1)];
            for bb = 1:nBlinks
                % <== DO THIS...
            end
        end
        
        % Compute Euclidean Error signal:
        errorX = eyeData.eyePosX{tt} - targPosX;
        errorY = eyeData.eyePosY{tt} - 0;
        eyeDataPro.errorEuclid{tt} = sqrt(errorX.^2 + errorY.^2);
        if plotYN && (mod(tt,plotEvery)==0) % visualisation check
            figure; hold on
            plot(errorX); plot(errorY)
            plot(eyeDataPro.errorEuclid{tt})
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
        eyeDataPro.filtEyeVelocity{tt} = [];
        
        % Compute eye velocity (saccades removed):
        % <== DO THIS...
        eyeDataPro.noSaccadeEyeVelocity{tt} = [];
        
        % Velocity error signal:
        % <== DO THIS...
        eyeDataPro.velocityError{tt} = [];
        
        % Classify samples (blink, fixation, saccade, smooth-pursuit):
        % <== DO THIS...
        eyeDataPro.eyeMovClassified{tt} = [];
        
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
        eyeDataPro.xcorr{tt} = [];
        
        % Lag-shifted error signal:
        % <== DO THIS...
        eyeDataPro.lagShiftedError{tt} = [];
        
        % Cleaned pupil signal (z-score?, norm by beginning?):
        % <== DO THIS...
        eyeDataPro.cleanedPupilSignal{tt} = [];
    end
    
    %% Export data for further Matlab processing:
    disp('... Exporting mat file ...')
    fname = [dataToPath, 'processedEyeData_s', num2str(sID), '.mat'];
    save(fname, 'eyeDataPro')
    
end

end