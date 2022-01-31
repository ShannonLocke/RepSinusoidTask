function [] = prep_trajectories(resDir)
% This script prepares the target trajectory information for archiving and
% further analysis. Specifically, interpolation of the target trajectory to 
% 1000 Hz sampling rate.
% 
% Repeated-Sinusoids task (Locke, Goettker, Gegenfurtner, & Mamassian, 2021).
%
% FILES GENERATED:
% 1) trajectory info for further processing in Matlab (trajectoryInfo.mat)
% 2) sinusoid components of the trajectories for OSF (trajectoryComponents.csv)
% 3) trajectory of targets for OSF (trajectoryTimeCourse.csv)
%
% Data: https://osf.io/j9asn/
% Code: https://github.com/ShannonLocke/RepSinusoidTask
%
% Created by SML Dec 2021
% License: CC-By Attribution 4.0 International

%% Preamble:

% Directory and file names:
dataFromPath = '../experimentFiles/';
dataToPath_matFiles = resDir;
dataToPath_osfFiles = [resDir 'forOSF/'];
fname = {'trajectories_training', 'trajectories_test'};
trainingYN = [true false];

%% Process trajectories:

for ii = 1:2

% Load trajectory data:
load([dataFromPath, fname{ii}], 'ampInfo', 'freqInfo', 't', 'targ', 'targVel')
[nSin,nTraj] = size(freqInfo);
nt = length(t);

% Get general info about each trajectory:
maxDev = max(abs(targ))';
tmpT_comp = table(repmat(trainingYN(ii),[nSin*nTraj,1]), ... % TraingingYN
    repelem((1:nTraj)',nSin), ... % trajectoryID
    freqInfo(:), ... % frequency (cycles per sec)
    ampInfo(:), ... % amplitude
    repelem(maxDev,nSin)); % maximum deviation
tmpT_comp.Properties.VariableNames = {'trainingYN', 'trajectoryID', ...
    'frequency', 'amplitude', 'maximumDeviationDeg'};

% Get specific trajectories:
% % trainingYN, trajectoryID, t, targPos
tmpT_traj = table(repmat(trainingYN(ii),[nt*nTraj,1]), ... % TraingingYN
    repelem((1:nTraj)',nt), ... % trajectoryID
    repmat(round(t,4),[nTraj,1]), ... % time in trial
    round(targ(:),4)); % target position
tmpT_traj.Properties.VariableNames = {'trainingYN', 'trajectoryID', ...
    'timeInTrialSec', 'TargPosDeg'};

% Update storage tables of OSF data:
if ii == 1
    T_comp = tmpT_comp;
    T_traj = tmpT_traj;
else
    T_comp = [T_comp; tmpT_comp];
    T_traj = [T_traj; tmpT_traj];
end

end

%% Get trajectory info for further processing (test trials only):

% Store existing info:
trajInfo.freqInfo = freqInfo;
trajInfo.ampInfo = ampInfo;
trajInfo.maxDev = maxDev;
trajInfo.t_60Hz = t+ (1/60)/2; % centre on middle of frame presentation time window
trajInfo.targ_60Hz = targ;
trajInfo.targVel_60Hz = targVel;

% Linear interpolation of trajectories for 1000Hz sampling rate:
trajInfo.t_1000Hz = (0:(1/1000):(6-1/1000))'; % 6 sec stimulus duration
nt = length(trajInfo.t_1000Hz);
trajInfo.targ_1000Hz = NaN([nt,nTraj]);
trajInfo.targVel_1000Hz = NaN([nt,nTraj]);
for tt = 1:nTraj
    intP = interp1(trajInfo.t_60Hz, targ(:,tt), trajInfo.t_1000Hz); % interpolate target position
    intV = interp1(trajInfo.t_60Hz, targVel(:,tt), trajInfo.t_1000Hz); % interpolate target velocity
    if tt == 1
        idx1 = find(diff(isnan(intP))==-1); % last beginning NaN
        idx2 = find(diff(isnan(intP))==1)+1; % first ending NaN
    end
    intP(1:idx1) = targ(1,tt);  % replace initial NaNs with first value
    intV(1:idx1) = targVel(1,tt);
    intP(idx2:end) = targ(end,tt); % replace final NaNs with last value
    intV(idx2:end) = targVel(end,tt);
    trajInfo.targ_1000Hz(:,tt) = intP; % store
    trajInfo.targVel_1000Hz(:,tt) = intV;
    
end

%% Save trajectory information:

% Main-task file for further Matlab processing:
fname = [dataToPath_matFiles, 'trajectoryInformation.mat'];
save(fname, 'trajInfo')

% Save .csv files for OSF:
fname = [dataToPath_osfFiles, 'trajectoryComponents.csv'];
writetable(T_comp,fname);
fname = [dataToPath_osfFiles, 'trajectoryTimeCourse.csv'];
writetable(T_traj,fname);

end