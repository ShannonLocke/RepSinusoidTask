% ---------------------
% Experiment Parameters
% ---------------------
%
% Created by SML October 2021

% Mode and condition phases;
expData.expDesign.modePhases = {'training','test'}; % experiment mode in each phase
expData.expDesign.condPhases = {'',''}; % here, just block number
expData.expDesign.nPhase = length(expData.expDesign.modePhases); % number of phases
expData.expDesign.phaseTracker = zeros([1,expData.expDesign.nPhase]); % set all to incomplete (0)
expData.expDesign.phaseTimeStamp = {}; % 
expData.expDesign.feedbackYN = {false,false}; % experiment mode in each phase
expData.expDesign.confYN = {false,true}; % experiment mode in each phase
expData.expDesign.nCalibrationRuns = [0 0]; % how many repeats for the blindspot calibration
nTestBlocks = sum(strcmp(expData.expDesign.modePhases,'test')); % number of test blocks

% Number of trajectories & repeats:
demoTraj = 5; 
demoRepeats = 1;
trainingTraj = 5;
trainingRepeats = 1; 
testTraj = 20;
testRepeats = 4; % 10; % test

% Make design matrices:
if experimenterYN
    switch setExperimenterMode
        case 'demo'
            conditions = {[1:demoTraj], [1]};
            expData.expDesign.designMat{1} = makeDesignMat(conditions,demoRepeats);
        case 'training'
            conditions = {[1:trainingTraj], [1]};
            expData.expDesign.designMat{1} = makeDesignMat(conditions,trainingRepeats);
        case 'test'
            conditions = {[1:testTraj], [-1,1]};
            expData.expDesign.designMat{1} = makeDesignMat(conditions,testRepeats/2);
    end
else
    conditions = {[1:trainingTraj], [1]};
    expData.expDesign.designMat{1} = makeDesignMat(conditions,trainingRepeats); % training
    conditions = {[1:testTraj], [-1,1]};
    for block = 1:nTestBlocks
        expData.expDesign.designMat{1+block} = makeDesignMat(conditions,testRepeats/2); % test
    end
end
expData.expDesign.designMat_Legend = {'No conditions'};

% Breaks, calibration, and phase progression:
modesWithBreak = {'test'}; % modes that are long enough for breaks
breaksEvery = 10; % how frequently to offer a short pause
forbidRepeatPhase = {}; % prevent from repeating
calibEvery = 200; % how frequently to calibrate the eye tracker

% Pre-allocate response storage:
expData.res.resp = {}; % confidence response
expData.res.RT = {}; % confidence reaction time