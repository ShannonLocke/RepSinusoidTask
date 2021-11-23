% ---------------------
% Experiment Parameters
% ---------------------
%
% Created by SML October 2021, updated November 2021

% Mode and condition phases;
expData.expDesign.modePhases = {'training','test','training','test'}; % experiment mode in each phase
expData.expDesign.condPhases = {'','','',''}; % nothing to note here
expData.expDesign.nPhase = length(expData.expDesign.modePhases); % number of phases
expData.expDesign.phaseTracker = zeros([1,expData.expDesign.nPhase]); % set all to incomplete (0)
expData.expDesign.phaseTimeStamp = {}; % 
expData.expDesign.feedbackYN = {false,false,false,false}; % experiment mode in each phase
expData.expDesign.confYN = {false,true,false,true}; % experiment mode in each phase
expData.expDesign.nCalibrationRuns = [0 0 0 0]; % how many repeats for the blindspot calibration

% Number of trajectories & repeats:
demoTraj = 5; 
demoRepeats = 1;
trainingTraj = 5;
trainingRepeats = 1; 
testTraj = 10;
testRepeats = 20;

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
    for ii = 1:expData.expDesign.nPhase
        conditions_training = {[1:trainingTraj], [1]};
        conditions_test = {[1:testTraj], [-1,1]};
        switch expData.expDesign.modePhases{ii}
            case 'training'
                expData.expDesign.designMat{ii} = makeDesignMat(conditions_training, trainingRepeats);
            case 'test'
                expData.expDesign.designMat{ii} = makeDesignMat(conditions_test, testRepeats/2); % test
        end
    end
end
expData.expDesign.designMat_Legend = {'Trajectory ID (1-10)', 'Direction (1: normal, -1: reversed)'};

% Breaks, calibration, and phase progression:
modesWithBreak = {'test'}; % modes that are long enough for breaks
breaksEvery = 10; % how frequently to offer a short pause
forbidRepeatPhase = {'test'}; % prevent from repeating
calibEvery = 250; % how frequently to re-calibrate the eye tracker

% Pre-allocate response storage:
expData.res.resp = {}; % confidence response
expData.res.RT = {}; % confidence reaction time