% run_ConfLeakEstExp: The run script for the Confidence-Leak &
% Estimation Experiment.
%
% Created by SML October 2021

% ----------------- %
% EXPERIMENT SET-UP
% ----------------- %

% Set fps manually if on laptop:
if strcmp(testLocation,'laptop')
    hardware.fps = 60;
end

% Get mode-specific parameters:
designMat = expData.expDesign.designMat{phase};
confYN = expData.expDesign.confYN{phase};

% Compute number of trials and stimulus timing:
nTrials = size(designMat,1);
nFrames = round(stimDur*hardware.fps);

% Get random trajectories and dot samples:
switch mode
    case 'demo'
        load('experimentFiles/trajectories_demo.mat')
    case 'training'
        load('experimentFiles/trajectories_training.mat')
    case 'test'
        load('experimentFiles/trajectories_test.mat')
end
sampleTrajectoriesDots;

% Keyboard key codes for responses:
leftKey = KbName('leftarrow');
rightKey = KbName('rightarrow');

% Set up stimulus and response matrices:
if confYN
    conf = NaN([nTrials,1]);
    RT = NaN([nTrials,1]);
end

% -------------- %
% EYELINK SET-UP
% -------------- %

% Initialise and calibrate Eyelink:
if useEyeLinkYN
    EyeLink_initialise;
end

% -------- %
% RUN TASK
% -------- %

% Start instructions:
expText_SPC; % get instructions text
nTextScreens = length(startTxt);
for tt = 1:nTextScreens
    quickDrawText(w,startTxt{tt},'keyPress',spaceKey);
end

% Trial loop:
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % prep for smoothed drawing
for trial = 1:nTrials
    
    % Prepare EyeLink:
    if useEyeLinkYN
        EyeLink_startRecord;
        disp('Recording started...')
    end
    
    % Play stimulus:
    for frame = 1:nFrames
        if halveSampleRate
            idx_dotRow = ceil(frame/2); % present every frame twice
        else
            idx_dotRow = frame; % present every frame once (as normal)
        end
        xvals = squeeze(dotPosX(idx_dotRow,:,trial)); % get dot x position
        yvals = squeeze(dotPosY(idx_dotRow,:,trial)); % get dot y position
        Screen('DrawDots', w, [xvals; yvals], dotSize_pix, dotCol, pos00, dotType);
        Screen('Flip', w);
    end
    
    % Record end of stimulus presentation:
    if useEyeLinkYN
        Eyelink('message','TRIAL_STIMSTOP');
    end
    
    % Blank grey screen:
    Screen('FillRect', w, bgCol)
    Screen('Flip', w); Screen('Flip', w);
    
    % Get confidence:
    if confYN
        
        % Fixation Screen:
        Screen('DrawDots', w, fixPos, fixSize_pix, fixCol, [0 0], 2);
        stimulusFinishTime = Screen('Flip', w); % record when stimulus finshed playing
        % INSERT GAZE CONTROL HERE
        
        % collect response:
        respKey = NaN;
        while ~(respKey == leftKey || respKey == rightKey)
            respKey = key_resp(-1);
        end
        
        % Record timing of response:
        RT(trial) = GetSecs - stimulusFinishTime; % reaction time
        if useEyeLinkYN
            Eyelink('message','TRIAL_RESP');
        end
        
        % save response:
        switch respKey
            case leftKey
                conf(trial) = -1; % left arrow, save as -1 (worse)
            case rightKey
                conf(trial) = 1; % right arrow, save as 1 (better)
        end
        
    end
    
    % Blank grey screen (ISI):
    Screen('FillRect', w, bgCol)
    Screen('Flip', w);
    WaitSecs(ISI/2);
    
    % End recording:
    if useEyeLinkYN
        EyeLink_stopRecord;
    end
    
    % Break screen:
    if any(strcmp(mode,modesWithBreak))
        if (mod(trial,breaksEvery) == 0) && (trial ~= nTrials) % time for a pause
            trialTxt = [num2str(round(100*(trial/nTrials))) '% of trials'];
            makeTxt = [breakTxt 'Your progress: ' trialTxt ', ' phaseTxt pressSpaceTxt];
            quickDrawText(w,makeTxt,'keyPress',spaceKey);
        end
    end
    
end

% --------------------------- %
% SAVE AND FINAL INSTRUCTIONS
% --------------------------- %

% Save results:
expData.hardware.testLocation{resprowidx,phase} = testLocation;
expData.res.tval{resprowidx,phase} = 0:(1/hardware.fps):nFrames;
expData.res.targX{resprowidx,phase} = targ / pixPerDeg;
expData.res.dotPosX{resprowidx,phase} = dotPosX / pixPerDeg;
expData.res.dotPosY{resprowidx,phase} = dotPosY / pixPerDeg;
if confYN
    expData.res.resp{resprowidx,phase} = conf;
    expData.res.RT{resprowidx,phase} = RT;
end

% Close and save eye-tracking data:
if useEyeLinkYN
    EyeLink_shutdown;
end

% Finish instructions:
nTextScreens = length(endTxt);
for tt = 1:nTextScreens
    quickDrawText(w,endTxt{tt},'keyPress',spaceKey);
end
