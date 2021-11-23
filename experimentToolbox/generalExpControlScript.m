% ----------------------------------------------- %
% GENERIC EXPERIMENT CONTROL SCRIPT - DO NOT EDIT!
% ----------------------------------------------- %
%
% Created by SML March 2021

% Confirm the participant has Psychtoolbox installed:
try
    versionPTB = PsychtoolboxVersion;
catch
    msg_PTB = ['You do not have Psychtoolbox installed or it is inaccessible. '...
        'The download files & instructions can be found here: http://psychtoolbox.org/download. '...
        'If you do not wish to install Psychtoolbox or are having any difficulties, '...
        'please contact the experimenter.'];
    ff = msgbox(msg_PTB,'','error');
    uiwait(ff);
    return
end

% Get participant info:
if experimenterYN % debug mode
    saveFile = ['results_' expName '_TESTRUN_' getTimeStamp '.mat'];
    participantInfo.sID = 999;
else % acutal experiment run, load subject
    [participantInfo, saveFile] = participantInfoCard(expName,pic_fieldNames,pic_prompt);
end

% Open screen if not already open with Sync-checks if possible:
try
    skipSyncChecksYN = false;
    [w, hardware] = openExpScreen(skipSyncChecksYN);
catch
    skipSyncChecksYN = true;
    [w, hardware] = openExpScreen(skipSyncChecksYN);
    warning('Sync-Checks had to be skipped.')
end
hardware.skipSyncChecksYN = skipSyncChecksYN;

% Keyboard set-up: silence input and find necessary key codes:
if ~experimenterYN; ListenChar(2); end % only switch off keyboard for participants
spaceKey = KbName('SPACE');
rKey = KbName('r');
qKey = KbName('q');

% Load existing results file or create new, updating for any new setups:
expDesignScript(); % get experiment info
try
    % Retrieve save file from previous testing session:
    load(saveFile,'expData');
    
    % Compare hardware to any previous versions:
    nSetups = length(expData.hardware);
    ff = fieldnames(hardware);
    for nn = 1:nSetups
        hardwareSameYN = true;
        for ii = 1:length(ff) % check each field for match against setup #nn
            if hardware.(ff{ii}) ~= expData.hardware(nn).(ff{ii})
                hardwareSameYN = false;
                break
            end
        end
        if hardwareSameYN; setupTracker = nn; break; end % match found, record setup number, move on
    end
    if ~hardwareSameYN % no match found
        expData.hardware(nSetups+1) = hardware;
        setupTracker = nSetups+1;
    end
    
catch
    % Create new save file:
    setupTracker = 1;
    expData.hardware = hardware; % add hardware information
    
    % Pre-allocate space for calibration results:
    if blindspotCalibrationYN % use blindspot calibration
        expData.calibration.pixPerDeg = NaN([1,expData.expDesign.nPhase]);
        expData.calibration.distToBlindspot = cell([1,expData.expDesign.nPhase]);
        expData.calibration.setup = NaN([1,expData.expDesign.nPhase]);
    else % use preset pixPerDeg value
        expData.calibration.pixPerDeg = preset_pixPerDeg;
        pixPerDeg = preset_pixPerDeg;
    end
    
    % Save the initial participant file:
    save(saveFile,'expData');
    
end

% Progress through experiment phases:
phase = find(expData.expDesign.phaseTracker==0, 1); % find next phase
continueTestingYN = true;
while continueTestingYN
    
    % Determine mode, condition:
    if experimenterYN % experimenter, enter mode (default: demo)
        if hardware.screenID == 0 % only one screen, default to demo
            mode = setExperimenterMode;
        else % using second screen, ask for mode
            phaseOK = false;
            while ~phaseOK
                mode = inputdlg('Mode:','',[1,10],{'demo'});
                mode = mode{1};
                if strcmp(mode,'demo') % demo setting
                    phaseOK = true;
                else  % not demo, find the correct phase number
                    phase = find(strcmp(mode,expData.expDesign.modePhases));
                    if isempty(phase) % mode input NOT recognised
                        warning('Mode not recognised. Try again.')
                    else
                        phaseOK = true;
                    end 
                end
            end
        end
        condition = '';
        expData.expDesign.nPhase = 1; % only do one phase
    else % participant, determine from file
        mode = expData.expDesign.modePhases{phase}; % find next mode
        condition = expData.expDesign.condPhases{phase}; % find next condition
    end
    
    % Announce beginning of phase:
    phaseTxt = ['phase ' num2str(phase) ' of ' num2str(expData.expDesign.nPhase) ' (' upper(mode) ')'];
    makeTxt= ['Starting ' phaseTxt '\n\n (press space to continue)'];
    quickDrawText(w,makeTxt,'keyPress',spaceKey);
    
    % Perform blindspot calibration and record set-up:
    if blindspotCalibrationYN % use blindspot calibration
        if strcmp(mode,'demo')
            nCalibrationRuns = 3;
        else
            nCalibrationRuns = expData.expDesign.nCalibrationRuns(phase);
        end
        [pixPerDeg, distToBlindspot] = blindspotCalibration(nCalibrationRuns,w,hardware);
        if ~experimenterYN
            expData.calibration.pixPerDeg(phase) = pixPerDeg;
            expData.calibration.distToBlindspot{phase} = distToBlindspot;
            expData.calibration.setup(phase) = setupTracker;
        end
    end
    
    % Find row index to save results:
    [nRow,nCol] = size(expData.res.resp); % get dimensions of results
    if nCol < phase % no entries for this phase
        resprowidx = 1; % make first entry
    else  % no entries for this phase
        resprowidx = nRow + 1; % default to adding a new row
        for ii = 1:nRow
            if isempty(expData.res.resp{ii,phase}) % empty spot found
                resprowidx = ii; % update row index to empty spot
                break
            end
        end
    end
    
    % Run next step in experiment:
    paramScript();
    taskScript();
    
    % Update the save file with the new results:
    expData.expDesign.phaseTracker(phase) = 1;
    expData.expDesign.phaseTimeStamp{end+1,phase} = getTimeStamp;
    save(saveFile,'expData');
    
    % Ask participant what to do next:
    if phase >= expData.expDesign.nPhase % end of experiment
        
        makeTxt = ['Experiment complete!\n\n Thank you for your participation :)\n\n'...
            'Please send all the files in the data folder to the experimenter.\n\n '...
            '(press space to close the screen)'];
        quickDrawText(w,makeTxt,'keyPress',spaceKey);
        if ~experimenterYN
            participantInfoCard(expName,pic_fieldNames,pic_prompt,participantInfo); % remove subject ID from participant info card
        end
        continueTestingYN = false;
        
    else % not complete, does the participant want continue, repeat, or exit?
        
        % Get participant's decision:
        if any(strcmp(mode,forbidRepeatPhase)) % not allowed to repeat this phase
            optionsTxt = ' Press space bar to continue or "q" to exit this session.';
            allowedKeys = [spaceKey qKey];
        else % repeats allowed
            optionsTxt = ' Press space bar to continue, "r" to repeat this phase,\n or "q" to exit this session.';
            allowedKeys = [spaceKey rKey qKey];
        end
        makeTxt = ['You have completed ' phaseTxt '\n\n' optionsTxt];
        respKey = quickDrawText(w,makeTxt,'keyPress',allowedKeys);
        
        % Prep for what's next:
        if respKey == spaceKey % CONTINUE experiment
            phase = phase + 1; % increment phase tracker
        elseif respKey == rKey % REPEAT this phase
            disp('You have selected to repeat this phase, yay!');
        elseif respKey == qKey % QUIT for now
            continueTestingYN = false; % exit while loop
        end
        
    end
    
end

% Shutdown:
sca;
ListenChar(0);
clear;