% START THE EXPERIMENT!
%
% Created by SML October 2021

try
           
    % Switchboard:
    experimenterYN = true; % toggle for a debug mode (true) or final verion mode (false)
    setExperimenterMode = 'demo'; % choose which mode to default to for the experimenter
    testLocation = 'laptop'; % Options: 'laptop', 'debugRoom', 'testRoom'
    switch testLocation
        case 'laptop'
            useEyeLinkYN = false; % use the Eyelink or not
            preset_pixPerDeg = 40; % 54; %  27; % pixels per degree
            halveSampleRate = false; % Default is 60 Hz display, halve if using 120 Hz display
        case 'debugRoom'
            useEyeLinkYN = true; % use the Eyelink or not
            preset_pixPerDeg = 43; % pixels per degree
            halveSampleRate = true; % Default is 60 Hz display, halve if using 120 Hz display
        case 'testRoom'
            useEyeLinkYN = true; % use the Eyelink or not
            preset_pixPerDeg = 80; % pixels per degree
            halveSampleRate = false; % Default is 60 Hz display, halve if using 120 Hz display
    end
    dummymode = 0; % Actually use the eye-tracker (0), or use the mouse instead (1)
    blindspotCalibrationYN = false; % toggle if you want to calibrate stimuli using blindspot calibration technique
    
    % Experiment and participant info:
    expName = 'SPC'; % experiment name to append to files
    pic_fieldNames = []; % use default participant info card field names
    pic_prompt = []; % use default participant info card prompts
    
    % Files and directories:
    addpath(genpath('experimentFiles/'))
    addpath(genpath('experimentToolbox/'))
    expDesignScript = @expDesign_SPC; % experimental design script
    paramScript = @param_SPC; % parameter specification script
    taskScript = @task_SPC; % task run script
    
    % Go!:
    generalExpControlScript;
    
    % File reorganisation (place .EDF files in subject data folder):
    dataFilePattern = ['*', num2str(participantInfo.sID), '.edf'];
    dataOnFile = dir(dataFilePattern);
    moveFilePath = ['data/s', num2str(participantInfo.sID), '/'];
    if ~isempty(dataOnFile) % data files found
        for ii = 1:length(dataOnFile) % EACH file
            movefile(dataOnFile(ii).name, moveFilePath);
        end
    end
    
catch e % error!!! Execute emergency shutdown...
    
    if useEyeLinkYN; Eyelink('Shutdown'); end
    sca;
    ListenChar(0);
    msg_PTB = ['Something went wrong. Please contact the experimenter. '... 
               'If possible, screenshot the error and send it to them.'];
    msgbox(msg_PTB,'','error');
    rethrow(e)
    
end