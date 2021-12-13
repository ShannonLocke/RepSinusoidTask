% START THE EXPERIMENT!
%
% Created by SML October 2021

try
           
    % Switchboard:
    experimenterYN = false; % toggle for a debug mode (true) or final verion mode (false)
    setExperimenterMode = 'demo'; % choose which mode to default to for the experimenter
    testLocation = 'laptop'; % Options: 'laptop', 'debugRoom', 'testRoom'
    switch testLocation
        case 'laptop'
            useEyeLinkYN = false; % use the Eyelink or not
            pixPerDeg = 40; % 54; %  27; % pixels per degree
            halveSampleRate = false; % Default is 60 Hz display, halve if using 120 Hz display
        case 'debugRoom'
            useEyeLinkYN = true; % use the Eyelink or not
            pixPerDeg = 43; % pixels per degree
            halveSampleRate = true; % Default is 60 Hz display, halve if using 120 Hz display
        case 'testRoom'
            useEyeLinkYN = true; % use the Eyelink or not
            pixPerDeg = 80; % pixels per degree
            halveSampleRate = false; % Default is 60 Hz display, halve if using 120 Hz display
    end
    dummymode = 0; % Actually use the eye-tracker (0), or use the mouse instead (1)
    
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
    
    % Clean up workspace:
    clear;
    
catch e % error!!! Execute emergency shutdown...
    
    sca;
    ListenChar(0);
    try; Eyelink('Shutdown'); end
    msg_PTB = ['Something went wrong. Please contact the experimenter. '... 
               'If possible, screenshot the error and send it to them.'];
    msgbox(msg_PTB,'','error');
    rethrow(e)
    
end
