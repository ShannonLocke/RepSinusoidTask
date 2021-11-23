% ---------------
% Experiment Text
% ---------------
%
% Created by SML October 2021

% Unchanging:
pressSpaceTxt = '\n\n (press space to continue)';
confTxt = 'Report CONFIDENCE: left arrow = worse than average, right arrow = better than average.';
breakIntroNoticeTxt = ['Short pauses will be offered every ' num2str(breaksEvery) ' trials.\n\n'...
                       'Try to wait until these pauses to take a longer rest because we \n'...
                       'are recording your response times in this experiment.' pressSpaceTxt];
blinkText = ['We are recording your eye movements during this experiment. Please try to \n'...
             'minimise blinking during the trials. It is best to blink when the calibration ring \n'...
             'is presented, before moving on to the next trial.'];
breakTxt = 'Take a pause :)\n\n';
                   
% Mode-specific:
switch mode
    case 'demo'
        startTxt = {};
        endTxt = {};
    case 'training'
        startTxt = {};
        endTxt = {};
    case 'test'
        startTxt = {breakIntroNoticeTxt, blinkText};
        endTxt = {};
end
