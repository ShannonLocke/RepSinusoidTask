% ---------------
% Experiment Text
% ---------------
%
% Created by SML October 2021, Updated 

% Unchanging:
pressSpaceTxt = '\n\n (press space to continue)';
welcomeTxt = ['Welcome to the metacognition for eye-tracking experiment.\n'...
              'In this experiment, you are going to track, with your eyes, a cloud of dots as it follows a random trajectory.\n'...
              'Your goal is to track the center of the dot cloud as closely as possible at all times.\n\n'...
              'Let us look at a few examples now...' pressSpaceTxt];
confInstructTxt = ['Now, you should have an idea of how well you can track these dot clouds.\n'...
                   'We can move on to the main experiment. Hooray!\n'...
                   'In the main experiment, you will have to judge your confidence in your tracking performance after each trial.\n'
                   'You are to report your confidence as follows:\n\n'];
confTxt = ['CONFIDENCE REPORT\n\n Estimate how well you think you can perform eye-tracking in this task.\n'...
           'Now, consider the tracking you just performed, from the first moment to the last.\n'...
           'Overall, do you think your tracking was...\n\n'... '
           'left arrow (<--) = worse than your average    or     right arrow (-->) = better than your average?'];
breakIntroNoticeTxt = ['Short pauses will be offered every ' num2str(breaksEvery) ' trials.\n\n'...
                       'Try to wait until these pauses to take a longer rest because we \n'...
                       'are recording your response times in this experiment.' pressSpaceTxt];
blinkText = ['Also, please note that we are recording your eye movements during this experiment.\n'...
             'Please try to minimise blinking during the trials. It is best to blink when\n'...
             'the calibration ring is presented, before moving on to the next trial.\n\n\'...
             'If you are having a hard time passing this calibration screen, let the experimenter know!' pressSpaceTxt];
breakTxt = 'Take a pause :)\n\n';
                   
% Mode-specific:
switch mode
    case 'demo'
        startTxt = {};
        endTxt = {};
    case 'training'
        startTxt = {welcomeTxt};
        endTxt = {[confInstructTxt, confTxt, pressSpaceTxt]};
    case 'test'
        startTxt = {breakIntroNoticeTxt, blinkText};
        endTxt = {};
end
