% Provide Eyelink with details about the graphics environment
% and perform some initializations. The information is returned
% in a structure that also contains useful defaults
% and control codes (e.g. tracker state bit and Eyelink key values).
el=EyelinkInitDefaults(w);

% Initialization of the connection with the Eyelink Gazetracker (exit if this fails):
if ~EyelinkInit(dummymode, 1)
        fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        return;
end

% Save EyeLink tracker version and tracking software version:
[expData.hardware.elTrackerVersion, expData.hardware.elSoftwareVersion] = Eyelink('GetTrackerVersion');

% Provide Eyelink with details about the graphics environment
% and perform some initializations. The information is returned
% in a structure that also contains useful defaults
% and control codes (e.g. tracker state bit and Eyelink key values).
% el=EyelinkInitDefaults(w);
% el.backgroundcolour = bgCol;
% el.foregroundcolour = [255 0 0];
% el.calibrationtargetsize = round(0.5 * pixPerDeg); % inner ring of annulus
% el.calibrationtargetwidth = round(0.6 * pixPerDeg); % outer ring of annulus
% el.targetbeep = 0;
% el.feedbackbeep = 0;
% EyelinkUpdateDefaults(el);

% Data to extract and file name (max length 8):
edfFilePath = cd; % [cd 'data/'];
switch mode
    case 'demo'
        edfFileName = 'D';
    case 'training'
        edfFileName = 'R';
    case 'test'
        edfFileName = ['T' num2str(phase-1) '_'];
end
edfFileName = [edfFileName num2str(participantInfo.sID)];
disp(['Eye-tracking results will be stored in the following file: ' edfFileName '.edf'])
Eyelink('command','file_sample_data=LEFT, RIGHT, GAZE, AREA');
Eyelink('OpenFile', edfFileName);
Eyelink('WaitForModeReady', 1000);

% Perform timing checks:
S1 = GetSecs;
TT = Eyelink('TrackerTime');
S2 = GetSecs;

% Saving timing variables, warn if call-time tolerance exceeded:
expData.TimeInfo.callTimeTolerance = 50; % Microseconds
expData.TimeInfo.initialSystemTime = S1;
expData.TimeInfo.initialTrackerTime = TT;
expData.TimeInfo.trackerTimeOffset =TT-S1;
call_time = (S2 - S1) * 1000 * 1000; % in microseconds
if call_time > expData.TimeInfo.callTimeTolerance
    fprintf('Warning: tracker timing calls took %2.2f, which exceeds the required tolerance of %2.2f microseconds!\n', call_time, expData.TimeInfo.callTimeTolerance);
end

% Update calibration settings:
% x = pos00(1);
% y = pos00(2);
% x_off = round(10 * pixPerDeg); % horizontal offset of calibration dots
% y_off = round(10 * pixPerDeg); % vertical offset of calibration dots
% calib = sprintf('%d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d',...
%     x,y,...
%     x,y-y_off,...
%     x,y+y_off,...
%     x-x_off,y,...
%     x+x_off,y,...
%     x-x_off,y-y_off,...
%     x+x_off,y-y_off,...
%     x-x_off,y+y_off,...
%     x+x_off,y+y_off);
% Eyelink('command','calibration_type = HV9');
% Eyelink('command','generate_default_targets = NO');
% Eyelink('command',sprintf('calibration_targets = %s',calib));
% Eyelink('command',sprintf('validation_targets = %s',calib));
% Eyelink('command','button_function 1 ''accept_target_fixation''');

% Calibrate the eye tracker:
result_el = EyelinkDoTrackerSetup(el);
if result_el == el.TERMINATE_KEY; return; end % abort experiment

% Start Block:
Eyelink('message', 'Block_Start');
Eyelink('WaitForModeReady', 500);

% WHY ARE WE NOT USING THE DEFAULT CALIBRATION TARGETS?