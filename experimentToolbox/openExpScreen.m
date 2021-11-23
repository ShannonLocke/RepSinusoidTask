function [w, hardware] = openExpScreen(skipSyncChecksYN,bgCol)
% This function is standard boilerplate code to open the screen and get
% some important properties.
% 
% Created by SML march 2021

% defaults:
if nargin < 2
    bgCol = 127;
    if nargin < 1
        skipSyncChecksYN = true;
    end
end

% Skip the sync checks if requested:
Screen('Preference', 'SkipSyncTests', double(skipSyncChecksYN));

% Set default settings for psychtoolbox
PsychDefaultSetup(2);

% Find the screen number, selecting the external screen, if available:
screens = Screen('Screens');
hardware.screenID = max(screens);

% Open gray screen window:
[w, hardware.screenRes] = Screen('OpenWindow', hardware.screenID, bgCol);
Priority(MaxPriority(w));
Screen('Flip', w); % do initial flip

% Screen center:
[hardware.sCenter(1), hardware.sCenter(2)] = RectCenter(hardware.screenRes);
hardware.screenRes = hardware.screenRes(3:4);

% Get frame rate:
hardware.fps = Screen('FrameRate',w); % nominal frame rate
for n = 1:10 
hardware.flipInterval(n) = Screen('GetFlipInterval', w); % Sample true flip interval
end
hardware.flipInterval = round(1000*mean(hardware.flipInterval))/1000; % average for final flip interval
if hardware.fps == 0 % no nominal frame rate set
    hardware.fps = round(1/hardware.flipInterval); % Calculate frame rate
end

% Record computer info and software versions:
hardware.computer = computer;
hardware.versionMatlab = version;
hardware.versionPTB = PsychtoolboxVersion;

% Mouse and keyboard:
HideCursor;
KbName('UnifyKeynames'); % unify key names

end
