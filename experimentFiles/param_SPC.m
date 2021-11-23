% --------------------------------
% Fixed Visual Stimulus Parameters
% --------------------------------
%
% Created by SML October 2021

% Screen coordinates:
pos00 = hardware.sCenter;

% Dot properties:
nDots = 4; % 2; % how many new dots drawn from the generating distribution
dotSize_deg = 0.25; % dot size in degrees
dotSize_pix = pixPerDeg * dotSize_deg; % dot size in pixels
dotType = 2; % circle with high quality anti-aliasing
dotCol = [255 255 255]; % white
cloudSD_deg = 2; % cloud spread in degrees
cloudSD_pix = pixPerDeg * cloudSD_deg; % cloud spread in degrees

% Stimulus timing (in sec):
stimDur = 6; % stimulus duration
ISI = 0.5; % interstimulus interval in seconds

% Fixation dot parameters:
fixPos = pos00;
fixSize_deg = 0.35;
fixSize_pix = pixPerDeg * fixSize_deg;
fixCol = [180 180 180]; % white 

% Misc display parameters:
bgCol = [127 127 127]; % background colour
lowerTextPos = {'center', 1.5 * pos00(2)}; % placement of text during task
