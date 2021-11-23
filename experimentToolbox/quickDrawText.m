function [respKey] = quickDrawText(win,txt,mode,modeOpt,screenOpt,loc,FS,colTxt,colBG)
% QUICKDRAWTEXT a shorthand for getting text to display on a screen. Uses
% the psychtoolbox commands. It is an updated version from quickPrintText
% in that it also collects responses and allows for other elements in the
% display.
%
% [] = quickPrintText(WIN,TXT,LOC,FS,COLTXT,COLBG,PAUSEYN,DISPDUR)
%
% WIN: screen ID.
% TXT: The string of text to be printed.
% MODE: 'keyPress' waits for any key press to clear the screen and will
%        report which key was pressed, and 'fixedTime' will display the
%        text for a fixed period of time specified by 'dispDur'.
% MODEOPT: If mode is 'fixedTime', modeOpt is the wait duration before
%          clearing the text. If mode is 'keyPress', it is a vector of
%          acceptable keypresses (e.g., [20, 44]). If you want to allow any
%          key, use [].
% SCREENOPT: Either 'clean', which sets the screen to clean-screen mode (i.e.,
%            only text presented on blank background. If set to 'busy',
%            other display elements can be drawn on the screen before the
%            text, and the text won't be automatically wiped at the end of
%            this script.
% LOC: location of text as a string (e.g. 'center')
% FS: Font size of text.
% COLTXT: colour of text. Either single value for greyscale, or RGB for
%         colour.
% COLBG: Colour of the background. Either single value for greyscale, or
%        RGB for colour.
%
% Created by SML March 2021

% Defaults
if nargin < 9
    colBG = [127 127 127]; % midgrey
    if nargin < 8
        colTxt = [255 255 255]; % white text
        if nargin < 7
            FS = 20; % 20 pt font size
            if nargin < 6
                loc = {'center', 'center'}; % horizontally and vertically centered text
                if nargin < 5
                    screenOpt = 'clean'; % use simple text screen
                    if nargin < 4
                        modeOpt = []; % fill in according to mode later
                        if nargin < 3
                            mode = 'keyPress'; % wait for any key press
                        end
                    end
                end
            end
        end
    end
end

% Completion if empty:
if isempty(mode); mode = 'keyPress'; end
if isempty(screenOpt); screenOpt = 'clean'; end
if isempty(loc); loc = {'center', 'center'}; end
if isempty(FS); FS = 20; end
if isempty(colTxt); colTxt = [255 255 255]; end
if isempty(colBG); colBG = [127 127 127]; end
if isempty(modeOpt)
    switch mode
        case 'keyPress'
            modeOpt = []; % any key press
        case 'fixedTime'
            modeOpt = 1; % 1 second wait
        otherwise
            warning('Mode input unknown. Defaulting to any key press...');
            mode = 'keyPress';
            modeOpt = [];
    end
end

% Check inputs:
assert((length(colBG)==1)||(length(colBG)==3),'Check that the colBG input is 1 or 3 valued.')
if length(colBG) == 1
    colBG = repmat(colBG,1,3);
end

% Determine if there are other elements:
switch screenOpt
    case 'clean'
        cleanScreenYN = true;
    case 'busy'
        cleanScreenYN = false;
end

% Print text to screen:
Screen('TextSize', win, FS);
if cleanScreenYN; Screen('FillRect', win, colBG); end
DrawFormattedText(win, txt, loc{1}, loc{2}, colTxt);
Screen('Flip', win);

% Pause until keypress or for specified time:
switch mode
    case 'keyPress'
        if isempty(modeOpt) % any key press will be accepted
            pause
        else % only allowed key presses accepted
            keyOK = false;
            while ~keyOK
                respKey = key_resp(-1);
                if any(respKey == modeOpt); keyOK = true; end % key accepted
            end
        end
    case 'fixedTime'
        respKey = [];
        WaitSecs(modeOpt);
end

% Clear screen:
if cleanScreenYN
    Screen('FillRect', win, colBG); % blank screen
    Screen('Flip', win);
end

end