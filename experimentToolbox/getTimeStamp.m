function [tStamp] = getTimeStamp
%GETTIMESTAMP Returns the current time in the format: yyyy-mm-dd_HH-MM.

tStamp = datestr(now, 31); % Current date and time in yyyy-mm-dd HH:MM:SS
tStamp = tStamp(1:end-3); % Chop seconds off
tStamp(11) = '_'; % Separate date and clock time
tStamp(14) = '-'; % Replace colon with dash

end

