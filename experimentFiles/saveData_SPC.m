% --------- %
% SAVE DATA
% --------- %
%
% Created by SML November 2021

% Save testing information:
expData.hardware.testLocation{resprowidx,phase} = testLocation;
expData.hardware.randomSeed{resprowidx,phase} = randomSeed;

% Save results:
expData.res.tval{resprowidx,phase} = 0:(1/hardware.fps):nFrames;
expData.res.targX{resprowidx,phase} = targ / pixPerDeg;
expData.res.dotPosX{resprowidx,phase} = dotPosX / pixPerDeg;
expData.res.dotPosY{resprowidx,phase} = dotPosY / pixPerDeg;
if confYN
    expData.res.resp{resprowidx,phase} = conf;
    expData.res.RT{resprowidx,phase} = RT;
end

% Close and save eye-tracking data:
if useEyeLinkYN
    EyeLink_shutdown;
end