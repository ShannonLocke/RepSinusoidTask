% Dot position randomly draw from generating distribution:
% --> (rows = frame, columns = dot)
temp_fps = 60;
temp_nFrames = stimDur * temp_fps;
targ = pixPerDeg * targ(:,designMat(:,1)); % Select trajectories from set, convert to pixels
targ = targ .* repmat(designMat(:,2)',[size(targ,1),1]); % Flip or not
mu_vals = repmat(reshape(targ,[temp_nFrames,1,nTrials]),[1,nDots,1]);
sig_vals = repmat(cloudSD_pix,[temp_nFrames,nDots,nTrials]);
dotPosX = normrnd(mu_vals, sig_vals); % Draw x locations
dotPosY = normrnd(0, sig_vals); % Draw y locations