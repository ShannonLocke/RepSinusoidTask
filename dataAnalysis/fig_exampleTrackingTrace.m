%% Get an example participant's data and trajectory info:
load('results/processed/processedEyeData_s555.mat','eyeDataPro')
load('results/trajectoryInformation.mat','trajInfo')
addpath('dataAnalysisToolbox')

%% Select example trial, get info for plotting:
tidx = 250; % trial index
traj = abs(eyeDataPro.trajectory(tidx));
dir = sign(eyeDataPro.trajectory(tidx)); 
t = eyeDataPro.t;
xEye = eyeDataPro.eyeX(:,tidx);
% vEye = abs(eyeDataPro.filtEyeVelocity(:,tidx)); % <== not done!
xTarg = dir * trajInfo.targ_1000Hz(:,traj);
vTarg = abs(trajInfo.targVel_1000Hz(:,traj));

% Smooth velocity signal here:
pOrder = 2;
fLength = 51;
dx = denoiseSG(xEye,1000,pOrder,fLength,0);
vEye = [dx(1:(end-1),:,2); 0];
vEye = abs(vEye);

%% Downsample for faster plotting:
skip = 10;
idx = [1, skip:skip:6000]; % keep first and then every N
t = t(idx);
xEye = xEye(idx);
vEye = vEye(idx);
xTarg = xTarg(idx);
vTarg = vTarg(idx);

%% Preview selected trace:
figure
subplot(2,1,1); hold on
plot(t,xTarg)
plot(t,xEye)
ylabel('Position (deg)')
subplot(2,1,2); hold on
plot(t,vTarg)
plot(t,vEye)
xlabel('Time in Trial (s)')
ylabel('Speed (deg/s)')

%% Export for VSS poster:
T = table(t, xEye, vEye, xTarg, vTarg);
fname = 'data_exampleTrackingTrace.txt';
writetable(T,fname,'Delimiter',' ');