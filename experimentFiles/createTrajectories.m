%
% Created by SML Oct 2020, edited for EyeMovementExp Oct 2021

% Setup details:
maxDev_mean = 12;
maxDev_sd = 1;
fps = 60;

% Trajectory parameters:
N = 20;
runDur = 6;
fundFreq = 0.05; % 1/(2*runDur); % fundamental frequency, will complete 1/2 period over the course of the run
harmFreq = fundFreq * (2:8); % select harmonic frequencies
nFreq = 5; % number of frequencies to sample, excluding the fundatmental frequency
maxSpeed = 30; % maximum speed of target, in degrees/sec
minFracVelSP = 0.75; % minimum fraction of time above 2 deg/s smooth pursuit threshold

% Pre-allocate:
t = 0 : (1/fps) : runDur-(1/fps); t = t'; % Time vector (in sec)
nt = length(t); % number of time steps
freqInfo = NaN([nFreq+1,N]);
ampInfo = NaN([nFreq+1,N]);
targ = zeros([nt,N]);
targSpeed = zeros([nt,N]);
targVel = zeros([nt,N]);

%%  Sample trajectories:
for nn = 1:N % EACH trajectory
    
    % Sample maximum deviation:
    maxDev = maxDev_mean + maxDev_sd * randn;
    
    % Get suitable trajectory:
    trajectoryOK = false;
    while ~trajectoryOK % Sample until required conditions are met
        
        % Draw frequency components and randomise amplitude:
        fidx = randperm(length(harmFreq),nFreq);
        xFreq = sort([fundFreq, harmFreq(fidx)],'ascend'); % sort, lowest to highest
        nF = nFreq + 1;
        xAmp = 2*rand([1,nF]) - 1;
        [~,aidx] = sort(abs(xAmp),'descend'); % sort by magnitude, highest to lowest
        xAmp = xAmp(aidx);% ./(1:nF);
        
        % Create time series:
        velX = zeros(length(t),1);
        for ii=1:nF
            targ(:,nn) = targ(:,nn) + xAmp(ii) * sin(2*pi*xFreq(ii)*t);
            % velX = velX + xAmp(ii) * 2*pi*xFreq(ii) * cos(2*pi*xFreq(ii)*t);
        end
        maxSignalVal = max(abs(targ(:,nn))); % most extreme peak/trough
        targ(:,nn) = targ(:,nn) / maxSignalVal; % normalise signal to (-1,1)
        targ(:,nn) = maxDev * targ(:,nn); % scale to tajectory limits, target horizontal position now in pixels
        
        % velX = velX * (maxDev / maxSignalVal); % apply same scaling to velocity, now in pixels/sec
        targVel(:,nn) = [0; diff(targ(:,nn))]/(1/fps);
        
        % Check speed (critical factor):
        targSpeed(:,nn) = abs(targVel(:,nn)); % target speed
        speedCheckOK_1 = all(targSpeed(:,nn)<=maxSpeed);
        speedCheckOK_2 = all(mean(targSpeed(:,nn)>=5) > minFracVelSP);
        if speedCheckOK_1 && speedCheckOK_2; trajectoryOK = true; end
        
    end
    
    % Store details of the selected trajectory:
    freqInfo(:,nn) = xFreq;
    ampInfo(:,nn) = xAmp;
    
end

%% Visualise %%

% Time series:
plotYN = true; % plot sampled trajectories, yes or no?
if plotYN
    set(0,'defaultAxesFontSize',18);
    for nn = 1:N
        figure
        subplot(2,1,1); hold on
        scatter(t,targ(:,nn),[],targSpeed(:,nn))
        plot([0 runDur],[0 0], 'k-')
        colorbar
        caxis([0,maxSpeed])
        ylim([-(maxDev_mean+2*maxDev_sd),(maxDev_mean+2*maxDev_sd)])
        xlabel('Time in Trial')
        ylabel('Horizontal Position (deg)')
        subplot(2,1,2); hold on
        scatter(t,targVel(:,nn),[],targSpeed(:,nn))
        colorbar
        caxis([0,maxSpeed])
        ylim([-maxSpeed,maxSpeed])
        plot([0 runDur],[0 0], 'k-')
        plot([0 runDur],[5 5], 'k--')
        plot([0 runDur],[-5 -5], 'k--')
        xlabel('Time in Trial')
        ylabel('Horizontal Velocity (deg/s)')
    end
end

save('trajectories_test.mat','freqInfo','ampInfo','targ','t','targVel')