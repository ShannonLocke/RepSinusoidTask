function [dx] = denoiseSG(X,sf,pOrder,fLength,plotYN)
%DENOISESG Summary of this function goes here
%
% Created by SML Nov 2018

% Apply filter, calculating first and second order derivatives:
[nt,nX] = size(X); % number of signals
dt = 1/sf; % calc timesteps from sampling frequency
[b,g] = sgolay(pOrder,fLength); % get Savitsky-Golay filter
dx = zeros(nt,nX,3);
for n = 1:nX
    for p = 0:2
        dx(:,n,p+1) = conv(X(:,n), factorial(p)/(-dt)^p * g(:,p+1), 'same');
    end
end
dx(end,:,2) = 0; % remove edge effect
dx(end,:,3) = 0; % remove edge effect
% dx(:,:,2) = dx(:,:,2)/sf; % rescale
% dx(:,:,3) = dx(:,:,3)/sf^2; 

% Visualise:
if plotYN == 1
    for n = 1:nX
        figure; hold on
        plot(X(:,n),'-k') % unfiltered position
        plot(dx(:,n,1),'-b') % filtered position
        plot(60*diff(X(:,n)),'-m') % unfiltered velocity
        plot(dx(:,n,2),'-r') % filtered velocity
        % plot(dx(:,n,3),'-g')
    end
end

end

