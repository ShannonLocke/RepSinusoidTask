function [] = H3_temporalMetacogSensitivity
% This script will perform the analysis for H3: The group-averaged temporal 
% metacognitive sensitivity curve shows a recency effect (i.e., greater 
% metacognitive sensitivity for later time bins). For more details see 
% Locke, Goettker, Gegenfurtner, & Mamassian (2021).
%
% FILES GENERATED:
% 1) Individual quantitle-quantile plots (*.pdf)
% 2) Individual-level metacognitive sensitivity results for OSF (*.csv)
% 3) Group-level metacognitive sensitivity results for OSF (*.csv)
%
% Data: https://osf.io/j9asn/
% Code: https://github.com/ShannonLocke/RepSinusoidTask
%
% Created by SML Jan 2022
% License: CC-By Attribution 4.0 International

%% Preamble:

% Reset seed to default for reproducability of bootstrapping CIs:
rng default

% Directories:
addpath('dataAnalysisToolbox')
dataFromPath = 'output_data/';
dataToPath_matFiles = 'output_data/hypothesis/';
dataToPath_osfFiles = 'output_data/forOSF/';
dataToPath_fig = 'output_figures/H3_temporalAUROCs/';

% Load data:
fname = [dataFromPath, 'summaryDataSPC.mat'];
load(fname, 'summaryData'); 
nSs = length(summaryData.sID);
sIDs = strsplit(num2str(summaryData.sID));

%% Compute temporal AUROCs:

% Preallocate:
AUROC = NaN([nSs,6]);
pLow = cell([nSs,6]);
pHigh = cell([nSs,6]);

% Get temporal AUROCs:
for nn = 1:nSs % EACH subject
    conf = summaryData.conf(:,nn);
    for tt = 1:6 % EACH 1s time bin
        RMSE = squeeze(summaryData.binnedRMSE(:,tt,nn));
        [AUROC(nn,tt), pLow{nn,tt}, pHigh{nn,tt}, ~] = getAUROC(RMSE,conf);
    end
end

%% Plot results:

% Individual plots:
fig = figure;
colorVals = [linspace(0,1,6); zeros([1,6]); linspace(1,0,6)]';
for nn = 1:nSs % EACH subject
    
    % Quantile-quantile plots:
    txt_sID =  num2str(summaryData.sID(nn));
    subplot(2,1,1); hold on
    plot([0,1], [0,1], 'k--', 'Linewidth', 1, 'HandleVisibility','off')
    for tt = 1:6
        plot(pLow{nn,tt}, pHigh{nn,tt}, 'Color', colorVals(tt,:),'Linewidth', 2)
    end
    title(['Participant #' txt_sID])
    xlabel('Cumulative P(RMSE | "worse")')
    ylabel('Cumulative P(RMSE | "better")')
    xticks(0:0.2:1)
    yticks(0:0.2:1)
    legend({'1 sec', '2 sec', '3 sec', '4 sec', '5 sec', '6 sec'}, ...
        'Location','bestOutside','Box','Off')
    axis square
    set(gca,'FontSize',16);
    set(gca,'linewidth',2);
    
    % Temporal plots:
    txt_sID =  num2str(summaryData.sID(nn));
    subplot(2,1,2); hold on
    plot([0,6], [0.5,0.5], 'k--', 'Linewidth', 1)
    plot(0.5:5.5,AUROC(nn,:),'k-','Linewidth', 2)
    for tt = 1:6
        plot(tt-0.5, AUROC(nn,tt), 'ko', 'MarkerFaceColor', colorVals(tt,:), 'MarkerSize', 10)
    end
    xlabel('Time in Trial (sec)')
    ylabel('Metacognitive Sensitivity')
    xlim([0,6])
    ylim([0.4,1])
    xticks(0:6)
    yticks(0:0.2:1)
    set(gca,'FontSize',16);
    set(gca,'linewidth',2);
    
    % Save figure
    fname = [dataToPath_fig 'H3_s' txt_sID];
    print(fig,fname,'-dpdf','-bestfit')
    
end

% All participant's temporal plots:
fig = figure; 
subplot(2,1,1); hold on 
plot([0,6], [0.5,0.5], 'k--', 'Linewidth', 1, 'HandleVisibility','off')
for nn = 1:nSs
    plot(0.5:5.5,AUROC(nn,:), 'Linewidth', 2)
end
title('Temporal Metacognitive Sensitivity Curves')
xlabel('Time in Trial (sec)')
ylabel('Metacognitive Sensitivity')
xlim([0,6])
ylim([0.4,1])
xticks(0:6)
yticks(0:0.2:1)
legend(sIDs,'Location','northWest','Box','Off')
set(gca,'FontSize',16);
set(gca,'linewidth',2);
subplot(2,1,2); hold on
plot([0,6], [0.5,0.5], 'k--', 'Linewidth', 1)
errorbar(0.5:5.5, mean(AUROC,1), std(AUROC,0,1)/sqrt(nSs), 'ro-', ...
    'Linewidth', 2, 'MarkerFaceColor', 'r', 'MarkerSize', 10)
title('Group Mean and SEM')
xlabel('Time in Trial (sec)')
ylabel('Metacognitive Sensitivity')
xlim([0,6])
ylim([0.4,1])
xticks(0:6)
yticks(0:0.2:1)
set(gca,'FontSize',16);
set(gca,'linewidth',2);
fname = [dataToPath_fig 'H3_all'];
print(fig,fname,'-dpdf','-bestfit')

%% Statistical test:

% <= FILL THIS IN!!!

%% Prep summary data for OSF:

% Individual results:
T = table(summaryData.sID, AUROC(:,1), AUROC(:,2), AUROC(:,3), AUROC(:,4), ...
    AUROC(:,5), AUROC(:,6));
T.Properties.VariableNames = {'subjectID', 'AUROCsec1', 'AUROCsec2', ...
    'AUROCsec3', 'AUROCsec4', 'AUROCsec5', 'AUROCsec6'};
fname = [dataToPath_osfFiles 'H3_temporalAUROCs_indivResults.csv'];
writetable(T,fname);

% Group results:
% <= FILL THIS IN!!!

end
