function [] = H3_temporalMetacogSensitivity(resDir)
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
dataFromPath = resDir;
dataToPath_osfFiles = [resDir 'forOSF/'];
dataToPath_fig = [resDir 'H3_temporalAUROCs/'];

% Load summary data:
fname = [dataFromPath, 'trialSummaryDataSPC.mat'];
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
    
    % Check if AUROC can be computed (there is both low and high conf reports):
    if all(summaryData.conf(:,nn)==1) || all(summaryData.conf(:,nn)==0)
        AUROC(nn,:) = nan([1,6]);
        for tt = 1:6 % EACH 1s time bin
            pLow{nn,tt} = NaN;
            pHigh{nn,tt} = NaN;
        end
        continue
    end
    
    % Compute AUROCs:
    idx_keep = summaryData.keepTrialYN(:,nn) & ~isnan(summaryData.RMSE(:,nn)); % exclude outlier/missing-data trials 
    N = sum(idx_keep); % number of trials
    conf = summaryData.conf(idx_keep,nn);
    for tt = 1:6 % EACH 1s time bin
        RMSE = squeeze(summaryData.binnedRMSE(idx_keep,tt,nn));
        [AUROC(nn,tt), pLow{nn,tt}, pHigh{nn,tt}, ~] = getAUROC(RMSE,conf);
    end
    
end

%% Plot individual results:

colorVals = [linspace(0,1,6); zeros([1,6]); linspace(1,0,6)]';
for nn = 1:nSs % EACH subject
    
    % Quantile-quantile plots:
    fig = figure;
    subplot(2,1,1); hold on
    plot([0,1], [0,1], 'k--', 'Linewidth', 1, 'HandleVisibility','off')
    for tt = 1:6
        plot(pLow{nn,tt}, pHigh{nn,tt}, 'Color', colorVals(tt,:),'Linewidth', 2)
    end
    title(['Participant #' sIDs{nn}])
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
    fname = [dataToPath_fig 'H3_s' sIDs{nn}];
    print(fig,fname,'-dpdf','-bestfit')
    
end

%% Summary statistics:

% Report summary statistics (all participants):
meanAUROC_all = nanmean(AUROC);
semAUROC_all = nanstd(AUROC)/sqrt(sum(~isnan(AUROC(:,1))));
disp('FOR ALL PARTICIPANTS...')
disp(['The mean AUROCs for 1-6s are ' num2str(meanAUROC_all,2)])
disp(['With SEMs of ' num2str(semAUROC_all,2)])

% Report summary statistics (exclude extreme bias):
passCheck = summaryData.passChecksYN';
meanAUROC_neb = nanmean(AUROC(passCheck,:));
semAUROC_neb = nanstd(AUROC(passCheck,:))/sqrt(sum(passCheck));
disp('EXCLUDING EXTREMELY BIASED PARTICIPANTS...')
disp(['The mean AUROCs for 1-6s are ' num2str(meanAUROC_neb,2)])
disp(['With SEMs of ' num2str(semAUROC_neb,2)])

%% Plot group results: 

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
ylim([0.2,1])
xticks(0:6)
yticks(0:0.2:1)
% legend(sIDs,'Location','northWest','Box','Off')
set(gca,'FontSize',16);
set(gca,'linewidth',2);
subplot(2,1,2); hold on
plot([0,6], [0.5,0.5], 'k--', 'Linewidth', 1)
errorbar(0.5:5.5, meanAUROC_all, semAUROC_all, 'ro-', ...
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

% Not extremely biased participant's only:
fig = figure; 
subplot(2,1,1); hold on 
plot([0,6], [0.5,0.5], 'k--', 'Linewidth', 1, 'HandleVisibility','off')
for nn = 1:nSs
    if ~passCheck(nn); continue; end
    plot(0.5:5.5,AUROC(nn,:), 'Linewidth', 2)
end
title('Temporal Metacognitive Sensitivity Curves')
xlabel('Time in Trial (sec)')
ylabel('Metacognitive Sensitivity')
xlim([0,6])
ylim([0.2,1])
xticks(0:6)
yticks(0:0.2:1)
% legend(sIDs,'Location','northWest','Box','Off')
set(gca,'FontSize',16);
set(gca,'linewidth',2);
subplot(2,1,2); hold on
plot([0,6], [0.5,0.5], 'k--', 'Linewidth', 1)
errorbar(0.5:5.5, meanAUROC_neb, semAUROC_neb, 'ro-', ...
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
fname = [dataToPath_fig 'H3_neb'];
print(fig,fname,'-dpdf','-bestfit')

%% Difference in AUROC, last 2 sec versus middle 2 sec:

splitHalfDiff_all = mean(AUROC(:,5:6),2) - mean(AUROC(:,3:4),2);
splitHalfDiff_neb = mean(AUROC(passCheck,5:6),2) - mean(AUROC(passCheck,3:4),2);

% All Participants:
fig = figure; subplot(2,1,1); hold on
plot([0,0], [0,10], 'k--', 'Linewidth', 1, 'HandleVisibility','off')
histogram(splitHalfDiff_all,'BinEdges',-0.5:0.05:0.5)
xlabel('Difference in AUROC (Last 2 sec - Middle 2 sec)')
ylabel('Frequency')
xlim([-0.5, 0.5])
axis square
set(gca,'FontSize',16);
set(gca,'linewidth',2);
title('All Participants')
sgtitle('Evidence for a Recency Effect','FontSize',18)

% Not extremely biased participant's only:
subplot(2,1,2); hold on
plot([0,0], [0,10], 'k--', 'Linewidth', 1, 'HandleVisibility','off')
histogram(splitHalfDiff_neb,'BinEdges',-0.5:0.05:0.5)
xlabel('Difference in AUROC (Last 2 sec - Middle 2 sec)')
ylabel('Frequency')
xlim([-0.5, 0.5])
axis square
set(gca,'FontSize',16);
set(gca,'linewidth',2);
title('Not Extreme-Bias Participants')
fname = [dataToPath_fig 'H3_all_histogram'];
print(fig,fname,'-dpdf','-bestfit')

%% Statistical test:

% All participants:
disp('FOR ALL PARTICIPANTS...')
disp('The t-test results are...')
[h_all,p_all,~,stats_all] = ttest(splitHalfDiff_all)
effectSize_all = stats_all.tstat/sqrt(sum(~isnan(splitHalfDiff_all)));
disp(['Effect size is ' num2str(effectSize_all,5)])

% Exclude extreme bias:
disp('EXCLUDING EXTREMELY BIASED PARTICIPANTS...')
disp('The t-test results are...')
[h_neb,p_neb,~,stats_neb] = ttest(splitHalfDiff_neb)
effectSize_neb = stats_neb.tstat/sqrt(sum(passCheck));
disp(['Effect size is ' num2str(effectSize_neb,5)])

%% Prep summary data for OSF:

% Individual results:
T = table(summaryData.sID', AUROC(:,1), AUROC(:,2), AUROC(:,3), AUROC(:,4), ...
    AUROC(:,5), AUROC(:,6));
T.Properties.VariableNames = {'subjectID', 'AUROCsec1', 'AUROCsec2', ...
    'AUROCsec3', 'AUROCsec4', 'AUROCsec5', 'AUROCsec6'};
fname = [dataToPath_osfFiles 'H3_temporalAUROCs_indivResults.csv'];
writetable(T,fname);

% Group results:
T = table([nanmean(splitHalfDiff_all); mean(splitHalfDiff_neb)], ...
    [nanstd(splitHalfDiff_all)/sqrt(sum(~isnan(splitHalfDiff_all))); std(splitHalfDiff_neb)/sqrt(sum(passCheck))], ...
    [stats_all.tstat; stats_neb.tstat], [stats_all.df; stats_neb.df], ...
    [p_all; p_neb], [h_all; h_neb], [effectSize_all; effectSize_neb], [0; 1]);
T.Properties.VariableNames = {'meanDiff', 'semDiff', 'tStat', 'df', ...
    'pVal', 'SigAboveChance', 'CohensD', 'ExcludeExtremeBias'};
fname = [dataToPath_osfFiles 'H3_temporalAUROCs_groupResults.csv'];
writetable(T,fname);

%% Investigate error mean and variance:

% Preallocate:
errTrace_mean = NaN([6000,nSs]);
errTrace_sd = NaN([6000,nSs]);

% Compute mean error curve:
for nn = 1:nSs % EACH subject
    
    % Get data:
    fname = [dataFromPath 'processed/processedEyeData_s' num2str(summaryData.sID(nn)) '.mat'];
    load(fname,'eyeDataPro')
    getError = eyeDataPro.errorEuclid;
    errTrace_mean(:,nn) = nanmean(getError,2);
    errTrace_sd(:,nn) = nanstd(getError,0,2);
    t = eyeDataPro.t;
    
    % Plot mean error curve of individuals:
    fig = figure; hold on
    plot(t, errTrace_mean(:,nn), 'k-', 'LineWidth', 2)
    plot(t, errTrace_mean(:,nn) - errTrace_sd(:,nn), 'r--', 'LineWidth', 2)
    plot(t, errTrace_mean(:,nn) + errTrace_sd(:,nn), 'r--', 'LineWidth', 2)
    title(['Participant #' sIDs{nn}])
    xlabel('Time in Trial (sec)')
    ylabel('Euclidean Error (deg)')
    xlim([0, 6])
    ylim([0, 6])
    set(gca,'FontSize',16);
    set(gca,'linewidth',2);
    fname = [resDir 'errorTraces/errorTrace_s' sIDs{nn}];
    print(fig,fname,'-dpdf','-bestfit')
    
end

% Plot mean error curve of group:
fig = figure; hold on
plot(t, errTrace_mean, 'k-', 'LineWidth', 1)
plot(t, mean(errTrace_mean,2), 'r-', 'LineWidth', 2)
title('Group Error Trace')
xlabel('Time in Trial (sec)')
ylabel('Euclidean Error (deg)')
xlim([0, 6])
ylim([0, 6])
set(gca,'FontSize',16);
set(gca,'linewidth',2);
fname = [resDir 'errorTraces/errorTrace_all'];
print(fig,fname,'-dpdf','-bestfit')

end
