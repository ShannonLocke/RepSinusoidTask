function [] = H2_repeatedTrajectories
% This script will perform the analysis for H4: metacognitive sensitivity 
% does not significantly differ between Session 1 and Session 2 (see Locke, 
% Goettker, Gegenfurtner, & Mamassian, 2021).
%
% FILES GENERATED:
% 1) Individual quantitle-quantile plots of both sessions (*.pdf)
% 2) Individual-level session effect results for OSF (*.csv)
% 3) Group-level session effect results for OSF (*.csv)
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
dataToPath_osfFiles = 'output_data/forOSF/';
dataToPath_fig = 'output_figures/H2_repeatedTrajectories/';

% Load data:
fname = [dataFromPath, 'summaryDataSPC.mat'];
load(fname, 'summaryData'); 
nSs = length(summaryData.sID);
sIDs = strsplit(num2str(summaryData.sID));
tidx = abs(summaryData.trajectory); % trajectory ID
N = length(unique(tidx(:,1))); % number of trajectories

%% Conf by trajectory:

% Preallocate:
conf = NaN([N,nSs]);
RMSE  = NaN([N,nSs]);
RMSE_byConf  = NaN([N,nSs,2]);

% Compute performance and proportion judged "better" per trajectory:
for nn = 1:nSs % EACH subject
    for tt = 1:N % EACH trajectory
        idx = tidx(:,nn) == tt;
        RMSE(tt,nn) = nanmean(summaryData.RMSE(idx,nn));
        conf(tt,nn) = nanmean(summaryData.conf(idx,nn)==1);
        idx = tidx(:,nn) == tt & summaryData.conf(:,nn) == 1; % "better"
        RMSE_byConf(tt,nn,1) = nanmean(summaryData.RMSE(idx,nn));
        idx = tidx(:,nn) == tt & summaryData.conf(:,nn) == -1; % "worse"
        RMSE_byConf(tt,nn,2) = nanmean(summaryData.RMSE(idx,nn));
    end
end
% <= REMOVE NANMEAN ONCE BLINK CORRECTION ADDED

% Determine trajectory difficulty ranking:
zRMSE = (RMSE - repmat(mean(RMSE),[N,1]))./repmat(std(RMSE),[N,1]);
zRMSE_group = mean(zRMSE,2);
[~,diffRank] = sort(zRMSE_group);

% Compute confidence difference:
diffRMSE = squeeze(RMSE_byConf(:,:,2)-RMSE_byConf(:,:,1));
avgDiffRMSE = mean(diffRMSE);
semDiffRMSE = std(diffRMSE)/sqrt(N);

%% Plot results:

% Individual results panel:
for nn = 1:nSs % EACH subject
    
    fig = figure;
    txt_sID =  num2str(summaryData.sID(nn));
    sgtitle(['Participant #' txt_sID],'FontSize',18)
    
    % RMSE by ranked trajectory:
    subplot(3,2,1); hold on
    bar(RMSE(diffRank,nn))
    title('Performance')
    xlabel('Ranked trajectory')
    ylabel('RMSE (deg)')
    xticks(1:1:10)
    yticks(0:0.5:3)
    xlim([0,11])
    ylim([0,3])
    set(gca,'FontSize',14);
    set(gca,'linewidth',2);
    
    % Proportion "better" scatter plot:
    subplot(3,2,3); hold on
    plot(zRMSE(:,nn),conf(:,nn),'o', 'Linewidth', 2)
    title('Confidence')
    xlabel('Normalised RMSE (a.u.)')
    ylabel('Proportion "better"')
    xticks(-3:3)
    yticks(0:0.2:1)
    xlim([-3,3])
    ylim([0,1])
    set(gca,'FontSize',14);
    set(gca,'linewidth',2);
    
    % Proportion "better":
    subplot(3,2,5); hold on
    bar(mean(summaryData.conf(:,nn)==1))
    plot([0,2],[0.5,0.5],'k--')
    title('Confidence Bias')
    ylabel('Proportion "better"')
    xticks(NaN)
    yticks(0:0.2:1)
    xlim([0,2])
    ylim([0,1])
    set(gca,'FontSize',14);
    set(gca,'linewidth',2);
    axis square
    
    % RMSE difference between "better" and "worse", split by confidence:
    selSplitRMSE = squeeze(RMSE_byConf(diffRank,nn,:));
    subplot(3,2,2); hold on
    bar(selSplitRMSE)
    title('Performance by Confidence')
    xlabel('Ranked trajectory')
    ylabel('RMSE (deg)')
    xticks(1:1:10)
    yticks(0:0.5:3.5)
    xlim([0,11])
    ylim([0,3.5])
    legend({'"better"','"worse"'},'Location','NorthWest','Box','Off')
    set(gca,'FontSize',14);
    set(gca,'linewidth',2);
    
    % RMSE difference between "better" and "worse":
    subplot(3,2,4); hold on
    bar(diffRMSE(diffRank,nn))
    title('Performance Difference')
    xlabel('Ranked trajectory')
    ylabel('\Delta RMSE (deg)')
    xticks(1:1:10)
    yticks(-0.4:0.2:0.4)
    xlim([0,11])
    ylim([-0.5,0.5])
    set(gca,'FontSize',14);
    set(gca,'linewidth',2);
    
    % Average RMSE difference between "better" and "worse":
    subplot(3,2,6); hold on
    bar(avgDiffRMSE(nn))
    er = errorbar(1, avgDiffRMSE(nn), semDiffRMSE(nn));
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    title('Avgerage Performance Diff.')
    xlabel('Ranked trajectory')
    ylabel('\Delta RMSE (deg)')
    xticks(NaN)
    yticks(-0.4:0.2:0.4)
    xlim([0,2])
    ylim([-0.5,0.5])
    set(gca,'FontSize',14);
    set(gca,'linewidth',2);
    axis square
    
    % Save:
    fname = [dataToPath_fig 'H2_s' txt_sID];
    print(fig,fname,'-dpdf','-bestfit')
    
end

%%

fig = figure;
sgtitle('All Participants','FontSize',18)

% RMSE by ranked trajectory:
subplot(2,2,1); hold on
barVals = mean(RMSE(diffRank,:),2);
semVals = std(RMSE(diffRank,:),0,2)/sqrt(nSs);
bar(barVals)
er = errorbar(1:N, barVals, semVals);
er.Color = [0 0 0];
er.LineStyle = 'none';
title('Performance')
xlabel('Ranked trajectory')
ylabel('RMSE (deg)')
xticks(1:1:10)
yticks(0:0.5:3)
xlim([0,11])
ylim([0,3])
set(gca,'FontSize',14);
set(gca,'linewidth',2);

% Proportion "better" scatter plot:
subplot(2,2,3); hold on
for nn = 1:nSs
    plot(zRMSE(:,nn),conf(:,nn),'o', 'Linewidth', 2)
end
title('Confidence')
xlabel('Normalised RMSE (a.u.)')
ylabel('Proportion "better"')
xticks(-3:3)
yticks(0:0.2:1)
xlim([-3,3])
ylim([0,1])
% legend(sIDs,'Location','BestOutside','Box','Off')
set(gca,'FontSize',14);
set(gca,'linewidth',2);

% RMSE difference between "better" and "worse":
subplot(2,2,2); hold on
barVals = mean(diffRMSE(diffRank,:),2);
semVals = std(diffRMSE(diffRank,:),0,2)/sqrt(nSs);
bar(barVals)
er = errorbar(1:N, barVals, semVals);
er.Color = [0 0 0];
er.LineStyle = 'none';
title('Performance Difference')
xlabel('Ranked trajectory')
ylabel('\Delta RMSE (deg)')
xticks(1:1:10)
yticks(-0.4:0.2:0.4)
xlim([0,11])
ylim([-0.5,0.5])
set(gca,'FontSize',14);
set(gca,'linewidth',2);

% Average RMSE difference between "better" and "worse":
subplot(2,2,4); hold on
barVals = mean(avgDiffRMSE);
semVals = std(avgDiffRMSE)/sqrt(nSs);
bar(barVals)
er = errorbar(1, barVals, semVals);
er.Color = [0 0 0];
er.LineStyle = 'none';
title('Avgerage Performance Diff.')
xlabel('Ranked trajectory')
ylabel('\Delta RMSE (deg)')
xticks(NaN)
yticks(-0.4:0.2:0.4)
xlim([0,2])
ylim([-0.5,0.5])
set(gca,'FontSize',14);
set(gca,'linewidth',2);
axis square

% Save:
fname = [dataToPath_fig 'H2_all' txt_sID];
print(fig,fname,'-dpdf','-bestfit')

%% Histogram of AUROC differences:
fig = figure; hold on
plot([0,0], [0,3], 'k--', 'Linewidth', 1, 'HandleVisibility','off')
histogram(sessionEffect,'BinEdges',-0.5:0.05:0.5)
title(['Session Effect in Sample Population'])
xlabel('Difference in AUROC (Session 2 - Session 1)')
ylabel('Frequency')
xlim([-0.5, 0.5])
axis square
set(gca,'FontSize',16);
set(gca,'linewidth',2);
fname = [dataToPath_fig 'H4_all'];
print(fig,fname,'-dpdf','-bestfit')

%% Statistical test:

disp('T-test results:...')
[h,p,ci,stats] = ttest(avgDiffRMSE)

%% Prep summary data for OSF:

% % Individual results:
% T = table(summaryData.sID, AUROC(:,1), AUROC(:,2));
% T.Properties.VariableNames = {'subjectID', 'session1AUROC', 'session2AUROC'};
% fname = [dataToPath_osfFiles 'H4_sessionEffect_indivResults.csv'];
% writetable(T,fname);
% 
% % Group results:
% T = table(mean(sessionEffect), std(sessionEffect)/sqrt(nSs), stats.tstat, stats.df, p, h);
% T.Properties.VariableNames = {'meanDiffAUROC', 'semDiffAUROC', 'tStat', 'df', 'pVal', 'SigAboveChance'};
% fname = [dataToPath_osfFiles 'H4_sessionEffect_groupResults.csv'];
% writetable(T,fname);

end
