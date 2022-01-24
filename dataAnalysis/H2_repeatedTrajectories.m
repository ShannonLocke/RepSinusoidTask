function [] = H2_repeatedTrajectories(resDir)
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
dataFromPath = resDir;
dataToPath_osfFiles = [resDir 'forOSF/'];
dataToPath_fig = [resDir 'H2_repeatedTrajectories/'];

% Load data:
fname = [dataFromPath, 'trialSummaryDataSPC.mat'];
load(fname, 'summaryData'); 
nSs = length(summaryData.sID);
sIDs = strsplit(num2str(summaryData.sID));
tidx = abs(summaryData.trajectory); % trajectory ID
N = length(unique(tidx(~isnan(tidx)))); % number of trajectories

%% Conf by trajectory:

% Preallocate:
conf = NaN([N,nSs]);
RMSE  = NaN([N,nSs]);
RMSE_byConf  = NaN([N,nSs,2]);

% Compute performance and proportion judged "better" per trajectory:
for nn = 1:nSs % EACH subject
    for tt = 1:N % EACH trajectory
        
        % Get overall statistics:
        idx = tidx(:,nn) == tt;
        getErr = summaryData.RMSE(idx,nn);
        getConf = summaryData.conf(idx,nn)==1;
        RMSE(tt,nn) = mean(getErr);
        conf(tt,nn) = mean(getConf);
        
        % Get statistics of median-split on performance:
        [~,midx] = sort(getErr); 
        midpt = length(midx)/2;
        idx = midx(1:midpt); % better 1/2 of trials
        conf_medSplit(tt,nn,1) = mean(getConf(idx)); 
        idx = midx((midpt+1):end); % worse 1/2 of trials
        conf_medSplit(tt,nn,2) = mean(getConf(idx)); 
        
    end
end

% Determine trajectory difficulty ranking:
zRMSE = (RMSE - repmat(mean(RMSE),[N,1]))./repmat(std(RMSE),[N,1]);
zRMSE_group = mean(zRMSE,2);
[~,diffRank] = sort(zRMSE_group);

% Compute confidence difference:
diffConf = squeeze(conf_medSplit(:,:,1)-conf_medSplit(:,:,2));
avgDiffConf = mean(diffConf);
semDiffConf = std(diffConf)/sqrt(N);

%% Plot results:

% Individual results panel:
for nn = 1:nSs % EACH subject
    
    fig = figure;
    sgtitle(['Participant #' sIDs{nn}],'FontSize',18)
    
    % RMSE by ranked trajectory:
    subplot(3,2,1); hold on
    bar(RMSE(diffRank,nn))
    title('Performance')
    xlabel('Ranked trajectory')
    ylabel('RMSE (deg)')
    yticks(0:0.5:3)
    xlim([0,(N+1)])
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
    
    % Proportion "better", split by performance:
    selSplitRMSE = squeeze(conf_medSplit(diffRank,nn,:));
    subplot(3,2,2); hold on
    bar(selSplitRMSE)
    title('Confidence by Performance')
    xlabel('Ranked trajectory')
    ylabel('Prop. "better"')
    yticks(0:0.2:1)
    xlim([0,(N+1)])
    ylim([0,1])
    legend({'better 1/2','worse 1/2'},'Location','NorthOutside','Box','Off','NumColumns',2)
    set(gca,'FontSize',14);
    set(gca,'linewidth',2);
    
    % Conf difference between better and worse halves:
    subplot(3,2,4); hold on
    bar(diffConf(diffRank,nn))
    title('Confidence Difference')
    xlabel('Ranked trajectory')
    ylabel('\Delta Prop. "better"')
    yticks(-0.4:0.2:0.4)
    xlim([0,(N+1)])
    ylim([-0.5,0.5])
    set(gca,'FontSize',14);
    set(gca,'linewidth',2);
    
    % Average RMSE difference between "better" and "worse":
    subplot(3,2,6); hold on
    bar(avgDiffConf(nn))
    er = errorbar(1, avgDiffConf(nn), semDiffConf(nn));
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    title('Avgerage Confidence Diff.')
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
    fname = [dataToPath_fig 'H2_s' sIDs{nn}];
    print(fig,fname,'-dpdf','-bestfit')
    
end

%% Group results figure:

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
yticks(0:0.5:3)
xlim([0,(N+1)])
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
barVals = mean(diffConf(diffRank,:),2);
semVals = std(diffConf(diffRank,:),0,2)/sqrt(nSs);
bar(barVals)
er = errorbar(1:N, barVals, semVals);
er.Color = [0 0 0];
er.LineStyle = 'none';
title('Confidence Difference')
xlabel('Ranked trajectory')
ylabel('\Delta RMSE (deg)')
yticks(-0.4:0.2:0.4)
xlim([0,(N+1)])
ylim([-0.5,0.5])
set(gca,'FontSize',14);
set(gca,'linewidth',2);

% Average RMSE difference between "better" and "worse":
subplot(2,2,4); hold on
barVals = mean(avgDiffConf);
semVals = std(avgDiffConf)/sqrt(nSs);
bar(barVals)
er = errorbar(1, barVals, semVals);
er.Color = [0 0 0];
er.LineStyle = 'none';
title('Avgerage Confidence Diff.')
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
fname = [dataToPath_fig 'H2_all'];
print(fig,fname,'-dpdf','-bestfit')

%% Statistical test:

disp('T-test results:...')
[h,p,ci,stats] = ttest(avgDiffConf)
effectSize = stats.tstat/sqrt(nSs);
disp(['Effect size is ' num2str(effectSize,5)])

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
