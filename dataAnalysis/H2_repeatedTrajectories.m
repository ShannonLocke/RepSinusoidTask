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

% Compute performance and proportion judged "better" per trajectory:
for nn = 1:nSs % EACH subject
    for tt = 1:N % EACH trajectory
        
        % Get overall statistics:
        idx = tidx(:,nn) == tt;
        getErr = summaryData.RMSE(idx,nn);
        getConf = summaryData.conf(idx,nn)==1;
        RMSE(tt,nn) = nanmean(getErr);
        conf(tt,nn) = mean(getConf);
        
        % Get statistics of median-split on performance:
        % (trials in which there is no recorded samples will be classified
        % as objectively worse half because sort places NaN values last).
        [~,midx] = sort(getErr); 
        midpt = length(midx)/2;
        idx = midx(1:midpt); % better 1/2 of trials
        conf_medSplit(tt,nn,1) = mean(getConf(idx)); 
        idx = midx((midpt+1):end); % worse 1/2 of trials
        conf_medSplit(tt,nn,2) = mean(getConf(idx)); 
        
    end
end

% Determine trajectory difficulty ranking:
zRMSE = (RMSE - repmat(nanmean(RMSE),[N,1]))./repmat(nanstd(RMSE),[N,1]);
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
    yticks(0:0.5:3.5)
    xlim([0,(N+1)])
    ylim([0,3.5])
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
    
    % Conf difference between better and worse 1/2:
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
    
    % Average confidence difference between better and worse 1/2:
    subplot(3,2,6); hold on
    bar(avgDiffConf(nn))
    er = errorbar(1, avgDiffConf(nn), semDiffConf(nn));
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    title('Avgerage Confidence Diff.')
    xlabel('Ranked trajectory')
    ylabel('\Delta Prop. "better"')
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

% Confidence difference between better and worse 1/2:
subplot(2,2,2); hold on
barVals = mean(diffConf(diffRank,:),2);
semVals = std(diffConf(diffRank,:),0,2)/sqrt(nSs);
bar(barVals)
er = errorbar(1:N, barVals, semVals);
er.Color = [0 0 0];
er.LineStyle = 'none';
title('Confidence Difference')
xlabel('Ranked trajectory')
ylabel('\Delta Prop. "better"')
yticks(-0.4:0.2:0.4)
xlim([0,(N+1)])
ylim([-0.5,0.5])
set(gca,'FontSize',14);
set(gca,'linewidth',2);

% Average confidence difference between better and worse 1/2:
subplot(2,2,4); hold on
barVals = mean(avgDiffConf);
semVals = std(avgDiffConf)/sqrt(nSs);
bar(barVals)
er = errorbar(1, barVals, semVals);
er.Color = [0 0 0];
er.LineStyle = 'none';
title('Avgerage Confidence Diff.')
xlabel('Ranked trajectory')
ylabel('\Delta Prop. "better"')
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

%% Summary statistics:

% Report summary statistics (all participants):
meanDiff_all = mean(avgDiffConf);
semDiff_all = std(avgDiffConf)/sqrt(length(avgDiffConf));
disp('FOR ALL PARTICIPANTS...')
disp(['The confidence difference mean & SEM is ' num2str(meanDiff_all,4) '+/-' num2str(semDiff_all,4)])

% Report summary statistics (exclude extreme bias):
passCheck = summaryData.passChecksYN';
meanDiff_neb = mean(avgDiffConf(passCheck));
semDiff_neb = std(avgDiffConf(passCheck))/sqrt(sum(passCheck));
disp('EXCLUDING EXTREMELY BIASED PARTICIPANTS...')
disp(['The confidence difference mean & SEM is ' num2str(meanDiff_neb,4) '+/-' num2str(semDiff_neb,4)])

%% Statistical test:

% All participants:
disp('FOR ALL PARTICIPANTS...')
disp('The t-test results:...')
[h_all,p_all,~,stats_all] = ttest(avgDiffConf)
effectSize_all = stats_all.tstat/sqrt(nSs);
disp(['Effect size is ' num2str(effectSize_all,5)])

% Exclude extreme bias:
disp('EXCLUDING EXTREMELY BIASED PARTICIPANTS...')
disp('The t-test results are...')
passCheck = summaryData.passChecksYN;
[h_neb,p_neb,~,stats_neb] = ttest(avgDiffConf(passCheck))
effectSize_neb = stats_neb.tstat/sqrt(sum(passCheck));
disp(['Effect size is ' num2str(effectSize_neb,5)])

%% Prep summary data for OSF:

% Individual results:
conf_betterRMSE = squeeze(conf_medSplit(:,:,1));
conf_worseRMSE = squeeze(conf_medSplit(:,:,2));
T = table(repelem(summaryData.sID',N), repmat((1:N)',[nSs,1]), RMSE(:), ...
    zRMSE(:), conf(:), conf_betterRMSE(:), conf_worseRMSE(:), diffConf(:));
T.Properties.VariableNames = {'subjectID', 'trajectoryID', 'RMSE', ...
    'normalisedRMSE', 'avgConf', 'medSplitConfBetter', ...
    'medSplitConfWorse', 'confDiff'};
fname = [dataToPath_osfFiles 'H2_repeatedTrajectories_indivResults.csv'];
writetable(T,fname);
 
% Group results:
T = table([meanDiff_all; meanDiff_neb], [semDiff_all; semDiff_neb], ...
    [stats_all.tstat; stats_neb.tstat], [stats_all.df; stats_neb.df], ...
    [p_all; p_neb], [h_all; h_neb], [effectSize_all; effectSize_neb], [0; 1]);
T.Properties.VariableNames = {'meanConfDiff', 'semConfDiff', 'tStat', ...
    'df', 'pVal', 'Sig', 'CohensD', 'ExcludeExtremeBias'};
fname = [dataToPath_osfFiles 'H2_repeatedTrajectories_groupResults.csv'];
writetable(T,fname);

%% Prep summary data for VSS poster:
y = 1:N; y = y(:);
barVals = mean(diffConf(diffRank,:),2);
semVals = std(diffConf(diffRank,:),0,2)/sqrt(nSs);
T = table(y, barVals, semVals);
fname = ['data_repeatedTrajectories.txt'];
writetable(T,fname,'Delimiter',' ')

end
