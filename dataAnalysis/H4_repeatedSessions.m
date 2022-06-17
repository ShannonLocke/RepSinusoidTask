function [] = H4_repeatedSessions(resDir)
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
dataToPath_fig = [resDir 'H4_sessionEffect/'];

% Load data:
fname = [dataFromPath, 'trialSummaryDataSPC.mat'];
load(fname, 'summaryData'); 
nSs = length(summaryData.sID);
sIDs = strsplit(num2str(summaryData.sID));
passCheck = summaryData.passChecksYN;

%% Compute AUROCs:

% Preallocate:
AUROC = NaN([nSs,2]);
pLow = cell([nSs,2]);
pHigh = cell([nSs,2]);

% Get AUROCs:
for nn = 1:nSs % EACH subject
    for ss = 1:2 % EACH session
        idx = summaryData.session(:,nn) == ss;
        RMSE = summaryData.RMSE(idx,nn);
        conf = summaryData.conf(idx,nn);
        [AUROC(nn,ss), pLow{nn,ss}, pHigh{nn,ss}, ~] = getAUROC(RMSE,conf);
    end
end
sessionEffect = diff(AUROC');

% Report summary statistics:
meanAUROC = num2str(nanmean(AUROC(:,1)),2);
semAUROC = num2str(nanstd(AUROC(:,1))/sqrt(sum(~isnan(AUROC(:,1)))),2);
display(['The AUROC mean & SEM for Session 1 is ' meanAUROC '+/-' semAUROC])
meanAUROC = num2str(nanmean(AUROC(:,2)),2);
semAUROC = num2str(nanstd(AUROC(:,2))/sqrt(sum(~isnan(AUROC(:,2)))),2);
display(['The AUROC mean & SEM for Session 2 is ' meanAUROC '+/-' semAUROC])

%% Plot results:

% Individual quantile-quantile plots:
for nn = 1:nSs % EACH subject
    txt_sID =  num2str(summaryData.sID(nn));
    fig = figure; hold on
    plot([0,1], [0,1], 'k--', 'Linewidth', 1, 'HandleVisibility','Off')
    plot(pLow{nn,1}, pHigh{nn,1}, 'b', 'Linewidth', 2)
    plot(pLow{nn,2}, pHigh{nn,2}, 'r', 'Linewidth', 2)
    title(['Participant #' txt_sID])
    xlabel('Cumulative P(RMSE | "worse")')
    ylabel('Cumulative P(RMSE | "better")')
    xticks(0:0.2:1)
    yticks(0:0.2:1)
    legend({'Session 1','Session 2'},'Location','NorthWest','Box','Off')
    axis square
    set(gca,'FontSize',16);
    set(gca,'linewidth',2);
    fname = [dataToPath_fig 'H4_s' txt_sID];
    print(fig,fname,'-dpdf','-bestfit')
end

% Histogram of AUROC differences (all participants):
fig = figure; hold on
plot([0,0], [0,3], 'k--', 'Linewidth', 1, 'HandleVisibility','off')
histogram(sessionEffect,'BinEdges',-0.5:0.05:0.5)
title(['Session Effect (All Particiants)'])
xlabel('Difference in AUROC (Session 2 - Session 1)')
ylabel('Frequency')
xlim([-0.5, 0.5])
axis square
set(gca,'FontSize',16);
set(gca,'linewidth',2);
fname = [dataToPath_fig 'H4_all'];
print(fig,fname,'-dpdf','-bestfit')

% Histogram of AUROC differences (all participants):
fig = figure; hold on
plot([0,0], [0,3], 'k--', 'Linewidth', 1, 'HandleVisibility','off')
histogram(sessionEffect(passCheck),'BinEdges',-0.5:0.05:0.5)
title(['Session Effect (Exclude Extremely Biased)'])
xlabel('Difference in AUROC (Session 2 - Session 1)')
ylabel('Frequency')
xlim([-0.5, 0.5])
axis square
set(gca,'FontSize',16);
set(gca,'linewidth',2);
fname = [dataToPath_fig 'H4_neb'];
print(fig,fname,'-dpdf','-bestfit')

%% Statistical tests:

% All participants:
disp('FOR ALL PARTICIPANTS...')
disp('The t-test results are...')
[h_all,p_all,~,stats_all] = ttest(sessionEffect)
effectSize_all = stats_all.tstat/sqrt(sum(~isnan(sessionEffect)));
disp(['Effect size is ' num2str(effectSize_all,5)])

% Exclude extreme bias:
disp('EXCLUDING EXTREMELY BIASED PARTICIPANTS...')
disp('The t-test results are...')
[h_neb,p_neb,~,stats_neb] = ttest(sessionEffect(passCheck))
effectSize_neb = stats_neb.tstat/sqrt(sum(passCheck));
disp(['Effect size is ' num2str(effectSize_neb,5)])

%% Summary statistics:

% Report summary statistics (all participants):
meanSE_all = nanmean(sessionEffect);
semSE_all = nanstd(sessionEffect)/sqrt(sum(~isnan(sessionEffect)));
disp('FOR ALL PARTICIPANTS...')
disp(['The session effect mean & SEM is ' num2str(meanSE_all,2) '+/-' num2str(semSE_all,2)])

% Report summary statistics (exclude extreme bias):
passCheck = summaryData.passChecksYN';
meanSE_neb = mean(sessionEffect(passCheck));
semSE_neb = std(sessionEffect(passCheck))/sqrt(sum(passCheck));
disp('EXCLUDING EXTREMELY BIASED PARTICIPANTS...')
disp(['The session effect mean & SEM is ' num2str(meanSE_neb,2) '+/-' num2str(semSE_neb,2)])

%% Prep summary data for OSF:

% Individual results:
T = table(summaryData.sID', AUROC(:,1), AUROC(:,2), ~passCheck);
T.Properties.VariableNames = {'subjectID', 'session1AUROC', 'session2AUROC', 'ExtremeBias'};
fname = [dataToPath_osfFiles 'H4_sessionEffect_indivResults.csv'];
writetable(T,fname);

% Group results:
T = table([meanSE_all; meanSE_neb], [semSE_all; semSE_neb], ...
    [stats_all.tstat; stats_neb.tstat], [stats_all.df; stats_neb.df], ...
    [p_all; p_neb], [h_all; h_neb], [effectSize_all; effectSize_neb], [0; 1]);
T.Properties.VariableNames = {'meanSessionEffect', 'semSessionEffect', ...
    'tStat', 'df', 'pVal', 'Sig', 'CohensD', 'ExcludeExtremeBias'};
fname = [dataToPath_osfFiles 'H4_sessionEffect_groupResults.csv'];
writetable(T,fname);

end
