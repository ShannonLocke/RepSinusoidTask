function [] = H1_metacogSensitivity
% This script will perform the analysis for H1: Participants have above 
% chance metacognitive sensitivity for sorting objectively better eye-
% tracking performance from objectively worse performance across the entire 
% set of repeated trajectories (see Locke, Goettker, Gegenfurtner, & 
% Mamassian, 2021).
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
dataToPath_osfFiles = 'output_data/forOSF/';
dataToPath_fig = 'output_figures/H1_metacogSensitivity/';

% Load data:
fname = [dataFromPath, 'summaryDataSPC.mat'];
load(fname, 'summaryData'); 
nSs = length(summaryData.sID);

% Bootstrap
iterBS = 1000; % number of bootstrap iterations
N = length(summaryData.conf); % number of trials

%% Compute AUROCs:

% Preallocate:
AUROC = NaN([nSs,1]);
CIs = NaN([nSs,2]);
pLow = cell([nSs,1]);
pHigh = cell([nSs,1]);

% Get AUROC:
for nn = 1:nSs % EACH subject
    
    % Point estimate:
    [AUROC(nn), pLow{nn}, pHigh{nn}, ~] = getAUROC(summaryData.RMSE(:,nn),summaryData.conf(:,nn));
    
    % Bootstrap 95% CIs:
    AUROC_BS = NaN([1,iterBS]);
    for ii = 1:iterBS
        idx = randi(N,[N,1]); % randomly sample indices with replacement
        AUROC_BS(ii) = getAUROC(summaryData.RMSE(idx,nn),summaryData.conf(idx,nn));
    end
    CIs(nn,:) = prctile(AUROC_BS,[2.5, 97.5]);
    
end

% Report summary statistics:
meanAUROC = num2str(mean(AUROC),2);
semAUROC = num2str(std(AUROC)/sqrt(nSs),2);
display(['The AUROC mean & SEM is ' meanAUROC '+/-' semAUROC])

%% Plot results:

% Individual quantile-quantile plots:
for nn = 1:nSs % EACH subject
    txt_sID =  num2str(summaryData.sID(nn));
    fig = figure; hold on
    plot([0,1], [0,1], 'k--', 'Linewidth', 1)
    plot(pLow{nn}, pHigh{nn}, 'r', 'Linewidth', 2)
    title(['Participant #' txt_sID ' (AUROC: ' num2str(AUROC(nn),2) ...
        ',  95% CI: ' num2str(CIs(nn,1),2) '-' num2str(CIs(nn,2),2) ')'])
    xlabel('Cumulative P(RMSE | "worse")')
    ylabel('Cumulative P(RMSE | "better")')
    xticks(0:0.2:1)
    yticks(0:0.2:1)
    axis square
    set(gca,'FontSize',16);
    set(gca,'linewidth',2);
    fname = [dataToPath_fig 'H1_s' txt_sID];
    print(fig,fname,'-dpdf','-bestfit')
end

% All participant's quantile-quantile plots:
fig = figure; hold on
plot([0,1], [0,1], 'k--', 'Linewidth', 1)
for nn = 1:nSs
    plot(pLow{nn}, pHigh{nn}, 'Linewidth', 1)
end
title(['All Participants (AUROC: ' meanAUROC '\pm' semAUROC ')'])
xlabel('Cumulative P(RMSE | "worse")')
ylabel('Cumulative P(RMSE | "better")')
xticks(0:0.2:1)
yticks(0:0.2:1)
axis square
set(gca,'FontSize',16);
set(gca,'linewidth',2);
fname = [dataToPath_fig 'H1_all'];
print(fig,fname,'-dpdf','-bestfit')

%% Statistical test:

disp('T-test results:...')
[h,p,ci,stats] = ttest(AUROC-0.5)

%% Prep summary data for OSF:

% Individual results:
T = table(summaryData.sID, AUROC, CIs(:,1), CIs(:,2), CIs(:,1)>0.5);
T.Properties.VariableNames = {'subjectID', 'AUROC', 'lower95CI', 'upper95CI', 'SigAboveChance'};
fname = [dataToPath_osfFiles 'H1_metacogSensitivity_indivResults.csv'];
writetable(T,fname);

% Group results:
T = table(mean(AUROC), std(AUROC)/sqrt(nSs), stats.tstat, stats.df, p, h);
T.Properties.VariableNames = {'meanAUROC', 'semAUROC', 'tStat', 'df', 'pVal', 'SigAboveChance'};
fname = [dataToPath_osfFiles 'H1_metacogSensitivity_groupResults.csv'];
writetable(T,fname);

end
