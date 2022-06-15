function [] = H1_metacogSensitivity(resDir)
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
dataFromPath = resDir;
dataToPath_osfFiles = [resDir 'forOSF/'];
dataToPath_fig = [resDir 'H1_metacogSensitivity/'];

% Load data:
fname = [dataFromPath, 'trialSummaryDataSPC.mat'];
load(fname, 'summaryData'); 
nSs = length(summaryData.sID);
sIDs = strsplit(num2str(summaryData.sID));

% Bootstrap
iterBS = 1000; % number of bootstrap iterations

%% Compute AUROCs:

% Preallocate:
AUROC = NaN([nSs,1]);
CIs = NaN([nSs,2]);
pLow = cell([nSs,1]);
pHigh = cell([nSs,1]);

% Get AUROC:
for nn = 1:nSs % EACH subject
    disp(['... Computing the AUROC of participant ' num2str(nn) ' of ' num2str(nSs)])
    
    % Check if AUROC can be computed (there is both low and high conf reports):
    if all(summaryData.conf(:,nn)==1) || all(summaryData.conf(:,nn)==0)
        AUROC(nn)= NaN;
        pLow{nn} = NaN;
        pHigh{nn} = NaN;
        CIs(nn,:) = [NaN, NaN];
        continue
    end
    
    % Point estimate:
    idx_keep = summaryData.keepTrialYN(:,nn); % exclude outlier trials
    N = sum(idx_keep); % number of trials
    [AUROC(nn), pLow{nn}, pHigh{nn}, ~] = getAUROC(summaryData.RMSE(idx_keep,nn),summaryData.conf(idx_keep,nn));
    
    % Bootstrap 95% CIs:
    AUROC_BS = NaN([1,iterBS]);
    for ii = 1:iterBS
        idx_rs = find(idx_keep); 
        getSamples = randi(N,[N,1]); % randomly sample indices with replacement
        idx_rs = idx_rs(getSamples); 
        AUROC_BS(ii) = getAUROC(summaryData.RMSE(idx_rs,nn),summaryData.conf(idx_rs,nn));
    end
    CIs(nn,:) = prctile(AUROC_BS,[2.5, 97.5]);
    
end

%% Plot individual results:

for nn = 1:nSs % EACH subject
    
    fig = figure;
    sgtitle(['Participant #' sIDs{nn} ' (AUROC: ' num2str(AUROC(nn),2) ...
        ',  95% CI: ' num2str(CIs(nn,1),2) '-' num2str(CIs(nn,2),2) ')'], ...
        'FontSize', 18)
    
    % Histograms of confidence split error:
    subplot(2,1,1); hold on
    idx_keep = summaryData.keepTrialYN(:,nn); % exclude outlier trials
    Err = summaryData.RMSE(:,nn);
    mE = mean(Err(idx_keep));
    plot(mE * [1,1], [0,30], 'k--', 'Linewidth', 1, 'HandleVisibility','off') % mean error
    cidx = idx_keep & summaryData.conf(:,nn) == 1; % high confidence trials
    histogram(Err(cidx), 'BinEdges', 0:0.1:5)
    cidx = idx_keep & summaryData.conf(:,nn) == -1; % low confidence trials
    histogram(Err(cidx), 'BinEdges', 0:0.1:5)
    title('Error Distributions')
    xlabel('Euclidean Error (deg)')
    ylabel('Frequency')
    xlim([0, 5])
    axis square
    legend({'High Conf','Low Conf'},'Location','NorthEast')
    set(gca,'FontSize',16);
    set(gca,'linewidth',2);
    
    % Individual quantile-quantile plots:
    subplot(2,1,2); hold on
    plot([0,1], [0,1], 'k--', 'Linewidth', 1)
    plot(pLow{nn}, pHigh{nn}, 'r', 'Linewidth', 2)
    title('Metacogntive Sensitivity')
    xlabel('Cumulative P(RMSE | "worse")')
    ylabel('Cumulative P(RMSE | "better")')
    xticks(0:0.2:1)
    yticks(0:0.2:1)
    axis square
    set(gca,'FontSize',16);
    set(gca,'linewidth',2);
    fname = [dataToPath_fig 'H1_s' sIDs{nn}];
    print(fig,fname,'-dpdf','-bestfit')
end

%% Summary statistics:

% Report summary statistics (all participants):
meanAUROC_all = nanmean(AUROC);
semAUROC_all = nanstd(AUROC)/sqrt(sum(~isnan(AUROC)));
disp('FOR ALL PARTICIPANTS...')
disp(['The AUROC mean & SEM is ' num2str(meanAUROC_all,2) '+/-' num2str(semAUROC_all,2)])

% Report summary statistics (exclude extreme bias):
passCheck = summaryData.passChecksYN';
meanAUROC_neb = mean(AUROC(passCheck));
semAUROC_neb = std(AUROC(passCheck))/sqrt(sum(passCheck));
disp('EXCLUDING EXTREMELY BIASED PARTICIPANTS...')
disp(['The AUROC mean & SEM is ' num2str(meanAUROC_neb,2) '+/-' num2str(semAUROC_neb,2)])

%% Plot group results:

% All participant's quantile-quantile plots:
fig = figure; hold on
plot([0,1], [0,1], 'k--', 'Linewidth', 1, 'HandleVisibility','off')
for nn = 1:nSs
    plot(pLow{nn}, pHigh{nn}, 'Linewidth', 1)
end
title(['All Participants (AUROC: ' num2str(meanAUROC_all,2) '\pm' num2str(semAUROC_all,2) ')'])
xlabel('Cumulative P(RMSE | "worse")')
ylabel('Cumulative P(RMSE | "better")')
xticks(0:0.2:1)
yticks(0:0.2:1)
% legend(sIDs,'Location','southEast','Box','Off')
axis square
set(gca,'FontSize',16);
set(gca,'linewidth',2);
fname = [dataToPath_fig 'H1_all'];
print(fig,fname,'-dpdf','-bestfit')

% Not extremely biased participant's quantile-quantile plots:
fig = figure; hold on
plot([0,1], [0,1], 'k--', 'Linewidth', 1, 'HandleVisibility','off')
for nn = 1:nSs
    if ~passCheck(nn); continue; end
    plot(pLow{nn}, pHigh{nn}, 'Linewidth', 1)
end
title(['Not Extreme-Bias Participants (AUROC: ' num2str(meanAUROC_neb,2) '\pm' num2str(semAUROC_neb,2) ')'])
xlabel('Cumulative P(RMSE | "worse")')
ylabel('Cumulative P(RMSE | "better")')
xticks(0:0.2:1)
yticks(0:0.2:1)
% legend(sIDs,'Location','southEast','Box','Off')
axis square
set(gca,'FontSize',16);
set(gca,'linewidth',2);
fname = [dataToPath_fig 'H1_neb'];
print(fig,fname,'-dpdf','-bestfit')

%% Statistical test:

% All participants:
disp('FOR ALL PARTICIPANTS...')
disp('The t-test results are...')
[h_all,p_all,~,stats_all] = ttest(AUROC-0.5)
effectSize_all = stats_all.tstat/sqrt(sum(~isnan(AUROC)));
disp(['Effect size is ' num2str(effectSize_all,5)])

% Exclude extreme bias:
disp('EXCLUDING EXTREMELY BIASED PARTICIPANTS...')
disp('The t-test results are...')
[h_neb,p_neb,~,stats_neb] = ttest(AUROC(passCheck)-0.5)
effectSize_neb = stats.tstat/sqrt(sum(passCheck));
disp(['Effect size is ' num2str(effectSize_neb,5)])

%% Prep summary data for OSF:

% Individual results:
T = table(summaryData.sID', AUROC, CIs(:,1), CIs(:,2), CIs(:,1)>0.5, ~passCheck);
T.Properties.VariableNames = {'subjectID', 'AUROC', 'lower95CI', 'upper95CI', ...
    'SigAboveChance', 'ExtremeBias'};
fname = [dataToPath_osfFiles 'H1_metacogSensitivity_indivResults.csv'];
writetable(T,fname);

% Group results:
T = table([meanAUROC_all; meanAUROC_neb], [semAUROC_all; semAUROC_neb], ...
    [stats_all.tstat; stats_neb.tstat], [stats_all.df; stats_neb.df], ...
    [p_all; p_neb], [h_all; h_neb], [effectSize_all; effectSize_neb], [0; 1]);
T.Properties.VariableNames = {'meanAUROC', 'semAUROC', 'tStat', 'df', ...
    'pVal', 'SigAboveChance', 'CohensD', 'ExcludeExtremeBias'};
fname = [dataToPath_osfFiles 'H1_metacogSensitivity_groupResults.csv'];
writetable(T,fname);

%% Prep summary data for VSS poster:

% Get subject curves:
maxLengthCell=max(cellfun('size',pLow,2)); % finding the longest vector in the cell array
curvesAUROC = NaN([maxLengthCell,2*nSs]); % pre allocate storage vector
colNames = {};
for ii = 1:length(pLow)
    idx_low = 2*(ii-1) + 1;
    idx_high = idx_low + 1;
    getLen = length(pLow{ii});
    curvesAUROC(1:getLen,idx_low) = pLow{ii}';
    curvesAUROC(1:getLen,idx_high) = pHigh{ii}';
    colNames{idx_low} = ['x' num2str(ii)];
    colNames{idx_high} = ['y' num2str(ii)];
end

% Get group average curve:
equiv_dPrime = norminv(mean(AUROC)) * sqrt(2); % equivalent d' for the mean AUROC
xvals = linspace(-3,3,maxLengthCell)'; % set sampling spacing to fit current matrix
curvesAUROC(:,idx_low+2) = normcdf(xvals-equiv_dPrime/2); % ``worse'' distribution cumulative value
curvesAUROC(:,idx_high+2) = normcdf(xvals+equiv_dPrime/2); % ``better'' distribution cumulative value
curvesAUROC([end-1 end],[end-1 end]) = [1, 1; 1.1, 0]; % replace final two values to create triangle shape for shading
colNames{idx_low+2} = 'xAll';
colNames{idx_high+2} = 'yAll';

% Get visuomotor curve (Locke et al., 2020):
equiv_dPrime = norminv(0.68) * sqrt(2); % equivalent d' for the mean AUROC
curvesAUROC(:,idx_low+4) = normcdf(xvals-equiv_dPrime/2); % ``worse'' distribution cumulative value
curvesAUROC(:,idx_high+4) = normcdf(xvals+equiv_dPrime/2); % ``better'' distribution cumulative value
curvesAUROC([end-1 end],[end-1 end]) = [1, 1; 1.1, 0]; % replace final two values to create triangle shape for shading
colNames{idx_low+4} = 'xVM';
colNames{idx_high+4} = 'yVM';

% Export data file:
T = array2table(curvesAUROC);
T.Properties.VariableNames = colNames;
fname = ['data_AUROCs.txt'];
writetable(T,fname,'Delimiter',' ')

end
