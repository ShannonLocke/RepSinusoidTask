function [AUROC, pLow, pHigh, eVals] = getAUROC(e,c)
% This function computes the Area Under the ROC for deriving metacognitive
% sensitivity. Inputs are:
%
% e:- error (e.g., RMSE)
% c:- confidence (1 = low, 0 = high or -1 = low, 1 = high)
%
% Created by SML.

% Convert confidence input to logical:
if min(c) == 0 % worse coded as 0:
    c = c == 0;
elseif min(c) == -1 % worse coded as -1
    c = c == -1;
end

% Remove NaN values which can mess with the calculation:
removeNaNs = isnan(e) || isnan(c);
if any(removeNaNs)
    warning([num2str(sum(removeNaNs)), ' NaN entries will be removed from the AUROC calculation'])
    e = e(~removeNaNs);
    c = c(~removeNaNs);
end

% Count number of low/high reports:
nLow = sum(c);
nHigh = sum(~c);

% Get bounds for sweep:
rval = 0.01;
lowBound = rval * round(min(e)/rval);
highBound = rval * round(max(e)/rval);
eVals = lowBound:rval:highBound;

% Get cumulatives for low and high:
e = repmat(e, 1, length(eVals));
eVals = repmat(eVals, size(e,1), 1);
belowCrit = e <= eVals;
pLow = sum(belowCrit(c,:))/nLow;
pHigh = sum(belowCrit(~c,:))/nHigh;

% Compute AUROC:
AUROC = trapz(pLow,pHigh);

end