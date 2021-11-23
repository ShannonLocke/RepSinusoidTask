function [desmat] = makeDesignMat(conditions,repeats)
% MAKEDESIGNMAT takes your conditions in cell format and creates a
% randomised design matrix with n repeats of each condition.
%
% Created by SML Dec 2014, Edited Jan 2015.

if nargin < 2
    repeats = 1;
end

% Check if conditions input is correctly formatted:
assert(iscell(conditions))
[i,j] = size(conditions);
assert(min(i,j)==1);
if j ~= 1 % make column cell array
    conditions = conditions';
end
ncond = size(conditions,1);


% -- Determine all unique conditions -- %

% First condition:
aa = conditions{1,1};
if isrow(aa)
    aa = aa';
end
desmat = aa;

% Add additional conditions:
for ii = 2:ncond
    
    % Next condition:
    aa = conditions{ii,1};
    if iscolumn(aa)
        aa = aa';
    end
    
    % Update design matrix  
    desmatL = size(desmat,1);
    aaL = length(aa);
    desmat = repmat(desmat,aaL,1);
    aa = repmat(aa,desmatL,1);
    desmat = [desmat aa(:)];
    
end

% --  -- %

% Add repeats and randomise:
desmat = repmat(desmat,repeats,1);
desmat = ShuffleRC(desmat);

end