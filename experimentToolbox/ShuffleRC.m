function [X] = ShuffleRC(X,dim,specrange)
% SHUFFLERC shuffles specific rows or columns of a given matrix or vector.
%
% [Y] = ShuffleRC(X,DIM,SPECRANGE) 
%
% X: the input 2D matrix or vector.
% DIM: specifies if rows or columns are to be shuffled (optional, default rows). 
% SPECRANGE: is the specific rows/columns to be shuffled (optional, default all).
% 
% Created by SML Dec 2014

% Defaults:
if nargin < 3
   specrange = [];
   if nargin < 2
       dim = 1;
   end
end

assert(dim==1||2) % valid dimension

% Transpose if column: 
if dim == 2
    X = X';
end

% Set specrange to be all rows if none selected:
if isempty(specrange)
    specrange = 1:size(X,1);
end

nShuffle = length(specrange); % number rows/col to be shuffled
newsort = randperm(nShuffle); % random order for shuffle
aa = X(specrange,:); % pull out section to be shuffled
aa = aa(newsort,:); % shuffle!
X(specrange,:) = aa; % put back in to matrix

% Transpose back if column: 
if dim == 2
    X = X';
end

end