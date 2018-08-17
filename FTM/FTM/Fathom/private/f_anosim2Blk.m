function r = f_anosim2Blk(dis,block,Rx,rank,permBlks)
% - utility function called by f_anosim2

% USAGE: r = f_anosim2Blk(dis,block,Rx,{rank},{permBlks});
%
% dis    = symmetric distance matrix
% block  = vector of integers (or chars) specifying BLOCKING factors
%          for rows/cols of distance matrix;
%          (e.g., block = [1 1 1 2 2 2 3 3])
%
% Rx     = vector of integers (or chars) specifying TREATMENT factors
%          for rows/cols of distance matrix;
%          (e.g., Rx = ['abcabcac'])
%
% rank   = optionally rank distances (default = 1)
% 
% permBlks = permute blocks for randomizations tests (default = 0)      
%
% r      = strength of treatment effect (averaged across all blocks)

% -----Dependencies:-----
% combnk.m STATISTICS Toolbox (could be replaced by choosenk)

% -----Author:-----
% by David L. Jones, Apr-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

if (nargin < 4), rank     = 1; end; % rank distances by default
if (nargin < 5), permBlks = 0; end; % no permutation by default

blockVar = unique(block);    % unique blocks
RxVar    = unique(Rx);       % unique treatments
noBlocks = length(blockVar); % # of blocks
noTreats = length(RxVar);    % # of treatments


for i = 1:noBlocks %%--Extract each block separately--%%
   xtract = zeros(noTreats,noTreats);  % initialize to max possible size
   index = find(block == blockVar(i)); % subset indices of block to extract
   
   % Extract subset of distance matrix corresponding to this block:
   sDis = dis(index,:);
   sDis = sDis(:,index);
   
   % Permute subset when called from a randomization test:
   if permBlks>0, sDis = f_shuffle(sDis); end; 
   
   % Find which treatments are PRESENT in this block: 
   pi = find(ismember(RxVar,Rx(index))>0);
   
   % Find which treatments are ABSENT in this block: 
   ai = find(ismember(RxVar,Rx(index))==0);
   
   % Use 'present indices' to correctly insert this subset into one
   % which is of full size when all treatments are present:
   xtract(pi,pi) = sDis;
   
   % Use 'absent indices' to specify missing values:
   xtract(ai,:) = NaN;
   xtract(:,ai) = NaN;
   
   % Make sure 0's are on diagonal:
   xtract(find(eye(noTreats))) = 0;
   
   % Unwrap as a column vector
   vectors(:,i) = f_unwrap(xtract);
end

% -----List of all pair-wise correlations:-----
pairs = combnk(1:noBlocks,2);    % list of all pairwise tests

for j = 1:size(pairs,1) % get correlation's between all possible pairs
   % Extract a pair of columns:
   thisPair = [vectors(:,pairs(j,1)) vectors(:,pairs(j,2))];
   
   % Keep only rows without NaN's:
   ki  = (isnan(thisPair)==0); % index of values to keep
   
   if size(ki,1)>0 % if there are any NaN's to exlcude:
      thisPair = thisPair(find(ki(:,1) .* ki(:,2)),:);
   end
   
   cMatrix(j) = f_corr(thisPair(:,1),thisPair(:,2),rank); % optionally rank
end;

% Average all pair-wise correlations:
r = mean(cMatrix);


