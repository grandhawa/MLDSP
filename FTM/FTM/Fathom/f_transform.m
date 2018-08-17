function Xt = f_transform(X,method)
% - several methods for data transformation
%
% Usage: Xt = f_transform(X,method);
%
% X      = input matrix (rows = observations, cols = variables)
% method = type of transformation
% 
% Xt     = transformed data
%
% ----- Methods: -----
% Biotic data  = 1: Square Root
%                2: Fourth Root
%                3: Log10(x+1)
%                4: Natural Log(x+1)
%                5: Log2(x+1)
%                6: Species Standardization
%            'hel': Hellinger transformation
%
% Abiotic data = 7: Normalize by col
%                8: Sum-of-Squares = 1 by col
%
% ----- Ecological Applications: -----
%"Transform" your BIOTIC raw data before calculating a dissimilarity matrix
% to give more weight to rarer species. Relative effect: squart root < fourth
% root < log.
%
%"Standardize" abundances by species (col-wise) in your BIOTIC raw data instead
% of "transforming" when using in a SPECIES ordination or cluster analysis.
%
%"Normalize" your ABIOTIC data col-wise (by variable) AFTER making any
% transformations. For each value of a variable (col), the mean is subtracted
% and divided by the standard deviation so the variance of each variable = 1.
% This makes all variables equally weighted & is especially useful for abiotic
% variables that are measured on different scales or in different units.
%
% SEE ALSO: f_normal, f_stnd

% -----Notes:-----
% When creating an ordination or cluster analysis of species (vs. samples),
% Clarke & Warwick (2001) suggest you use "Species Standardization", remove the
% rarer taxa, then create a dissimilarity matrix from the "untransformed" data.

% -----References:-----
% Clarke, K. R. and R. M. Warwick. 2001. Change in Marine Communities: An
% Approach to Statistical Analysis and Interpretation. 2nd edition. PRIMER-E,
% Plymouth, UK 172 p.

% -----Author:-----
% by David L. Jones, April-2001
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% 03-Apr-2002: fixed 'divide by zero' error in normalize
% 18-Apr-2002: improved error message
% Oct-2004: added more documentation
% Mar-2011: format changed so row = obs, cols = variables to maintain
%           consistency with the rest of the toolbox; removed for loops;
%           now support Hellinger's transformation
% Jun-2012: updated documentation

switch method
   case 1
      Xt = X.^0.5;       % square-root
   case 2
      Xt = X.^0.25;      % fourth-root
   case 3
      Xt = log10(X + 1); % Log(x+1)
   case 4
      Xt = log(X + 1);   % Ln(x+1)
   case 5
      Xt = log2(X + 1);  % Log2(x+1)
      
   case 6 % Standardize species abundances by col sum:
      Xt = (X ./ repmat(sum(X),size(X,1),1) ) * 100;
   
   case 7 % 'Normalize abiotic data' (= standardize, zscores):
      Xt = f_stnd(X);
   
   case 8 % Normalize so sum-of-squares = 1
      Xt = X ./ repmat(sqrt(sum(X.^2)),size(X,1),1);
   
   case 'hel' % Hellinger transform (eq. 13 in Legendre & Gallagher, 1998):
      nc        = size(X,2);        % get # columns
      Xt        = zeros(size(X));   % preallocate
      rowSum    = sum(X,2);
      idx       = find(rowSum > 0); % skip these to prevent divide-by-zero error
      Xt(idx,:) = sqrt(X(idx,:) ./ repmat(rowSum(idx),1,nc)); % Hellinger-transformed
      
   otherwise
      error('Unknown transformation unknown!');
end

% Final check for imaginary numbers:
if sum(~isreal(Xt(:)))>0
   disp('');
   disp('--------------------------------------------------');
   disp('Transformed data contain imaginary numbers!')
   disp('Did you take the square-root of a negative number?')
   disp('--------------------------------------------------');
   disp('');
end





















