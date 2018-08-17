function Xn = f_normal(x,method)
% - normalize values of X
%
% Usage: Xn = f_normal(x,'method');
%
% method = type of transformation to use
% Xn     = normalized data
%
% ----- Methods: -----
% Method  = '01': Presence/Absence
%           '2' : Square Root
%           '3' : Cube Root
%           '4' : Fourth Root
%        'log2' : Log2(x+1)
%    'log2_alt' : Log2(x) + 1, preserve 0's
%          'ln' : Ln(x+1)
%      'ln_alt' : Ln(x) + 1, preserve 0's
%         'log' : Log10(x+1)
%     'log_alt' : Log10(x) + 1, preserve 0's
%        'ang'  : arcsin(sqrt(x))
%        'sum1' : sum = 1
%
% SEE ALSO: f_center, f_stnd, f_ranging, f_transform

% -----Notes:-----
% This function is used to transform (normalize) the values in X by 12 methods.
% The angular ('ang') transformation is appropriate for percentages and
% proportions.
%
% Relative effect: '2' < '3' < '4' < log2 < ln < 'log' < 'ang' < '01'.
%
% Clarke (1993) recommended normalizing species abundance data for use in
% ordinations and cluster analyses to give more weight to less abundant
% species.
%
% Palmer (1993) strongly recommended log-transformation of explanatory
% variables in canonical analyses when, for biological reasons, small
% changes in LOW levels of a predictor has more of an impact than small
% changes in HIGH levels.

% -----References:-----
% Anderson, M. J., K. E. Ellingsen, and B. H. McArdle. 2006. Multivariate
%   dispersion as a measure of beta diversity. Ecology Letters 9: 683-693.
% Clarke, K. R. 1993. Non-parametric multivariate analyses of changes
%   in community structure. Aust. J. Ecol. 18: 117-143.%
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam. xv + 853 pp.
% Palmer, M. W. 1993. Putting things in even better order: the advantages
%   of canonical correspondence analysis. Ecology 74(8): 2215-2230.

% -----Author:-----
% by David L. Jones, April-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Oct-2003: updated doc
% Jan-2004: added 'sum1'
% May-2008: updated doc
% May-2009: '01' now returnes double vs. logical data
% Feb-2010: added 'alt' methods

% Convert to lower case:
method = lower(method);

switch method
   case {1,'01'} % Presence/Absence
      Xn = double(x>0);
   
   case {2,'2'} % Square-root transform
      Xn = x.^0.5;
      
   case {3, '3'} % Cubic-root transform
      Xn = x.^(1/3);
      
   case {4,'4'} % Fourth-root transform
      Xn = x.^0.25;
      
   case {'log2'} % Log2 transform
      Xn = log2(x + 1);
      
   case {'log2_alt'} % Alternate Log2 transform
      Xn      = zeros(size(x));     % preallocate
      idx     = find(x~=0);         % get indices to non-zero elements
      Xn(idx) = log2( x(idx) ) + 1; % transform non-zero elements
      
   case {'ln'} % Natural Log (Ln) transform
      Xn = log(x + 1);
      
      case {'ln_alt'} % Alternate Natural Log (Ln) transform
      Xn      = zeros(size(x));    % preallocate
      idx     = find(x~=0);        % get indices to non-zero elements
      Xn(idx) = log( x(idx) ) + 1; % transform non-zero elements
            
   case {'log'} % Log10 transfom
      Xn = log10(x + 1);
      
   case {'log_alt'} % Alternate Log10 transfom
      Xn      = zeros(size(x));      % preallocate
      idx     = find(x~=0);          % get indices to non-zero elements
      Xn(idx) = log10( x(idx) ) + 1; % transform non-zero elements
      
   case {'ang'}
      Xn = asin(x.^0.5); % Angular or Arcsin
      
   case {'sum1'}; % Sum of all elements = 1
      Xn = x ./ repmat(sum(x),size(x,1),1);
      
   otherwise  
      error('Unknown transformation method !');
end

