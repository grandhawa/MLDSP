function pA = f_adjustP(p,method)
% - adjust p-values for multiple comparison tests
%
% USAGE: pA = f_adjustP(p,method);
%
% p      = column vector of p-values obtained from multiple testing
% method = 'bon' (= Bonferroni), 'prog' (= progressive Bonferroni), 'ds' (=
%           Dunn-Sidak), 'holm' (= Holm's), or 'none' (= no correction applied)
%
% pA = column vector of adjusted p-values

% -----Notes:-----
% The Bonferroni method is the most conservative, while the Holmes
% correction is the most powerful for non-independent tests.
%
% NaN's are ignored by this function

% -----References:-----
% Borcard, D., F. Gillet, and P. Legendre. 2011. Numerical Ecology with R.
%  Springer, NY.
% Legendre, P. & L. Legendre. 2012. Numerical ecology. 3rd English ed.
%   Elsevier Science BV, Amsterdam.

% -----Author:-----
% by David L. Jones, Apr-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jul-2012: fixed formula for Dunn-Sidak
% Oct-2012: updated documentation
% Jan-2014: now supports progressive Bonferroni correction

method = lower(method);   % force lower case
p      = p(:);            % force as col vector
pOld   = p;               % keep the original
idx    = find(~isnan(p)); % get index to non-NaN values
p      = p(idx);          % extract non-NaN values
n      = numel(p);        % # tests
k      = repmat(n,n,1);   % simultaneous correction

switch method
   case 'bon' % Bonferroni Method (Legendre & Legendre, 2012: p. 23)
      pA = p.*k;
      
   case 'prog' % Progressive Bonferroni Method (Legendre & Legendre, 2012: p.745)
      pA = p .* (1:n)';
      
   case 'ds'  % Dunn-Sidak Method (Legendre & Legendre, 2012: p. 23)
      pA = 1 - (1 - p).^n;
      
   case 'holm' % Holm's Method (Legendre & Legendre, 2012: p. 23):
      [pS,idxS,idxU] = f_sort(p);     % sort in increasing order of p values
      pA             = pS.*(n:-1:1)'; % progressively relax the correction factor
      
      % Make sure adjusted p-value isn't smaller than the previous:
      for i = 2:n
         if pA(i) < pA(i-1)
            pA(i) = pA(i-1);
         end
      end
      
      % Restore original sort order:
      pA = pA(idxU);
      
   case 'none' % No correction applied
      pA = p;
      
   otherwise
      error('Unknown METHOD!')
end

% Set values > 1 to 1:
pA(pA>1) = 1;

% Add NaN's back into result:
pOld(idx) = pA;
pA        = pOld; % rename for output
