function out = f_outlier(prox,Y,plt)
% - calculate outlier measures for a symmetric proximity (similarity) matrix
%
% USAGE: out = f_outlier(prox,Y,plt);
%
% prox = symmetric similarity matrix
% Y    = column vector of integers specifying class membership
%        (default = all belong to the same class)
% plt  = create plot    (default = 0)
%
% out  = outlier measures for each observation
%
% SEE ALSO: f_grpOutlier

% -----References:-----
% This code follows that used by outlier.R from the RandomForest package for R.
% 
% Breiman, L., and A. Cutler. 2003. Manual on setting up, using, and
%   understanding Random Forests v4.0. Technical Report. 
%   ftp://ftp.stat.berkeley.edu/pub/users/breiman/Using_random_forests_v4.0.pdf

% -----Notes:-----
% For each observation, an outlying measure is computed as n / sum(squared
% proximity), which is then normalized by subtracting the median and divided by
% the median absolute deviation within each class.
% 
% Values > 10 suggest an observation may be an outlier (Breiman & Cutler, 2003)

% -----Author:-----
% by David L. Jones, Sep-2009
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Oct-2010: replaced unique with f_unique; the inverse is taken slightly
%           different than the R code (following Breiman & Cutler); negative
%           values are converted to 0 (following Breiman &* Cutler).

% -----Check input & set defaults:-----
if (nargin < 2), Y   =  ones(size(prox,1),1); end % default only 1 class
if (nargin < 3), plt =  0; end                    % default no plot

% Check if symmetric similarity matrix:
if (f_issymSim(prox,1) == 0)
   error('Input PROX must be square symmetric similarity matrix');
end

% Replace empty inputs with defaults:
if isempty(Y), Y =  ones(size(prox,1),1); end % default only 1 class

if (size(Y,1) ~= size(prox,1)), error('Y & DIS must have same # rows!'); end

% Check proximities:
if (min(prox(:))<0) || (max(prox(:))>1)
   error('Proximities must range from 0-1!')
end
% -------------------------------------

% Recode Y to begin at 1:
Y = f_recode(Y);

grps   = f_unique(Y);         % get list of groups, unsorted
noGrps = size(grps(:),1);     % get # groups
n      = size(Y,1);           % get # observations
out    = repmat(NaN,n,1);     % preallocate

for i=1:noGrps
   % Extract subset corresponding to each group:
   idx     = find((Y==grps(i))==1);  % index to obs in this group
   subProx = prox(idx,idx);          % extract corresponding proximities
   ns      = size(subProx,1);        % get size of subset
   
   % Sum of squared proximities:
   SS = sum(subProx.^2)';
       
   % Prevent divide by zeros:
   SS(SS==0)=1;
   
   % Take the inverse:
   %    SS = repmat(n,ns,1) ./ SS;
   SS = repmat(1,ns,1) ./ SS; % Beiman & Cutler, 2003
     
   % Median & MAD:
   med_SS = median(SS); % median
   mad_SS = f_mad(SS);  % median absolute deviation
   
   % Outlierness:
   out(idx) = (SS - med_SS) / mad_SS;
end

% Set negative values to 0 (Beiman & Cutler, 2003):
out(out<0) = 0;


% -----Plot Oultier Measure:-----
if (plt>0)
   figure('Name','Outlierness');
   set(gcf,'color','w'); % set bg color to white
   box on;
   hold on;
      
   X   = [1:n]';    % number the observations
   null = zeros(n,1);
   
   for i = 1:noGrps % plot classes with different colors
      idx  = find((Y==grps(i))==1); % index to obs in this group
      
      plot([X(idx) X(idx)]',[null(idx) out(idx)]','LineStyle','-',...
         'LineWidth',1,'color',f_rgb(i));
   end
   
   xlabel('\bfObservation');
   ylabel('\bfMagnitude');
end

