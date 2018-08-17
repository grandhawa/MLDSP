function res = f_pnnSm(X,grp,vec,iter,plt)
% - determine optimal smoothing factor for a PNN
%
% USAGE: res = f_pnnSm(X,grp,vec,iter,plt);
%
% X    = matrix of training data (rows = obs, cols = variables) 
% grp  = vector of integers specifying group membership of X
% vec  = vector containing a range of values for Smoothing Factor
% iter = # iterations for bootstrap resampling        (default = 50)
% plt  = optionally plot results                      (default = 1)
%
% res  = [smoothing_factor error_rate]
%
% SEE ALSO: f_pnn, f_pnnCV, f_pnn632

% -----Notes:-----
% This function is used to determine the optimal smoothing factor to use
% for a Probabilistic Neural Network (PNN).

% -----Dependencies:-----
% This program requires the Matlab NNET Toolbox.

% -----References:-----
% Matlab NNET Toolbox reference &
% news://comp.ai.neural-nets

% -----Author:-----
% by David L. Jones, Feb-2004
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin<4), iter  = 50;  end; % default 50 bootstraps
if (nargin<5), plt   = 1; end;   % create plot by default

grp = grp(:); % force column vector

if size(X,1) ~= size(grp,1)
   error('X & GRP must have same # of rows!');
end
% -------------------------------------

vec   = vec(:);        % force column vector
nVec  = size(vec,1);   % get # of spread values
res   = zeros(nVec,2); % initialize

% Notify user:
if plt>0   
   fprintf('\nCalculating PNN Error Rate for %d smoothing factors... \n',nVec);
end

% -----Compute error rate:-----
for i=1:nVec
   fprintf(' -> examining factor %d of %d \n',i,nVec);
   temp      = f_pnn632(X,grp,vec(i),iter,0);
   res(i,1)  = vec(i);
   res(i,2)  = temp.err_632plus;
end
% -----------------------------------

% Select optimal smoothing value:
idx = find(res(:,2) == min(res(:,2)));
sm  = res(idx(1),1); % if tied, only take 1 sm value

% Plot output:
if (plt>0)
   figure;
   plot(res(:,1),res(:,2),'-bo');
   
   title('\bfPNN Optimal Smoothing Factor');
   smFixed = sprintf('%2.2f',sm);
   xText   = ['Smoothing Factor (optimum = ' num2Str(smFixed) ')'];
   xlabel(xText);
   ylabel('Error Rate');
   
   % Set figure options:
   box on;
   grid on;   
   fprintf('\n Optimal smoothing factor = %2.2f \n\n',sm);
end



