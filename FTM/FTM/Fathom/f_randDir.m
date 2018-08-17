function  R = f_rndDir(A,n)
% - draw random numbers from a Dirichlet distribution
% 
% USAGE: R = f_rndDir(A,n);
% 
% A = alpha parameters of a Dirichlet distribution
% n = # random numbers to return
% 
% R = random numbers sampled from the specified Dirichlet distribution

% -----References:-----
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/47175

% -----Author:-----
% by David L. Jones, 
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

A  = A(:)';                            % force as a column vector
nc = numel(A);                         % get # of alpha parameters
R  = gamrnd(repmat(A,n,1),ones(n,nc)); % sample from a gamma distribution
R  = R ./ repmat(sum(R,2),1,nc);       % normalize so rows sum to 1
