function R = f_randWH(n,plt,seed)
% - Wichmann & Hill's (2006) good pseudo-random number generator
%
% USAGE: R = f_randWH(n,plt,seed)
%
% n    = # random numbers required
% plt  = create a plot to visually assess uniformity               (default = 0)
% seed = optional vector of 4 random seeds,                (default = rand(4,1))
%        e.g., seed = [0.81472 0.90579 0.12699 0.91338] 
% 
% R = pseudo-random uniform numbers

% -----References:-----
% Wichmann, B. A. and I. D. Hill. 2006. Generating good pseudo-random numbers.
% Computational Statistics & Data Analysis 51:1614-1622.
%
% - modified (and corrected) after whrandom.m from:
% https://mailman.cae.wisc.edu/pipermail/help-octave/2012-June/052340.html

% -----Notes:-----
% The function returns uniformly distributed pseudo-random numbers using an
% algorithm that has passed the 'Big Crush' test. 

% -----Author:-----
% by David L. Jones, Mar-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 2), plt  = 0;  end % default don't create a plot
if (nargin < 3), seed = []; end % default random seed

% Check random seed:
if ~isempty(seed)
   if (numel(seed)~=4)
      error('SEED must specify only 4 values!');
   end
      
   if (any(seed>=1) || any(seed<=0))
      error('Values of SEED must range between 0-1!');
   end
else
   % Get default random seeds:
   % rng('shuffle'); % re-initialize random number generator
   seed = rand(4,1);
end
% -------------------------------------

% Parse seed: 
ix   = seed(1);
iy   = seed(2);
iz   = seed(3);
it   = seed(4);

% Define constants:
p1 = 2147483579;
p2 = 2147483543;
p3 = 2147483423;
p4 = 2147483123;

% Note: the loop is necessary to update each ix value, etc.
R  = zeros(n,1); % preallocate
for i=1:n
   ix   = mod(11600*ix,p1);
   iy   = mod(47003*iy,p2);
   iz   = mod(23000*iz,p3);
   it   = mod(33000*it,p4);
   R(i) = mod(ix/p1+iy/p2+iz/p3+it/p4,1);
end

% -----Create a plot:-----
if (plt>0)
   figure; set(gcf,'color','w');
   hist(R,100);
   h = findobj(gca,'Type','patch');
   set(h,'FaceColor','k','EdgeColor','none')
   title('Wichman & Hill Pseudo-Random Numbers')
   xlabel('Random variate');
   ylabel('Frequency');
end
