function f_bubble(x,y,z,mx,sq,srt)
% - create B&W bubble plots
%
% USAGE: f_bubble(x,y,z,mx,sq,srt)
% 
% x   = x-coordinates 
% y   = y-coordinates 
% z   = symbol size at each x,y coordinate (pos values produce black symbols, neg
%       vaues produce white symbols)
% mx  = maximum symbol size
% sq  = scale symbols according to the sqrt of z              (default = 0)
% srt = sort data so larger bubbles don't obscure smaller one (default = 0)
% 
% SEE ALSO: m_bubble

% -----Author:-----
% by David L. Jones, Apr-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Sep-2008: added 'srt' option

% -----Check input & set defaults:-----
if (nargin < 5), sq  = 0; end; % default don't use sqrt of z
if (nargin < 6), srt = 0; end; % default don't sort data

if (isscalar(mx)==0)
   error('MX should be a scalar specifying maximum symbol size')
end
% -------------------------------------

if srt>0 % Sort the data descending:
   [nul,des] = sort(z,1,'descend');
   x = x(des);
   y = y(des);
   z = z(des);
end

n = size(x,1);

if sq>0
   idx    = z<0;     % get logical indices to neg numbers
   s      = sqrt(f_ranging(abs(z),2)) * mx;
   s(idx) = -s(idx); % make negative again
else
   s = f_ranging(z,2) * mx;
end

hold on;

for i=1:n

   if s(i)<0 % Negative values:
      plot(x(i),y(i),'o','MarkerEdgeColor','k','MarkerFaceColor','w',...
         'MarkerSize', abs(s(i)));
   else      % Positive values:
      plot(x(i),y(i),'o','MarkerEdgeColor','w','MarkerFaceColor','k',...
         'MarkerSize', abs(s(i)));
   end
end


