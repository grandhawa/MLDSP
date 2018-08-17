function m_bubble(lon,lat,z,mx,sq,srt)
% - M_Map compatible version of f_bubble
%
% USAGE: m_bubble(lon,lat,z,mx,sq,srt);
%
% lon  = longitude
% lat  = latitude
% z    = symbol size at each x,y coordinate (pos values produce black symbols,
%        neg vaues produce white symbols)
% mx   = maximum symbol size
% mn   = minimum symbol size                                   (default = 0)
% sq   = scale symbols according to the sqrt of z              (default = 0)
% srt  = sort data so larger bubbles don't obscure smaller one (default = 0)
%
% SEE ALSO: f_bubble

% -----Notes:-----
% This function should be called WITHIN a script that sets the projection/grid
% of an m_map figure.
% 
% Zero or missing values (NaN's) are plotted as red X's.

% -----Author:-----
% by David L. Jones, Sep-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 5), sq  = 0; end; % default don't use sqrt of z
if (nargin < 6), srt = 0; end; % default don't sort data

if (isscalar(mx)==0)
   error('MX should be a scalar specifying maximum symbol size')
end
% -------------------------------------

if srt>0 % Sort the data descending:
   [nul,des] = sort(z,1,'descend');
   lon = lon(des);
   lat = lat(des);
   z   = z(des);
end

n = size(lon,1);

if sq>0
   idx    = z<0;     % get logical indices to neg numbers
   s      = sqrt(f_ranging(abs(z),2)) * mx;
   s(idx) = -s(idx); % make negative again
else
   s = f_ranging(z,2) * mx;
end

% -----This doesn't quite work with NEG values:-----
% if mn>0 % Set minimum symbol size
%    idx          = find(s<0 & abs(s)<mn); % get logical indices to neg numbers
%    s(abs(s)<mn) = mn;      % increase symbol size to minimum
%    s(idx)       = -s(idx); % make negative again
% end

hold on;


for i=1:n

   if s(i)<0 % Negative values:
      m_line(lon(i),lat(i),'marker','o','MarkerEdgeColor','k',...
         'MarkerFaceColor','w','MarkerSize',abs(s(i)));

   elseif isnan(s(i)) || s(i)==0 % Zero/Missing values:
      m_line(lon(i),lat(i),'marker','x','markersize', 4,...
         'MarkerFaceColor', 'none', 'MarkerEdgeColor','r','LineStyle','none');

   else % Positive values:
      m_line(lon(i),lat(i),'marker','o','MarkerEdgeColor','w',...
         'MarkerFaceColor','k','MarkerSize', abs(s(i)));
   end
end
