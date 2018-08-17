function hdl = f_plotError(x,y,e,type,width)
% - add error bars to current figure
%
% USAGE: hdl = f_plotError(x,y,e,type,width);
%
% x,y  = coordinates of midpoint of error bar
% e    = error which is value to be added/subtracted
%       to x,y midpoint
% type  = vertical (= 1, default) or horizontal (= 2) error bars
% width = width of error bar (default = 1)
%
% hdl = figure handle to errorbars

% -----Notes:-----
% This function is used to add error bars to the current figure by simply
% plotting lines separated by NaN's. If E is a column vector then
% SYMMETRICAL error bars are produced; NON-SYMMETRICAL error bars, such as
% those specifying min,max values, etc. are produced if E is a 2 column
% matrix specifying [lower upper].

% -----Author:-----
% by David L. Jones, May-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% June-2007: added support for horizontal error bars
% Dec-2012: updated documentation; now plots horizontal error bars correctly 

% -----Check input:-----
if (nargin <4),type  = 1; end; % default vertical error bars
if (nargin <5),width = 1; end; % default width = 1

if (size(e,2)<2)
   e = [e e];
end

n = size(x,1);

if (n ~= size(y,1)) || (n ~= size(e,1))
   error('X,Y, & E must have same # of rows !')
end
% ----------------------


% -----Create error bars:-----
idx = 0; % initialize counter

switch type
   case 1 % Vertical bars
      for i=1:n
         idx     = idx + 1;
         xx(idx) = x(i);
         yy(idx) = y(i) - e(i,1);
         
         idx     = idx + 1;
         xx(idx) = x(i);
         yy(idx) = y(i) + e(i,2);
         
         idx     = idx + 1;
         xx(idx) = NaN;
         yy(idx) = NaN;
      end
      
   case 2 % Horizontal bars
      for i=1:n
         idx     = idx + 1;
         xx(idx) = x(i) - e(i,1);
         yy(idx) = y(i);
         
         idx     = idx + 1;
         xx(idx) = x(i) + e(i,2);
         yy(idx) = y(i);
         
         idx     = idx + 1;
         xx(idx) = NaN;
         yy(idx) = NaN;
      end
   otherwise
      error('TYPE must be 1 or 2!')
end
% -----Plot error bars:-----
hold on;
hdl = plot(xx,yy,'-k','LineWidth',width,'Color',[1 1 1]*0.70);



