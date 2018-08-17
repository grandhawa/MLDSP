function idx = f_brush2idx
% - create index to rows of data currently selected by 'brushing'
% 
% USAGE: idx = f_brush2idx
% 
% idx = index to rows of brushed data
% 
% SEE ALSO: brush

% -----Author:-----
% by David L. Jones, Jul-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

try
   idx = find(get(get(gca,'Children'),'BrushData')==1);
   idx = idx';
catch
   error('Please BRUSH some data points first!')
end



