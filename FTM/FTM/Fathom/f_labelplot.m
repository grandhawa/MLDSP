function f_labelplot(crds,labels,fntcolor,fntsize)
% - create 2-d label plot
%
% Usage: f_labelplot(crds,labels,{fntcolor},{fntsize});
%
% ----- Input: -----
% crds   = matrix of 2-d coordinates (rows = sites; cols = dimensions)
%
% labels = cell array of vector labels (if empty, autocreate)
%           e.g.,labels = {'sal' 'tmp' 'elev'};
%
% fntcolor = optional font color (default = 'b');
% fntsize  = optional font size  (default = 8);

% -----Author:-----
% by David L. Jones, Aug-2001
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Mar-2003: changed figure properties so labels can be plotted
%           over existing plot

% ----- Check input and set default values: -----
if (size(crds,1) == size(labels,1))<1
   error('The # of rows (sites) in LABELS & CRDS must be equal!');
end;

if (nargin < 2)
   error('At least 2 input variables are required!');
end;

if (nargin < 3), fntcolor = 'b'; end;
if (nargin < 4), fntsize  = 8; end;

% if labels are not cell arrays, try forcing them:
if iscell(labels)<1, labels = num2cell(labels); end;


nloop = size(crds,1); % get number of points

% figure;
hold on;

for i = 1:nloop
   plot(crds(i,1),crds(i,2),'w.');
   h = text(crds(i,1),crds(i,2),labels(i)); % plot labels
   set(h,'HorizontalAlignment','center','FontSize',fntsize,'Color',fntcolor);
end;

% axis equal;
hold off;
% box on;