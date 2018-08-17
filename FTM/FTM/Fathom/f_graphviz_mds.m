function f_graphviz_mds(dis,fname,labels)
% - create a Graphviz DOT file to perform multidimensional scaling (MDS)
%
% USAGE: f_graphviz_mds(dis,'fname',labels);
%
% dis    = symmetric distance matrix
% fname  = filename to export to, WITHOUT the extension
% labels = cell array of labels; if empty, autocreate
%            e.g., labels = {'A' 'B' 'C' 'D' 'E' 'F'};
%
% SEE ALSO: f_graphviz_mst

% -----Notes:----
% You will need 'Graphviz' properly installed on your system to render the
% '*.dot' file created by this function. This application is cross-platform and
% freely available from: http://www.graphviz.org
% 
% If the text labels don't fit inside the nodes, you can manually edit the *.dot
% file and increase the size of the 'height' parameter.

% -----Author:-----
% by David L. Jones, Jun-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), labels = []; end % default set labels it empty

% Autocreate labels:
if isempty(labels)
   labels = f_num2cell(1:size(dis,1));
end

% Make sure labels are of compatible size:
if numel(labels) ~= size(dis,1)
   error('Size mismatch between DIS & LABELS!')
end
% -------------------------------------

% Set name of output file:
fnameOUT = [fname '.dot'];

% % Don't overwrite existing file:
% if exist(fnameOUT,'file')
%    error('Destination file already exists!')
% end

% Create DOT file:
header = 'graph G {\n';
footer = '}\n';

% -----Print results to an external file:-----
fid = fopen(fnameOUT,'w'); % open file for writing

% Print file creation comments:
fprintf(fid,'// ==================================================\n');
fprintf(fid,'// This file was created by ''f_graphviz_mds''\n');
fprintf(fid,'// which is part of the ''FATHOM TOOLBOX FOR MATLAB''\n');
fprintf(fid,'// by David L. Jones\n');
fprintf(fid,'// ==================================================\n\n');

% Header:
fprintf(fid,header);

% Options:
fprintf(fid,'\n');
fprintf(fid,'\tgraph [model=mds,center=true,overlap=scale,size=7];\n');
fprintf(fid,'\tnode [shape=circle,style=filled,height=0.3,fixedsize=true,fontname=Helvetica,fontsize=12];\n');
fprintf(fid,'\n');

n = size(dis,2); % size of distance matrix

% Get indices of lower tridiagonal elements:
[r,c] = find(triu(ones(n,n))==0);

% Unwrap distances, sort ascending:
[vec,key] = sort(f_unwrap(dis,0),'ascend');

% Get list of corresponding row/col indices:
list = [r(key) c(key)];

% Include these as INVISIBLE edges:
fprintf(fid,'// Invisible edges:\n');
n = size(list,1); % get # edges
for i = 1:n
   fprintf(fid,'\t%s -- %s [len=%s,style="invis"];\n',...
      labels{list(i,1)},labels{list(i,2)},num2str(vec(i)));
end
fprintf(fid,'\n');


% ---------------------------------------------

% Footer:
fprintf(fid,footer);

fclose(fid);
fprintf('\nDOT-formated data saved in file: %s\n\n',fnameOUT);
% --------------------------------------------