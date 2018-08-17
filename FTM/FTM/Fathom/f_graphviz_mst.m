function f_graphviz_mst(dis,fname,method)
% - export a minimal spanning tree to Graphviz DOT format
%
% USAGE: f_graphviz_mst(dis,'fname');
%
% dis   = symmetric distance matrix
% fname = filename to export to, WITHOUT the extension
% method = as 'all' or 'mst'
%
% SEE ALSO: f_graphviz_mds, f_graphviz_neato

% -----Notes:----
% You will need 'Graphviz' properly installed on your system to render the
% '*.dot' file created by this function. This application is cross-platform and
% freely available from: http://www.graphviz.org
%
% Invisible edges are optionally included to get a more accurate depiction of
% the distance between nodes that are not connected by the MST.

% -----Author:-----
% by David L. Jones, Jun-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), method = 'all'; end % default include all pair-wise distances

% Force lower case:
method = lower(method);
% -------------------------------------

% Set name of output file:
fnameOUT = [fname '.dot'];

% % Don't overwrite existing file:
% if exist(fnameOUT,'file')
%    error('Destination file already exists!')
% end

% Create Minimal Spanning Tree:
[len,con] = f_mst_mex(dis);

% Create DOT file:
header = 'graph G {\n';
footer = '}\n';

% -----Print results to an external file:-----
fid = fopen(fnameOUT,'w'); % open file for writing

% Print file creation comments:
fprintf(fid,'// ==================================================\n');
fprintf(fid,'// This file was created by ''f_graphviz_mst''\n');
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

% Edges of MST to plot:
fprintf(fid,'// Edges to plot:\n');
n = size(con,1); % get # connections
for i = 1:n
   fprintf(fid,'\t%s -- %s [len=%s];\n',num2str(con(i,1)),num2str(con(i,2)),...
      num2str(len(i)));
end
fprintf(fid,'\n');

% -----Optionally include invisible edges:-----
switch method
   case 'all'
      n = size(dis,2); % size of distance matrix
      
      % Get indices of lower tridiagonal elements:
      [r,c] = find(triu(ones(n,n))==0);
      
      % Unwrap distances, sort ascending:
      [vec,key] = sort(f_unwrap(dis,0),'ascend');
      
      % Get list of corresponding row/col indices:
      list = [r(key) c(key)];
      
      % Get list of edges NOT part of MST (swap cols of list to match MST):
      temp = setdiff([list(:,2) list(:,1) vec],[con len],'rows');
      list = temp(:,1:2);
      vec  = temp(:,3);
      
      % Include these as INVISIBLE edges:
      fprintf(fid,'// Invisible edges:\n');
      n = size(list,1); % get # edges
      for i = 1:n
         fprintf(fid,'\t%s -- %s [len=%s,style="invis"];\n',num2str(list(i,1)),num2str(list(i,2)),...
            num2str(vec(i)));
      end
      fprintf(fid,'\n');
   case 'mst'
      % do nothing
   otherwise
      error('Method must be ''all'' or ''mst''!')
end

% ---------------------------------------------

% Footer:
fprintf(fid,footer);

fclose(fid);
fprintf('\nDOT-formated data saved in file: %s\n\n',fnameOUT);
% --------------------------------------------