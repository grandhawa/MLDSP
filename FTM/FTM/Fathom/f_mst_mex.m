function  [len,con] = f_mst_mex(dis,method)
% - create a minimal spanning tree (MEX version)
%
% USAGE: [len,con] = f_mst_mex(dis);
%
% dis    = symmetric distance matrix
% method = 1: use MatlabBGL Toolbox      (default)
%          2: use Bioinformatics Toolbox
% 
% len  = length of edges
% con  = corresponding pairs of objects connected by triangle edges
% 
% SEE ALSO: f_mst

% -----Notes:-----
% This function creates a minimal spanning tree from a symmetric dissimilarity
% matrix. Since it uses sparse matrices and calls an external Matlab MEX
% routine, it is faster than the F_MST function included in the Fathom Toolbox.
% It calls either the MST routine in David Gleich's 'MatlabBGL: a Matlab Graph
% Library' or the GRAPHMINSPANTREE routine in the Matlab Bioinformatics Toolbox.
% 
% The former is available from: 
% https://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/
% 
% A compiled version of MatlabBGL for 64-bit OSX is available from:
% http://dgleich.wordpress.com/2010/07/08/matlabbgl-osx-64-bit/
% http://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/old/matlab_bgl_4.0_osx64.zip

% -----Author:-----
% by David L. Jones, May-2014
%
% This file is part of the 'FATHOM Toolbox for Matlab' and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 2), method = 1; end % default don't send output to display
% -------------------------------------

% Create MST using Kruskal's algorithm:
switch method
   case 1 % MatlabBGL Toolbox:
      if (nargout > 1)
         [i,j,len] = find(triu(mst(sparse(dis))));
         con       = [i j];
      else
         [~,~,len] = find(triu(mst(sparse(dis))));
      end
      
   case 2 % Bioinformatics Toolbox:
      if (nargout > 1)
         [i,j,len] = find(graphminspantree(sparse(dis),'Method', 'Kruskal'));
         con       = [i j];
      else
         [~,~,len] = find(graphminspantree(sparse(dis),'Method', 'Kruskal'));
      end
   otherwise
      error('Unknown method!')
end
