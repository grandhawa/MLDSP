function [evects,evals] = f_eig(x)
% - eigenanalysis of a square matrix
%
% USAGE: [evects,evals] = f_eig(x);
%
% x      = square input matrix
% evects = sorted eigenvectors
% evals  = sorted eigenvalues
%
% SEE ALSO: f_svd, f_pca, f_pcoa

% -----Notes:-----
% This function is used to perform an Eigenanalysis decomposition of a square
% input matrix. The Eigenvalues and corresponding Eigenvectors are sorted
% descending. The imaginary components of any complex values produced are
% discarded. 

% -----Author:-----
% by David L. Jones, Apr-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Mar-2010: removed unnecessary brackets

% Eigenanalysis:
[U,D] = eig(x); 
evals = diag(D);

% Discard imaginary components:
U     = real(U);
evals = real(evals);

% Sort by eigenvalues, descending:
sorted = flipud(sortrows([U' evals],size(U',1)+1));
evects = sorted(:,(1:end-1))';
evals  = sorted(:,end);
