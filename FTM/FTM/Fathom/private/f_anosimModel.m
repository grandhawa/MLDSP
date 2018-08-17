function model = f_anosimModel(grps)
% - utility function called by f_anosim

% This function creates a specially scaled
% Model Matix to perform ANOSIM via a Mantel Test
%
% grps = row vector designating group membership
%        e.g., grps = [1 1 2 2 3 3 2];

% -----References:-----
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam. pp.552;561-562

% -----Author:-----
% by David L. Jones, Mar-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

noObj = length(grps);
dist  = f_euclid(grps); 

z = logical(dist==0) .* logical(eye(noObj) == 0);
z = f_unwrap(z);

noZeros = length(find(z==0)); % # of between group values
noOnes  = length(find(z==1)); % # of within group values

betweenVar = (1/noZeros)/(noObj*(noObj-1)/4);
withinVar  = -1*(1/noOnes)/(noObj*(noObj-1)/4);

[r1,c1] = find(z==0); % indices for "between"
[r2,c2] = find(z==1); % indices for "within"

z(r1,c1) = betweenVar;
z(r2,c2) = withinVar;

model = f_rewrap(z); % actual model matrix
