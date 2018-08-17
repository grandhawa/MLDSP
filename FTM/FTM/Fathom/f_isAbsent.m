function sol = f_isAbsent(A,B)
% - determine which rows of A are absent from B
% 
% USAGE: sol = f_isAbsent(A,B)
% 
% A = column vector of positive integers
% B = column vector of positive integers
% 
% sol = boolean (logical) indicating 1 = true, 0 = false
% 
% SEE ALSO: f_isPresent, ismember, setdiff

% -----Notes:-----
% This function is an altermative to the Matlab 'ismember' command that may run
% faster if you're only working with two sets of positive integers (such as row
% or column indices).

% -----References:-----
% modified after 'fastintersect' by Jonathan Epperl (17-Oct-2012)
% http://www.mathworks.com/matlabcentral/answers/51102-ismember-function-too-slow

% -----Author:-----
% by David L. Jones, Jan-2014
%
% This file is part of the 'FATHOM Toolbox for Matlab' and
% is released under the GNU General Public License, version 2.

P    = nan(max([max(A);max(B)]),1); % initialize
P(A) = 1;                           % mark values present in A
P(B) = 0;                           % mark values absent in B 
sol  = P(A);                        % true if present in A & absent in B

% Use this to return the actual values:
% C = A(logical(P(A)));




