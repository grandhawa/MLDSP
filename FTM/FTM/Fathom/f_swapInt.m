function [x,y] = f_swapInt(x,y)
% - swap 2 integer values
% 
% USAGE: [xx,yy] = f_swapInt(x,y);
% 
% x,y   = input scalars
% xx,yy = swapped output

% -----References:-----
% C code to swap integers x and y x = x ^ y ; y = x ^ y ; x = x ^ y ;
% See: http://blogs.mathworks.com/loren/2006/10/25/cute-tricks-in-matlab-adapted-from
% -other-languages/

% -----Notes:-----
% This is a port of C-code macros from the randomForest package for R:
% #define swapInt(a, b) ((a ^= b), (b ^= a), (a ^= b))             // from reg_RF.cpp
% #define swapInt2(a, b) ((a = a ^ b), (b = b ^ a), (a = a ^ b))   // equivalent form

% -----Author:-----
% by David L. Jones, Jan-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

x = bitxor(x,y);
y = bitxor(x,y);
x = bitxor(x,y);
