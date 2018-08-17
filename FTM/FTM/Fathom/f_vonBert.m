function Lhat = f_vonBert(para,t)
% - von Bertalanffy growth function
%
% USAGE: Lhat = f_vonBert([La k t0],t);
%
% para = [La k t0];
% La   = theoretical asymptotic length
% k    = growth coefficient
% t0   = theoretical age when length = 0
% t    = age
%
% Lhat = predicted length at age t
%
% SEE ALSO: f_vonBertModel, f_vonBertAge

% -----Author:-----
% by David L. Jones, Jun-2006
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% force column vector
t  = t(:);

La = para(1);
k  = para(2);
t0 = para(3);

% Von Bertalanffy growth equation:
Lhat = La*(1 - exp(-k*(t-t0)));
