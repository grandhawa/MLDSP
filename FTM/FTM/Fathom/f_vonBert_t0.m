function t0 = f_vonBert_t0(La,k)
% - calculate t0 parameter for the von Bertalanffy growth function
%
% USAGE: t0 = f_vonBert_t0(La,k);
%
% La   = theoretical asymptotic length (in millimeters)
% k    = growth coefficient
% 
% t0   = theoretical age when length = 0
%
% SEE ALSO: f_vonBert

% -----References:-----
% Pauly, 1980 (in Martinez-Andrade, 2003)
 
% -----Author:-----
% by David L. Jones, Jul-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% convert lengths from mm to cm:
La = La/10;

log_t0 = -0.3922 - 0.2752*log10(La) - 1.038*log10(k);
t0     = -(10^(log_t0)); % inverse log
