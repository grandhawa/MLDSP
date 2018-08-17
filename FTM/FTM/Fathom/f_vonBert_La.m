function La = f_vonBert_La(Lmax)
% - calculate La parameter for the von Bertalanffy growth function
%
% USAGE: La = f_vonBert_La(Lmax);
%
% Lmax = maximum observed total length (in millimeters) 
% La   = theoretical asymptotic length
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
Lmax = Lmax/10;

log_La = 0.044 + 0.9841*log10(Lmax);
La     = 10^log_La; % inverse log

% Convert lengths back to mm:
La = La*10;
