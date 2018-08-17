function age = f_vonBertAge(length,La,k,t0,plt)
% - age a fish of known length with the von Bertalanffy growth equation
%
% USAGE: age = f_vonBertAge(length,La,k,t0,plt);
%
% length = length of fish
% La     = theoretical asymptotic length
% k      = growth coefficient
% t0     = theoretical age when length = 0
% plt    = create plot (default = 0)
%
% age = estimated ages
%
% SEE ALSO: f_vonBert

% -----References:-----
% http://www.fao.org/docrep/X5685E/x5685e03.htm

% -----Author:-----
% by David L. Jones, Jun-2006
% 
% This file is part of the FATHOM Toolbox for Matlab and is
% released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 5), plt =  0; end; % no plot by default
% -------------------------------------

age = (1/k)*log(La./(La-length))+ t0; % eq. 3.8

% Make plot:
if (plt==1)
	figure;
   set(gcf,'color','w'); % set background color to white
	hold on;
	ta   = f_vonBertAge(La-0.01,La,k,t0,0); % get age @ length La
	yFit = f_vonBert([La k t0],t0:ta/100:ta);
	plot(age,length,'bo');
	plot(t0:ta/100:ta,yFit,'r-');
	xlabel('Age');
	ylabel('Length');
	title('von Bertalanffy Growth Model');
	grid on;
end



