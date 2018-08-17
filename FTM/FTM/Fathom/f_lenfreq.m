function res = f_lenfreq(y,x,plotflag)
% - length-frequency plot (adjusted for sampling effort)
%
% Usage: res = f_lenfreq(y,x,{plotflag});
% ----------------------
% y  = input matrix; col 1 = fish_length, col 2 = vol_filtered (effort)
% x  = vector of bins for frequency histogram
% plotflag = optionally create plot (default = 1)
%
% res = row vector of length frequencies
%
% Plots length-frequencies from MOCNESS samples, adjusted
% for sampling effort using 'vol_filtered" as a weighting term

% -----Author:-----
% by David L. Jones, Oct-2001
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

if (nargin < 3), plotflag = 1; end; % make plot by default

y = sortrows(y,1); % sort input data by lengths

w    = y(:,2);    % extract vol_filtered (weights)
maxW = max(w);    % get maximum vol_filtered
ww   = maxW ./ w; % divide each by maximum_vol

y       = y(:,1);    % extract lengths
n       = hist(y,x); % get bin counts
noBins  = size(n,2); % get number of bins

res     = zeros(1,noBins); % initialize row vector
counter = 1;               % initialize counter

for i = 1:noBins
   if (n(i) ~= 0)
      res(i) = sum(ww(counter:(counter + (n(i)-1))));
      counter = counter + n(i); % increment counter
   end;
end;

res = (res ./ sum(res)) * 100; % convert to percent frequency

%-----Plot histogram:-----
if plotflag>0
   bar(x,res,0.9);
   xlabel('Length (mm)');
   ylabel('Frequency (%)');
   h = findobj(gca,'Type','patch');
   set(h,'FaceColor',[0 0 0],'EdgeColor','k');
   
   set(gca,'FontSize', 8);
%  set(gca,'TickDir','out');
end;

