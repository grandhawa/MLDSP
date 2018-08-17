function prob = f_permtest(x,y,iter,test,plotflag);
% - two sample permutation test of means
%
% USAGE: prob = f_permtest(x,y,{iter},{test},{plotflag});
%
% x,y      = column vectors of input data
% iter     = # of randomization iterations (default = 1000)
% test     = 0: (default) two-tailed test whether mean of x = y 
%            1: one-tailed test whether mean of x > y 
% plotflag = show histogram plot of randomized probabilities (default = 0)
%
% prob     = randomized probability that null Ho: is true
%
% SEE ALSO: f_npManova

% -----References:-----
% Manly, B. F. J. 1997. Randomization, Bootstrap and Monte Carlo Methods in Biology.
% Second edition. Chapman and Hall, New York. 399 pp.
%
% Resampling Stats. http://www.resample.com

% -----Author:-----
% by David L. Jones, Oct-2001

% 14-Mar-2002: corrected formula for p-value

if (nargin < 3), iter     = 1000; end; % set default # of resampling iterations
if (nargin < 4), test     = 0; end; % default testing of equal means
if (nargin < 5), plotflag = 0; end; % by default don't plot histogram

x = x(:); y = y(:); % force into column vectors

meanX = mean(x); meanY = mean(y);
sizeX = size(x,1); sizeY = size(y,1);

if test >0
   tStat = (meanX - meanY); % test statistic for Ho: not mean(x) > mean(y)
else
   tStat = abs(meanX - meanY); % test statistic for Ho: mean(x) = mean(y)
end

z = [x' y']'; % combine into 1 variable
sizeZ = size(z,1);

randStat = zeros(iter,1); % preallocate results array

%-----Permutation Test:-----
rand('state',sum(100*clock)); % set random generator to new state
if test > 0 %% One-Tailed Test
   for i = 1:(iter-1) % resampling without replacement
      zz = f_shuffle(z);    % randomly permute the order of z
      xx = zz(1:sizeX);     % extract xx
      yy = zz(sizeX+1:end); % extract yy
      randStat(i) = (mean(xx) - mean(yy)); % keep list of randomized statistic
   end;
else %% Two-Tailed Test
   for i = 1:(iter-1) % resampling without replacement
      zz = f_shuffle(z);    % randomly permute the order of z
      xx = zz(1:sizeX);     % extract xx
      yy = zz(sizeX+1:end); % extract yy
      randStat(i) = abs(mean(xx) - mean(yy)); % keep list of randomized statistic
   end;
end;
j    = find(randStat >= tStat); % get randomized values >= to test statistic
prob = (length(j)+1)./iter;     % count those vales & convert to probability 
%-------------------------------

%-----Plot histogram:-----
if plotflag > 0 
   figure;
   [n,xout] = hist(randStat); % create bins
   bar(xout,(n/iter)*100);
   xlabel('Value'); 
   ylabel('Freqency (%)');
   h = findobj(gca,'Type','patch');
   set(h,'FaceColor',[0 0 0.75],'EdgeColor','b')
   hold on;
   plot(tStat,0,'r*'); % plot test statistic from input data
   hold off;
   title(['Randomized probability = ' num2str(prob) ' (' num2str(iter) ' iterations)']);
end;

