% Use Johnson Curve to calculate probability densities of Greenwood's Statistic
% by David L. Jones, Apr-2014

% NOTE: this example requires the JOHNSON CURVE TOOLBOX FOR MATLAB:
% http://www.marine.usf.edu/user/djones/jctm/jctm.html
% http://www.mathworks.com/matlabcentral/fileexchange/46123-johnson-curve-toolbox

% Clear the workspace:
clz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    N = 4                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Fit a Johnson Curve:
n    = 4; % sample size
para = f_greenwood_par(n); % parameters for a Johnson Curve
jsn  = f_johnson_M(para(1),para(2),para(3),para(4));

% Get cumuative probability densities:
G = 0.294941;
CD = f_johnson_cdf(G,jsn.coef,jsn.type)
% CD = 0.40921;

% Get percentage points:
Y = f_johnson_inv(0.4,jsn.coef,jsn.type)
% Y = 0.29337;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    N = 9                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Fit a Johnson Curve:
n    = 9;                  % sample size
para = f_greenwood_par(n); % parameters for a Johnson Curve
jsn  = f_johnson_M(para(1),para(2),para(3),para(4))
% jsn = 
% 
%     coef: [6.4076 1.8175 0.10399 2.3764]
%     type: 'SB'

% Get percentage points:
P = [0.01 0.05 0.1:0.1:0.9 0.95 0.99]';
G = f_johnson_inv(P,jsn.coef,'SB')
% 
% G =
% 
%       0.12328
%       0.13196
%       0.13805
%       0.14721
%       0.15528
%       0.16332
%       0.17194
%       0.18177
%       0.19381
%       0.21017
%       0.23761
%       0.26518
%       0.33149


% -> Hill's method failed to converge on an SB solution for n = 5-10, although
%    my code did fit an SB curve.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    N = 9                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Fit a Johnson Curve:
n    = 12;                 % sample size
para = f_greenwood_par(n); % parameters for a Johnson Curve
jsn  = f_johnson_M(para(1),para(2),para(3),para(4))
% jsn = 
% 
%     coef: [11.967 2.0284 0.081652 19.863]
%     type: 'SB'
% 
% - we can fit an SB curve up to n = 12, afterwards it fits an SU curve;
%   adjusting the default tolerance values might be useful here.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    N = 9                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Fit a Johnson Curve:
n    = 100;                % sample size
para = f_greenwood_par(n); % parameters for a Johnson Curve
jsn  = f_johnson_M(para(1),para(2),para(3),para(4))
% jsn = 
% 
%     coef: [-2.0192 2.485 0.016512 0.0031552]
%     type: 'SU'



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Reproduce Burrows's (1979) Table 1 using Johnson Curves:           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
P = [0.01 0.05 0.1:0.1:0.9 0.95 0.99]';
n = (4:10);
G = nan(numel(P),numel(n));
for i=1:numel(n)
   par    = f_greenwood_par(n(i));                    % get parameters
   jsn    = f_johnson_M(par(1),par(2),par(3),par(4)); % fit Johnson Curve
   G(:,i) = f_johnson_inv(P,jsn.coef,jsn.type);
end

[[NaN;P] [n;G]]

% NaN   4         5         6         7         8         9         10
% 0.01   0.21893   0.18827   0.16576   0.14841   0.13459   0.12328   0.11384
% 0.05   0.23305   0.20101   0.17724   0.15881   0.14406   0.13196   0.12183
%  0.1   0.24384   0.21049   0.18562   0.16628   0.15078   0.13805    0.1274
%  0.2   0.26104    0.2253   0.19853   0.17768   0.16095   0.14721   0.13571
%  0.3   0.27696   0.23878   0.21016   0.18786   0.16997   0.15528   0.14299
%  0.4   0.29337   0.25254   0.22193    0.1981   0.17899   0.16332   0.15022
%  0.5   0.31139   0.26754    0.2347   0.20915    0.1887   0.17194   0.15795
%  0.6   0.33232   0.28488   0.24939   0.22184    0.1998   0.18177   0.16674
%  0.7   0.35827   0.30631   0.26752   0.23744   0.21342   0.19381   0.17748
%  0.8   0.39371   0.33561   0.29227   0.25872   0.23198   0.21017   0.19205
%  0.9   0.45247   0.38457   0.33378   0.29444   0.26312   0.23761   0.21647
% 0.95   0.50943   0.43286    0.3751   0.33018   0.29435   0.26518   0.24099
% 0.99   0.63396   0.54284   0.47147   0.41479   0.36902   0.33149   0.30025



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Reproduce Currie's (1981) Table 1 using Johnson Curves:           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
P = [0.01 0.05 0.1:0.1:0.9 0.95 0.99]';
n = (11:15);
G = nan(numel(P),numel(n));
for i=1:numel(n)
   par    = f_greenwood_par(n(i));                    % get parameters
   jsn    = f_johnson_M(par(1),par(2),par(3),par(4)); % fit Johnson Curve
   G(:,i) = f_johnson_inv(P,jsn.coef,jsn.type);
end

[[NaN;P] [n;G]]
% 
%           NaN           11           12           13           14           15
%          0.01      0.10582      0.09892     0.092799     0.087231     0.082344
%          0.05      0.11322      0.10581     0.099372     0.093655     0.088592
%           0.1      0.11834      0.11054      0.10379      0.09784      0.09256
%           0.2      0.12594      0.11752      0.11024      0.10386     0.098194
%           0.3      0.13256      0.12358      0.11581        0.109      0.10297
%           0.4      0.13911      0.12956      0.12128      0.11403      0.10761
%           0.5       0.1461      0.13592      0.12708      0.11935      0.11252
%           0.6      0.15402      0.14312      0.13363      0.12535      0.11804
%           0.7      0.16369      0.15187       0.1416      0.13263      0.12473
%           0.8      0.17677      0.16371      0.15237      0.14247      0.13376
%           0.9      0.19867      0.18348      0.17035      0.15889      0.14884
%          0.95      0.22066      0.20333      0.18841      0.17541      0.16401
%          0.99      0.27394      0.25151      0.23234      0.21564        0.201



P = [0.01 0.05 0.1:0.1:0.9 0.95 0.99]';
n = (16:20);
G = nan(numel(P),numel(n));
for i=1:numel(n)
   par    = f_greenwood_par(n(i));                    % get parameters
   jsn    = f_johnson_M(par(1),par(2),par(3),par(4)); % fit Johnson Curve
   G(:,i) = f_johnson_inv(P,jsn.coef,jsn.type);
end

[[NaN;P] [n;G]]
% 
%           NaN           16           17           18           19           20
%          0.01     0.078019      0.07416     0.070695     0.067564     0.064721
%          0.05     0.084075     0.080017     0.076351     0.073021     0.069983
%           0.1     0.087841     0.083596     0.079756     0.076266     0.073078
%           0.2      0.09313     0.088574     0.084453     0.080708     0.077288
%           0.3     0.097576     0.092732     0.088354     0.084378      0.08075
%           0.4      0.10189     0.096753     0.092114     0.087905     0.084068
%           0.5      0.10643      0.10097     0.096051     0.091591     0.087529
%           0.6      0.11153      0.10571      0.10046     0.095716     0.091397
%           0.7      0.11771      0.11144      0.10579      0.10069     0.096059
%           0.8      0.12604      0.11915      0.11297      0.10739      0.10232
%           0.9      0.13995      0.13203      0.12493      0.11854      0.11275
%          0.95      0.15393      0.14498      0.13696      0.12975      0.12323
%          0.99      0.18808       0.1766      0.16634      0.15713      0.14881





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Use Johnson Curves to estimate probabilities of Greenwood's Statistic:    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% -> n = 4-12 are Johnson SB curves:
% 
% Set range of sample sizes:
n        = (4:12)';
nr       = numel(n);
coef     = nan(nr,4); % preallocate
type{nr} = '';

for i=1:nr
   par       = f_greenwood_par(n(i));                    % get parameters
   jsn       = f_johnson_M(par(1),par(2),par(3),par(4)); % fit Johnson Curve
   coef(i,:) = jsn.coef;                                 % collect coefficients
   type{i}   = jsn.type;                                 % collect Johnson type
end

% Use these SB coefficients for n = 4-12:
coef_SB = [
   2.857488476105308   1.334637520818872   0.197015971011635   1.087464532474250
   3.522093754720662   1.460931647408669   0.166210944052008   1.230475737966411
   4.177842384852712   1.565877316633673   0.144185681219991   1.394959462944271
   4.852761052854307   1.657526273501108   0.127567669618079   1.605984235748621
   5.579920434113547   1.740461420219146   0.114529013281242   1.904677541939498
   6.407643210354968   1.817481029197412   0.103991339524473   2.376373947416354
   7.425757210809234   1.890390822910679   0.095275447176263   3.247479159509095
   8.864862720200591   1.960410397327272   0.087933057769960   5.410180450026539
  11.967478262704688   2.028395819139190   0.081652024553840  19.862931573264103];



% -> n = 13+ are Johnson SU curves:
% 
% Set range of sample sizes:
tic;
n        = (13:1000)';
nr       = numel(n);
coef     = nan(nr,4); % preallocate
type{nr} = '';

for i=1:nr
   par       = f_greenwood_par(n(i));                    % get parameters
   jsn       = f_johnson_M(par(1),par(2),par(3),par(4)); % fit Johnson Curve
   coef(i,:) = jsn.coef;                                 % collect coefficients
   type{i}   = jsn.type;                                 % collect Johnson type
end
toc;
[n coef]

% Elapsed time is 0.155006 seconds.
% -> this took a fraction of a second to fit nearly 1000 Johnson Curves, so
% there's little reason not to just calculate these on the fly.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Plot PDF of Greenwood's Statistic:                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
n = 1; 
f_greenwood_plt(n,'PDF');
axisVar = axis;
axis([0 axisVar(2:4)]);
f_pdf('PDF_Greenwood_Stat_01');

% 
n = 2; 
f_greenwood_plt(n,'PDF');
axisVar = axis;
axis([0 axisVar(2:4)]);
f_pdf('PDF_Greenwood_Stat_02');

% 
n = 3:20;
figure; set(gcf,'color','w');
hold on;
for i=1:numel(n)
   [G,den] = f_greenwood_plt(n(i),'PDF',0);
   plot(G,den,'k-');
end
title('PDF');
xlabel('Greenwood''s Statistic');
ylabel('Probability Density');
xTxt = sprintf('  (n=%d-%d)',min(n),max(n));
xlabel(['Greenwood''s Statistic' xTxt]);
grid on; box on;
clear G den;
f_pdf('PDF_Greenwood_Stat_03-20');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Plot CDF of Greenwood's Statistic:                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
n = 1:20;
figure; set(gcf,'color','w');
hold on;
for i=1:numel(n)
   [G,den] = f_greenwood_plt(n(i),'CDF',0);
   plot(G,den,'k-');
end
title('CDF');
xlabel('Greenwood''s Statistic');
ylabel('Cumulative Probability Density');
xTxt = sprintf('  (n=%d-%d)',min(n),max(n));
xlabel(['Greenwood''s Statistic' xTxt]);
grid on; box on;
axis([0 1 0 1]);
f_pdf('CDF_Greenwood_Stat_01-20');
