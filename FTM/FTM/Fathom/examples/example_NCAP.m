% Example of Nonlinear Canonical Analysis of Principal Coordinates (NCAP)
% 
% by David L. Jones, Feb-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/holdfasts.mat'
% These data consist of numbers of 351 taxa of marine invertebrates found on the
% holdfasts of 80 subtidal kelp off the NE coast of New Zealand. From: Anderson,
% M. J., S. D. Connell, B. M. Gillanders, C. E. Diebel, W. M. Blom, J. E.
% Saunders, & T. J. Landers. 2005. Relationships between taxonomic resolution
% and spatial scales of multivariate variation. Journal of Animal Ecology 74:
% 636-646.

% Variables:
% bio     = biotic data (80 obs, 351 species)
% env     = environmental data (80 obs, 3 variables)
% env_txt = corresponding variable labels {'density' 'depth' 'Volume'}

% Load data:
load holdfasts.mat

% Forth-root transformed response data:
bio_4 = f_normal(bio,'4');

% Setup explanatory variable (= Volume):
X     = env(:,3);
X_txt = env_txt(3);

% Find optimal value of m:
resR2 = f_ncapOptimal(bio_4,'bc',X,'von','R2',25,1);
% 
% > Evaluating values of m from 1 to 25...
% > The optimal value of m (between 1 & 25) = 7

% Evaluate m according to RDA statistic:
resRDA = f_ncapOptimal(bio_4,'bc',X,'von','RDA',25,1);

% Perform NCAP using optimal m:
m      = 5;
result = f_ncap(bio_4,'bc',X,'von','R2',1,m);
% 
% -----Finding nonlinear coefficients for m = 5:-----
% stat = 0.676828185   SSres =        0.104440022    b = 0.002018410
% stat = 0.676841648   SSres =        0.104431320    b = 0.002019410
% stat = 0.682874601   SSres =        0.100568519    b = 0.026023693
% stat = 0.682869078   SSres =        0.100572022    b = 0.026024693
% stat = 0.688187748   SSres =        0.097226880    b = 0.025055568
% stat = 0.688182299   SSres =        0.097230279    b = 0.025056568
% stat = 0.698350120   SSres =        0.090992650    b = 0.023155676
% stat = 0.698344889   SSres =        0.090995806    b = 0.023156676
% stat = 0.716190195   SSres =        0.080548005    b = 0.019502795
% stat = 0.716185753   SSres =        0.080550527    b = 0.019503795
% stat = 0.735968278   SSres =        0.069712750    b = 0.012798004
% stat = 0.735967396   SSres =        0.069713216    b = 0.012799004
% stat = 0.725914638   SSres =        0.075122786    b = 0.007366652
% stat = 0.725919747   SSres =        0.075119985    b = 0.007367652
% stat = 0.736371056   SSres =        0.069500220    b = 0.011440166
% stat = 0.736371376   SSres =        0.069500051    b = 0.011441166
% stat = 0.696842344   SSres =        0.091904564    b = 0.023442807
% stat = 0.696837073   SSres =        0.091907760    b = 0.023443807
% stat = 0.733491386   SSres =        0.071026841    b = 0.014440826
% stat = 0.733489295   SSres =        0.071027956    b = 0.014441826
% stat = 0.736349696   SSres =        0.069511483    b = 0.012190331
% stat = 0.736349328   SSres =        0.069511677    b = 0.012191331
% stat = 0.736414423   SSres =        0.069477356    b = 0.011627707
% stat = 0.736414565   SSres =        0.069477281    b = 0.011628707
% stat = 0.736424855   SSres =        0.069471857    b = 0.011815248
% stat = 0.736424824   SSres =        0.069471873    b = 0.011816248
% stat = 0.736425295   SSres =        0.069471625    b = 0.011794863
% stat = 0.736425283   SSres =        0.069471632    b = 0.011795863
% stat = 0.736425354   SSres =        0.069471594    b = 0.011774495
% stat = 0.736425360   SSres =        0.069471591    b = 0.011775495
% stat = 0.736424648   SSres =        0.069471966    b = 0.011821380
% stat = 0.736424611   SSres =        0.069471985    b = 0.011822380
% stat = 0.736425367   SSres =        0.069471587    b = 0.011786216
% stat = 0.736425362   SSres =        0.069471590    b = 0.011787216
% stat = 0.736425368   SSres =        0.069471587    b = 0.011785898
% stat = 0.736425364   SSres =        0.069471589    b = 0.011786898
% 
% First 44 axes of Q explains 99.79% of Y

% Plot gradient:
XX    = [0:300]';
g     = 1 - exp(-XX*result.b); % fitted gradient
model = f_mregress(result.Qm,1 - exp(-X*result.b),0,0,0);
figure;
set(gcf,'color','w'); % set background color to white
hold on;
plot(XX,g,'k-'); % fitted gradient
plot(X,model.yfit,'bo');
xlabel('X');
ylabel('Gradient');
box on;


% Get bootstrapped 95% CI's:
boot = f_ncap(bio_4,'bc',X,'von','R2',0,m,0,200)
% Bootstrapping the data 200 times...
% 
% boot = 
% 
%          Q: [80x44 double]
%     Qevals: [44x1 double]
%      Qexpl: [44x2 double]
%         Qm: [80x5 double]
%          m: 5
%          b: 0.011786
%       stat: 0.73643
%       type: 'R2'
%        AIC: -446.42
%     lin_R2: 0.64643
%      lin_b: 0.0020184
%          F: 0.34146
%          p: [1x1 struct]
%         CI: [0.010002 0.033574]
