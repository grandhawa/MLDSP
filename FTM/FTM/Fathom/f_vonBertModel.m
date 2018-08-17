function result = f_vonBertModel(t,L,Lm,La_i,k_i,t0_i,iter,tol)
% - fit a von Bertalanffy growth model using nonlinear regression 
%
% USAGE: result = f_vonBertModel(age,length,Lm,La_i,k_i,t0_i,iter,tol)
%
% age    = observed t
% length = observed L
% Lm     = maximum total length (cm) reported for this species  (default = 0)
% La_i   = initial estimate of La, set to '0' to use the default; if Lm is
%          provided then La_i is derived from it empirically. 
% k_i    = initial estimate of k, set to '0' to use the default (= random number)
% t0_i   = initial estimate of t0,set to '0' to use the default (= eps) 
% iter   = number of iterations for nonlinear fitting           (default = 100)
% tol    = tolerance to terminate fit of estimated parameters   (default = 1e-8)
%
% result = structure with the following fields:
% .La  = theoretical asymptotic length
% .k   = growth coefficient
% .t0  = theoretical age at length = 0
% .ta  = theoretical age at length = La
% .L0  = theoretical length at birth
% .obs = observed lengths
% .fit = fitted lengths
% .res = residuals
% .R2  = R-squared
%
% SEE ALSO: f_vonBert, f_vonBertAge

% -----References:-----
% http://www.fao.org/docrep/X5685E/x5685e03.htm
% 
% Froese, R. and C. Binohlan. 2000. Empirical relationships to estimate
%  asymptotic length, length at first maturity and length at maximum yield per
%  recruit in fishes, with a simple method to evaluate length frequency data. J.
%  Fish Biol. 56: 758-773.

% -----Notes:-----
% Starting estimates of LA, K, and T0 are needed as La_i, k_i, and t0_i in order
% to initialize the model before running the nonlinear fitting routine. You can
% provide some or all of these on the command line; setting one or more of them
% = 0 (or not specifying them at all) will result in default values being used:
% LA = maximum observed length; K = a random number ranging from 0.1-1.0; T0 =
% eps. If the maximum reported total length of the species in question is known
% (LM), it can be used to empirically estimate an initial value of LA, provided
% both LM and your observed lengths are total lengths in cm (Froese & Binohlan, 
% 2000); setting LM = 0 (the default) will skip this.

% -----Author:-----
% by David L. Jones, Jun-2006
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Nov-2008: replaced call to f_vonBertModelInit with a simplier way of
%           initializing the model; updated documentation  

% -----Check input & set defaults:-----
if (nargin < 3), Lm   = 0;    end; % no Lm by default
if (nargin < 4), La_i = 0;    end; % default to max length
if (nargin < 5), k_i  = 0;    end; % default to a random number
if (nargin < 6), t0_i = 0;    end; % default to eps
if (nargin < 7), iter = 100;  end; % 100 iterations by default
if (nargin < 8), tol  = 1e-8; end; % default tolerance

L = L(:); % force col vector
t = t(:); 
if (size(L,1) ~= size(t,1)), error('Size mismatch b/n L & t'); end;
% -------------------------------------

% -----Setup default initial parameter estimates:-----
if (La_i==0), La_i = max(L(:));                end;
if (k_i ==0), k_i  = f_randRange(1,100,1)/100; end;
if (t0_i==0), t0_i = eps;                      end;

% Estimate La empirically using Lm:
if (Lm>0), La_i = 10^(0.044 + 0.9841*log10(Lm)); end

init = [La_i k_i t0_i];
% -----------------------------------------------------

% Use nonlinear regression to estimate model parameters:
opt         = statset('MaxIter',iter,'Display','iter','TolFun',tol);
[out,resid] = nlinfit(t,L,'f_vonBert',init,opt); 

% Wrap results up into a structure:
result.La  = out(1);
result.k   = out(2);
result.t0  = out(3);
result.ta  = f_vonBertAge(result.La-0.01,result.La,result.k,result.t0,0);
result.Lo  = result.La*(1 - exp(result.k*result.t0)); % length at birth
result.obs = L;                                       % observed
result.fit = -1*(resid - L);                          % fitted 
result.res = resid;                                   % residuals
result.R2  = f_corr(result.fit,result.obs)^2;         % R-squared


% -----Make plots:-----
figure;
set(gcf,'color','w'); % set background color to white
yFit = f_vonBert([result.La result.k result.t0],result.t0:result.ta/100:result.ta);
plot(t,result.obs,'bo'); hold on; plot(result.t0:result.ta/100:result.ta,yFit,'r-');
xlabel('\bfAge');
ylabel('\bfLength');
title('\bfvon Bertalanffy Growth Model');
grid on;
