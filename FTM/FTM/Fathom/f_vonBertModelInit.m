function result = f_vonBertModelInit(t,L)
% - initial parameter estimate of a von Bertalanffy growth model
%
% USAGE: result = f_vonBertModelInit(t,L)
%
% t = observed age at length L
% L = observed length at age t
%
% La = theoretical asymptotic length
% k  = growth coefficient
%
% SEE ALSO: f_vonBert, f_vonBertModel

% -----References:-----
% http://www.fao.org/docrep/X5685E/x5685e03.htm

% -----Author:-----
% by David L. Jones, Jun-2006
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
t = t(:); % force col vector
L = L(:); % force col vector

if (size(L,1) ~= size(t,1)), error('Size mismatch b/n L & t'); end;
% -------------------------------------

% Sort the data:
sortVar = sortrows([L t]);
L       = sortVar(:,1);
t       = sortVar(:,2);

% Instantaneous growth rate:
y = [L(2:end) - L(1:(end-1))]./[t(2:end) - t(1:(end-1))];

% Mean length:
x = mean([L(1:(end-1)) L(2:end)],2);

figure
plot(x,y,'k.');

% Linear regression:
model = f_mregress(x(:),y(:),0,0,0);
k     = -1*model.b(2);              % slope of regression
La    = -1*[model.b(1)/model.b(2)]; % x-intercept

% Wrap results up into a structure:
result.model = model;
result.k     = k;
result.La    = La;
result.x     = x;
result.y     = y;

% Can't find a good estimate of t0, so will use t0=eps for now.
%
% estimate t0 for each age, keep the mean; there are other methods in the
% paper referenced for estimating t0, this is probably not a very accurate
% one.
%t0 = mean(t + [(1/k) * real(log([(La-L)/La]))]);
%
% z     = real(log((La - L)/La));      % eq. 3.5
% model = f_mregress(t(:),z(:),0,0,0);
% t0    = -1 * [model.b(1)/model.b(2)];       % x-intercept
% 
% plot(t,z,'ro');
% 
% t0
% model.b
