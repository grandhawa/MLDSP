% Examples of spike removal from time series data
% 
% by David L. Jones, Apr-2011
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/Olsson_Newell.mat'

% The file 'Olsson_Newell.mat' contains two sets of processed plant time
% series data (Y and W) used in example 6.8 from Olsson, G., and B. Newell.
% 1999. Wastewater treatment systems: modelling, diagnosis and control. IWA
% Publishing, London. pp.140-141.

% Load data file:
load Olsson_Newell.mat;

% Add X variable:
X = [1:size(Y,1)]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           ROSNER TEST:          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----Y:-----
[Ys,S,res] = f_rosner(Y,1,0.5,0.1,2); % spike removal:
idx        = find(S==1);              % get index to spikes

figure;
subplot(2,1,1); 
box on;
hold on;
h(1) = plot(X,Y,'b-');                % original data
h(2) = plot(X,Ys,'r-');               % spikes eliminated
legend(h,'Original','Filtered');
title('ROSNER TEST');

subplot(2,1,2); 
box on;
hold on;
h(1) = plot(X,res,'b-');              % residuals
h(2) = plot(X(idx),res(idx),'ro');    % identify spikes
legend(h,'Residuals','Spikes');

% -----W:-----
[Ws,S,res] = f_rosner(W,1,0.5,0.1,2); % spike removal
idx        = find(S==1);              % get index to spikes

figure;
subplot(2,1,1); 
box on;
hold on;
h(1) = plot(X,W,'b-');                % original data
h(2) = plot(X,Ws,'r-');               % spikes eliminated
legend(h,'Original','Filtered');
title('ROSNER TEST');

subplot(2,1,2); 
box on;
hold on;
h(1) = plot(X,res,'b-');              % residuals
h(2) = plot(X(idx),res(idx),'ro');    % identify spikes
legend(h,'Residuals','Spikes');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           GRUBBS TEST:          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----Y:-----
[Ys,S] = f_grubbs(Y,1,0.05,5); % spike removal
idx    = find(S==1);           % get index to spikes

figure;
box on;
hold on;
h(1) = plot(X,Y,'b-');         % original data
h(2) = plot(X,Ys,'r-');        % spikes eliminated
legend(h,'Original','Filtered');
title('GRUBBS TEST');


% -----W:-----
[Ws,S] = f_grubbs(W,1,0.05,5);  % spike removal:
idx    = find(S==1);            % get index to spikes

figure;
box on;
hold on;
h(1) = plot(X,W,'b-');           % original data
h(2) = plot(X,Ws,'r-');          % spikes eliminated
legend(h,'Original','Filtered');
title('GRUBBS TEST');

% Note: the Grubbs Test does not appear to be very effective in detecting
% local spikes, as illustrated in the above example involving 'W'

