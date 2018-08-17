% Examples of creating Moran's Eigenvector Maps (PCNM's and MEM's)
% 
% by David L. Jones, Apr-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/Griffith_Peres-Neto_2006.mat'
% File: '.../examples/Dray_etal_2006.mat'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Example from Figure 1 of Griffith & Peres-Neto, 2006:            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Load the data file:
load Griffith_Peres-Neto_2006.mat

% Create a Minimum Spanning Tree:
mst = f_mst(f_euclid(xy'),xy);

% Show the corresponding binary connectivity matrix ('B'):
f_table(mst.B,'%d','s');
%  
% 0 0 0 0 1 0 0
% 0 0 1 0 1 1 0
% 0 1 0 0 0 0 0
% 0 0 0 0 1 0 0
% 1 1 0 1 0 0 0
% 0 1 0 0 0 0 1
% 0 0 0 0 0 1 0

% Create a weighting matrix ('A') based on the truncated distance matrix:
% (note: Dray et al. 2006 use full vs. truncated distances at this step)
A = f_dis2sim(mst.tDis,2,2);

% Create Eigenvector Maps, use 1000 iterations, don't keep neg eigenvalues:
MEM = f_eigenMaps(mst,A,1000,0);

% Show the (truncated) spatial weighting matrix ('W') used to create MEM's:
f_table(MEM.W,'%0.2f','s');
%  
% 0.00 0.00 0.00 0.00 0.96 0.00 0.00
% 0.00 0.00 0.98 0.00 0.94 0.94 0.00
% 0.00 0.98 0.00 0.00 0.00 0.00 0.00
% 0.00 0.00 0.00 0.00 0.99 0.00 0.00
% 0.96 0.94 0.00 0.99 0.00 0.00 0.00
% 0.00 0.94 0.00 0.00 0.00 0.00 0.98
% 0.00 0.00 0.00 0.00 0.00 0.98 0.00

% Show positive eigenvectors:
f_table(MEM.evects,'%+0.2f','s')
%  
% -0.37 -0.19  0.72
%  0.16  0.49 -0.01
%  0.15  0.63  0.01
% -0.38 -0.19 -0.70
% -0.48 -0.02 -0.01
%  0.50 -0.21 -0.01
%  0.43 -0.50 -0.00

% Show positive eigenvalues:
f_table(MEM.evals,'%+0.2f','s')
%  
%  1.18
%  0.61
%  0.00
% 
% -> Note: the eigenvalues provided in Fig. 1 of Griffith & Peres-Neto (2006)
% are apparently wrong, as the values above are also returned by Peres-Neto's
% EigenvectorMaps.m Matlab function.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Example from Figure 1b of Dray et al. (2006):                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Load the data file:
load Dray_etal_2006.mat

% Create PCNM's:
pcnm = f_pcnm(dat);

% Show complete distances:
f_table(pcnm.dis,'%+0.2f','s')
% 
%  0.00  2.00  5.00  3.16  3.61
%  2.00  0.00  3.00  1.41  2.24
%  5.00  3.00  0.00  2.24  2.83
%  3.16  1.41  2.24  0.00  3.00
%  3.61  2.24  2.83  3.00  0.00

% Show truncated distances:
f_table(pcnm.tDis,'%+0.2f','s')
% 
%  0.00  2.00  8.94  8.94  8.94
%  2.00  0.00  8.94  1.41  2.24
%  8.94  8.94  0.00  2.24  8.94
%  8.94  1.41  2.24  0.00  8.94
%  8.94  2.24  8.94  8.94  0.00

% Calculate complete/truncated similarities:
pcnm.sim  = f_dis2sim(pcnm.dis,2,2);
pcnm.tSim = f_dis2sim(pcnm.tDis,2,2);

% Show complete similarities:
f_table(pcnm.sim,'%+0.2f','s')
%  
%  1.00  0.84  0.00  0.60  0.48
%  0.84  1.00  0.64  0.92  0.80
%  0.00  0.64  1.00  0.80  0.68
%  0.60  0.92  0.80  1.00  0.64
%  0.48  0.80  0.68  0.64  1.00

% Show truncated similarities:
f_table(pcnm.tSim,'%+0.2f','s')
% 
%  1.00  0.95  0.00  0.00  0.00
%  0.95  1.00  0.00  0.97  0.94
%  0.00  0.00  1.00  0.94  0.00
%  0.00  0.97  0.94  1.00  0.00
%  0.00  0.94  0.00  0.00  1.00

% Plot Minimum Spanning Tree:
mst = f_mst(f_euclid(dat'),dat);
f_plotNeigh(mst,2,'k');
axis([-1 6 -1 5]);
grid on;
title('MST of truncated distances');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Example from Figure 2 of Dray et al. (2006):                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Generate 100 irregularly spaced points along a transect:
A.dat = [sort(rand(100,1)*100) zeros(100,1)];
B.dat = [sort(rand(100,1)*100) zeros(100,1)];

% Calculate PCNM's:
A.pcnm = f_pcnm(A.dat);
B.pcnm = f_pcnm(B.dat);

% Correlation between PCNM's:
pcnmR_1   = f_corr(A.pcnm.evects(:,1), B.pcnm.evects(:,1));
pcnmR_15  = f_corr(A.pcnm.evects(:,15), B.pcnm.evects(:,15));
% 
fprintf(['\n\n -> r for PCNM  1: A vs. B = ' num2str(pcnmR_1) '\n'])
fprintf(['\n\n -> r for PCNM 15: A vs. B = ' num2str(pcnmR_15) '\n'])
% -> r for PCNM  1: A vs. B = -0.15314
% -> r for PCNM 15: A vs. B = -0.15182
% 
% NOTE: A & B are randomly generated, so your values of 'r' will vary with
%       each run!

% Create a MST to use as the binary connectivity matrix ('B') : 
A.mst = f_mst(f_dis(A.dat,'euc'),A.dat);
B.mst = f_mst(f_dis(B.dat,'euc'),B.dat);

% Create a weighting matrix ('A'):
A.A = f_dis2sim(A.mst.dis,1);
B.A = f_dis2sim(B.mst.dis,1);

% Create Eigenvector Maps, use 1000 iterations, keep neg eigenvalues:
A.MEM = f_eigenMaps(A.mst,A.A,1000,1);
B.MEM = f_eigenMaps(B.mst,B.A,1000,1);
 
% Correlation between MEM's:
memR_1   = f_corr(A.MEM.evects(:,1), B.MEM.evects(:,1));
memR_15  = f_corr(A.MEM.evects(:,15), B.MEM.evects(:,15));
%  
fprintf(['\n\n -> r for MEM  1: A vs. B = ' num2str(memR_1) '\n'])
fprintf(['\n\n -> r for MEM 15: A vs. B = ' num2str(memR_15) '\n'])
% 
% -> r for MEM  1: A vs. B = 0.89395
% -> r for MEM 15: A vs. B = 0.99279

%-----Make some plots:-----
figure; set(gcf,'color','w');

% Irregularly spaced data:
subplot(5,2,1)
plot(A.dat(:,1),A.dat(:,2),'k.');
axisVar = axis;
axis([min(A.dat(:,1))-1 max(A.dat(:,1))+1 axisVar(3:4)]);
title(['Sample A']);

subplot(5,2,2)
plot(B.dat(:,1),B.dat(:,2),'k.');
axisVar = axis;
axis([min(B.dat(:,1))-1 max(B.dat(:,1))+1 axisVar(3:4)]);
title(['Sample B']);

% PCNM's:
subplot(5,2,3)
plot(A.pcnm.evects(:,1),'k-')
title(['A: PCNM 1'])

subplot(5,2,4)
plot(B.pcnm.evects(:,1),'k-')
title(['B: PCNM 1'])

subplot(5,2,5)
plot(A.pcnm.evects(:,15),'k-')
title(['A: PCNM 15'])

subplot(5,2,6)
plot(B.pcnm.evects(:,15),'k-')
title(['B: PCNM 15'])

% MEM's:
subplot(5,2,7)
plot(A.MEM.evects(:,1),'k-')
title(['A: MEM 1'])

subplot(5,2,8)
plot(B.MEM.evects(:,1),'k-')
title(['B: MEM 1'])

subplot(5,2,9)
plot(A.MEM.evects(:,15),'k-')
title(['A: MEM 15'])

subplot(5,2,10)
plot(B.MEM.evects(:,15),'k-')
title(['B: MEM 15'])

% Clean up:
clear axisVar;
