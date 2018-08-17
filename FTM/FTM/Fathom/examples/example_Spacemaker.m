% Example of using Moran's eigenvector maps (MEM's) in spatial ecology
% 
% by David L. Jones, May-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/spacemaker.mat'
% File: '.../examples/oribatid_mites.mat''

% -----Notes:-----
% These examples follows that presented in 'Dray, S. 2006. Moran's eigenvectors
% of spatial weighting matrices in R', which is included in the 'spacemakeR for
% R' documentation. Available from:
% http://biomserv.univ-lyon1.fr/~dray/software.php 

% Load the data file:
load spacemaker.mat

% Create a Euclidean distance matrix:
xyir.dis = f_dis(xyir.dat,'euc');

% Show the distance matrix:
f_table(xyir.dis,'%0.3f','s')
% 
% 0.000 0.640 0.218 0.166 0.561 0.540 0.403 0.230 0.562 0.518
% 0.640 0.000 0.424 0.483 0.417 0.383 0.788 0.550 0.455 0.286
% 0.218 0.424 0.000 0.062 0.398 0.368 0.496 0.192 0.411 0.354
% 0.166 0.483 0.062 0.000 0.414 0.388 0.490 0.150 0.422 0.411
% 0.561 0.417 0.398 0.414 0.000 0.038 0.895 0.349 0.038 0.589
% 0.540 0.383 0.368 0.388 0.038 0.000 0.864 0.334 0.073 0.551
% 0.403 0.788 0.496 0.490 0.895 0.864 0.000 0.616 0.907 0.534
% 0.230 0.550 0.192 0.150 0.349 0.334 0.616 0.000 0.343 0.541
% 0.562 0.455 0.411 0.422 0.038 0.073 0.907 0.343 0.000 0.620
% 0.518 0.286 0.354 0.411 0.589 0.551 0.534 0.541 0.620 0.000
 
% Create a Gabriel Graph to use as the binary connectivity matrix ('B'):
xyir.gab = f_gabriel(xyir.dat);

% Create a weighting matrix ('A'):
xyir.A = f_dis2sim(xyir.dis,1);

% Create Eigenvector Maps, use 1000 iterations, keep neg eigenvalues:
xyir.MEM = f_eigenMaps(xyir.gab,xyir.A,1000,1);

% Show the spatial weighting matrix used to create MEM's:
f_table( xyir.MEM.W,'%0.3f','s')
% 
% 0.000 0.000 0.000 0.817 0.000 0.000 0.556 0.000 0.000 0.000
% 0.000 0.000 0.533 0.000 0.000 0.578 0.000 0.000 0.000 0.685
% 0.000 0.533 0.000 0.932 0.000 0.594 0.000 0.000 0.000 0.610
% 0.817 0.000 0.932 0.000 0.000 0.000 0.000 0.835 0.000 0.000
% 0.000 0.000 0.000 0.000 0.000 0.958 0.000 0.000 0.958 0.000
% 0.000 0.578 0.594 0.000 0.958 0.000 0.000 0.631 0.000 0.000
% 0.556 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.412
% 0.000 0.000 0.000 0.835 0.000 0.631 0.000 0.000 0.000 0.000
% 0.000 0.000 0.000 0.000 0.958 0.000 0.000 0.000 0.000 0.000
% 0.000 0.685 0.610 0.000 0.000 0.000 0.412 0.000 0.000 0.000

% Plot magnitude of eigenvalues:
figure;
set(gcf,'color','w'); % set bg color to white
bar(xyir.MEM.evals);
title('\bfEigenvalues of spatial weighting matrix (W)')

% Plot the eigenvectors in geographic space:
figure;
n = size(xyir.MEM.evects,2);
for i = 1:n
   subplot(3,3,i)
   f_plotNeigh(xyir.gab,0,'k')
   f_bubble(xyir.dat(:,1),xyir.dat(:,2),xyir.MEM.evects(:,i),10,1);
   title(['Eigenvector ' num2str(i) ' (' num2str(f_round(xyir.MEM.evals(i),3)) ')'])
end
clear n i;

% Show Moran's I and corresponding p-values:
f_table([xyir.MEM.MC xyir.MEM.p],'%+0.5f','s')
% 
%  MC:      p:
%  0.70159  0.00100
%  0.52255  0.00600
%  0.41162  0.01300
%  0.02725  0.27900
% -0.13668  0.45700
% -0.28734  0.25300
% -0.49867  0.05200
% -0.72968  0.00500
% -1.01065  0.00100

% Show the congruence among the Eigenvalues and Moran's I:
figure;
set(gcf,'color','w'); % set bg color to white
hold on; box on;
axis([-2 1.5 -1.25 1]);
axisVar = axis;
slope   = size(xyir.dat,1)/sum(xyir.MEM.W(:));
plot([axisVar(1);axisVar(2)],[axisVar(1);axisVar(2)]*slope,'-','LineWidth',1.5,...
   'Color', [1 1 1]*0.85);
plot(xyir.MEM.evals,xyir.MEM.MC,'o','MarkerEdgeColor','none',...
   'MarkerFaceColor','k','MarkerSize',13);
txt = cellstr(num2str([1:size(xyir.MEM.evals,1)]')); % create some labels:
text(xyir.MEM.evals,xyir.MEM.MC,txt,'HorizontalAlignment','center',...
   'Color','w','FontSize',10,'FontWeight','bold')
title([ 'Correlation (R^2) = ' num2str((f_corr(xyir.MEM.evals,xyir.MEM.MC))^2) ])
xlabel('Eigenvalues')
ylabel('Moran''s I')
axis equal;
f_origin('v',':'); % mark division b/n +/- eigenvalues
clear slope axisVar txt;

% Plot significant eigenvectors in geographic space:
idx = find(xyir.MEM.p<=0.05);
n = size(idx,1);
figure;
for i = 1:n
   subplot(3,2,i)
   f_plotNeigh(xyir.gab,0,'k')
   f_bubble(xyir.dat(:,1),xyir.dat(:,2),xyir.MEM.evects(:,idx(i)),10,1);
   title(['ev ' num2str(idx(i)) ' (' num2str(f_round(xyir.MEM.evals(idx(i)),3)) ')'])
end
clear idx n i;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Oribatid Mite example from Dray's 'spacemakeR for R' manual:         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Load data:
load oribatid_mites.mat

% Hellinger transform the species data:
bio.H = f_hellinger(bio.num);

% Perform preliminary RDA on X,Y coordinates to identify linear trends:
result_XY = f_rda(bio.H,[env.x env.y],[0],1000,1);

% Residuals serve as the de-trended data:
bio.dt = result_XY.res;

% Create a Delaunay Triangulation based neighbor graph:
tri = f_delaunay([env.x env.y]);

% Show the binary connectivity matrix ('B'):
f_table(tri.B,'%d','s')

% Plot connectivity graph:
f_plotNeigh(tri,1,'r');
axis([-0.25 2.75 0 10])
title('\bfDelaunay Triangulation');

% Create a weighting matrix ('A') that doesn't change 'B':
A = ones(size(tri.B));

% Create Eigenvector Maps, use 1000 iterations, KEEP neg eigenvalues:
MEM = f_eigenMaps(tri,A,1000,1);

% AIC-based stepwise selection of MEM's:
model_b = f_rdaAIC(bio.dt,MEM.evects,0,0);

% Show the results (i.e., the optimal model):
% model_b = 
%        RSS: [69x1 double]
%         R2: [69x1 double]
%      R2adj: [69x1 double]
%        AIC: [69x1 double]
%        var: {1x69 cell}
%        idx: [3 2 4 7 1 5 16 53]
%     minAIC: -95.365
%       null: -87.471
%          X: [70x8 double]
% 
% -> Note: the example in the 'spacemakeR' manual shows slightly different
% eigenvectors were selected (i.e., 3 2 7 1 16). This is due to 1) slight
% differences in the values of the eigenvectors returned by 'R' vs. Matlab,
% resulting in differences in AIC, and to 2) differences in precision and
% rounding between the two platforms. The tutorial in 'exampleAIC.m'
% demonstrates that 'f_rdaAIC' and 'ortho.AIC' return the same results when
% provided identical input. The examples in 'exampleMatlab_v_R.m' demonstrates
% there are only VERY slight differences in the ouptput generated by the two
% platforms.

% -----Test a range of concave down weighting functions:-----
% 
% Initialize:
clear A;
e = [2:10]'; % range of exponents

% Create neighbor graphs and corresponding weighting matrices:
for i = 1:size(e,1)
   N{i} = tri;                       % neighbor graph
   A{i} = f_dis2sim(tri.dis,2,e(i)); % weighting matrix ('A')
end
clear i e;

% Select optimal spatial weighting matrix, keep negative eigenvalues:
model_f2 = f_selectW(bio.dt,N,A,1);
% model_f2 = 
%        RSS: [69x1 double]
%         R2: [69x1 double]
%      R2adj: [69x1 double]
%        AIC: [69x1 double]
%        var: {1x69 cell}
%        idx: [3 2 7 1 16 4 5 53]
%     minAIC: -95.653
%       null: -87.471
%          X: [70x8 double]
% 
% -> Note: the optimal model, among the 9 compared, was the first one considered
% with e = 2 and 8 MEM's.



% -----Evaluate a range of distance-based neighbor matrices-----:
% 
% Get minimun distance criterion based on a Minimun Spanning Tree:
dis  = f_euclid([env.x env.y]');
mst  = f_mst(dis,[env.x env.y]);
dMin = max(mst.len);

% Get maximum distance criterion based on a Variogram:
variogram = f_variogram(bio.dt,[env.x env.y],20,1000,1);
% Examine variogram plot, choose highest distance class that is still
% significant (but not more than half the maximum distance)
dMax = variogram.d(9);

% Construct 10 distance-based Nearest Neighbor graphs:
dRange = linspace(dMin,dMax,10)';
% 
clear A;
for i=1:10
   dnn{i} = f_dnn(dis,[0 dRange(i)],[env.x env.y]);
   A{i}   = ones(size(dnn{i}.B)); % weighting matrix ('A') that doesn't change 'B'
end
clear i;

% Select optimal spatial weighting matrix, keep negative eigenvalues:
[model_dnn,idx] = f_selectW(bio.dt,dnn,A,1);
% model_dnn = 
%        RSS: [69x1 double]
%         R2: [69x1 double]
%      R2adj: [69x1 double]
%        AIC: [69x1 double]
%        var: {1x69 cell}
%        idx: [3 2 1 6]
%     minAIC: -100.02
%       null: -87.471
%          X: [70x4 double]

% Show distance class associated with optimal model:
dRange(idx)
% ans = 2.0368




% ==============================
% The following is NOT working yet, Dray is not comparing 90 W's:
% ==============================


% -----Similar approach using a spatial weighting function:-----
% Construct Nearest Neighbor graphs based on 10 distance classes, for each
% distance class evaluate each of 9 exponents defining a weighting function:
clear A dnn;
e    = [2:10]'; % range of exponents
k    = 0;
list = [];
for i=1:10
   for j=1:size(e,1)
      k      = k+1;
      dnn{k} = f_dnn(dis,[0 dRange(i)],[env.x env.y]);
      A{k}   = f_dis2sim(dnn{k}.dis,2,e(j)); % weighting matrix ('A')
      list   = [ list; [i j] ];              % keep track of parameters
   end
end

% Select optimal spatial weighting matrix, keep negative eigenvalues:
[model_dnn_f2,idx,dave] = f_selectW(bio.dt,dnn,A,1);
