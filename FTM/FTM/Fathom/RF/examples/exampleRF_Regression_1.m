% Example Random Forest Regression 1
% 
% -----Author:-----
% by David L. Jones, Oct-2011
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% File: '.../examples/diabetes.mat'
% File: '.../examples/biscayne.mat'

%%%%%%%%%%%%%%%%%%%%%%%
%      DIABETES:      %
%%%%%%%%%%%%%%%%%%%%%%%
% 
% Load diabetes data:
load diabetes.mat X Y;

% Fit regression:
rf_dia = f_RFreg(X,Y,[],[],1,0,'stnd',1);

% Create diagnostic plots:
f_RFregPlot(rf_dia,1,1,1);
% -> variables 3 & 9 seem to be the most important

% Save PDF's:
f_pdf('RF_regression_diabetes_1');
f_pdf('RF_regression_diabetes_2');

%%%%%%%%%%%%%%%%%%%%%%%
%      Gray Snapper:  %
%%%%%%%%%%%%%%%%%%%%%%% 
% 
% File 'gray_snapper.mat' contains data on gray snapper collected from visual
% surveys of mangrove prop root habitats in Biscayne Bay, FL.

% Load data:
load gray_snapper.mat Y X X_txt;

% Variables:
% Y = # gray snapper per 60 m^2
% X = matrix of predictors
% 
% X_txt  = cell array of predictor labels as follows:
% tra    = transect number
% str    = strata (1=Key, 2=Mainland)
% sea    = season (1=dry, 2=wet)
% yr     = YYYY
% tem    = water temperature (degrees C)
% do     = dissolved O2 (mg/L)
% sal    = salinity (PPT)
% dep    = average depth (cm)
% E      = UTM Easting
% N      = UTM Northing
% jdate  = julian date
% FW     = distance (km) from nearest FW canal


% Create a regression model using a Random Forest :
rf_gray = f_RFreg(X,Y,[],[],1,0,'stnd',1,X_txt);

rf_null = f_RFreg(ones(size(Y,1),1),Y,[],[],1,0,'stnd',1,{'null'});




% Create diagnostic plots to assess the relative importance of each variable
% to the regression:
f_RFregPlot(rf_gray,1,1,1);

% What are the most important variables influencing the abundance/distribution
% of gray snapper in Biscayne Bay?

X_txt([8 9 10 12])
% ans = 'dep'    'E'    'N'    'FW'

% Save PDF's
f_pdf('RF_regression_gray_1');
f_pdf('RF_regression_gray_2')
