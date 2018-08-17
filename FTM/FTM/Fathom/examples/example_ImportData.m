% Example of Importing and Parsing Data
% 
% by David L. Jones, Jun-2012
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% File: '.../examples/boston_housing.csv'

% Note: this is also a good way to import data from an XLS spreadsheet that
% has column labels you want to retain. Simply export it from Excel as a
% 'comma-separated-values' (*.csv) file, then use the f_importCSV function to
% load it into the Matlab workspace.
% 

% Import Boston Housing data:
raw = f_importCSV('boston_housing.csv',1);
% 
% Variables:
% medv    = median value of owner-occupied homes in $1000's
% crim    = per capita crime rate by town
% zn      = proportion of residential land zoned for lots over 25,000 sq.ft.
% indus   = proportion of non-retail business acres per town
% chas    = Charles River dummy variable (= 1 if tract bounds river; 0 otherwise)
% nox     = nitric oxides concentration (parts per 10 million)
% rm      = average number of rooms per dwelling
% age     = proportion of owner-occupied units built prior to 1940
% dis     = weighted distances to five Boston employment centres
% rad     = index of accessibility to radial highways
% tax     = full-value property-tax rate per $10,000
% ptratio = pupil-teacher ratio by town
% b       = 1000(Bk - 0.63)^2 where Bk is the proportion of blacks by town
% lstat   = % lower status of the population
% lon     = longitude
% lat     = latitude
% 
% These are Boston house-price data from: Harrison, D. and Rubinfeld, D. L.
% 1978. Hedonic prices and the demand for clean air. J. Environ. Economics
% & Management 5: 81-102.

% Optionally parse the data as separate fields in a structure:
boston = f_parse(raw.dat,raw.txt);

 % Convert lat/lon to UTM coordinates
[boston.E,boston.N] = f_deg2utm(boston.lat,boston.lon);
