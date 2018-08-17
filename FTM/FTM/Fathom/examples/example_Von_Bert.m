% Example of using the von Bertalanffy Growth Model
% 
% by David L. Jones, 2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Note: data and procedures used for development and testing of tthe F_VONBERT
% functions were from:
%
% http://www.fao.org/docrep/X5685E/x5685e03.htm
%
% FAO31: worksheet 3.1
% herring: section 3.3 exercise 2

% -----Tutorial:-----
% load data set:
load fishGrowth.mat

% Create growth model for herring data:
result = f_vonBertModel(herring(:,1),herring(:,2));
%  
%                                      Norm of         Norm of
%    Iteration             SSE        Gradient           Step 
%   -----------------------------------------------------------
%            0         147.225
%            1          5.3479         570.981         1.63908
%            2         1.23425         103.226        0.869825
%            3         1.08069         10.0662        0.263687
%            4          1.0795        0.145617       0.0307697
%            5          1.0795      0.00121037      0.00242069
%            6          1.0795     1.17146e-05     0.000216098
%            7          1.0795     2.38207e-07     1.98455e-05
% Iterations terminated: relative change in SSE less than OPTIONS.TolFun


% Create some data for 10 fish whose lengths are known, but not age:
fishLen = [15 32 24 29 23 19 26 15 24 32]';

% Estimate fish age using previously derived vonBert model:
age = f_vonBertAge(fishLen,result.La,result.k,result.t0,1);

% Examine data:
[fishLen age]
