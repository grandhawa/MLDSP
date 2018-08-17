function [p,t,net] = f_nnMLP(X,grp,hid)
% - feedforward multilayer perceptron neural network for classification
%
%
% USAGE: [p,t,net] = f_nnMLP(X,grp,hid);
%
% X   = matrix of training data   (rows = obs, cols = variables) 
% grp = vector of integers specifying group membership of X
% hid = number of hidden neurons  (default = 5)
%
% p   = network input
% t   = network target
% net = neural network
%
% SEE ALSO: f_nnMLP632, f_pnn

% -----Notes:-----
% This function is used to create a two-layer feedforward, backprogagating
% neural network in order to solve classification problems.
%
% LOGSIG must be the transfer function for the output as the target is uses
% dummy variables (0/1) to code for groups membership. The hidden layer
% could use LOGSIG as well, but initial tests indicate TANSIG works better.
%
% This net consists of 2 layers using the DOTPROD weight function, NETSUM
% net input function, and the specified transfer functions. The first layer
% has weights coming from the input.  Each subsequent layer has a weight
% coming from the previous layer.  All layers have biases.  The last layer
% is the network output. Each layer's weights and biases are initialized
% with INITNW. Training is done with the specified training function.
% MSE is the performance measure.

% -----Author:-----
% by Dave Jones, Jan-2004
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), hid  = 5; end; % default 5 hidden neurons

grp = grp(:); % force column vector

if size(X,1) ~= size(grp,1)
   error('X & GRP must have same # of rows!');
end

if (f_isScalar(hid)==0)
   error('HID must be a scalar!')
end
% -------------------------------------

uGrp  = unique(grp);  % unique groups
noGrp = size(uGrp,1); % # of groups

p   = [f_stnd(X)]';             % input vectors
t   = [f_designMatrix(grp,0)]'; % target vectors

% -----Create Neural Network:-----
%
% net = newff(inputrange,[hid,O],{transFnc1,transFnc2},trainFnc,learnFnc,performFnc)
%
% inputrange = a matrix of ranges the inputs take
% hid        = # neurons in hidden layer
% O          = # neurons in output layer (depends on your target)
% transFnc1  = transfer function for hidden layer
% transFnc2  = transfer functio fof output layer
% trainFnc   = training function
% performFnc = performance function

%net = newff(minmax(p),[hid,noGrp],{'tansig','logsig'},'traincgb','learngdm','mse');

net = newff(minmax(p),[hid,noGrp],{'tansig','logsig'},'trainlm','learngdm','mse');

% Train parameters:
net.trainParam.epochs   = 1000;  % max. times complete data used in training
net.trainParam.show     = NaN;   % time between status reports
net.trainParam.lr       = 0.01;  % learning rate (was 0.01)
net.trainParam.goal     = 0.001; % performance goal (was 0.001)
