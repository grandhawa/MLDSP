function model = f_RFreg(X,Y,nTree,mTry,imp,sim,stnd,verb,X_txt,opt)
% - fit a regression via Random Forest
% 
% USAGE: model = f_RFreg(X,Y,nTree,mTry,imp,sim,stnd,verb,X_txt,opt);
% 
% X     = predictors                                (row = obs, col = variables)
% Y     = column vector specifying response variable
% nTree = # random trees to grow             (default = 1000 if set to empty [])
% mTry  = # random variable subsets            (default used if set to empty [])
%                           (default = max(floor(D/3),1), where D = # predictor)
% imp   = compute variable importance measures                     (default = 0)
% sim   = return proximity matrix                                  (default = 0)
% stnd  = use 'raw', standardized ('stnd'), or centered('ctr') variables for X
%                                                              (default = 'raw')
% verb  = verbose listing of input specifications                  (default = 0)
% X_txt = cell array of variable labels                   (if empty, autocreate)
%         e.g., X_txt = {'sal' 'tmp' 'elev'};
% opt   = structure of extra options                        (see details inside)
% 
% model = structure of Random Forest with the following fields:
% .nTree      = # random trees grown
% .mTry       = size of random variable subsets
% .X          = transformed predictors
% .X_txt      = cell array of predictor labels
% .stnd       = transformation applied to predictors
% .Y          = observed response variable
% .Y_hat      = fitted values of response variable
% .prox       = symmetric proximity (similarity) matrix based on frequency
%               obs occur in the same node
% .localImp   = importance measures for each observation
% .meanAcc    = mean decrease in accuracy (col: variable)
% .meanAccSD  = standard error of meanAcc
% .zScore     = z-score of meanAcc
% .zP         = parametric significance of z-score
% .meanMSE    = mean decrease in node impurity (measured by MSE)
% .mse        = mean squared error with each iteration (= RSS/n)
% .rsq        = pseudo R-squared with each iteration
% .inbag      = in bag
% .oob_times  = # times obs are out-of-bag (& used for OOB error estimates)
% .nPerm      = # times OOB data were permuted per tree for importance measures
% .biasCorr   = bias correction
% .coef       = cooefficients for bias correction
% .lDau       = left daughter
% .rDau       = right daughter
% .nodestatus = node status
% .nrnodes    = 
% .upper      = 
% .avnode     = 
% .mbest      = 
% .treeSize   = tree size
% .type       = type of Random Forest (= 'regression')
% 
% SEE ALSO: f_RFregPredict, f_RFregPlot, f_RFclass


% ----- Author(s): -----
%**************************************************************
%* mex interface to Andy Liaw et al.'s C code (used in R package randomForest)
%* Added by Abhishek Jaiantilal ( abhishek.jaiantilal@colorado.edu )
%* License: GPLv2
%* Version: 0.02
%
% Calls Classification Random Forest
% A wrapper matlab file that calls the mex file
% This does training given the data and labels
% Documentation copied from R-packages pdf
% http://cran.r-project.org/web/packages/randomForest/randomForest.pdf
% Tutorial on getting this working in tutorial_ClassRF.m
%**************************************************************
% 
% modified by David L. Jones, Feb-2010
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB'.

% edits by David L. Jones:
% - edited documentation so easier to read
% renamed 'extra_otions' to 'opt', replaced '|' with '||', replaced 'DEBUG_ON'
% with 'verb', lots of re-organization of the code, TRUE/FALSE replaced with
% 1/0, renamed 'oobprox' to 'oob_prox', 'doProx' to 'sim',

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                         OPTIONS                                     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opt = strucuture of extra options with the following fields (TRUE=1, FALSE=0):
%
% .replace = bootstrap sampling uses replacement                   (default = 1)
%
% .sampsize = Size(s) of sample to draw; (default = N if bootstrap uses
%             replacement, else default = 0.632*N)
%
% .nodesize = Minimum size of terminal nodes. Setting this number larger causes
%             smaller trees to be grown (and thus take less time). Default = 5
%             for regression (5)
%
% .localImp = compute casewise importance measures                 (default = 0)
%             (if = 1, then over-rides 'imp')
%
% .oob_prox = calculate proximities only on 'out-of-bag' data      (default = 1)
%
% .do_trace   = If set to 1, gives verbose output as Random Forest is run. If set
%             to some integer, then running output is printed for every do_trace 
%             trees.                                               (default = 0)
%
% .keep_inbag = should an n by nTree matrix be returned that keeps track of
%               which samples are 'in-bag' in which trees (but not how many
%               times, if sampling with replacement)               (default = 0)
% 
% .biasCor  = performs bias correction for regression if = 1; experimental so
%             use at your own risk                                 (default = 0)
% 
% .nPerm    = number of times the 'out-of-bag' data are permuted per tree for
%             assessing variable importance. Number larger than 1 gives slightly
%             more stable estimate, but not very effective (default = 1)
% ------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                 DOCUMENTATION by David L. Jones                     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nTree: at least 1000 trees is recommended if requesting variable importance
%   measures and/or proximities. Breiman recommends nTree=5000 when there are
%   many variables and you want stable importance measures.
% 
% For the importance measures, higher values equate to greate importance. If
% imp=0, column 1 of impout = a vector of NaN's
% ------------------------------------------------------------------------------

% Oct-2011: localImp set to imp by default; added error message when IMP=1
%           & SIM=1; renamed 'jprint' to 'do_trace'; added '.type' field in
%           output

% -----Check input & set defaults:-----
if (nargin < 3), nTree = [];    end % use default # trees
if (nargin < 4), mTry  = [];    end % use default size of variable subsets
if (nargin < 5), imp   = 0;     end % default don't compute importance measures
if (nargin < 6), sim   = 0;     end % default don't return proximity matrix
if (nargin < 7), stnd  = 'raw'; end % default use raw predictors in X
if (nargin < 8), verb  = 0;     end % default don't create plots/listings
if (nargin < 9), X_txt = [];    end % use default variable labels

% Check size of X:
N = size(X,1); % # observations
D = size(X,2); % # variables
if (N==0), error('X has 0 rows!'); end

% Replace empty inputs with defaults:
if isempty(nTree), nTree = 1000;              end
if isempty(mTry),  mTry  = max(floor(D/3),1); end

% Check mTry:
if (mTry<1 || mTry>D), error('MTRY must range from 1 to D!'); end
mTry  = max(floor(D/3),1); 

% Check for missing data:
if any(isnan(X)); error('NaNs in X'); end
if any(isnan(Y)); error('NaNs in Y'); end

% Check for compatible sizes:
if (size(X,1) ~= size(Y,1)), error('Size mismatch b/n X & Y'); end
if size(Y,2)~=1, error('Y should be a column vector!'); end

% Check whether this is a classification (vs. regression) problem:
if length(unique(Y))<=5,
   warning('Do you want regression? there are just 5 or less unique response values!');
end

% Check imp + sim:
if (imp>0) && (sim>0)
   disp('Due to a bug in the original FORTRAN, proximities calculated when')
   disp('IMP=1 are incorrect:')
   error('Try again with IMP=0 and SIM=1!');
end
% -------------------------------------

% -----Check/generate variable labels:-----
if isempty(X_txt)
   X_txt = cellstr(num2str((1:D)'))'; % default X labels
end

% if labels are not cell arrays, try forcing them:
if iscell(X_txt)<1, X_txt = num2cell(X_txt); end;

X_txt = X_txt(:); % force as row vector
if (D == size(X_txt,1))<1
   error('Size mismatch b/n X & X_txt!');
end
% -------------------------------------------

% -----Standardize/Center data:-----
switch stnd
   case 'raw' % use raw data for X
      stndTxt = '(Using raw data for X) \n';
   case 'stnd' % standardize variables in X
      X = f_stnd(X);
      stndTxt = '(Using standardized data for X) \n';
   case 'ctr' % center variables in X
      X = f_center(X);
      stndTxt = '(Using centered data for X) \n';
   otherwise
      error('Invalid option for STND!');
end

if (verb>0)
   % List input parameters:
   fprintf('\nGenerating %d trees with %d random subsets of variables\n',nTree,mTry);
end

% -----Parse extra options:-----
if exist('opt','var')
   if isfield(opt,'replace');    replace     = opt.replace;    end
   if isfield(opt,'sampsize');   sampsize    = opt.sampsize;   end
   if isfield(opt,'nodesize');   nodesize    = opt.nodesize;   end
   if isfield(opt,'localImp');   localImp    = opt.localImp;   end
   if isfield(opt,'nPerm');      nPerm       = opt.nPerm;      end
   if isfield(opt,'oob_prox');   oob_prox    = opt.oob_prox;   end
   if isfield(opt,'do_trace');   do_trace    = opt.do_trace;   end
   if isfield(opt,'biasCorr');   biasCorr    = opt.biasCorr;   end
   if isfield(opt,'keep_inbag'); keep_inbag  = opt.keep_inbag; end
end


% -----Set defaults if not already set:-----
if ~exist('replace','var');  replace = 1; end
if ~exist('sampsize','var');
   if (replace)
      sampsize = N;             % bootstrap with replacement
   else
      sampsize = ceil(0.632*N); % resampling without replacement
   end
end
if ~exist('nodesize','var');   nodesize   = 5;   end % regression = 5
if ~exist('localImp','var');   localImp   = imp; end
if ~exist('oob_prox','var'),   oob_prox   = sim; end
if ~exist('nPerm','var');      nPerm      = 1;   end
if ~exist('do_trace','var');   do_trace   = 0;   end
if ~exist('biasCorr','var');   biasCorr   = 0;   end
if ~exist('keep_inbag','var'); keep_inbag = 0;   end

% -----Categories:----- [??? DLJ]
% now handle categories. Problem is that categories in R are more enhanced. In
% this i ask the user to specify the column/features to consider as categories,
% 1 if all the values are real values else specify the number of categories here
if exist ('opt','var') && isfield(opt,'categories')
   cat = opt.categories;
else
   cat = ones(1,D);
end

maxcat = max(cat);
if maxcat>32
   error('Can not handle categorical predictors with more than 32 categories');
end
% ----------------------

% When localImp = 1, imp must = 1:
if (localImp>0), imp = 1; end

if imp
   if (nPerm<1)
      nPerm = int32(1);
   else
      nPerm = int32(nPerm);
   end
end

%  Make sure nTree is sufficient:
if ( nTree<1000 && (imp==1 || sim==1) )
   error('Robust importances &/or proximities require at least 1000 trees!')
end

% -----Stratified data:-----
Stratify = (length(sampsize)>1);
if (~Stratify && sampsize>N)
   error('Sampsize too large')
end

if Stratify
   error('Sampsize should be of length one')
end

% -----Verbose reporting of inputs:-----
if (verb>0)
   % Print the parameters that i am sending in:
   fprintf('=======================================================\n');
   fprintf('=====     INPUT PARAMETERS FOR RANDOM FOREST:     =====\n');
   fprintf('=======================================================\n');
   fprintf('nTree          = %d\n',nTree);
   fprintf('mTry           = %d\n',mTry);
   fprintf('X: # obs = %d, # var''s = %d \n',N,D);
   fprintf('size(cat)     = %d\n',size(cat));
   fprintf('maxcat         = %d\n',maxcat);
   fprintf('size(sampsize) = %d\n',size(sampsize));
   fprintf('sampsize[0]    = %d\n',sampsize(1));
   fprintf('Stratify       = %d\n',Stratify);
   fprintf('imp            = %d\n',imp);
   fprintf('sim            = %d\n',sim);
   fprintf('oob_prox       = %d\n',oob_prox);
   fprintf('nodesize       = %-5.4f\n',nodesize);
   fprintf('replace        = %d\n',replace);
   fprintf('=======================================================\n');
   fprintf(stndTxt);
end

% -----Options:-----
Options = [imp,localImp,nPerm];

% =================================
%     Call *.mexmaci64 MEX file:
% =================================
[lDaughter, rDaughter, nodestatus, nrnodes, upper, avnode, mbest, treeSize, Y_hat,...
   mse, impOut, impMat, impSD, prox, coef, oob_times, inbag]...
   = mexRF_train(X', Y, nTree, mTry, sampsize, nodesize, int32(Options), int32(cat),...
   int32(maxcat), int32(do_trace), int32(sim), int32(oob_prox), int32(biasCorr),...
   keep_inbag, replace);

% Cleanup:
clear mexRF_train

% Done in R file so doing it too:
Y_hat(oob_times==0)=NaN;

% Extract importance measures:
if (imp>0)
   meanAcc = impOut(:,1); % mean decrease in accuracy
   meanMSE = impOut(:,2); % mean decrease in MSE
else
   meanAcc = NaN;
   meanMSE = impOut(:)'; % if imp = 0, meanMSE is still returned
end

% Extract Standard Errors & calculate Z-scores:
if (imp>0)
   meanAccSD  = impSD;
   denom      = meanAccSD;         % make a copy
   denom(find(denom == 0)) = eps ; % prevent divide by 0 errors:
   zScore = meanAcc./denom;        % raw importance divided by its standard error
   zP     = normcdf(-abs(zScore),0,1); % after Matlab's ztest.m
   % zP  = tcdf(-abs(zScore), N-1);
else
   meanAccSD  = NaN;
   zScore     = NaN;
   zP         = NaN;
end

% Local importance:
if (localImp==0)
   impMat = NaN;
end

% -----Wrap results up into a structure:-----
model.nTree        = nTree; % # random trees grown
model.mTry         = mTry;  % size of random variable subsets
model.X            = X;     % transformed predictors
model.X_txt        = X_txt; % cell array of variable labels
model.stnd         = stnd;  % transformation applied to predictors
model.Y            = Y;     % observed response variable
model.Y_hat        = Y_hat; % fitted values of response variable
% 
if (sim>0)
   model.prox      = prox;   % symmetric proximity (similarity) matrix, based on
else                         % the frequency that pairs of observations occur in
   model.prox = [];          % the same terminal nodes
end
% 
model.localImp     = impMat';     % importance measures for each obs
model.meanAcc      = meanAcc*100; % mean decrease in accuracy    (col: variable)
model.meanAccSD    = meanAccSD;   % standard error
model.zScore       = zScore;      % z-score of meanAcc
model.zP           = zP;          % parametric significance of z-score
model.meanMSE      = meanMSE;     % mean decrease in node impurity (measured by MSE)
model.mse          = mse;         % mean squared error with each iteration (= RSS/n)
model.rsq          = 1 - mse / (var(Y) * (N-1) / N); % pseudo R-squared with each iteration
model.inbag        = inbag;     % in bag
model.oob_times    = double(oob_times); % # times obs are out-of-bag (& used for OOB error estimates)
model.nPerm        = nPerm;     % # times OOB data are permuted per tree for imp measures 
if (biasCorr>0)
   model.coef      = coef;      % coefficients for bias correction
else
   model.coef      = NaN;
end
% 
model.lDau         = lDaughter;
model.rDau         = rDaughter;
model.nodestatus   = nodestatus;
model.nrnodes      = nrnodes;
model.upper        = upper;
model.avnode       = avnode;
model.mbest        = mbest;
model.treeSize     = treeSize;
model.type         = 'regression';
