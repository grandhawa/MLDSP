function [model,err] = f_RFclass(X,Y,nTree,mTry,imp,sim,stnd,verb,X_txt,opt)
% - train a Random Forest classifier
%
% USAGE: [model,err] = f_RFclass(X,Y,nTree,mTry,imp,sim,'stnd',verb,X_txt,opt);
%
% X     = input data                              (rows = obs, cols = variables)
% Y     = column vector of integers specifying class membership
%                                          (Y = [] forces unsupervised learning)
% nTree = # random trees to grow             (default = 1000 if set to empty [])
% mTry  = # random variable subsets          (default =  sqrt # variables if [])
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
% .nTree        = # random trees grown
% .mTry         = size of random variable subsets
% .X            = transformed predictors
% .X_txt        = cell array of variable labels
% .stnd         = transformation applied to predictors in X
% .Y_old        = original class labels
% .Y            = new labels starting at 1
% .Y_hat        = predicted class for each observation
% .Y_txt        = cell array of class labels
% .nClass       = # classes
% .votes        = # votes each obs recieved by class
% .margin       = proportion of correct votes minus max incorrect
% .prox         = symmetric proximity (similarity) matrix based on frequency
%                 obs occur in the same node
% .localImp     = importance measures for each observation
% .classAcc     = class-specific mean decrease in accuracy (row:class,col:variable)
% .classAccSD   = standard error
% .meanAcc      = mean decrease in accuracy                          (col:variable)
% .meanAccSD    = standard error of meanAcc
% .zScore       = z-score of meanAcc
% .zP           = parametric significance of z-score
% .gIndex       = mean decrease in node impurity (= Gini Index)      (col:variable)
% .errtr        = error with each iteration
% .inbag        = in bag
% .oob_times    = # times obs are out-of-bag (& used for OOB error estimates)
% .classwt      = prior probabilities for each class
% .cutoff       =
% .nrnodes      =
% .xbestsplit   =
% .treemap      =
% .nodestatus   =
% .nodeclass    =
% .bestvar      =
% .ndbigtree    =
% .type         = type of Random Forest (= 'classification')
%
% err = structure of classification error rates with the following fields:
% .tot  = total error
% .grp  = total error for each class
% .conf = confusion matrix
%
% SEE ALSO: f_RFclassPredict, f_outlier, f_RFplot, f_RFvis, f_RFreg


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

% modified by David L. Jones, Aug-2009
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB'

% -----Notes:-----
% Options eliminated:
% corr_bias  : which happens only for regression ommitted
% norm_votes : always set to return total votes for each class

% edits by David L. Jones:
% - edited documentation so easier to read
% renamed 'extra_otions' to 'opt', replaced '|' with '||', replaced 'DEBUG_ON'
% with 'verb', lots of re-organization of the code, TRUE/FALSE replaced with
% 1/0, new code for labelling classes in Y, renamed nclass to nClass, create
% diagnostic plot, standardize/center X, calculate error rates/confusion matrix,
% calculate 'importance' by default, split 'importance' measures into
% corresponding variables, plot measures of importance.

% -----NOTES:----- (by David L. Jones)
% Used standardized variables (stnd = 1) to create a RF for DISCRIMINATION and
% assessing variable importance; use centered data (stnd = 2) to create a RF for
% CLASSIFYING unknowns.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                         OPTIONS                                     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opt = strucuture of extra options with the following fields (TRUE=1, FALSE=0):
%
% .replace = bootstrap sampling with replacement                   (default = 1)
%
% .classwt = priors of classes. Here the function first gets the labels in
%   ascending order and assumes the priors are given in the same order. So if
%   the class labels are [-1 1 2] and classwt is [0.1 2 3] then there is a 1-1
%   correspondence. (ascending order of class labels). Once this is set the freq
%   of labels in train data also affects.
%
% .cutoff = vector of length equal to number of classes; The 'winning' class for
%   an observation is the one with the maximum ratio of proportion of votes to
%   cutoff. Default is 1/k where k is the number of classes (i.e., majority vote
%   wins).
%
% .strata = (not yet stable in code) variable that is used for stratified
%   sampling. I don't yet know how this works. Disabled by default
%
% .sampsize = Size(s) of sample to draw. For classification, if sampsize is a
%   vector of the length the number of strata, then sampling is stratified by
%   strata, and the elements of sampsize indicate the numbers to be drawn from
%   the strata.
%
% .nodesize = Minimum size of terminal nodes. Setting this number larger causes
%   smaller trees to be grown (and thus take less time). Note that the default
%   values are different for classification (1) and regression (5).
%
% .localImp   = compute casewise importance measures               (default = 0)
%               (if = 1, then over-rides 'imp')
%
% .oob_prox  = calculate proximities only on 'out-of-bag' data     (default = 1)
%
% .do_trace  = If set to TRUE, give a more verbose output as randomForest is
%   run. If set to some integer, then running output is printed for every
%   do_trace trees.
%
% .keep_inbag = should an n by nTree matrix be returned that keeps track of
%   which samples are 'in-bag' in which trees (but not how many times, if
%   sampling with replacement)
% ------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                 DOCUMENTATION by David L. Jones                     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nTree: at least 1000 trees is recommended if requesting variable importance
%   measures and/or proximities. Breiman recommends nTree=5000 when there are
%   many variables and you want stable importance measures.
%
% meanAcc: this is 100 * the change in margins averaged over all classes; higher
%   values indicate a variable is more important. Variables can be ranked
%   according to the z-score of the meanAcc (zScore) and this, along with zP,
%   can be used to select a subset of the most important variables. Note that
%   Darby & Reissland (1981) used z-scores in variable selection, not as a test
%   statistic for a significance test, but rather as a numerical guide to
%   selecting variables. Darby, S.C., Reissland, J.A. (1981) "Low levels of
%   ionizing radiation and cancer ? are we underestimating the risk? (with
%   discussion)". Journal of the Royal Statistical Society, Series A, 144(3),
%   298-331.
% ------------------------------------------------------------------------------

% Oct-2010: added error message when IMP=1 & SIM=1
% Aug-2011: replaced 'unique' with 'f_unique'
% Oct-2011: localImp set to imp by default; added '.type' field in output
% May-2012: set errProp=0 when idx is empty
% May-2013: output now returns normalized votes

% -----Check input & set defaults:-----
if (nargin < 2), Y     = [];    end % default unsupervised learning
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
if isempty(Y)
   Y        = ones(N,1);
   addclass = 1; % flag indicating unsupervised learning
else
   addclass = 0;
end
if isempty(nTree), nTree = 1000;           end
if isempty(mTry),  mTry  = floor(sqrt(D)); end

% Check mTry:
if (mTry<1 || mTry>D), error('MTRY must range from 1 to D!'); end
mTry = max(1,min(D,round(mTry)));

% Check for missing data:
if any(isnan(X)); error('NaNs in X'); end
if any(isnan(Y)); error('NaNs in Y'); end

% Check for compatible sizes:
if (size(X,1) ~= size(Y,1)), error('Size mismatch b/n X & Y'); end
if size(Y,2)~=1, error('Y should be a column vector!'); end

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

% Unsupervised learning:
if (addclass==1)
   Y   = [Y;Y+1];
   X   = [X;X];
   sim = 1;
end

% Create new class labels:
Y_txt    = cellstr(num2str(f_unique(Y))); % create class labels for f_RFclassPlot
Y_old    = Y;                             % keep copy of old labels
Y        = f_recode(Y);                   % recode groups to start at 1
grps_old = f_unique(Y_old);               % unique groups (old)
grps     = f_unique(Y);                   % unique groups
nClass   = size(grps(:),1);               % # of classes

% Need at least 2 classes (even with unsupervised learning):
if (nClass<2), error('Y must specify at least 2 classes!'); end

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
   fprintf('for %d classes...\n\n',nClass);
end


% -----Parse extra options:-----
if exist('opt','var')
   if isfield(opt,'replace'),    replace    = opt.replace;    end
   if isfield(opt,'classwt'),    classwt    = opt.classwt;    end
   if isfield(opt,'cutoff'),     cutoff     = opt.cutoff;     end
   if isfield(opt,'strata'),     strata     = opt.strata;     end
   if isfield(opt,'sampsize'),   sampsize   = opt.sampsize;   end
   if isfield(opt,'nodesize'),   nodesize   = opt.nodesize;   end
   if isfield(opt,'localImp'),   localImp   = opt.localImp;   end
   if isfield(opt,'oob_prox'),   oob_prox   = opt.oob_prox;   end
   if isfield(opt,'do_trace'),   do_trace   = opt.do_trace;   end
   if isfield(opt,'keep_inbag'), keep_inbag = opt.keep_inbag; end
   if isfield(opt,'coef'),       coef       = opt.coef(:);    end % [added by DLJ]
   if isfield(opt,'cut'),        cut        = opt.cut;        end % [added by DLJ]
   
end
keep_forest = 1; % always save the trees :)

% -----Set defaults for options not already specified:-----
if ~exist('replace','var'),  replace  = 1;  end
if ~exist('sampsize','var');
   if (replace)
      sampsize = N;             % bootstrap with replacement
   else
      sampsize = ceil(0.632*N); % bootstrap without replacement
   end;
end
if ~exist('nodesize','var'),   nodesize    = 1;   end % classification = 1
if ~exist('localImp','var'),   localImp    = imp; end
if ~exist('oob_prox','var'),   oob_prox    = sim; end
if ~exist('do_trace','var'),   do_trace    = 0;   end
if ~exist('keep_inbag','var'), keep_inbag  = 0;   end
if ~exist('coef','var'),       coef        = ones(nClass,1); end
if ~exist('cut','var'),        cut         = NaN; end

% Check input for coef:
if length(coef)~=nClass, error('Size mismatch b/n COEF and nClass!'); end

% When localImp = 1, imp must = 1:
if (localImp>0), imp = 1; end

%  Make sure nTree is sufficient:
if ( nTree<1000 && (imp==1 || sim==1) )
   error('Robust importances &/or proximities require at least 1000 trees!')
end

switch cut
   case 'prop'% set cutoff proportional to class size:
      cutoff = zeros(nClass,1); % preallocate
      for i = 1:nClass
         cutoff(i) = sum(Y==grps(i))/N;
      end
                 
      % Adjust proportions according to coefficients (higher values give less weight
      % to smaller classes):
      cutoff = (cutoff .* coef) ./ sum(cutoff .* coef);
      cutoff = f_round(cutoff,4); % round off so sum=1;
end

% -----Categories:----- [??? DLJ]
% now handle categories. Problem is that categories in R are more
% enhanced. In this i ask the user to specify the column/features to
% consider as categories, 1 if all the values are real values else
% specify the number of categories here
if exist ('opt','var') && isfield(opt,'categories')
   ncat = opt.categories;
else
   ncat = ones(1,D);
end

maxcat = max(ncat);
if maxcat>32
   error('Can not handle categorical predictors with more than 32 categories');
end
% ----------------------


% -----Cutoff & Class Weights:-----
% classRF - line 88 in randomForest.default.R
if ~exist('cutoff','var')
   cutoff = ones(1,nClass)* (1/nClass);
else
   if sum(cutoff)>1 || sum(cutoff)<0 || length(find(cutoff<=0))>0 || length(cutoff)~=nClass
      error('Incorrect cutoff specified');
   end
end
if ~exist('classwt','var')
   classwt = ones(1,nClass);
   ipi     = 0;
else
   if length(classwt)~=nClass
      error('Length of classwt not equal to the number of groups')
   end
   if ~isempty(find(classwt<=0))
      error('classwt must be positive');
   end
   ipi = 1;
end


if (addclass>0)
   nsample = 2*N; % unsupervised learning
else
   nsample = N;
end

% -----Stratified data:-----
Stratify = (length(sampsize)>1);
if (~Stratify && sampsize>N)
   error('Sampsize too large')
end

if Stratify
   if ~exist('strata','var')
      strata = Y;
   end
   nsum = sum(sampsize);
   if ( ~isempty(find(sampsize<=0)) || nsum==0)
      error('Bad sampsize specification');
   end
else
   nsum = sampsize;
end

if Stratify
   strata = int32(strata);
else
   strata = int32(1);
end


% -----Options:-----
Options = int32([addclass,imp,localImp,sim,oob_prox,do_trace,...
   keep_forest,replace,Stratify,keep_inbag]);

% Verbose reporting of inputs:
if (verb>0)
   % Print the parameters that i am sending in:
   fprintf('=======================================================\n');
   fprintf('=====     INPUT PARAMETERS FOR RANDOM FOREST:     =====\n');
   fprintf('=======================================================\n');
   fprintf('nTree          = %d\n',nTree);
   fprintf('mTry           = %d\n',mTry);
   fprintf('X: # obs = %d, # var''s = %d \n',N,D);
   fprintf('# Y classes    = %d \n',nClass);
   fprintf('addclass       = %d\n',addclass);
   fprintf('size(ncat)     = %d\n',size(ncat));
   fprintf('maxcat         = %d\n',maxcat);
   fprintf('size(sampsize) = %d\n',size(sampsize));
   fprintf('sampsize[0]    = %d\n',sampsize(1));
   fprintf('stratify       = %d\n',Stratify);
   fprintf('imp            = %d\n',imp);
   fprintf('sim            = %d\n',sim);
   fprintf('oob_prox       = %d\n',oob_prox);
   fprintf('strata         = %d\n',strata);
   fprintf('ipi            = %d\n',ipi);
   fprintf('classwt        = %-5.4f\n',classwt);
   fprintf('cutoff         = %-5.4f\n',cutoff);
   %fprintf('coef           = %-5.4f\n',coef);
   fprintf('nodesize       = %-5.4f\n',nodesize);
   fprintf('=======================================================\n');
   fprintf(stndTxt);
end

% -----Call mex file:-----
[nrnodes,nTree,xbestsplit,classwt,cutoff,treemap,nodestatus,nodeclass,bestvar,...
   ndbigtree,mTry,Y_hat,counttr,prox,impMat,impOut,impSD,errTr,inbag]...
   = mexClassRF_train(X',int32(Y),nClass,nTree,mTry,int32(ncat),...
   int32(maxcat),int32(sampsize),strata,Options,int32(ipi),classwt,cutoff,...
   int32(nodesize),int32(nsum),int32(N),int32(D),int32(nsample));
% ------------------------

% Format output:
errTr = errTr';

% Format predicted values:
Y_hat = double(Y_hat);

% Translate labels to original format (new -> old):
YY = repmat(NaN,size(Y_hat)); % preallocate
for i = 1:nClass
   YY(Y_hat==grps(i)) = grps_old(i);
end
%
Y_hat = YY; % replace with new labels


% Extract importance measures:
if (imp>0)
   classAcc  = impOut(:,1:nClass)'; % class specific mean decrease in accuracy
   meanAcc   = impOut(:,nClass+1)'; % mean decrease in accuracy
   gIndex    = impOut(:,nClass+2)'; % mean decrease in Gini index
else
   classAcc  = NaN;
   meanAcc   = NaN;
   gIndex    = impOut(:)';          % if importance = 0, gIndex still returned
end

% Extract Standard Errors & calculate Z-scores:
if (imp>0)
   classAccSD = impSD(:,1:nClass)';
   meanAccSD  = impSD(:,nClass+1)';
   denom      = meanAccSD;         % make a copy
   denom(find(denom == 0)) = eps ; % prevent divide by 0 errors:
   zScore = meanAcc./denom;        % raw importance divided by its standard error
   zP    = normcdf(-abs(zScore),0,1); % after Matlab's ztest.m
   % zP  = tcdf(-abs(zScore), N-1);
else
   classAccSD = NaN;
   meanAccSD  = NaN;
   zScore     = NaN;
   zP         = NaN;
end

% Local importance:
if (localImp==0)
   impMat = NaN;
end

% Votes & Margin of Predictions:
votes = double(counttr');
nVotes = votes./repmat(sum(votes,2),1,size(votes,2)); % normalize by row
if (addclass==0)
   margin = repmat(NaN,N,1); % preallocate
   for i=1:N
      corrProp  = nVotes(i,Y(i));        % correct proportion
      idx       = find((nVotes(i,:) ~= corrProp)==1);
      if ~isempty(idx) % DLJ: fix for when all are correct?
         errProp   = max(nVotes(i,idx));    % maximum of incorrect proportions
      else
         errProp = 0;
      end
      margin(i) = corrProp - errProp(1); % difference
   end
else
   margin = NaN; % unsupervised learning
end

% -----Wrap results up into a structure:-----
model.nTree       = nTree;   % # random trees grown
model.mTry        = mTry;    % size of random variable subsets
model.X           = X;       % transformed predictors
model.X_txt       = X_txt;   % cell array of variable labels
model.stnd        = stnd;    % transformation applied to predictors
model.Y_old       = Y_old;   % original class labels
model.Y           = Y;       % new labels starting at 1
model.Y_hat       = Y_hat;   % predicted class for each observation
model.Y_txt       = Y_txt;   % cell array of class labels for plotting
model.nClass      = nClass;  % # classes
model.nVotes      = nVotes;  % normalized votes each obs received by class
model.margin      = margin;  % proportion of correct votes minus max incorrect;
%                              +/- margins = correct/incorrect classifications
if (sim>0)
   model.prox = prox;        % symmetric proximity (similarity) matrix, based on
else                         % the frequency that pairs of observations occur in
   model.prox = NaN;         % the same terminal nodes
end
model.localImp   = impMat';    % importance measures for each obs
model.classAcc   = classAcc;   % class-specific mean decrease in accuracy
model.classAccSD = classAccSD; % standard error
model.meanAcc     = meanAcc;   % mean decrease in accuracy
model.meanAccSD   = meanAccSD; % standard error
model.zScore      = zScore;    % z-score of meanAcc
model.zP          = zP;        % parametric significance of z-score
model.gIndex      = gIndex;    % mean decrease in node impurity (= Gini Index)
model.errtr       = errTr;     % error with each iteration
model.inbag       = inbag;
model.oob_times   = sum(counttr)'; % times obs are out-of-bag (& used for OOB error estimates)
%
model.classwt     = classwt;       % prior probabilities for each class
model.cutoff      = cutoff;
model.nrnodes     = nrnodes;
model.xbestsplit  = xbestsplit;
model.treemap     = treemap;
model.nodestatus  = nodestatus;
model.nodeclass   = nodeclass;
model.bestvar     = bestvar;
model.ndbigtree   = ndbigtree;
model.type        = 'classification';
% -------------------------------------------

% Clean up:
clear mexClassRF_train

% -----Trim results for Unsupervised Learning:-----
if (addclass>0)
   model.X      = model.X(1:N,:);
   model.Y_old  = model.Y_old(1:N,:);
   model.Y      = model.Y(1:N,:);
   model.Y_hat  = model.Y_hat(1:N,:);
   model.Y_txt  = {'1'};
   model.nClass = 1;
end
% -------------------------------------------------

% Classification error rates:
if (addclass==0)
   err = f_errRate(Y_old,Y_hat);
else
   err = NaN; % unsupervised learning
end

% -----Send output to display:-----
if (verb>0 && addclass==0)
   fprintf('Ave. error rate (after %d iterations) = %-5.4f',nTree,errTr(end,1));
   
   fprintf('\n==================================================\n');
   fprintf('       RANDOM FOREST ''Majority Rules'' \n'           );
   fprintf('  Internal Cross-validation Classification Success: \n');
   fprintf('--------------------------------------------------\n' );
   
   fprintf('Class        Corrrect  \n');
   for j=1:nClass
      
      fprintf('%s %d %s %10.2f %s \n',['  '],grps_old(j),['     '],[1-err.grp(j)]*100,['%']);
   end
   
   fprintf('\n\n');
   fprintf('Total Correct  = %4.2f %s \n',(1-err.tot)*100,['%']);
   fprintf('Total Error    = %4.2f %s \n',err.tot*100,['%']);
   fprintf('Prior prob     = class size \n');
   
   fprintf('\n--------------------------------------------------\n' );
   fprintf('     Confusion Matrix (%s): \n',['%'])
   hdr = sprintf('%6.0f ',grps_old(:)');
   fprintf(['class: ' hdr]);
   fprintf('\n')
   for j=1:nClass
      txt = [sprintf('%6.0f ',grps_old(j)) sprintf('%6.2f ',err.conf(j,:)*100)];
      fprintf(txt);
      fprintf('\n')
   end
   fprintf('==================================================\n\n');
end

