function anova = f_rdaDB_manova(yDis,X,iter,verb)
% - nonparametric (permutation-based) MANOVA via db-RDA
%
% USAGE: anova = f_rdaDB_manova(yDis,X,iter,verb);
%
% yDis  = square symmetric distance matrix derived from response variables
% X     = matrix of integers specifying group membership for objects in yDis
% iter  = # iterations for permutation test  (default = 0)
% verb  = optionally send results to display (default = 1)
%
% SEE ALSO: f_rdaDB, f_npManova, f_rda

% -----Author:-----
% by David L. Jones, Oct-2012
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Notes:-----
% This function performs MANOVA using db-RDA by re-coding the matrix of group
% membership (X) as orthogonal contrast codes. For 2 or more factors, each
% factor is tested separately via a partial RDA analysis by setting it as the
% explanatory variable with the remaining factors serving as covariables.

n   = size(yDis,1); % # observations
ncX = size(X,2);    % # factors

% -----Check input & set defaults:-----
if (nargin < 3), iter = 0; end % default iterations for permutation test
if (nargin < 4), verb = 1; end % send output to display by default

if n ~= size(X,1), error('yDis & X need same # of rows'); end

if (ncX>3)
   error('Only up to 3 factor analyses are currently supported');
end

% Check imput if this isn't a permutation run:
if (f_issymdis(yDis) == 0)
   error('Input yDis must be a square symmetric distance matrix!');
end

% Determine if data has replication:
if (size(unique(X,'rows'),1) < n), rep = 1; else rep = 0; end;
% -------------------------------------

%===============================================================
% Determine which model to use:
if (ncX > 1)
   switch ncX
      case 2
         fprintf('\n========================================\n');
         fprintf('  Please specify the ANOVA model \n');
         fprintf('  for 2-way ANOVA factors 1 & 2:  \n');
         fprintf('-----------------------------------------\n');
         fprintf('1 & 2 are fixed                    [21] \n');
         fprintf('1 is fixed or random, 2 is random  [22] \n');
         fprintf('1 is fixed or random, 2 is nested  [23] \n\n');

      case 3
         fprintf('\n===========================================================\n');
         fprintf('  Please specify the ANOVA model    \n');
         fprintf('  for 3-way ANOVA factors 1, 2, & 3: \n');
         fprintf('-----------------------------------------------------------\n');
         fprintf('All factors fixed                                     [31] \n');
         fprintf('1 & 2 are fixed, 3 is random                          [32] \n');
         fprintf('1 is fixed or random, 2 & 3 are random                [33] \n\n');
         %
         fprintf('1 & 2 fixed,     3 nested in 1        (Cross-Nested)  [34] \n');
         fprintf('1 &/or 2 random, 3 nested in 1        (Cross-Nested)  [35] \n');
         fprintf('3 nested in 2 nested in 1             (Fully Nested)  [36] \n\n');
      otherwise
         error('Too many ANOVA factors...reduce to 3 or less \n');
   end

   % Get user to specify which ANOVA model to use:
   inputText = ['Select model...[O will cancel]' ['\n']];
   model     = input(inputText);
   if (model==0),error('User selected 0 to cancel and exit \n');end


   % -----Check user-selected ANOVA model:-----
   % Determine if model is nested:
   nestedModels =[23 34:36];
   if ((find(model == nestedModels))>0)
      nest = 1;
   else
      nest = 0;
   end

   if (rep<1) && (nest==1)
      error('Nested designs require REPLICATION in your data')
   end

   % Check input against list of current models:
   if (ncX == 2)
      supportedModels = 21:23;
   else
      supportedModels = 31:36;
   end
   if (sum(find(model == supportedModels))==0)
      error([num2str(model) ' is an invalid ANOVA Model spec...try again.']);
   end
   % ------------------------------------------
end
%===============================================================

switch ncX
   case 1 % 1-Way ANOVA:
      result = f_rdaDB_manova1(yDis,X,0,iter); % 1-way ANOVA via db-RDA
      anova  = sub_format(result,1);           % format as an ANOVA table
      
   case 2 % 2-way ANOVA:
      
      if (model == 21) || (model == 22) 
         anova = f_rdaDB_manova2(yDis,X(:,1),X(:,2),iter,model);
         
         
      elseif (model == 23) % NESTED:
         A        = f_xMatrix(X(:,1),1);        % factor 1
         Bn       = f_xMatrix(X(:,1:2),1,1);    % factor 2 as nested
         
         % This produces the wrong F-stat (see: f_npManova2n.m):
         result.A = f_rda(Y,A,[Bn],iter,0,0,0);
         
         result.B = f_rda(Y,Bn,[A],iter,0,0,0);
         anova    = sub_format(result,-2);      % format as nested anova table
      else
         error('That''s currently not supported')
      end
      
   otherwise
      error('Only 1 & 2-way ANOVA''s are supported');
end



% -----Send output to display:-----
if (verb>0)
   headerCell = {'Source','df','SS','MS','F','p'}; % Set up column labels
   resultCell = f_struct2flat(anova); % flatten struct to cell array
   resultCell = [headerCell' resultCell']'; % concat

   fprintf('\n==================================================\n');
   fprintf('    Permutation-based MANOVA via db-RDA:\n');
   fprintf('--------------------------------------------------\n');
   disp(resultCell);
   fprintf('\n %20s %10d \n',['# iterations = '],iter);
   fprintf('--------------------------------------------------\n');
   fprintf('(Note: NaNs are placeholders for the ANOVA table)\n\n');
   if (rep>0)
      fprintf('(Data with replication) \n')
   else
      fprintf('(Data has NO replication) \n')
   end
end
% ---------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%       SUBFUNCTION     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function anova = sub_format(result,fmt)
% - format db-RDA results as an ANOVA table
switch fmt
   case 1 % 1-way ANOVA:
      anova(1).so = {'factor 1'};
      anova(2).so = {'residual'};
      anova(3).so = {'total'};

      anova(1).df = result.dfR;
      anova(2).df = result.dfE;
      anova(3).df = result.dfT;

      anova(1).SS = result.SSr;
      anova(2).SS = result.SSe;
      anova(3).SS = result.SSt;

      anova(1).MS = result.MSr;
      anova(2).MS = result.MSe;
      anova(3).MS = NaN;

      anova(1).F = result.F;
      anova(2).F = NaN;
      anova(3).F = NaN;

      anova(1).p = result.p;
      anova(2).p = NaN;
      anova(3).p = NaN;

   case 2 % 2-way ANOVA:
      anova(1).so = {'factor 1'};
      anova(2).so = {'factor 2'};
      anova(3).so = {'factor 1x2'};
      anova(4).so = {'residual'};
      anova(5).so = {'total'};

      anova(1).df = result.A.dfR;  % df factor 1
      anova(2).df = result.B.dfR;  % df factor 2
      anova(3).df = result.AB.dfR; % df interaction
      anova(4).df = result.A.dfE;  % df error
      anova(5).df = result.A.dfT;  % df total

      anova(1).SS = result.A.SSr;  % SS factor 1
      anova(2).SS = result.B.SSr;  % SS factor 2
      anova(3).SS = result.AB.SSr; % SS interaction
      anova(4).SS = result.A.SSe;  % SS error
      anova(5).SS = result.A.SSt;  % SS total

      anova(1).MS = result.A.MSr;  % MS factor 1
      anova(2).MS = result.B.MSr;  % MS factor 2
      anova(3).MS = result.AB.MSr; % MS interaction
      anova(4).MS = result.A.MSe;  % MS error
      anova(5).MS = NaN;

      anova(1).F  = result.A.F;    % F-ratio factor 1
      anova(2).F  = result.B.F;    % F-ratio factor 2
      anova(3).F  = result.AB.F;   % F-ratio interaction
      anova(4).F  = NaN;
      anova(5).F  = NaN;
      
      anova(1).p = result.A.p;     % p-value factor 1
      anova(2).p = result.B.p;     % p-value factor 2
      anova(3).p = result.AB.p;    % p-value interaction
      anova(4).p = NaN;
      anova(5).p = NaN;

   case -2 % NESTED 2-way ANOVA:
      anova(1).so = {'factor 1'};
      anova(2).so = {'factor 2'};
      anova(3).so = {'residual'};
      anova(4).so = {'total'};

      anova(1).df = result.A.dfR;  % df factor 1
      anova(2).df = result.B.dfR;  % df factor 2
      anova(3).df = result.A.dfE;  % df error
      anova(4).df = result.A.dfT;  % df total

      anova(1).SS = result.A.SSr;  % SS factor 1
      anova(2).SS = result.B.SSr;  % SS factor 2
      anova(3).SS = result.A.SSe;  % SS error
      anova(4).SS = result.A.SSt;  % SS total

      anova(1).MS = result.A.MSr;  % MS factor 1
      anova(2).MS = result.B.MSr;  % MS factor 2
      anova(3).MS = result.A.MSe;  % MS error
      anova(4).MS = NaN;

      anova(1).F  = result.A.F;    % F-ratio factor 1
      anova(2).F  = result.B.F;    % F-ratio factor 2
      anova(3).F  = NaN;
      anova(4).F  = NaN;
      
      anova(1).p = result.A.p;     % p-value factor 1
      anova(2).p = result.B.p;     % p-value factor 2
      anova(3).p = NaN;
      anova(4).p = NaN;
   
    
      
      % case 3 % 3-way ANOVA:
      
   otherwise
      error('Unsupported format')
end
