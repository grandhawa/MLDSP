function result = f_npManova(yDis,x,iter,verb,model)
% - nonparametric (permutation-based) MANOVA for ANY distance matrix
%
% USAGE: result = f_npManova(yDis,x,iter,verb,model)
%
% yDis  = square symmetric dissimilarity matrix derived from response variables
% x     = matrix of integers specifying group membership for objects in yDis
%         (column-wise)
% iter  = # iterations for permutation test  (default = 0)
% verb  = optionally send results to display (default = 1)
% model = optionally bypass user-prompt to specify model
%
% result.so = source of variation
% result.df = degrees of freedom
% result.SS = sum of squares
% result.MS = mean square
% result.F  = F-statistics
% result.p  = permutation-based significance probabilities
%
% See also: f_npManovaPW, f_rdaAnova, f_anosim, f_anosim2, f_distlm

% -----References:-----
% Anderson, M. J. 2000. NPMANOVA: a FORTRAN computer program for non-parametric
%  multivariate analysis of variance (for any two-factor ANOVA design) using
%  permutation tests. Dept. of Statistics, University of Auckland.
%  (http://www.stat.auckland.ac.nz/PEOPLE/marti/)
% Anderson, M. J. 2001. A new method for non-parametric multivariate
%   analysis of variance. Austral Ecology 26: 32-46.
% Anderson, M. J. 2002. DISTML v.2: a FORTRAN computer program to calculate a
%   distance-based multivariate analysis for a linear model. Dept. of Statistics
%   University of Auckland. (http://www.stat.auckland.ac.nz/PEOPLE/marti/)
% McArdle, B. H. and M. J. Anderson. 2001. Fitting multivariate models to
%   community data: a comment on distance-based redundancy analysis. Ecology
%   290-297.
% Sokal, R. R. and F. J. Rohlf. 1995. Biometry - The principles and
%   practice of statistics in bioligical research. 3rd ed. W. H.
%   Freeman, New York. xix + 887 pp.
% Underwood, A. J. 1981. Techniques of analysis of variance in experimental
%   marine biology and ecology. Oceanogr. Mar. Biol. Ann. Rev. 19: 513-605.
% Zar, J. H. 1999. Biostatistical analysis. 4th ed. Prentice Hall, Upper Saddle
%   River, NJ.

% -----Notes:-----
% Special care must be taken when coding levels of nested factors, i.e.
% treatment levels of a nested factor must not be REPEATED across different
% levels of the main factor. For example:
%
% main factor   = [1 1 1 2 2 2 3 3 3]'
% nested factor = [1 2 3 1 2 3 1 2 3]' < don't
% nested factor = [1 2 3 4 5 6 7 8 9]' < do
%
% More examples are in the Appendix of the User's Manual.
%
% This program takes a regression-approach to ANOVA using the General Linear
% Model (GLM) and constructs F-ratios using the UNRESTRICTED form of the model.
% The F-ratios used for each type of test are provided in the Appendix of the
% User's Manual. Some of these differ somewhat from those used in textbook
% examples, especially for balanced, mixed-model designs, but are the same as
% those used in most computer programs that use GLM (e.g., SAS and MINITAB).

% -----Author:-----
% by David L. Jones, May-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Jul-2002: improved display of results
% Aug-2002: results are now returned as a structure, but displayed
%           in the console as a 'flattened' cell array
% Nov-2002: added support for fixed, random, and nested factors
% Jan-2008: changed '&' to '&&' and '|' to '||', added support to specify ANOVA
%           model on the command line
% Mar-2008: removed RANK option

if (nargin < 3), iter   = 0; end; % default iterations for permutation test
if (nargin < 4), verb   = 1; end; % send output to display by default
if (nargin < 5), model = -1; end; % prompt user to specify ANOVA model

n = size(yDis,1);

% -----Check input:-----
if n ~= size(x,1), error('yDis & X need same # of rows'); end;

if (f_issymdis(yDis) == 0)
   error('Input yDIS must be a square symmetric distance matrix');
end;

if (size(x,2)>3)
   error('Only up to 3 factor analyses are currently supported');
end

% Determine if data has replication:
if (size(unique(x,'rows'),1) < n), rep = 1; else rep = 0; end;


%===============================================================
noFactors = (size(x,2)); % get # of input factors

% Determine which model to use:
if (noFactors > 1)
   switch noFactors
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

   % Get user to specify which ANOVA model to use if not already specified:
   if (model == -1)
      inputText = ['Select model...[O will cancel]' ['\n']];
      model     = input(inputText);
      if (model==0),error('User selected 0 to cancel and exit \n');end;
   end

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
   if (noFactors == 2)
      supportedModels =[21:23];
   else
      supportedModels =[31:36];
   end
   if (sum(find(model == supportedModels))==0)
      error([num2str(model) ' is an invalid ANOVA Model spec...try again.']);
   end
   % ------------------------------------------
end
%===============================================================

A   = -0.5*(yDis.^2);
I   = eye(n,n);
uno = ones(n,1);
G   = (I-(1/n)*(uno*uno'))*A*(I-(1/n)*(uno*uno')); % Gower's centered matrix

ncX = size(x,2); % # factors

switch ncX
   case 1; % 1-Way MANOVA:

      fac = f_designMatrix(x,1);      % create ANOVA design matrix
      fac = [ones(n,1) fac];          % add intercept
      [Q1,R1] = qr(fac,0); H= Q1*Q1'; % compute Hat-matrix via QR
      m = size(fac,2);                % # of parameters in design matrix

      % just keep the result:
      [null1,null2,null3,null4,null5,result] = f_npManova1(n,m,H,I,G,iter,verb);


   case 2 % 2-WAY ORTHOGONAL MANOVA:

      %-----Main effects (factors 1 & 2):-----
      for i=1:2
         fac{i} = f_designMatrix(x(:,i),1);     % create ANOVA design matrix
         fac{i} = [ones(n,1) fac{i}];           % add intercept
         [Q1,R1] = qr(fac{i},0); H{i} = Q1*Q1'; % compute Hat-matrix via QR
         m(i) = size(fac{i},2);                 % # of parameters in design matrix
      end

      %-----Interaction term (1x2):-----
      [ignore,fac{3}] = f_designMatrix(x,1);    % create ANOVA design matrix
      fac{3} = [ones(n,1) fac{3}];              % add intercept
      [Q1,R1] = qr(fac{3},0); H{3} = Q1*Q1';    % compute Hat-matrix via QR
      m(3) = size(fac{3},2);                    % # of parameters in design matrix

      if (rep<1) || (model == 23)
         if (rep<1) %-----1,2:----- (no replication)
            result = f_npManova2n(n,m,H,I,G,iter,0,model); % SS.resid = interaction (no replication)
         else %-----1,2:-----       (1 main, 2 nested)
            nNestLevels = size(unique(x(:,2)),1); % # of levels of nested factor
            result = f_npManova2n(n,m,H,I,G,iter,1,model,nNestLevels);
         end
      else %-----1,2,1x2:----       (replication)
         result = f_npManova2(n,m,H,I,G,iter,model);
      end


   case 3 % 3-WAY ORTHOGONAL MANOVA:

      %-----Main effects (factors 1-3):-----
      for i=1:3
         fac{i} = f_designMatrix(x(:,i),1);     % create ANOVA design matrix
         fac{i} = [ones(n,1) fac{i}];           % add intercept
         [Q1,R1] = qr(fac{i},0); H{i} = Q1*Q1'; % compute Hat-matrix via QR
         m(i) = size(fac{i},2);                 % # of parameters in design matrix
      end;


      %-----Interaction terms (1x2,1x3,2x3,1x2x3):-----
      [ignore,iTerms] = f_designMatrix(x,1);    % design matrices specifying interactions
      for i=4:7
         fac{i} = [ones(n,1) iTerms{i-3}];      % add intercept
         [Q1,R1] = qr(fac{i},0); H{i} = Q1*Q1'; % compute Hat-matrix via QR
         m(i) = size(fac{i},2);                 % # of parameters in design matrix
      end;

      if (rep<1) || (nest == 1)
         if (model == 34)|| (model == 35)  % CROSS-NESTED Model
            nNestLevels = size(unique(x(:,3)),1); % # of levels of nested factor
            result = f_npManova3Nest1(n,m,H,I,G,iter,1,model,nNestLevels);

         elseif (model == 36) % FULLY NESTED Model
            nNestLevels(1) = size(unique(x(:,2)),1); % # of levels of 1st nested factor
            nNestLevels(2) = size(unique(x(:,3)),1); % # of levels of 2nd nested factor
            result = f_npManova3Nest2(n,m,H,I,G,iter,1,model,nNestLevels);

         else % NO REPLICATION
            %-----1,2,3,1st-order interactions: (no replication)
            result = f_npManova3noRep(n,m,H,I,G,iter,model); % SS.resid = 2nd-order interaction
         end

      else
         %-----1,2,3,1st-order,2nd-order interactions: (replication)
         result = f_npManova3(n,m,H,I,G,iter,model);
      end

   otherwise
      error('Only 1, 2, & 3-way ANOVA''s are supported');
end



% -----Send output to display:-----
if (verb>0)
   headerCell = {'Source','df','SS','MS','F','p'}; % Set up column labels
   resultCell = f_struct2flat(result); % flatten struct to cell array
   resultCell = [headerCell' resultCell']'; % concat

   fprintf('\n==================================================\n');
   fprintf('    Nonparametric (Permutation-based) MANOVA:\n');
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

