function result = f_npManova2n(n,m,H,I,G,iter,rep,model,nNestLevels)
% - utility function called by f_npManova
% - 2-way MANOVA with no replication or nested factors.
% 
% n          = # rows/colums in distance matrix
% m          = array of # of parameters for each factor
% H          = cell array of hat matrix for each factor
% I          = I matrix
% G          = Gower's centered matrix
% iter       = # iterations for permutation test
% rep        = data with replication (=1) or without (=0)
% model      = ANOVA model specifying factor types
% nNestLevels = # of levels of nested factor

% -----Notes:-----
% This function is similar to f_npManova2, but is for data NO REPLICATION or
% NESTED factors. With no replication, there is no SS.residual term so the
% interaction term is used instead. As a consequence you cannot test for an
% interaction. Use 'Tukey's test for non-additivity' to test for interactions. 
% With nested factors there is no interaction term and replication is REQUIRED.

% -----Author:-----
% by David L. Jones, 
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Aug-2002: formatted output as a structure
% Nov-2002: added support for fixed, random, and nested factors
% Mar-2008: removed loop for 'Factors 1 & 2'; added trace

% Specify sources of variation:
result(1).so = {'factor 1'};
result(2).so = {'factor 2'};
result(3).so = {'residual'};
result(4).so = {'total'};

if (iter<1),
   iter = 1; % do at least once
else
   Fstat1 = zeros(iter,1); % preallocate results array
   Fstat2 = zeros(iter,1);
   fprintf('\nPermuting the data %d times...\n',iter-1);
end

for k = 1:iter

   if (k==1);
      Gvar = G;              % use observed G for 1st iteration
   else
      Gvar = f_shuffle(G,2); % permute square symmetric matrix
   end


   % Factors 1 & 2:
   [df(1),SS(1),MS(1)] = f_npManova1(n,m(1),H{1},I,Gvar,0);
   [df(2),SS(2),MS(2)] = f_npManova1(n,m(2),H{2},I,Gvar,0);


   if (rep>0) % WITH REPLICATION (factor 2 is nested, so no interaction)

      % Adjust degrees of freedom: (Sokal & Rohlf, 1995:293):
      df(2).among = nNestLevels - (df(1).among+1);
      residual_df = n - nNestLevels;
      df(1).total = n-1;

      % Adjust SS & MS of nested factors: (Sokal & Rohlf, 1995:276-277):
      SS(2).among = SS(2).among - SS(1).among;
      MS(2).among = SS(2).among/df(2).among;

      % Adjust error terms:
      residual_SS = SS(1).total - SS(1).among - SS(2).among;
      residual_MS = residual_SS/residual_df;

   else % NO REPLICATION (interaction term used as the error term):

      % Interaction term:
      SS(3).among = trace(H{3}*Gvar*H{3});
      residual_SS = SS(3).among - SS(1).among - SS(2).among; % SS interaction
      residual_df = df(1).among * df(2).among;               % df interaction
      residual_MS = residual_SS/residual_df;                 % MS interaction
   end

   %=======================================================

   % Compute pseudo-F's:

   if (rep<1) % No replication
      Fstat1(k) = MS(1).among/residual_MS;
      Fstat2(k) = MS(2).among/residual_MS;

   else	% Factor 1 fixed or random, factor 2 nested within 1 [WITH REPLICATION]
      Fstat1(k) = MS(1).among/MS(2).among;
      Fstat2(k) = MS(2).among/residual_MS;

   end
   %=======================================================

   % collect observed values:
   if (k==1)
      result(1).df = df(1).among;
      result(2).df = df(2).among;
      result(3).df = residual_df;
      result(4).df = df(1).total;

      result(1).SS = SS(1).among;
      result(2).SS = SS(2).among;
      result(3).SS = residual_SS;
      result(4).SS = SS(1).total;

      result(1).MS = MS(1).among;
      result(2).MS = MS(2).among;
      result(3).MS = residual_MS;
      result(4).MS = NaN;

      result(1).F  = Fstat1(1);
      result(2).F  = Fstat2(1);
      result(3).F  = NaN;
      result(4).F  = NaN;
   end
end

if (iter==1)
   result(1).p = NaN;
   result(2).p = NaN;
   result(3).p = NaN;
   result(4).p = NaN;

else
   % Get randomized stats >= to observed statistic:
   j1 = find(Fstat1(2:end) >= result(1).F);
   j2 = find(Fstat2(2:end) >= result(2).F);

   % Count values & convert to probability:
   result(1).p = (length(j1)+1)./(iter);
   result(2).p = (length(j2)+1)./(iter);
   result(3).p = NaN;
   result(4).p = NaN;
end

