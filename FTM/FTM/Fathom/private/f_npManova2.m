function result = f_npManova2(n,m,H,I,G,iter,model)
% - utility function called by f_npManova (2-way MANOVA with interaction)
% 
% n     = # rows/colums in distance matrix
% m     = array of # parameters for each factor
% H     = cell array of hat matrix for each factor
% I     = I matrix
% G     = Gower's centered matrix
% iter  = # iterations for permutation test
% model = specifies factor type (fixed = 1, random = 0)

% -----Author:-----
% by David L. Jones
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Aug-2002: formatted results as a structure
% Nov-2002: added support for fixed and random factors
% May-2003: added abs() of F-ratios
% Mar-2008: removed loops for 'Factors 1 & 2'
% Nov-2013: changed SS(3).among to avoid negative SS interaction (now follows
%           Sokal & Rolf, p.325 Box 11.1 eq. 7) 

% Specify sources of variation:
result(1).so = {'factor 1'};
result(2).so = {'factor 2'};
result(3).so = {'factor 1x2'};
result(4).so = {'residual'};
result(5).so = {'total'};

if (iter<1),
   iter = 1; % do at least once
else   
   Fstat1 = zeros(iter,1); % preallocate results array
   Fstat2 = zeros(iter,1);
   Fstat3 = zeros(iter,1);
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

   
   % Interaction term:
   [ignore,SS(3)] = f_npManova1(n,m(3),H{3},I,Gvar,0);
   % SS(3).among  = SS(3).total - (SS(1).among + SS(2).among + SS(3).resid); % SS interaction
   SS(3).among    = SS(3).among - SS(1).among - SS(2).among; % SS interaction
   % Handle negative SS interaction terms:
   if SS(3).among < 0
      SS(3).among = 0;
   end
   df(3).among    = df(1).among * df(2).among;               % df interaction
   MS(3).among    = SS(3).among/df(3).among;                 % MS interaction
   
   % Global residual df, SS, & MS:
   residual_df = (n-1) - sum([df(1).among df(2).among df(3).among]); 
   residual_SS = SS(1).total - sum([SS(1).among SS(2).among SS(3).among]); 
   residual_MS = residual_SS/residual_df;
   
   % Compute pseudo-F's:
   switch model
	  
   case 21 % Both factors fixed:
	  Fstat1(k) = abs(MS(1).among/residual_MS);
	  Fstat2(k) = abs(MS(2).among/residual_MS); 
	  Fstat3(k) = abs(MS(3).among/residual_MS);
	  
   case 22 % Factor 1 fixed or random, factor 2 random	  
	  Fstat1(k) = abs(MS(1).among/MS(3).among);
	  Fstat2(k) = abs(MS(2).among/MS(3).among); 
	  Fstat3(k) = abs(MS(3).among/residual_MS);
	  
   otherwise
	  error('Unsupported ANOVA factor specification');
   end
      
   % Collect observed values:
   if (k==1)
	  result(1).df = df(1).among;
	  result(2).df = df(2).among;
	  result(3).df = df(3).among;
	  result(4).df = residual_df;
	  result(5).df = (n-1);
	  
	  result(1).SS = SS(1).among;
	  result(2).SS = SS(2).among;
	  result(3).SS = SS(3).among;
	  result(4).SS = residual_SS;
	  result(5).SS = SS(1).total;
	  
	  result(1).MS = MS(1).among;
	  result(2).MS = MS(2).among;
	  result(3).MS = MS(3).among;
	  result(4).MS = residual_MS;
	  result(5).MS = NaN;
	  
	  result(1).F  = Fstat1(1);
	  result(2).F  = Fstat2(1);
	  result(3).F  = Fstat3(1);
	  result(4).F  = NaN;
	  result(5).F  = NaN;
	  
   end
end

if (iter==1)
   result(1).p = NaN;
   result(2).p = NaN;
   result(3).p = NaN;
   result(4).p = NaN;
   result(5).p = NaN;
else
   % get randomized stats >= to observed statistic:
   j1 = find(Fstat1(2:end) >= result(1).F);
   j2 = find(Fstat2(2:end) >= result(2).F);
   j3 = find(Fstat3(2:end) >= result(3).F);
   
   % count values & convert to probability:
   result(1).p = (length(j1)+1)./(iter);
   result(2).p = (length(j2)+1)./(iter);
   result(3).p = (length(j3)+1)./(iter);
   result(4).p = NaN;
   result(5).p = NaN;
end
