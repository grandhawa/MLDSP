function result = f_npManova3Nest2(n,m,H,I,G,iter,rep,model,nNestLevels)
% - utility function called by f_npManova
% - 3-way MANOVA: factor 3 nested in factor 2 nested in factor 1
% 
% n = # rows/colums in distance matrix
% m = array of # parameters for each factor
% H = cell array of hat matrix for each factor
% I = I matrix
% G = Gower's centered matrix
% iter = # iterations for permutation test
% rep = data with replication (1) or without (0)
% model = specifies ANOVA design
% nNestLevels = # of levels of nested factors

% -----Notes:-----
% This function is similar to f_npManova3, but is for completely
% nested designs. Replication within each nested factor is
% REQUIRED.

% -----Author:-----
% by David L. Jones, Nov-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Specify sources of variation:
result(1).so = {'factor 1'};
result(2).so = {'factor 2'};
result(3).so = {'factor 3'};
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
		Gvar = G; % use observed G for 1st iteration   
	else
		Gvar = f_shuffle(G,2); % permute square symmetric matrix
	end
	
	for i=1:3 
		% Factors 1 thru 3:
		[df(i),SS(i),MS(i)] = f_npManova1(n,m(i),H{i},I,Gvar,0);
	end;
	
	% adjust degrees of freedom: (Sokal & Rohlf, 1995:293)
	df(2).among = nNestLevels(1) - (df(1).among+1);
	df(3).among = nNestLevels(2) - (nNestLevels(1));
	residual_df = n - nNestLevels(2);
	df(1).total = n-1;
	
	% adjust SS & MS of nested factors: (Sokal & Rohlf, 1995:276-277)
	SS(2).among = SS(2).among - SS(1).among;
	MS(2).among = SS(2).among/df(2).among;
	SS(3).among = SS(3).among - SS(1).among - SS(2).among;
	MS(3).among = SS(3).among/df(3).among;
	
	% adjust error terms:
	residual_SS = SS(1).total - (SS(1).among + SS(2).among + SS(3).among);
	residual_MS = residual_SS/residual_df;
	
	%-----Compute F ratio's:-----
	switch model 
	case 36 % 3 nested in 2 nested in 1:
		Fstat1(k) = MS(1).among/MS(2).among;
		Fstat2(k) = MS(2).among/MS(3).among;
		Fstat3(k) = MS(3).among/residual_MS;
		
	otherwise
		error('Unsupported ANOVA factor specification');
	end
	
		%-----Collect observed values:-----
	if (k==1)
		result(1).df = df(1).among;
		result(2).df = df(2).among;
		result(3).df = df(3).among;
		result(4).df = residual_df;
		result(5).df = df(1).total;
			
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


