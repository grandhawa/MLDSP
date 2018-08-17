function result = f_npManova3Nest1(n,m,H,I,G,iter,rep,model,nNestLevels)
% - utility function called by f_npManova
% - 3-way MANOVA: Crossed-Nested design
% 
% n = # rows/colums in distance matrix
% m = array of # parameters for each factor
% H = cell array of hat matrix for each factor
% I = I matrix
% G = Gower's centered matrix
% iter = # iterations for permutation test
% rep = data with replication (=1) or without (=0)
% model = specifies ANOVA design
% nNestLevels = # of levels of nested factors

% -----Notes:-----
% This function is similar to f_npManova3, but is for crossed-
% nested designs (A & B are main factors, C is nested within A).
% Replication within the nested factor is REQUIRED. Note that
% there is no 3-way interaction or interaction b/n AxC.

% -----Author:-----
% by David L. Jones, Nov-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Specify sources of variation:
result(1).so = {'factor 1'};
result(2).so = {'factor 2'};
result(3).so = {'factor 3'};
result(4).so = {'factor 1x2'};
result(5).so = {'factor 2x3'};
result(6).so = {'residual'};
result(7).so = {'total'};


if (iter<1),
	iter = 1; % do at least once
else   
	Fstat1 = zeros(iter,1); % preallocate results array
	Fstat2 = zeros(iter,1);
	Fstat3 = zeros(iter,1);
	Fstat4 = zeros(iter,1);
	Fstat5 = zeros(iter,1);
	
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
	
	% Interaction 1x2:
	[ignore,SS(4)] = f_npManova1(n,m(4),H{4},I,Gvar,0);
	SS(4).among = SS(4).among - SS(1).among - SS(2).among; % SS interaction
	df(4).among = df(1).among * df(2).among;               % df interaction
	MS(4).among = SS(4).among/df(4).among;                 % MS interaction
	
	% Adjust degrees of freedom (Sokal & Rohf, 1995:293):
	df(3).among = nNestLevels - (df(1).among+1);
	df(5).among = df(2).among * df(3).among; % df 2x3
	df(1).total = n-1;
		
	% Adjust SS & MS:
	SS(3).among = SS(3).among - SS(1).among;
	MS(3).among = SS(3).among/df(3).among;
	
	% Interaction 2x3:
	[ignore,SS(6)] = f_npManova1(n,m(6),H{6},I,Gvar,0);
	SS(5).among = SS(6).among - sum([SS(1).among SS(2).among SS(3).among SS(4).among]); % SS interaction
	MS(5).among = SS(5).among/df(5).among;
	
	% Residual terms:
	residual_df = df(1).total - sum([df(1).among df(2).among df(3).among df(4).among df(5).among]);
	residual_SS = SS(1).total - sum([SS(1).among SS(2).among SS(3).among SS(4).among SS(5).among]);
	residual_MS = residual_SS/residual_df;
			
	%-----Compute F ratio's:-----
	switch model
		
	case 34 % Cross-Nested design (1 & 2 FIXED)
		Fstat1(k) = MS(1).among/MS(3).among;
		Fstat2(k) = MS(2).among/MS(5).among;
		Fstat3(k) = MS(3).among/MS(5).among;
		Fstat4(k) = MS(4).among/MS(5).among;
		Fstat5(k) = MS(5).among/residual_MS;		
		
	case 35 % Cross-Nested design (1 &/or 2 RANDOM)
		Fstat1(k) = MS(1).among/(MS(3).among + MS(4).among - MS(5).among) ;
		Fstat2(k) = MS(2).among/MS(4).among;
		Fstat3(k) = MS(3).among/MS(5).among;
		Fstat4(k) = MS(4).among/MS(5).among;
		Fstat5(k) = MS(5).among/residual_MS;
	otherwise
		error('Unsupported ANOVA model specification');
		
	end
	
	%-----Collect observed values:-----
	if (k==1)
		result(1).df = df(1).among;
		result(2).df = df(2).among;
		result(3).df = df(3).among;
		result(4).df = df(4).among;
		result(5).df = df(5).among;
		result(6).df = residual_df;
		result(7).df = (n-1);
		
		result(1).SS = SS(1).among;
		result(2).SS = SS(2).among;
		result(3).SS = SS(3).among;
		result(4).SS = SS(4).among;
		result(5).SS = SS(5).among;
		result(6).SS = residual_SS;
		result(7).SS = SS(1).total;
		
		result(1).MS = MS(1).among;
		result(2).MS = MS(2).among;
		result(3).MS = MS(3).among;
		result(4).MS = MS(4).among;
		result(5).MS = MS(5).among;
		result(6).MS = residual_MS;
		result(7).MS = NaN;
		
		result(1).F  = Fstat1(1);
		result(2).F  = Fstat2(1);
		result(3).F  = Fstat3(1);
		result(4).F  = Fstat4(1);
		result(5).F  = Fstat5(1);
		result(6).F  = NaN;
		result(7).F  = NaN;
	end
end

if (iter==1)
	result(1).p = NaN;
	result(2).p = NaN;
	result(3).p = NaN;
	result(4).p = NaN;
	result(5).p = NaN;
	result(6).p = NaN;
	result(7).p = NaN;
else
	% get randomized stats >= to observed statistic:
	j1 = find(Fstat1(2:end) >= result(1).F);
	j2 = find(Fstat2(2:end) >= result(2).F);
	j3 = find(Fstat3(2:end) >= result(3).F);
	j4 = find(Fstat4(2:end) >= result(4).F);
	j5 = find(Fstat5(2:end) >= result(5).F);
	
	% count values & convert to probability:
	result(1).p = (length(j1)+1)./(iter);
	result(2).p = (length(j2)+1)./(iter);
	result(3).p = (length(j3)+1)./(iter);
	result(4).p = (length(j4)+1)./(iter);
	result(5).p = (length(j5)+1)./(iter);
	result(6).p = NaN;
	result(7).p = NaN;
end

