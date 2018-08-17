function result = f_npManova3noRep(n,m,H,I,G,iter,model)
% - utility function called by f_npManova
% - 3-way MANOVA without replication (no 2nd order interaction term)
% 
% n = # rows/colums in distance matrix
% m = array of # parameters for each factor
% H = cell array of hat matrix for each factor
% I = I matrix
% G = Gower's centered matrix
% iter = # iterations for permutation test
% model = specifies ANOVA design

% -----Notes:-----
% This function is similar to f_npManova3, but is for data with no
% replication. In this case there is no SS.residual term so the 2nd-order
% interaction term is used instead. As a consequence you cannot test for an
% interaction. Use 'Tukey's test for non-additivity' to test interactions.

% -----Author:-----
% by David L. Jones, 
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Aug-2002: formatted results as a structure
% Nov-2002: added input variable model

% Specify sources of variation:
result(1).so = {'factor 1'};
result(2).so = {'factor 2'};
result(3).so = {'factor 3'};
result(4).so = {'factor 1x2'};
result(5).so = {'factor 1x3'};
result(6).so = {'factor 2x3'};
result(7).so = {'residual'};
result(8).so = {'total'};


if (iter<1),
	iter = 1; % do at least once
else   
	Fstat1 = zeros(iter,1); % preallocate results array
	Fstat2 = zeros(iter,1);
	Fstat3 = zeros(iter,1);
	Fstat4 = zeros(iter,1);
	Fstat5 = zeros(iter,1);
	Fstat6 = zeros(iter,1);
	
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
	
	%-----Interaction terms-----:
	% 1x2:
	[ignore,SS(4)] = f_npManova1(n,m(4),H{4},I,Gvar,0);
	SS(4).among = SS(4).among - SS(1).among - SS(2).among; % SS interaction
	df(4).among = df(1).among * df(2).among;               % df interaction
	MS(4).among = SS(4).among/df(4).among;                 % MS interaction
	
	% 1x3:
	[ignore,SS(5)] = f_npManova1(n,m(5),H{5},I,Gvar,0);
	SS(5).among = SS(5).among - SS(1).among - SS(3).among; % SS interaction
	df(5).among = df(1).among * df(3).among;               % df interaction
	MS(5).among = SS(5).among/df(5).among;                 % MS interaction
	
	% 2x3:
	[ignore,SS(6)] = f_npManova1(n,m(6),H{6},I,Gvar,0);
	SS(6).among = SS(6).among - SS(2).among - SS(3).among; % SS interaction
	df(6).among = df(2).among * df(3).among;               % df interaction
	MS(6).among = SS(6).among/df(6).among;                 % MS interaction
	
	
	% 1x2x3: (used as the error term)
	SS(7).among = sum(diag(H{7}*Gvar*H{7})); % Sum-of-Squares:
	residual_SS = SS(7).among - SS(1).among - SS(2).among - SS(3).among - SS(4).among - SS(5).among - SS(6).among; % SS interaction
	residual_df = df(1).among * df(2).among * df(3).among; % df interaction
	residual_MS = residual_SS/residual_df;                 % MS interaction
	
	%-----Compute F ratio's:-----
	switch model
		
	case 31 % all factors fixed:
		Fstat1(k) = MS(1).among/residual_MS;
		Fstat2(k) = MS(2).among/residual_MS;
		Fstat3(k) = MS(3).among/residual_MS;
		Fstat4(k) = MS(4).among/residual_MS;
		Fstat5(k) = MS(5).among/residual_MS;
		Fstat6(k) = MS(6).among/residual_MS;
		
	case 32 % factors 1 & 2 fixed, 3 random:
		Fstat1(k) = MS(1).among/MS(5).among;
		Fstat2(k) = MS(2).among/MS(6).among;
		Fstat3(k) = MS(3).among/(MS(5).among + MS(6).among - residual_MS);
		Fstat4(k) = MS(4).among/residual_MS;
		Fstat5(k) = MS(5).among/residual_MS;
		Fstat6(k) = MS(6).among/residual_MS;
		
	case 33 % factor 1 fixed or random, 2 & 3 random:
	  Fstat1(k) = MS(1).among/(MS(4).among + MS(5).among - residual_MS);
	  Fstat2(k) = MS(2).among/(MS(4).among + MS(6).among - residual_MS);
	  Fstat3(k) = MS(3).among/(MS(5).among + MS(6).among - residual_MS);
	  Fstat4(k) = MS(4).among/residual_MS;
	  Fstat5(k) = MS(5).among/residual_MS;
	  Fstat6(k) = MS(6).among/residual_MS;
			
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
		result(6).df = df(6).among;
		result(7).df = residual_df;
		result(8).df = df(1).total;
		
		result(1).SS = SS(1).among;
		result(2).SS = SS(2).among;
		result(3).SS = SS(3).among;
		result(4).SS = SS(4).among;
		result(5).SS = SS(5).among;
		result(6).SS = SS(6).among;
		result(7).SS = residual_SS;
		result(8).SS = SS(1).total;
		
		result(1).MS = MS(1).among;
		result(2).MS = MS(2).among;
		result(3).MS = MS(3).among;
		result(4).MS = MS(4).among;
		result(5).MS = MS(5).among;
		result(6).MS = MS(6).among;
		result(7).MS = residual_MS;
		result(8).MS = NaN;
		
		result(1).F  = Fstat1(1);
		result(2).F  = Fstat2(1);
		result(3).F  = Fstat3(1);
		result(4).F  = Fstat4(1);
		result(5).F  = Fstat5(1);
		result(6).F  = Fstat6(1);
		
		result(7).F  = NaN;
		result(8).F  = NaN;
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
	result(8).p = NaN;
else
	% get randomized stats >= to observed statistic:
	j1 = find(Fstat1(2:end) >= result(1).F);
	j2 = find(Fstat2(2:end) >= result(2).F);
	j3 = find(Fstat3(2:end) >= result(3).F);
	j4 = find(Fstat4(2:end) >= result(4).F);
	j5 = find(Fstat5(2:end) >= result(5).F);
	j6 = find(Fstat6(2:end) >= result(6).F);
	
	% count values & convert to probability:
	result(1).p = (length(j1)+1)./(iter);
	result(2).p = (length(j2)+1)./(iter);
	result(3).p = (length(j3)+1)./(iter);
	result(4).p = (length(j4)+1)./(iter);
	result(5).p = (length(j5)+1)./(iter);
	result(6).p = (length(j6)+1)./(iter);
	result(7).p = NaN;
	result(8).p = NaN;
end

