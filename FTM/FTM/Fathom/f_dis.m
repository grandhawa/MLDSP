function dis = f_dis(X,method,adj,minkR,nanFlag)
% - create symmetric dissimilarity (or distance) matrix
%
% USAGE: dis = f_dis(X,'method',adj,minkR);
%
%  X     = input matrix (row = observations, cols = variables)
% method = dissimilarity coefficient to use
% adj    = adjust joint absences by adding dummy species           (default = 0)
% minkR  = value of r for Minkowski's metric                       (default = 1)
%
% dis    = symmetric dissimilarity matrix
%
% ----- Methods: -----
%  'ave': Average distance                                             (=  'd2')
%   'bc': Bray-Curtis dissimilarity
%  'bin': Binomial deviance dissimilarity (scale-invariant)
% 'binU': Binomial deviance dissimilarity (unscaled)
%  'can': Canberra metric                                              (= 'd10')
% 'can2': Canberra dissimilarity
%  'chi': Chi-square metric                                            (= 'd15')
% 'chi2': Chi-square distance                                          (= 'd16')
%  'chJ': Chao's abundance-based Jaccard dissimilarity
% 'chJB': Bias-corrected version of 'chJ'
%  'chS': Chao's abundance-based Sorensen dissimilarity
%  'cho': Orloci's Chord distance                                      (=  'd3')
%  'cod': Coefficient of Divergence + exclude double-zeros             (=  d11')
%  'cor': Pearson Correlation dissimilarity
% 'cor2': Spearman Rank Correlation dissimilarity
%  'cy':  Chao's dissimilarity for count data
% 'czek': Czekanowski distance + exclude double-zeros                  (=  'd8')
%  'euc': Euclidean distance (default)                                 (=  'd1')
% 'euc2': Euclidean distance + exclude double-zeros
%  'esq': Euclidean distance squared
% 'eucs': Euclidean similarity
%  'geo': Geodesic metric                                              (=  'd4')
%  'gow': Gower's metric                                               (= 's15')
% 'gow2': Gower's metric + exclude double-zero                         (= 's19')
% 'gow3': Gower's metric + exclude double-zeros (log2 transform)
% 'gow4': Gower's metric + exclude double-zeros (log10 transform)
%  'hel': Hellinger dissimilarity                                      (= 'd17')
%  'ioa': Whittaker's Index of Association dissimilarity               (=  'd9')
%  'jac': Jaccard dissimilarity                                        (=  's7')
%  'kul': Kulczynski quantitative dissimilarity                        (= 's18')
%  'man': Manhattan (City-block) metric                                (=  'd7')
% 'man2': Manhattan (City-block) metric + exclude double-zeros         (=  'd8')
% 'mink': Minkowski's metric                                           (=  'd6')
%  'mor': Morisita's Index of Overlap for count data
% 'morH': Morisita-Horn dissimilarity
%  'och': Ochiai quantitative dissimilarity
%  'sor': Sorensen's dissimilarity                                     (=  's8')
%  'wat': Watson's nonmetric coefficient dissimilarity                 (= 'd13')
%   's1': Simple Matching dissimilarity
%   's2': Rogers & Tanimoto dissimilarity
%   's3': Sokal & Sneath dissimilarity
%   's4': Sokal & Sneath dissimilarity
%   's5': Sokal & Sneath dissimilarity
%   's6': Sokal & Sneath dissimilarity
%   's9': Jaccard dissimilarity variant
%  's10': Sokal & Sneath dissimilarity
%  's11': Russel & Rao dissimilarity
%  's13': Binary Kulczynski dissimilarity
%  's14': Binary Ochiai dissimilarity
%  's26': Faith dissimilarity
%  'spe': Species Profiles distance
%
% SEE ALSO: f_latlong

% -----Notes:-----
% The Kulczynski dissimilarity coefficient is well suited for raw species
% abundance data and was determined by Faith et al. (1987) to be one of the more
% robust measures of ecological distance. Legendre & Legendre (1998) suggest the
% similarity version of this coefficient is usually monotonic with the Steinhaus
% similiarity metric (= Bray-Curtis dissimilarity).
%
% For the Canberra metric/dissimilarity, a given difference between rare species
% contributes more to the distance than the same distance between more abundant
% species, which contrasts the behavior of the Bray-Curtis where abundant and
% rare species contribute equally (Legendre & Legendre, 1998:p.283,287).

% -----Testing:-----
% The following coefficients give the same results as PERMDISP for the 'tikus'
% data:
% BC, CAN, CHI2, CHO, GOW, GOW2, HEL, JAC, KUL, MAN
%
% The following coefficients give the same results as PRIMER 6 for the 'ekofisk'
% data:
% GOW3, GOW4, BIN, BINU, CY
%
% The following coefficients give the same results as the "vegan package
% for R" for the "varespec" data set:
% - original data: MORH
% - data were converted to integers using ceil: MOR, CHJB, CY
%
% The following coefficients give the same results as the "EstimateS"
% computer program for the "seedbank" data set:
% BC, MORH, JAC, SOR, CHJ, CHS, CHJB, CHSB

% -----References:-----
% Anderson, M. J., K. E. Ellingsen, and B. H. McArdle. 2006. Multivariate
%   dispersion as a measure of beta diversity. Ecology Letters 9: 683-693.
% Anderson, M. J., and R. B. Millar. 2004. Spatial variation and effects of
%   habitat on temperate reef fish assemblages in northeastern New Zealand.
%   J. Exp. Mar. Biol. Ecol. 305: 191?221.
% Anderson, M. J. and A. A. Thompson. 2004. Multivariate control charts for
%   ecological and environmental monitoring. Ecological Applications 14:
%   1921-1935.
% Chao, A., R. L. Chazdon, R. K. Colwell, and T.-J. Shen. 2006.
%   Abundance-based similarity indices and their estimation when there are
%   unseen species in samples. Biometrics 62: 361-371.
% Clarke, K. R., P. J. Somerfield, and M. G. Chapman. 2006. On resemblance
%   measures for ecological studies, including taxonomic dissimilarities and a
%   zero-adjusted Bray-Curtis coefficient for denuded assemblages. J. Exp. Mar.
%   Biol. Ecol. 330: 55-80.
% Faith, D. P., P. R. Minchin, and L. Belbin. 1987. Compositional dissimilarity
%   as a robust measure of ecological distance. Vegetatio 69: 57-68.
% Elmore, K. L., and M. B. Richman. 2001. Euclidean distance as a
%   similarity metric for principal component analysis. Monthly Weather
%   Review 129: 540-549.
% Legendre, P. and E. D. Gallagher. 2001. Ecologically meaningful
%   transformations for ordination of species data. Oecologia 129: 271-280.
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam.
%
% The indices for Morisita and Horn-Morisita were ported to Matlab
% following the formulas after those in vegdist.c from the "vegan package
% for R" and:
% "http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/vegdist.html" and the

% -----Author:-----
% by David L. Jones, Feb-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Feb-2011: added 'eucs'
% Mar-2011: edited documentation
% Jul-2012: included final error check for NaN's
% Aug-2012: edited documentation
% Nov-2012: updated "remove cols of zero's"; added support for Morisita,
%           Horn-Morisita, Binomial Deviance, and Chao's abundance-based Jaccard.
% Dec-2012: added chS, chJB, chSB
% Apr-2013: updated documentation
% Jan-2014: added internal nanFlag

% -----Check input & set defaults:-----
if (nargin < 2), method  = 'euc'; end % default Euclidean distance
if (nargin < 3), adj     =  0;    end % default no adjustment via dummy species
if (nargin < 4), minkR   =  1;    end % default Minkowski's r = 1
if (nargin < 5), nanFlag =  1;    end % internal flag for final check for NaN's

% if (sum(find(isnan(X(:)==1))>1)), error('Input matrix contains NaN''s!'); end
n = size(X,1); % get # rows, # variables
% if (n<2), error('Input matrix requires > 1 row!'); end

% Convert to lower case:
method = lower(method);

% Remove cols with only 0's (these species are NEVER present):
switch method
   case {'euc','euc2','esq','eucs'} % skip for these methods
   otherwise
      idx = find(sum(X==0) == n);
      X(:,idx) = [];
end

% For some indices, check input is only count data:
switch method
   case {'mor','chjb','chsb','cy'}
      if ~(isequal(X,floor(X))) || sum(sum(X<0)>0)
         error([method ' requires count data (positive whole numbers)!'])
      end
end
% -------------------------------------

% Adjust for joint-absences by adding a dummy species (Clarke et al., 2006):
if (adj>0)
   X = [ones(n,1) X];
end

D   = repmat(NaN,n-1,1); % preallocate
idx = [0 0];             % initialize

switch method
   
   % Average Distance (eq. 7.36 in Legendre & Legendre, 1998):
   case {'ave','d2'}
      p = size(X,2); % get # of variables
      for i = 1:n-1 % repeat for each row except the last
         blk    = numel(i+1:n);              % get size of this block of distances
         r      = repmat(X(i,:),blk,1);      % extract this row, replicate
         R      = X(i+1:n,:);                % extract remaining rows
         idx    = idx(end)+1:idx(end)+blk;   % get index to this block of distances
         D(idx) = sqrt( ( sum( (r-R).^2 ,2) )/p ); % sum by row prodces a column
      end
      
      
      % Bray-Curtis Dissimilarity (eq. 7.57 in Legendre & Legendre, 1998):
   case {'bc'}
      for i = 1:n-1 % repeat for each row except the last
         blk    = numel(i+1:n);                 % get size of this block of distances
         r      = repmat(X(i,:),blk,1);         % extract this row, replicate
         R      = X(i+1:n,:);                   % extract remaining rows
         idx    = idx(end)+1:idx(end)+blk;      % get index to this block of distances
         D(idx) = sum(abs(r-R),2)./ sum(r+R,2); % sum by row produces a column
      end
      
      % Binomial Deviance Dissimilarity (p.199 in Anderson & Millar, 2004):
      % (after from 'veg_millar' in 'vegdist.c' from 'vegan for R')
   case {'bin','binu'}
      for i = 1:n-1 % repeat for each row except the last
         blk    = numel(i+1:n);                 % get size of this block of distances
         r      = repmat(X(i,:),blk,1);         % extract this row, replicate
         R      = X(i+1:n,:);                   % extract remaining rows
         idx    = idx(end)+1:idx(end)+blk;      % get index to this block of distances
         
         % Exclude double-zeros:
         binZ    = ( (r>0)+(R>0) )==0; % binary matrix with 1 indicating double-zeros
         idxZ    = find(binZ==1);      % indices to elements with double-zeros
         r(idxZ) = NaN;                % replace double-zeros with NaN's
         R(idxZ) = NaN;
         
         % Set up coefficients:
         nk        = (r+R);            % combined total:
         c_r       = (log(r)-log(nk)); % coefficient for r
         c_r(r<=0) = 0;                % replace with 0 when r = 0
         c_R       = (log(R)-log(nk)); % coefficient for R
         c_R(R<=0) = 0;                % replace with 0 when R = 0
         
         if isequal(method,'bin')
            % Scale-invariant:
            D(idx) = nansum( ( (r.*c_r) + (R.*c_R) + (nk*abs(log(0.5))) ) ./ nk ,2); % Nan's ignored in calculations
         else
            % Unscaled:
            D(idx) = nansum( ( (r.*c_r) + (R.*c_R) + (nk*abs(log(0.5))) ) ,2); % Nan's ignored in calculations
         end
      end
      
      % Canberra Meric/Dissimilarity (eq. 7.49/7.50 in Legendre & Legendre, 1998):
   case {'can','d10','can2'}
      for i = 1:n-1 % repeat for each row except the last
         blk   = numel(i+1:n);            % get size of this block of distances
         r     = repmat(X(i,:),blk,1);    % extract this row, replicate
         R     = X(i+1:n,:);              % extract remaining rows
         idx   = idx(end)+1:idx(end)+blk; % get index to this block of distances
         
         % Exlude double-zeros:
         binZ    = ( (r>0)+(R>0) )==0; % binary matrix with 1 indicating double-zeros
         idxZ    = find(binZ==1);      % indices to elements with double-zeros
         r(idxZ) = NaN;                % replace double-zeros with NaN's
         R(idxZ) = NaN;
         
         D(idx) = nansum( abs(r-R)./(r+R) ,2); % NaN's ignored in calculations
         
         % Canberra Dissimilarity (scaled from 0-1):
         if isequal(method,'can2')
            nz     = sum( (binZ==0),2); % get # non-zero variables, sum by row produces a column
            D(idx) = D(idx)./nz;        % scale from 0-1
         end
      end
      
      % Chao's abundance-based indices:
   case {'chj','chs'}
      % Note presence = 1, absence = 0:
      % site X = 1 1 1 0 0 1 0;
      % site Y = 0 0 1 1 1 0 0;
      %      I = 1 1 0 1 1 1 0
      %
      %      S = 0 0 1 0 0 0 0
      %
      % 1) the sum of X.*Y gives the number of species present in both sites
      % 2) abs(X-Y) provides a vector (I) where 1 = only 1 of the 2
      %    species were present and 0 = both were present
      % 3) the sum of X.*I gives the number of species only present in site X
      % 4) the sum of Y.*I gives the number of species only present in site Y
      % 5) the number of 0's from the sum of X+Y gives the number of species
      %    absent from both sites (see Sorenson coefficient)
      
      % Transform data to proportional abundances:
      nc          = size(X,2);        % get # columns
      rowSum      = sum(X,2);         % row sums
      X           = X ./ repmat(rowSum,1,nc);
      X(isnan(X)) = 0;                % handle row sums of zero
      % ------------------------------------------
      
      for i = 1:n-1 % repeat for each row except the last
         B   = double(X>0);              % convert to binary (presence/absence)
         blk = numel(i+1:n);             % get size of this block of distances
         
         % Binary values:
         Br = repmat(B(i,:),blk,1);      % extract this row, replicate
         BR = B(i+1:n,:);                % extract remaining rows
         S  = Br.*BR;                    % Boolean indicating both species present
         
         % Abundance values:
         r = repmat(X(i,:),blk,1);       % extract this row, replicate
         R = X(i+1:n,:);                 % extract remaining rows
         
         % Chao's disimilarity:
         idx  = idx(end)+1:idx(end)+blk; % get index to this block of distances
         U    = sum(r.*S,2);             % both species present
         V    = sum(R.*S,2);             % both species present
         % chao = ( (U.*V) ./ (U + V - U.*V) ); % eq. 5
         
         % page 363 in Chao et al., 2006:
         a   = U.*V;
         b   = U .* (1-V);
         c   = V .* (1-U);
         
         switch method
            case 'chj' % Chao's abundanced-based Jaccard
               sim = ( a ./ (a + b + c) );
            case 'chs' % Chao's abundanced-based Sorsenson
               sim = ( 2*a ./ (2*a + b + c) );
         end
         
         % Handle NaN's resulting from 'zero divided by zero':
         sim(isnan(sim)) = 0;
         
         % Convert similarity to dissimiarity:
         D(idx) = 1 - sim;
      end
      
      % Bias-corrected versions of Chao's abundance-based indices:
   case {'chjb', 'chsb'}
      nc = size(X,2); % get # columns
      for i = 1:n-1 % repeat for each row except the last
         B   = double(X>0);             % convert to binary (presence/absence)
         blk = numel(i+1:n);            % get size of this block of distances
         idx = idx(end)+1:idx(end)+blk; % get index to this block of distances
         
         % Abundance values:
         r = repmat(X(i,:),blk,1);      % extract this row, replicate
         R = X(i+1:n,:);                % extract remaining rows
         
         % Binary values:
         Br   = repmat(B(i,:),blk,1);   % extract this row, replicate
         BR   = B(i+1:n,:);             % extract remaining rows
         S    = Br.*BR;                 % Boolean: both species present
         
         Bfp1 = Br.*(R==1);             % Boolean: X is present, Y = 1
         Bfp2 = Br.*(R==2);             % Boolean: X is present, Y = 2
         Bf1p = BR.*(r==1);             % Boolean: Y is present, X = 1,
         Bf2p = BR.*(r==2);             % Boolean: Y is present, X = 2
         
         fp1 = sum(Bfp1,2);             % # shared species when Y = 1
         fp2 = sum(Bfp2,2);             % # shared species when Y = 2
         f1p = sum(Bf1p,2);             % # shared species when X = 1
         f2p = sum(Bf2p,2);             % # shared species when X = 2
         
         % Prevent divide-by-zero errors (p. 365 in Chao et al., 2006):
         fp2(fp2==0) = 1;
         f2p(f2p==0) = 1;
         
         % Get sample totals (row sums):
         N = sum(r,2);
         M = sum(R,2);
         
         % Convert data to proportional abundances:
         r_p = r ./ repmat(N,1,nc);
         R_p = R ./ repmat(M,1,nc);
         
         % Handle row sums of zero:
         r_p(isnan(r_p)) = 0;
         R_p(isnan(R_p)) = 0;
         
         % Bias correction (eq. 1 & 2 in Chao et al., 2006):
         U_1 = sum((r_p).*S,2); % sum of X when both are present
         V_1 = sum((R_p).*S,2); % sum of Y when both are present
         
         U_2 = ((M-1) ./ M) .* (fp1 ./ (2*fp2));
         V_2 = ((N-1) ./ N) .* (f1p ./ (2*f2p));
         
         U_3 = sum((r_p) .* Bfp1,2); % sum of X when X is present, Y = 1
         V_3 = sum((R_p) .* Bf1p,2); % sum of Y when Y is present, X = 1
         
         U_hat = U_1 + (U_2 .* U_3);
         V_hat = V_1 + (V_2 .* V_3);
         
         % Adjust data for highly overlapped communities:
         U_hat(U_hat>1) = 1;
         V_hat(V_hat>1) = 1;
         
         % page 363 in Chao et al., 2006:
         a = U_hat .* V_hat;
         b = U_hat .* (1-V_hat);
         c = V_hat .* (1-U_hat);
         
         switch method
            case 'chjb' % Chao's abundanced-based Jaccard
               sim = ( a ./ (a + b + c) );
            case 'chsb' % Chao's abundanced-based Sorsenson
               sim = ( 2*a ./ (2*a + b + c) );
         end
         
         % Handle NaN's resulting from 'zero divided by zero':
         sim(isnan(sim)) = 0;
         
         % Convert similarity to dissimiarity:
         D(idx) = 1 - sim;
      end
      
      
      % Chi-square Metric/Distance (eq. 7 & 9 in Legendre & Gallagher, 1998):
   case {'chi','d15','chi2','d16'}
      nc        = size(X,2);        % get # columns
      rowSum    = sum(X,2);
      colSum    = sum(X);
      
      % Chi-square-transformed data:
      CS = X ./ ( repmat(rowSum,1,nc) .* repmat(colSum.^0.5,n,1) );
      
      % Chi-square Distance:
      if isequal(method,'chi2') || isequal(method,'d16')
         totSum = sum(X(:))^0.5;
         CS      = CS * totSum; % eq. 9
      end
      
      
      % Chord Distance (eq. 3 in Legendre & Gallagher, 1998):
   case {'cho','d3','geo','d4'}
      nc        = size(X,2);        % get # columns
      CD        = zeros(n,nc);      % preallocate
      rowSum    = (sum(X.^2,2)).^0.5;
      idx       = find(rowSum > 0); % skip these to prevent divide-by-zero error
      CD(idx,:) = X(idx,:) ./ repmat(rowSum(idx),1,nc);
      
      
      % Coefficient of Divergence (eq. 7.51 in Legendre & Legendre, 1998):
   case {'cod','d11'}
      for i = 1:n-1 % repeat for each row except the last
         blk   = numel(i+1:n);            % get size of this block of distances
         r     = repmat(X(i,:),blk,1);    % extract this row, replicate
         R     = X(i+1:n,:);              % extract remaining rows
         idx   = idx(end)+1:idx(end)+blk; % get index to this block of distances
         
         % Exlude double-zeros:
         binZ    = ( (r>0)+(R>0) )==0; % binary matrix with 1 indicating double-zeros
         nz      = sum( (binZ==0),2);  % get # non-zero variables, sum by row produces a column
         idxZ    = find(binZ==1);      % indices to elements with double-zeros
         r(idxZ) = NaN;                % replace double-zeros with NaN's
         R(idxZ) = NaN;
         
         D(idx) = sqrt((nansum( ((r-R)./(r+R)).^2 , 2))./nz); % NaN's ignored in calculations
      end
      
   % Chao's CY dissimilarity (Anderson & Thompson, 2004; vegdist.c from vegan):   
   case {'cy'}
      for i = 1:n-1 % repeat for each row except the last
         blk = numel(i+1:n);                 % get size of this block of distances
         r   = repmat(X(i,:),blk,1);         % extract this row, replicate
         R   = X(i+1:n,:);                   % extract remaining rows
         idx = idx(end)+1:idx(end)+blk;      % get index to this block of distances
         
         % Exclude double-zeros:
         binZ    = ( (r>0)+(R>0) )==0; % binary matrix with 1 indicating double-zeros
         idxZ    = find(binZ==1);      % indices to elements with double-zeros
         r(idxZ) = NaN;                % replace double-zeros with NaN's
         R(idxZ) = NaN;
         
         % Get # species involved in each pair of sites;
         BrR = (r>0)+(R>0);     % add presence/absence
         S   = sum((BrR>0) ,2); % sum binary version
         
         % Replace remaining zeros with 0.1:
         r(r==0) = 0.1; % <- this value should be user defined (Anderson & Thompson, 2004)
         R(R==0) = 0.1;
         
         % Set up terms:
         t1 = (r+R).^(-1);
         t2 = (r+R) * log(0.5);
         t3 = r .* log( R ./ (r+R) );
         t4 = R .* log( r ./ (r+R) );
         
         % CY dissimilarity:
         disCY = (nansum( t1 .* (t2 - t3 - t4), 2)) ./ S; % Nan's ignored in calculations
         
         % Replace negative values with 0:
         disCY(disCY<0) = 0;
         
         % Return dissimilarities:
         D(idx) = disCY;
      end
      
      
      % Pearson Correlation Dissimilarity (Clarke et al, 2006 p.75):
   case {'cor'}
      % Note: the upper vs. lower tri-diagonals returned by 'corrcoef' are not
      % always symmetrical to 16 digits, hence the unwrap is required
      D = ( 1 - f_unwrap(corrcoef(X'),0) )/2; % convert to distance, range 0-1
      
      
      % Spearman Rank Correlation Dissimilarity (Clarke et al, 2006 p.75)
   case {'cor2'}
      D = ( 1 - f_unwrap(corrcoef(f_ranks(X')),0) )/2; % convert to distance, range 0-1
      
      
      % Euclidean Distance (eq. 7.34/7.35 in Legendre & Legendre, 1998):
   case {'euc','esq','eucs'}
      for i = 1:n-1 % repeat for each row except the last
         blk    = numel(i+1:n);              % get size of this block of distances
         r      = repmat(X(i,:),blk,1);      % extract this row, replicate
         R      = X(i+1:n,:);                % extract remaining rows
         idx    = idx(end)+1:idx(end)+blk;   % get index to this block of distances
         
         if isequal(method,'euc')
            D(idx) = sqrt( sum( (r-R).^2 ,2) ); % sum by row prodces a column
         else % 'esq'
            D(idx) = sum( (r-R).^2 ,2);         % sum by row prodces a column
         end
      end
      
      % Euclidean Distance + exclude double-zeros:
   case {'euc2'}
      for i = 1:n-1 % repeat for each row except the last
         blk    = numel(i+1:n);              % get size of this block of distances
         r      = repmat(X(i,:),blk,1);      % extract this row, replicate
         R      = X(i+1:n,:);                % extract remaining rows
         idx    = idx(end)+1:idx(end)+blk;   % get index to this block of distances
         
         % Exlude double-zeros (same effect as eq. 8 in Anderson et al., 2006):
         binZ    = ( (r>0)+(R>0) )==0; % binary matrix with 1 indicating double-zeros
         idxZ    = find(binZ==1);      % indices to elements with double-zeros
         r(idxZ) = NaN;                % replace double-zeros with NaN's
         R(idxZ) = NaN;
         nz      = sum( (binZ==0),2);  % get # non-zero variables, sum by row produces a column
         
         D(idx) = (sqrt( nansum( (r-R).^2 ,2) ))./nz; % NaN's ignored, sum by row produces a column
      end
      
      
      % Gower's Coefficient (eq. 12 in Clarke et al., 2006):
   case {'gow','s15'}
      S = max(X) - min(X); % get range (span) of each descriptor
      p = size(X,2);       % get # of variables
      for i = 1:n-1 % repeat for each row except the last
         blk    = numel(i+1:n);            % get size of this block of distances
         r      = repmat(X(i,:),blk,1);    % extract this row, replicate
         R      = X(i+1:n,:);              % extract remaining rows
         idx    = idx(end)+1:idx(end)+blk; % get index to this block of distances
         D(idx) = (sum( abs(r-R) ./ repmat(S,blk,1) ,2))/p; % sum by row prodces a column
      end
      
      
      % Gower's Coefficient + exclude double-zeros:
   case {'gow2','s19'}
      S = max(X) - min(X);                 % get range (span) of each descriptor
      for i = 1:n-1 % repeat for each row except the last
         blk    = numel(i+1:n);            % get size of this block of distances
         r      = repmat(X(i,:),blk,1);    % extract this row, replicate
         R      = X(i+1:n,:);              % extract remaining rows
         idx    = idx(end)+1:idx(end)+blk; % get index to this block of distances
         
         % Exlude double-zeros (same effect as eq. 4 in Anderson et al., 2006):
         binZ    = ( (r>0)+(R>0) )==0; % binary matrix with 1 indicating double-zeros
         idxZ    = find(binZ==1);      % indices to elements with double-zeros
         r(idxZ) = NaN;                % replace double-zeros with NaN's
         R(idxZ) = NaN;
         nz      = sum( (binZ==0),2); % get # non-zero variables, sum by row produces a column
         
         D(idx) = (nansum( abs(r-R) ./ repmat(S,blk,1) ,2))./nz; % sum by row prodces a column
      end
      
      
      % Gower's Coefficient Modified (eq. 6 in Anderson et al., 2006):
   case {'gow3','gow4'}
      % Transformm data:
      if isequal(method,'gow3')
         X = f_normal(X,'log2_alt');
      else
         X = f_normal(X,'log_alt');
      end
      
      for i = 1:n-1 % repeat for each row except the last
         blk    = numel(i+1:n);            % get size of this block of distances
         r      = repmat(X(i,:),blk,1);    % extract this row, replicate
         R      = X(i+1:n,:);              % extract remaining rows
         idx    = idx(end)+1:idx(end)+blk; % get index to this block of distances
         
         % Exlude double-zeros (same effect as eq. 4 in Anderson et al., 2006):
         binZ    = ( (r>0)+(R>0) )==0; % binary matrix with 1 indicating double-zeros
         idxZ    = find(binZ==1);      % indices to elements with double-zeros
         r(idxZ) = NaN;                % replace double-zeros with NaN's
         R(idxZ) = NaN;
         nz      = sum( (binZ==0),2); % get # non-zero variables, sum by row produces a column
         
         D(idx) = (nansum( abs(r-R) ,2))./nz; % sum by row prodces a column
      end
      
      
      % Hellinger Dissimilarity (eq. 13 in Legendre & Gallagher, 1998):
   case {'hel','d17'}
      nc       = size(X,2);        % get # columns
      H        = zeros(n,nc);      % preallocate
      rowSum   = sum(X,2);
      idx      = find(rowSum > 0); % skip these to prevent divide-by-zero error
      H(idx,:) = sqrt(X(idx,:) ./ repmat(rowSum(idx),1,nc)); % Hellinger-transformed
      
      
      % Whittaker's Index of Association (eq. 7.45 in Legendre & Legendre, 1998):
   case {'ioa','d9'}
      p = size(X,2);
      for i = 1:n-1 % repeat for each row except the last
         blk    = numel(i+1:n);            % get size of this block of distances
         r      = repmat(X(i,:),blk,1);    % extract this row, replicate
         R      = X(i+1:n,:);              % extract remaining rows
         idx    = idx(end)+1:idx(end)+blk; % get index to this block of distances
         D(idx) = ( sum( abs( r./repmat(sum(r,2),1,p) - R./repmat(sum(R,2),1,p) ), 2) )/2; % sum by row prodces a column
      end
      
      
      % Kulczynski Dissimilarity (eq. 7.25 in Legendre & Legendre, 1998):
   case {'kul','s18'}
      for i = 1:n-1 % repeat for each row except the last
         blk    = numel(i+1:n);            % get size of this block of distances
         r      = repmat(X(i,:),blk,1);    % extract this row, replicate
         R      = X(i+1:n,:);              % extract remaining rows
         idx    = idx(end)+1:idx(end)+blk; % get index to this block of distances
         D(idx) = 1 - (( (sum((min(r,R)),2) ./ sum(r,2)) + (sum((min(r,R)),2) ./ sum(R,2)) )/2);
      end
      
      
      % Manhattan Distance (eq. 7.45 in Legendre & Legendre, 1998):
   case {'man','d7'}
      for i = 1:n-1 % repeat for each row except the last
         blk    = numel(i+1:n);            % get size of this block of distances
         r      = repmat(X(i,:),blk,1);    % extract this row, replicate
         R      = X(i+1:n,:);              % extract remaining rows
         idx    = idx(end)+1:idx(end)+blk; % get index to this block of distances
         D(idx) = sum(abs(r-R),2);         % sum by row prodces a column
      end
      
      
      % Manhattan Distance + exclude double zeros:
      % (= Czekanowski Distance + exclude double zeros)
   case {'man2','czek','d8'}
      for i = 1:n-1 % repeat for each row except the last
         blk    = numel(i+1:n);            % get size of this block of distances
         r      = repmat(X(i,:),blk,1);    % extract this row, replicate
         R      = X(i+1:n,:);              % extract remaining rows
         idx    = idx(end)+1:idx(end)+blk; % get index to this block of distances
         
         % Exlude double-zeros (same effect as eq. 7 in Anderson et al., 2006):
         binZ    = ( (r>0)+(R>0) )==0; % binary matrix with 1 indicating double-zeros
         idxZ    = find(binZ==1);      % indices to elements with double-zeros
         r(idxZ) = NaN;                % replace double-zeros with NaN's
         R(idxZ) = NaN;
         nz      = sum( (binZ==0),2);  % get # non-zero variables, sum by row produces a column
         
         D(idx) = (nansum(abs(r-R),2))./nz; % ignore NaN's in calculations, sum by row prodces a column
      end
      
      
      % Minkowskis Metric (eq. 7.44 in Legendre & Legendre, 1998):
   case {'mink','d6'}
      for i = 1:n-1 % repeat for each row except the last
         blk    = numel(i+1:n);            % get size of this block of distances
         r      = repmat(X(i,:),blk,1);    % extract this row, replicate
         R      = X(i+1:n,:);              % extract remaining rows
         idx    = idx(end)+1:idx(end)+blk; % get index to this block of distances
         D(idx) = (sum( (abs(r-R)).^minkR ,2)).^(1/minkR); % sum by row prodces a column
      end
      
      % Morisita's Index of Overlap (after vegan's vegdist.c):
   case {'mor'}
      for i = 1:n-1 % repeat for each row except the last
         blk      = numel(i+1:n);                 % get size of this block of distances
         r        = repmat(X(i,:),blk,1);         % extract this row, replicate
         R        = X(i+1:n,:);                   % extract remaining rows
         idx      = idx(end)+1:idx(end)+blk;      % get index to this block of distances
         lambda_r = sum(r.*(r-1),2) ./ sum(r,2) ./ (sum(r,2)-1);
         lambda_R = sum(R.*(R-1),2) ./ sum(R,2) ./ (sum(R,2)-1);
         sim      = 2*sum(r.*R,2) ./ ( lambda_r + lambda_R ) ./ sum(r,2) ./ sum(R,2) ; % sum by row produces a column
         
         % Handle NaN's that arise from 'divide zero by zero':
         sim(isnan(sim)) = 0;
         
         % Convert similarity to distance:
         D(idx)   = 1 - sim;
      end
      
      % Morisita-Horn Dissimilarity (after vegan for R documentation):
      % (see also eq. 5 in Chao et al., 2006)
   case {'morh'}
      for i = 1:n-1 % repeat for each row except the last
         blk      = numel(i+1:n);                 % get size of this block of distances
         r        = repmat(X(i,:),blk,1);         % extract this row, replicate
         R        = X(i+1:n,:);                   % extract remaining rows
         idx      = idx(end)+1:idx(end)+blk;      % get index to this block of distances
         lambda_r = sum(r.^2,2) ./ sum(r,2).^2;
         lambda_R = sum(R.^2,2) ./ sum(R,2).^2;
         sim      = 2*sum((r.*R),2) ./ ( (lambda_r + lambda_R) .* (sum(r,2) .* sum(R,2)) ); % sum by row produces a column
         
         % Handle NaN's that arise from 'divide zero by zero':
         sim(isnan(sim)) = 0;
         
         % Convert similarity to distance:
         D(idx)   = 1 - sim;
      end
      
      % Ochiai Quantitative Dissimilarity (eq. 17 in Clarke et al., 2006)
   case {'och'}
      % BC     = W / [(A+B)/2]         -> denom is arithmetic mean
      % KUL    = W / 1/[(1/A + 1/B)/2] -> denom is harmonic mean
      % OCHIAI = W / sqrt(A*B)         -> denom is geomeric mean
      for i = 1:n-1 % repeat for each row except the last
         blk    = numel(i+1:n);            % get size of this block of distances
         r      = repmat(X(i,:),blk,1);    % extract this row, replicate
         R      = X(i+1:n,:);              % extract remaining rows
         idx    = idx(end)+1:idx(end)+blk; % get index to this block of distances
         D(idx) =  1 - ( sum((min(r,R)),2) ./ sqrt(sum(r,2) .* sum(R,2)) );
      end
      
      
      % Presence/Absence Metrics:
   case {'s1','s2','s3','s4','s5','s6','s7','jac','s8','sor','s9','s10','s11',...
         's13','s14','s26','wat','d13'}
      for i = 1:n-1 % repeat for each row except the last
         % Note presence = 1, absence = 0:
         % site X = 1 1 1 0 0 1 0;
         % site Y = 0 0 1 1 1 0 0;
         %      I = 1 1 0 1 1 1 0
         %
         % 1) the sum of X.*Y gives the number of species present in both sites
         % 2) abs(X-Y) provides a vector (I) where 1 = only 1 of the 2
         %    species were present and 0 = both were present
         % 3) the sum of X.*I gives the number of species only present in site X
         % 4) the sum of Y.*I gives the number of species only present in site Y
         % 5) the number of 0's from the sum of X+Y gives the number of species
         %    absent from both sites (see Sorenson coefficient)
         X   = double(X>0);             % convert to binary (presence/absence)
         blk = numel(i+1:n);            % get size of this block of distances
         r   = repmat(X(i,:),blk,1);    % extract this row, replicate
         R   = X(i+1:n,:);              % extract remaining rows
         idx = idx(end)+1:idx(end)+blk; % get index to this block of distances
         a   = sum(r.*R,2);             % both species present
         b   = sum(abs(r-R).*r,2);      % only first species present
         c   = sum(abs(r-R).*R,2);      % only second species present
         d   = sum(double((r+R)==0),2); % both species absent
         
         switch method
            case {'s1'} % Simple Matching Dissimilarity:   (eq.7.1 in L&L, 1998)
               D(idx) = 1 - ( (a + d) ./ (a + b + c + d) );
            case {'s2'} % Rogers & Tanimoto Dissimilarity:  (eq.7.2 in L&L, 1998)
               D(idx) = 1 - ( (a + d) ./ (a + 2*b + 2*c + d) );
            case {'s3'} % Sokal & Sneath, 1963:            (eq.7.3 in L&L, 1998)
               D(idx) = 1 - ( (2*a + 2*d) ./ (2*a + b + c + 2*d) );
            case {'s4'} % Sokal & Sneath, 1963:            (eq.7.4 in L&L, 1998)
               D(idx) = (a + d) ./ (b + c);
            case {'s5'} % Sokal & Sneath, 1963:            (eq.7.5 in L&L, 1998)
               D(idx) = 1 - ( ( (a./(a+b) + a./(a+c) + d./(b+d) + d./(c+d)) )/4 );
            case {'s6'} % Sokal & Sneath, 1963:            (eq.7.6 in L&L, 1998)
               D(idx) = 1 - ( (a./sqrt((a+b).*(a+c))) .* (d./sqrt((b+d).*(c+d))) );
            case {'s7','jac'} % Jaccard Dissimilarity: (eq.7.10/7.58 in L&L, 1998)
               D(idx) = 1 - ( a ./ (a + b + c) );
            case {'s8','sor'} % Sorensen Dissimilarity:(eq.7.11/7.56 in L&L, 1998)
               D(idx) = 1 - ( 2*a ./ (2*a + b + c) );
            case {'s9'} % Jaccard variant:                (eq.7.12 in L&L, 1998)
               D(idx) = 1 - ( 3*a ./ (3*a + b + c) );
            case {'s10'} % Sokal & Sneath, 1963:          (eq.7.13 in L&L, 1998)
               D(idx) = 1 - ( a ./ (a + 2*b + 2*c) );
            case {'s11'} % Russell & Rao, 1940:           (eq.7.14 in L&L, 1998)
               D(idx) = 1 - ( a ./ (a + b + c + d) );
            case {'s13'} % Binary Kulczynski:             (eq.7.16 in L&L, 1998)
               D(idx) = 1 - ( ( a./(a+b) + a./(a+c) )/2 );
            case {'s14'} % Ochiai Dissimilarity:          (eq.7.17 in L&L, 1998)
               D(idx) = 1 - ( a ./ sqrt((a + b).*(a + c)) );
            case {'s26'} % Faith Dissimilarity:           (eq.7.18 in L&L, 1998)
               D(idx) = 1 - ( (a + d/2) ./ (a + b + c + d) );
            case {'wat','d13'} % Watson's Nonmetric Coefficient: (eq.7.56 in L&L, 1998)
               D(idx) = (b + c) ./ (2*a + b + c);
         end
      end
      % Convert similarity to dissimilarity:
      if isequal(method,'s4')
         maxVar = max(D);   % range is 0-inf
         D = maxVar(1) - D; % convert similarity to distance
      end
      
      
      % Species Profiles Distance (eq. 11 in Legendre & Gallagher, 1998):
   case {'spe'}
      nc        = size(X,2);        % get # columns
      SP        = zeros(n,nc);      % preallocate
      rowSum    = sum(X,2);
      idx       = find(rowSum > 0); % skip these to prevent divide-by-zero error
      SP(idx,:) = X(idx,:) ./ repmat(rowSum(idx),1,nc); % Species Profiles-transformed
      
   otherwise
      error('Unknown dissimilarity coefficient!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Symmetric distance matrix:        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch method
   case {'chi','d15','chi2','d16'}
      dis = f_dis(CS,'euc'); % Euclidean distance of Chi-square-transformed data
      
   case {'cho','d3','geo','d4'}
      dis = f_dis(CD,'euc'); % Euclidean distance of Chord-transformed data
      % Geodesic metric: (eq. 7.39 in L&L, 1998)
      if isequal(method,'geo') || isequal(method,'d4')
         dis = f_rewrap(acos(1 - ((f_unwrap(dis).^2)/2))); % convert cord to geodesic
      end
      
   case {'hel','d17'}
      dis = f_dis(H,'euc');  % Euclidean distance of Hellinger-transformed data
      
   case {'spe'}
      dis = f_dis(SP,'euc'); % Euclidean distance of SP-transformed data
      
   case {'eucs'} % Euclidean Similarity (Elmore & Richman, 2001):
      dis = 1 - f_rewrap(D/max(D));
      
   otherwise
      dis = f_rewrap(D);     % wrap up into a symmetric distance matrix
end


if (nanFlag>0)
   % -----Final error check:-----
   if sum(isnan(dis(:)))>0
      error('Some distances are NaN! Remove rows that contain only zeros?');
   end
else
   % Replace NaN's with 1's
   idx      = find(isnan(dis));
   dis(idx) = 1;
end
