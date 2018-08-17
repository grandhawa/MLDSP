function [err,PP] = f_cdaCV(x,grp,method,verb)
% - leave-one-out cross validation for Canonical Discriminant Analysis
%
% USAGE: [err,PP] = f_cdaCV(x,grp,method,verb);
%
% x   = input data (rows = objects, cols = variables)
% grp = column vector of integers specifying group membership
%
% method = method of classification:
%          1 : linear     
%          2 : quadratic
%          3 : mahalanobis
%          4 : centroid              (default)
%          5 : spatial median
%
% verb  = optionally display results (default = 1)
% 
% err = structure of results having the following fields:
%  .tot  = total error
%  .grp  = total error for each group
%  .conf = confusion matrix (proportion of CLS classified as GRP)
% 
% PP = posterior probabilities
%
% SEE ALSO: f_cda, f_cda632, classify

% -----Notes:-----
% This function is used to diagnose a Canonical Discriminant Analysis using
% the Leave-One-Out method of Cross-Validation. For methods 1-3, most of
% the work is done by the CLASSIFY function in the Matlab Statistics
% Toolbox. Method 4 ('centroid') follows Anderson (2002) and Anderson &
% Willis (2003).

% -----References:-----
% Anderson, M. J. 2002. CAP: a FORTRAN program for canonical analysis of
%  principal coordinates. Dept. of Statistics University of Auckland.
%  Available from: http://www.stat.auckland.ac.nz/PEOPLE/marti/
% Anderson, M. J. & T. J. Willis. 2003. Canonical analysis of principal
%  coordinates: a useful method of constrained ordination for ecology.
%  Ecology 84(2): 511-525.
% White, J. W. and B. I. Ruttenberg. 2007. Discriminant function analysis in
%  marine ecology: some oversights and their solutions. Mar. Ecol. Prog. Ser.
%  329: 301-305.

% -----Author:-----
% by David L. Jones, Dec-2003
%
% This file is part of the FATHOM Toolbox for Matlab and is released under
% the GNU General Public License, version 2.

% modified after f_capValid, written April-2003
% Feb-2004: re-writted to call f_errRate, now calculates confusion matrix,
%           fixed setting of defaults
% Oct-2010: added support for method = 4 ('centroid'), adapted from code
%           originally in f_capValid and f_cdaValid; changed default method to
%           4.
% Feb-2011: replaced 'unique' with 'f_unique'; use non-informative priors
%           following recommendations of Wright & Ruttenberg; return
%           posterior probabilities; added support for spatial median

% -----Set defaults & check input:-----
if (nargin < 3), method = 4; end % 'centroid' method by default
if (nargin < 4), verb   = 1; end % send output to display by default

switch method
   case 1, type = 'linear';
   case 2, type = 'quadratic';
   case 3, type = 'mahalanobis';
   case 4, type = 'centroid';
   case 5, type = 'spatial median';
   otherwise, error('Classification method must be 1, 2, 3, or 4!');
end
% -----------------------

n      = size(x,1);               % # of rows (observations)
uGrp   = f_unique(grp);           % unique groups, unsorted
noGrp  = size(uGrp,1);            % # of groups
priors = repmat(1/noGrp,1,noGrp); % non-informative PRIORS

cls    = repmat(NaN,n,1);         % preallocate
PP     = repmat(NaN,n,noGrp);

for omit = 1:n
   idx       = 1:n; % index of all observations
   idx(omit) = [];  % leave one observation out
   
   if (~isequal(type,'centroid')) && (~isequal(type,'spatial median'))
      [cls_loo,null,PP_loo] = classify(x(omit,:),x(idx,:),grp(idx),type,priors);
      cls(omit)  = cls_loo;
      PP(omit,:) = PP_loo;
   
   else
      % Use Centroid or Spatial Median:
      if (method==4), sm = 0; elseif (method==5), sm = 1; end
      
      % Perform CDA with 1 obs omitted, unwrap structure of results:
      result_loo = f_cda(x(idx,:),grp(idx),1,0,0,sm);
      Cvects     = result_loo.Cvects;
      centroids  = result_loo.centroids;
      
      % Center omitted obs using means of others:
      loo_ctr = x(omit,:) - mean(x(idx,:));
      
      % Project omitted obs in canonical space:
      loo_scores = loo_ctr*Cvects;
      
      % Find group centroid closest to omitted obs and classify accordingly:
      D = sum((repmat(loo_scores,noGrp,1) - centroids).^2,2)'; % squared distances
      L                = 1 - (D / max(D(:)));  % convert distance to likelihood
      PP(omit,:)       = L / sum(L);           % posterior probabilities
      [null,cls(omit)] = max(PP(omit,:),[],2); % classify omitted obs
   end
end

% Classification error rates:
err = f_errRate(grp,cls);

% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   fprintf('            F_CDA CROSS-VALIDATION\n'                 );
   fprintf('            Classification Success: \n'               );
   fprintf('--------------------------------------------------\n' );
   
   fprintf('Group        Corrrect  \n');
   for j=1:noGrp
      fprintf('%s %d %s %10.1f %s \n',['  '],uGrp(j),['     '],[1-err.grp(j)]*100,['%']);
   end
   
   fprintf('\n\n');
   fprintf('Total Correct  = %4.2f %s \n',(1-err.tot)*100,['%']);
   fprintf('Total Error    = % 4.2f %s \n',err.tot*100,['%']);
   fprintf('Class. method  = %s \n',type);
   if method<4, fprintf('Priors = non-informative\n'); end
   fprintf('\n--------------------------------------------------\n' );
   fprintf('     Confusion Matrix (%s): \n',['%'])
   hdr = sprintf('%-6.0f ',uGrp(:)');
   fprintf(['group: ' hdr]);
   fprintf('\n')
   for j=1:noGrp
      txt = [sprintf('%6.0f ',uGrp(j)) sprintf('%6.1f ',err.conf(j,:)*100)];
      fprintf(txt);
      fprintf('\n')
   end
   fprintf('\n==================================================\n\n');
end
