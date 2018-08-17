function result = f_ranks(X)
% - ranks data in X (column-wise) with averaging of ties
%
% USAGE: result = f_ranks(X)
%
% X  = column vector, matrix, or symmetric distance matrix
% Xr = ranked data

% -----References:-----
% McCune, B., J. B. Grace, and D. L. Urban. 2002. Analysis of Ecological
% Communities. MjM Software Design. Gleneden Beach, Oregon. iv + 300 pp.

% -----Author:-----
% originally Ranks.m by Lutz Duembgen, 23.02.1999
%
% modified by David L. Jones, 2001
%
% This file is part of the FATHOM Toolbox for Matlab.

% 26-Mar-02: added special handling of symmetric distance matrices
% Apr-2003:  improved checking for symmetric distance matrices; preserve 0
%            values
% Feb-2010: overhauled coded to add column-wise support for matrices

% Symmetric distance matrix requires special handling:
if (f_issymdis(X) == 1)
   result = f_rewrap(f_ranks(f_unwrap(X)));
else
   [n,nc] = size(X);
   result = zeros(n,nc); % preallocate
   
   for i=1:nc % repeat for each column separately
      Xc       = X(:,i); % extract column
      Xr       = zeros(n,1); % preallocate
      [hv,ar] = sort(Xc);
      a       = 1;
      
      for b = 2:n
         if hv(b) > hv(a)
            hv(a:b-1) = (a+b-1)/2;
            a = b;
         end
      end
      
      hv(a:n) = (a+n)/2;
      Xr(ar)  = hv;
      
      % Preserve 0 values:
      if (sum(logical(Xc(:) == 0))>0);
         minVar = min(Xr(:));
         Xr(:)  = Xr(:) - minVar(1);
      end
      result(:,i) = Xr; % add this column to result
   end
end



