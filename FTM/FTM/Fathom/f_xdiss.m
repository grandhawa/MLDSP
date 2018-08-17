function xdist = f_xdiss(dist,dcrit,usemin,eps,bigd)
% - calculate Extended Dissimilarities from a symmetric distance matrix
%
% USAGE: xdist = f_xdiss(dist,dcrit,usemin,eps,bigd);
%
% dist   = symmetric distance matrix
%
% dcrit  = critical distance threshold (default = 1)
%          all distances greater than or equal to dcrit
%          will be recalculated as extended dissimilarities
%
% usemin = use the minimum extended dissimilarity
%          (~shortest path); 1=true (default), 0=false
%
% eps    = epsilon: precision value (defaut = 0.00001)
%
% bigd   = maximum extended dissimilarity allowed
%
% xdist  = symmetric extended dissimilarity matrix
%
% SEE ALSO: f_dis

% -----Reference:-----
% De'ath, G. 1999. Extended dissimilarity: a method of robust
% estimation of ecological distances from high beta diversity data.
% Plant Ecology 144(2):191-199.

% -----Authors:-----
% original C code "xdiss.c" by Glenn De'ath<glenn.death@jcu.edu.au>
% obtained from: http://esa.sdsc.edu/Archive/E080-013/pcurve.zip
%
% ported to Matlab by David L. Jones, April-2001
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% 28-Mar-02: check dist, make call to f_rewrap

% -----Check input and set default values:------
if (f_issymdis(dist) == 0)
   error('Input DIST must be square symmetric distance matrix');
end;

if (nargin < 2), dcrit    = 1;      end
if (nargin < 3), usemin   = 1;      end
if (nargin < 4), eps      = 0.0001; end 
if (nargin < 5), bigd     = 10000;  end
% ------------------------------------------------
big   = bigd;            % maximum dissimilarity
n     = size(dist,1);    % get dimensions of input distance matrix
dist  = f_unwrap(dist);  % extract lower tridiagonal as a column vector
ndist = n*(n-1)/2;       % number of dissimilarity values
nbig  = 0;               % initialize counter

for i=0:(ndist-1)
   if (dist(i+1)>dcrit-eps) % see if dissimilarity is greater than critical distance 
      dist(i+1) =  -1;      % flag dissimilarities > critical distance
      nbig = nbig + 1;      % count # of dissimilarities > critical distance
   end
end

nloop = 0;            % initialize WHILE loop counter
fprintf('Loop %1d : Number of critical distances = %5d \n', nloop, nbig);  

while (nbig>0)
   nbigold=nbig;
   nloop = nloop + 1; % increment WHILE loop counter
   
   if (nloop==1)      % first iteration thru WHILE loop
      maxd = 0.0;     % added to prevent MEX compile error [DJ]
      adj=dcrit;
   else               % all subsequent iterations thru WHILE loop
      if (nloop>1); adj=maxd; end;
   end;
   
   maxd = 0.0;
   mind = 1e+009;
   
   for j=1:(n-1)       
      for i=0:(j-1)
         pos = n*i-i*(i+1)/2+j-i-1;
         if (dist(pos+1) < 0) % do the following for dissimilarities flagged above
            mink  = 1e+009;
            sum   = 0.0;
            count = 0;
            
            for k=0:(n-1)
               if ((k~=i) & (k~=j))
                  if (k>i)
                     ki = n*i-i*(i+1)/2+k-i-1;
                  else 
                     ki = n*k-k*(k+1)/2+i-k-1;
                  end;
                  if (k>j)
                     kj = n*j-j*(j+1)/2+k-j-1;
                  else 
                     kj = n*k-k*(k+1)/2+j-k-1;
                  end;
                  dki=dist(ki+1);
                  dkj=dist(kj+1);
                  if ((dki>=0) & (dkj>=0) & (dki<big) & (dkj<big))
                     if ((dki+dkj) < mink); mink=dki+dkj; end;
                     sum = sum + (dki+dkj);
                     count = count+1;
                  end;
               end;
            end;
            
            if (count>0)
               if (usemin>0)  % usemin = True
                  dist(pos+1) = mink + big;
               else           % usemin = False
                  dist(pos+1) = sum/count + big;
               end;
               if (dist(pos+1)>maxd); maxd = dist(pos+1); end;
               if (dist(pos+1)<mind); mind = dist(pos+1); end;
            end
         end
      end
   end
   
   maxd=maxd-big;
   mind=mind-big;
   fprintf('Max and min distances = %7.2f %7.2f \n', maxd, mind);
   
   
   nbig = 0;
   for i=0:(ndist-1)
      if (dist(i+1)<0)
         nbig = nbig + 1; % count again the number of distances flagged
      end;
   end;
   fprintf('Loop %1d : Number of critical distances = %5d \n', nloop, nbig);
   
   if (nbigold==nbig) % ----- See if data are disjoint -----
      nbig=0;         % ----- exits WHILE loop when nbig=0 -----
      fprintf('Data disjoint...Ordinate groups separately! \n');  
   end;
   
   for i=0:(ndist-1)
      if (dist(i+1) > big)
         if (usemin>0)
            dist(i+1) = dist(i+1) - big;
         else
            dist(i+1) = dist(i+1) - big - mind + adj + eps;
         end;
      end;
   end;
   
   maxd = maxd - mind + adj;
   
end; % end the WHILE loop

xdist = f_rewrap(dist); % turn back into symmetric distance matrix

