function D = f_latlong(X)
% - compute terrestrial distance matrix
%
% Usage: km = f_latlong(X);
% 
% -----INPUT:-----
% X = matrix with degrees of longitude (1st column) and degrees of latitudes (2nd column)
%
% -----OUTPUT:-----
% - symmetric pairwise distance matrix (in kilometers)

% -----Authors:-----
% ==============================================================================
% Copyright (c) 1997 B. Planque - Sir Alister Hardy Foundation for Ocean Science
% bp@wpo.nerc.ac.uk
% Permission is granted to modify and re-distribute this code
% in any manner as long as this notice is preserved.
% All standard disclaimers apply.
% ==============================================================================
% slightly modified after Planque's distance.m
% by David L. Jones, Aug-2001
%
% This file is part of the FATHOM Toolbox for Matlab.

[n,p]=size(X);
D=zeros(n);
% TERRESTRIAL DISTANCE
D=zeros(n);
R=6370;				% Earth Radius (in kilometers)
r=(360/(2*pi));	% rapport degres:radians
for i=1:n-1
   for j=i+1:n
      lon1=X(i,1);lon2=X(j,1);lat1=X(i,2);lat2=X(j,2);
      D(i,j)=R*(acos((sin(lat1/r)*sin(lat2/r))+(cos(lat1/r)*cos(lat2/r)*cos((lon2-lon1)/r))));
      D(j,i)=D(i,j);
   end;
end;

   
   