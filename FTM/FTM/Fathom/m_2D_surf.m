function [h] = m_2D_surf(long,lat,data)
%  - draw a 2D surface on an M_Map map
%    M_2D_SURF(LONG,LAT,DATA,...) draws a 2D surface on a map. Behavior
%    is the same as for SURF(X,Y,Z) with viewpoint above and interpolated shading
%    except that LONG and LAT vectors or matrices must be specified.
%
%    [H] = M_2D_surf(...) returns a handle H to the SURF object.
%
%    See also SURF

% Brian Farrelly, brian.Farrelly@nho.hydro.com, 3 June 1999
% based on m_contour by Rich Pawlowicz (rich@ocgy.ubc.ca) 17/Jan/1998

% This file is provided with the FATHOM Toolbox for Matlab.


global MAP_PROJECTION

% Have to have initialized a map first

if isempty(MAP_PROJECTION),
   disp('No Map Projection initialized - call M_PROJ first!');
   return;
end;


if min(size(long))==1 & min(size(lat))==1,
   if size(data) ==  [length(lat) length(long)]
      [long,lat]=meshgrid(long,lat);
   elseif size(data) ==  [length(long) length(lat)]
      [lat,long]=meshgrid(lat,long);
   else
      disp('Data size not in agreement with longitude/latitude vectors');
      return;
   end;
end;   

if size(lat)~=size(long)
   disp('Latitude size not in agreement with longitude size');
   return;
end

if size(lat)==size(data')
   data=data';
end

if size(lat)~=size(data)
   disp('Data size not in agreement with longitude/latitude matrices')
end

[X,Y]=m_ll2xy(long,lat,'clip','on');

i=isnan(X);      % For these we set the *data* to NaN...
data(i)=NaN;

% And then recompute positions without clipping. THis
% is necessary otherwise contouring fails (X/Y with NaN
% is a no-no. 
if any(i(:)), [X,Y]=m_ll2xy(long,lat,'clip','off'); end;  

Z=zeros(size(data));

if any(~i(:)),
   h=surf(X,Y,Z,data);
   shading interp;
   view(0,90);
   set(h,'tag','m_2D_surf');
else
   h=[];
end;