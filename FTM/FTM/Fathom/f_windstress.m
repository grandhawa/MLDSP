function tau = f_windstress(mag,z)
% - wind stress in dynes/cm^2
%
% Usage: tau = f_windstress(mag,z);
%
% mag = wind speed in m/s
% z     = height above sea surface of wind sensor in meters
%         (default = 10)
%
% tau  = wind stress (dynes/cm^2)
%
% See also: f_windCman, f_vecUV, f_vecPlot, f_ekmanDepth

% ----- Notes: -----
% This function is used to calculate wind stress (tau) given wind speed
% (in m/s). Wind stress increases with the square of the wind speed.
%
% For conversions: 1 dynes/cm^2 = 0.1 pascals (N/m^2)

% -----References:-----
% Kourafalou, V. H., T. N. Lee, and L.-Y. Oey. 1996. The fate of river discharge
% on the continental shelf 2. Transport of coastel low-salinity waters under
% realistic wind and tidal forcing. J. Geophys. Res. 101: 3435-3455.
%
% Lee, T. N. and E. Williams. 1999. Mean distribution and seasonal variability
% of coastal currents and temperature in the Florida Keys with implications
% for larval recruitment. Bull. Mar. Sci. 64: 35-56.

% -----Author:-----
% by David L. Jones, Sep-2001
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Oct-2003: re-corrected units for density of air

% ----- Check input & set defaults: -----
if (nargin < 2), z = 10; end; % default no correction for instrument height

if size(z) ~= [1 1]
	error('Z must be a scalar');
end

if (size(mag,2)>1)
   error('MAG must be a column vector');
end
% ---------------------------------------

Pa = 1.225;  % density of air at sealevel in kg/m^2
Cd = 0.0015; % drag coefficient

mag_10 = mag .* (10/z)^0.1;      % calculate wind speed at standard 10 m height
tau    = (Pa*Cd) .* (mag_10).^2; % tau in N/m^2
tau    = tau .* 10;              % convert to dynes/cm^2

