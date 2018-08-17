function f_beep(n)
% - play a custom sound file
% 
% USAGE: f_beep(n);
% 
% n = # times to play sound file (default = 1)
% 
% SEE ALSO: beep, wavread

% -----Author:-----
% by David L. Jones 
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 1), n = 1; end % default play sound file once

% Specify location of sound file on your system:
fname = '~/work/mWork/Fathom/resources/sounds/Sosumi.wav';

if ~exist(fname,'file')
   disp(['Cannot find: ' fname '...'])
   error('Edit the FNAME variable in F_BEEP!')
end
% -------------------------------------

% Load audio data:
[Y, Fs] = wavread(fname);

% Play audio n times:
for i=1:n
   sound(Y,Fs);
end

