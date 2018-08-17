%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Standard Reference Material for LA-ICP-MS:   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by David L. Jones
% 
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Data compiled in file 'SRM.xls', then imported into Malab file 'SRM.mat'

%%%%%%%%%%%%%%%%%%%%%%%%
%    GP-4 carbonate:   %
%%%%%%%%%%%%%%%%%%%%%%%%
% 
% gp4 = USGS 'internal' SRM used mostly for cross validation
%  .txt = cell array of isotope labels
%  .ppm = parts per million
% 
% -> Data from certificate of analysis that ships with the standard; values for
%    Ca provided via email from Stephen A Wilson <swilson@usgs.gov>


%%%%%%%%%%%%%%%%%%%%%%%%
%       MACS-3:        %
%%%%%%%%%%%%%%%%%%%%%%%%
% 
% macs3 = Microanalytical Carbonate Standard (MACS-3)
%  .txt = cell array of isotope labels
%  .ppm = parts per million
% 
% -> Data from certificate of analysis that ships with the standard
% 
% -----Check Fe:-----
% Define molecular wts:
Fe_mw = 55.845; 
O_mw  = 15.9994;
% 
% Define percent of oxides in matrix (Overall average, Table 3, Pearce et al, 1997):
Fe2O3 = 1.60;
% 
% Convert percent of oxide to ppm of element:
% Note: ppm/10000 = percent
Fe_ppm = (Fe2O3*( (2*Fe_mw) /(2*Fe_mw + 3*O_mw)))*10000; % = 11191 ppm

%%%%%%%%%%%%%%%%%%%%%%%%
%       NIST 610:      %
%%%%%%%%%%%%%%%%%%%%%%%%
% 
% nist610 = NBS-610 concentrations from Pearce et al. 1997
%  .txt = cell array of isotope labels
%  .ppm = parts per million
% 
% -> Data for trace elements from 'Preferred average' column of Table 8 in
% Pearce et al., 1997
% 
% -----Data for elements comprising the matrix are as follows:-----
% 
% Define molecular wts:
Ca_mw = 40.078; 
Na_mw = 22.9897;
Al_mw = 26.9815;
Si_mw = 28.0855;
O_mw  = 15.9994;
% 
% Define percent of oxides in matrix (Overall average, Table 1, Pearce et al, 1997):
CaO   = 11.450;
Na2O  = 13.352;
Al2O3 = 2.039;
SiO2  = 69.975;
% 
% Convert percent of oxide to ppm of element:
% Note: ppm/10000 = percent
Ca_ppm = (CaO*(Ca_mw/(Ca_mw + O_mw)))*10000;             % =  81832.09 ppm
Na_ppm = (Na2O*( (2*Na_mw) /(2*Na_mw + O_mw)))*10000;    % =  99052.73 ppm
Al_ppm = (Al2O3*( (2*Al_mw) /(2*Al_mw + 3*O_mw)))*10000; % =  10791.41 ppm
Si_ppm = (SiO2 * (Si_mw/(Si_mw + 2*O_mw)))*10000;        % = 327087.59 ppm
% 
% -> reverse engineer AMS's percents:
(85263.71/10000)/(Ca_mw/(Ca_mw + O_mw))
(97190.21/10000)/( (2*Na_mw) /(2*Na_mw + O_mw))
(10053.47/10000)/( (2*Al_mw) /(2*Al_mw + 3*O_mw))
(334048.00/10000)/(Si_mw/(Si_mw + 2*O_mw))
% 
% Define SD of percent of oxides (for Pepita):
CaO_sd   = 0.231;
Na2O_sd  = 0.681;
Al2O3_sd = 0.157;
SiO2_sd  = 0.391;
% 
% Convert SD (in percent of oxide) to SD (in ppm of element):
Ca_sd = (CaO_sd*(Ca_mw/(Ca_mw + O_mw)))*10000;             % = 1650.94 ppm
Na_sd = (Na2O_sd*( (2*Na_mw) /(2*Na_mw + O_mw)))*10000;    % = 5052.05 ppm
Al_sd = (Al2O3_sd*( (2*Al_mw) /(2*Al_mw + 3*O_mw)))*10000; % =  830.92 ppm
Si_sd = (SiO2_sd * (Si_mw/(Si_mw + 2*O_mw)))*10000;        % = 1827.67 ppm



%%%%%%%%%%%%%%%%%%%%%%%%
%       NIST 612:      %
%%%%%%%%%%%%%%%%%%%%%%%%
% 
% nist612 = NBS-612 concentrations from Pearce et al. 1997
%  .txt = cell array of isotope labels
%  .ppm = parts per million
% 
% -> Data for trace elements from 'Preferred average' column of Table 9 in
% Pearce et al., 1997; Data from Mg from Gao et al., 2002
%
% -----Data for elements comprising the matrix are as follows:-----
% 
% Define molecular wts:
Ca_mw = 40.078; 
Na_mw = 22.9897;
Al_mw = 26.9815;
Si_mw = 28.0855;
O_mw  = 15.9994;
% 
% Define percent of oxides in matrix (Overall average, Table 2, Pearce et al, 1997):
CaO   = 11.93;
Na2O  = 13.98;
Al2O3 = 2.11;
SiO2  = 71.90;
% 
% Convert percent of oxide to ppm of element:
% Note: ppm/10000 = percent
Ca_ppm = (CaO*(Ca_mw/(Ca_mw + O_mw)))*10000;             % =   85262.61 ppm
Na_ppm = (Na2O*( (2*Na_mw) /(2*Na_mw + O_mw)))*10000;    % =  103711.59 ppm
Al_ppm = (Al2O3*( (2*Al_mw) /(2*Al_mw + 3*O_mw)))*10000; % =   11167.18 ppm
Si_ppm = (SiO2 * (Si_mw/(Si_mw + 2*O_mw)))*10000;        % =  336085.72 ppm
% 
% Define SD of percent of oxides (for Pepita):
CaO_sd   = 0.22;
Na2O_sd  = 0.56;
Al2O3_sd = 0.16;
SiO2_sd  = 0.96;
% 
% Convert SD (in percent of oxide) to SD (in ppm of element):
Ca_sd = (CaO_sd*(Ca_mw/(Ca_mw + O_mw)))*10000;             % = 1572.32 ppm
Na_sd = (Na2O_sd*( (2*Na_mw) /(2*Na_mw + O_mw)))*10000;    % = 4154.40 ppm
Al_sd = (Al2O3_sd*( (2*Al_mw) /(2*Al_mw + 3*O_mw)))*10000; % =  846.80 ppm
Si_sd = (SiO2_sd * (Si_mw/(Si_mw + 2*O_mw)))*10000;        % = 4487.38 ppm


%%%%%%%%%%%%%%%%%%%%%%%%
%       NIST 614:      %
%%%%%%%%%%%%%%%%%%%%%%%%
% 
% nist614 = NIST 614 concentrations from Gao et al., 2002
%  .txt = cell array of isotope labels
%  .ppm = parts per million
% 
% -> Data for trace elements from Table 3 in Gao et al., 2002
% 
% Define molecular wts:
Ca_mw = 40.078; 
Si_mw = 28.0855;
O_mw  = 15.9994;
% 
% Define percent of oxides in matrix (Table 10 in Jochum et al., 2011):
CaO   = 11.90;
SiO2  = 72.1;
% 
% Convert percent of oxide to ppm of element:
% Note: ppm/10000 = percent
Ca_ppm = (CaO*(Ca_mw/(Ca_mw + O_mw)))*10000;             % =   85048.00 ppm
Si_ppm = (SiO2 * (Si_mw/(Si_mw + 2*O_mw)))*10000;        % =  337020.58 ppm
% 
% Define SD of percent of oxides (for Pepita):
CaO_sd   = 0.20;
SiO2_sd  = 0.90;
% 
% Convert SD (in percent of oxide) to SD (in ppm of element):
Ca_sd = (CaO_sd*(Ca_mw/(Ca_mw + O_mw)))*10000;             % = 1429.40 ppm
Si_sd = (SiO2_sd * (Si_mw/(Si_mw + 2*O_mw)))*10000;        % = 4206.90 ppm


%%%%%%%%%%%%%%%%%%%%%%%%
%     NIST GLASSES:    %
%%%%%%%%%%%%%%%%%%%%%%%%
% To determine the spatial resolution of our Laser Ablation instrument, we
% follow the method of Sanbord & Telmer (2003) and consider all 3 NIST
% glasses (610, 612, & 614) to have 12% CaO.
% 
% Define molecular wts:
Ca_mw = 40.078; 
O_mw  = 15.9994;
% 
% Define percent of oxides in matrix:
CaO   = 12.0;
% 
% Convert percent of oxide to ppm of element:
% Note: ppm/10000 = percent
Ca_ppm = (CaO*(Ca_mw/(Ca_mw + O_mw)))*10000;             % = 85763 ppm
% 
% Setup Internal Standard (IS):
IS.txt = {'Ca'};
IS.ppm = 85763;  % NIST Glasses
