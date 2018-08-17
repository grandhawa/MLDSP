function f_exportR(x,x_row,x_col,fname)
% - export data for import into R
%
% USAGE: f_exportR(x,x_row,x_col,'fname');
%
% -----Input:-----
% x     = matrix of numerical values
% x_row = cell array of row labels
% x_col = cell array of column labels
% fname = output filename
%
% -----R commands for import/export:-----
% read.table("fname.txt",header=TRUE,row.names=1,sep=",")
%
% write.table(x,file="fname.txt",sep=",",quote=FALSE)
%
% SEE ALSO: f_import, f_export

% -----Notes:-----
% This script is used to exchange data frames with R. Data exported from R
% using write.table() can be imported using F_IMPORT. That data can be
% subsequently exported from Matlab using F_EXPORTR and imported into R using
% read.table().
%
% This function requires TBLWRITE from the Statistics Toolbox.

% -----Author:-----
% by David L. Jones, Sep-2003
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

tblwrite(x,char(x_col),char(x_row),fname,',');
