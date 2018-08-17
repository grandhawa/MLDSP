function f_RFplotPP(votes,txt)
% - plot posterior probabilities from f_RFclassPredict votes
% 
% USAGE: f_RFplotPP(votes)
% 
% votes = normalized weights obtained from f_RFclassPredict
% txt   = cell array of plot labels
% 
% SEE ALSO: f_RFclassPredict

% -----Author:-----
% by David L. Jones, Aug-2012
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% May-2013: add white space 

nc = size(votes,2); % get # of columns

figure;
for i = 1:nc
   subplot(1,nc,i)
   hist(votes(:,i),0.05:0.1:0.95);
   
   % Customize plot:
   set(gca,'TickDir','out');
   title(txt{i})
   xlabel(' '); % add some white space
   ylabel(' '); % add some white space
      
   % Adjust fill color:
   h = findobj(gca,'Type','patch');
   set(h,'FaceColor','k','EdgeColor','k')
   
   axisVar = axis;
   axis([axisVar(1) axisVar(2) 0 45])
   daspect([1 40 1])
   if i>1
      set(gca,'YTickLabel',[]);
   end
end
