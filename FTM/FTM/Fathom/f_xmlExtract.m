function result = f_xmlExtract(xml,opt,nodeName)
% - extract node, name, attribute, or tag from XML tree
%
% USAGE: result = f_xmlExtract(xml,opt,'nodeName');
%
% xml    = xml tree as a char array
% opt: 1 = full node
%      2 = full node (self-closing)
%      3 = full node with dummy
%      4 = name part
%      5 = attributes part
%      6 = attributes last
%     
% nodeName = name of node to extract (for opt = 1-3)
% 
% result   = structure with the following fields
%   .match = cell array containing the text of each match
%   .idx   = array of [start end] indices of each match

% by D. Jones, July-2007

% -----Notes & References:-----
% The regex's used here were posted by "marielle" to 'forums.runrev.com. She
% notes the advantage of using [\s\S]+ is that it will capture any character,
% including an end of line one, while .+ doesn't if the right tag isn't used.   

% -----To do:-----
% add support for dummy chunk

% Jan-2008: updated documentation


% -----Check input and set defaults:-----
if (nargin < 2); opt = 999; end;


% Setup pattern to match:
switch opt
   case 1 % Full node:
      patn = ['(?s)(\s*<' nodeName '(\s+[^>]*)?>[\s\S]+?</' nodeName '>\s*)'];
   
   case 2 % Full node (self-closing):
      patn = ['(?s)(\s*<' nodeName '(\s+[^>]*)? />\s*)'];
   
   case 3 % Full node with dummy:
      patn = ['(?s)(<' nodeName' '[^>]*>[\s\S]*' chunkVar...
         '[\s\S]*</' nodeName '>)'];

   case 4 % Name part:
      patn = ['<(\w+)'];
      
   case 5 % Attributes part:
   patn = ['<[^\s]+\s+([\s\S]*?)(\s+/)?>'];
      
   case 6 % Attributes last:
      patn = ['(\s*(\w+)=([^=]+?)$)'];
   
   otherwise
   error('Method unknown!');
end


% Match regular expression:
[q1,q2,q3] = regexp(xml,patn,'match','start','end');


% Wrap results into a structure:
result.match = q1';           % cell array containing the text of each match
result.idx   = [q2(:) q3(:)]; % array of [start end] indices of each match



