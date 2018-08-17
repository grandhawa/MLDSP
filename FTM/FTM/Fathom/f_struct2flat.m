function out = f_struct2flat(d,fname)
% - turn a structure into a flat table & optionally export
%
% USAGE out = f_struct2flat(struct,'fname')
%
% struct = input structure
% fname  = optional filename for export
% out    = cell array of flattened structure
%
% SEE ALSO: f_extractFields

% STRUCT2FLAT turns a structure into a flat table
% also produces a CSV (comma sep. text file)
% for importing int Excel, Access, or any other DB
%
% It supports fields that are different lengths
% It removes pesky nested cells
% It can easily be modified to output relational
% database tables E-mail me to find out how
% IT'S NOT FANCY, BUT IT'S USEFUL

% -----Authors:---
% This is a slight modification of Michael Robbins'
% struct2flat.m <michael.robbins@us.cibc.com>
% Original code posted to news://comp.soft-sys.mat &
% http://www.mathworks.com/matlabcentral/fileexchange/
%
% slightly modified by David L. Jones, Aug-2002
%
%
% This file is part of the FATHOM Toolbox for Matlab.

% -----Check input:-----[DJ]
if (isstruct(d)<1)
   error('A structure is required as input');
end

col_width=30;
fn=fieldnames(d);
out=[];
for record=1:length(d)
   [l,w]=size(out);
   if l>0 recnum=l+1; else recnum=1; end;
   for col=1:length(fn)
      temp=eval(sprintf('[d(%d).%s]', ...
         record,fn{col}));
      if iscell(temp)&~iscellstr(temp)
         temp2=temp{:};
      else
         temp2=temp;
      end;
      temprecnum=recnum;
      for linenum=1:length(temp2)
         if iscellstr(temp2(linenum))
            out{temprecnum,col}= ...
               temp2{linenum};
         else
            temp3=temp2(linenum);
            while iscell(temp3) 
               temp3=[temp3{:}]; 
            end;
            if ~isempty(temp3)
               out{temprecnum,col}=temp3;
            end;
         end; % if iscellstr
         temprecnum=temprecnum+1;
      end; % for linenum
   end; % for col
end; % for record

if (nargin >1) % optionally export if given a filename [DJ]
   %%%%%%%%%%%%%
   %% CSV FILE
   %%%%%%%%%%%%%
   [l,w]=size(out);
   fid=fopen(fname,'w');

   temp=sprintf('%s,',fn{:});
   fprintf(fid,'%s\n',temp(1:end-1));
   %%%%%%%% MAKE FORMAT STRING
   for j=1:w
      temp=[out{:,j}];
      if isnumeric([out{:,j}])
         if all([out{:,j}]-floor([out{:,j}])==0)
            f{j}='%d';
         else
            f{j}='%f';
         end;
      else
         f{j}='%s';
      end;
   end;
   
   %%%%%%% PRINT THE FILE
   if 0
      % I'd like this to work, but it doesn't
      aaa=sprintf('%s,',f{:});
      
      bbb=out';
      fprintf(fid,[aaa(1:end-1) '\n'],bbb{:});
   else
      for i=1:l
         if i>1
            for j=1:w
               if isempty(out{i,j})
                  out{i,j}=out{i-1,j};
               end;
               fprintf(fid,f{j},out{i,j});
               if j<w
                  fprintf(fid,',');
               else
                  fprintf(fid,'\n');
               end;
            end;
         end;
      end;
   end;
   fclose(fid);        
end
