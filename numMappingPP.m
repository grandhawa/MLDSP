function [ numSeq ] = numMappingPP( sq )
%function for Purine/Pyramidine representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gurjit Singh Randhawa  %
% Department of Computer Science,%
% Western University, Canada     %
% email: grandha8@uwo.ca         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    len = length(sq);  
    numSeq = zeros(1,len,'double');
    for K = 1:len
       t = sq(K);
       if(strcmpi(t,'A'))
           numSeq(K) = -1;
       elseif(strcmpi(t,'C'))
           numSeq(K) = 1;
       elseif(strcmpi(t,'G'))
           numSeq(K) = -1; 
       elseif(strcmpi(t,'T'))
           numSeq(K) = 1;
       end           
   end   
end

