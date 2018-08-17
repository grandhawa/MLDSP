function [ numSeq ] = numMappingAtomic( sq )
%function for Atomic representation
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
       if(t=='A')
           numSeq(K) = 70;
       elseif(t=='C')
           numSeq(K) = 58;
       elseif(t=='G')
           numSeq(K) = 78; 
       elseif(t=='T')
           numSeq(K) = 66;
       end           
   end   
end

