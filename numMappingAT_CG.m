function [ numSeq ] = numMappingAT_CG( sq )
%function for PairedNumeric representation
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
           numSeq(K) = 1;
       elseif(t=='C')
           numSeq(K) = -1;
       elseif(t=='G')
           numSeq(K) = -1; 
       elseif(t=='T')
           numSeq(K) = 1;
       end           
   end   
end

