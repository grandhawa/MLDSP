function [ numSeq ] = numMappingEIIP( sq )
%function for EIIP representation
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
           numSeq(K) = 0.1260;
       elseif(t=='C')
           numSeq(K) = 0.1340;
       elseif(t=='G')
           numSeq(K) = 0.0806; 
       elseif(t=='T')
           numSeq(K) = 0.1335;
       end           
   end   
end

