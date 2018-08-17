function [ numSeq ] = numMappingInt( sq )
%function for Integer representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gurjit Singh Randhawa  %
% Department of Computer Science,%
% Western University, Canada     %
% email: grandha8@uwo.ca         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dob = {'T', 'C' , 'A', 'G' };
    len = length(sq);  
    numSeq = zeros(1,len,'double');
    
    for K = 1:len
       t = sq(K);
       tp = find(strcmpi(dob, t))-1;
       numSeq(K) = tp;        
   end       
end

