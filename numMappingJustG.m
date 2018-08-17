function [ numSeq ] = numMappingJustG( sq)
%function for JustG representation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gurjit Singh Randhawa  %
% Department of Computer Science,%
% Western University, Canada     %
% email: grandha8@uwo.ca         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a = 'G';
    len = length(sq);  
    numSeq = zeros(1,len,'double');
    for K = 1:len
       t = sq(K);
       if(strcmpi(t,a))
           numSeq(K) = 1;
       else
           numSeq(K) = 0;
       end         
   end   
end



