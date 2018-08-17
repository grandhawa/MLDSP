function [ result ] = numMappingRandom3( s )
%function for random among 3 representations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gurjit Singh Randhawa  %
% Department of Computer Science,%
% Western University, Canada     %
% email: grandha8@uwo.ca         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a = 3;
    r = randi([1 a],1,1);
    if(r==1)
         result = numMappingPP(s);    
    elseif(r==2)
         result = numMappingReal(s);    
    elseif(r==3)
         result = numMappingJustA(s);   
    end 
end
