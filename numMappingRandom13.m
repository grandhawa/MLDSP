function [ result ] = numMappingRandom13( s )
%%function for random among 13 representations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gurjit Singh Randhawa  %
% Department of Computer Science,%
% Western University, Canada     %
% email: grandha8@uwo.ca         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a=13;
    r = randi([1 a],1,1);
    if(r==1)
         result = numMappingPP(s);    
    elseif(r==2)
         result = numMappingReal(s);    
    elseif(r==3)
        result = numMappingJustA(s);   
    elseif(r==4)
       result = numMappingAtomic(s); 
    elseif(r==5)
       result = numMappingEIIP(s);  
    elseif(r==6)
        result = numMappingInt(s); 
    elseif(r==7)
        result = numMappingAT_CG(s); 
    elseif(r==8)
        result = numMappingDoublet(s); 
    elseif(r==9)
       result = numMappingCodons(s); 
    elseif(r==10)
       result = numMappingIntN(s);  
    elseif(r==11)
       result = numMappingJustC(s); 
    elseif(r==12)
       result = numMappingJustG(s);  
    elseif(r==13)
       result = numMappingJustT(s); 
    end 
end
