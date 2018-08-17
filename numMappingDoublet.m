function [ numSeq ] = numMappingDoublet( sq )
%function for Nearest-neighbor based doublet representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gurjit Singh Randhawa  %
% Department of Computer Science,%
% Western University, Canada     %
% email: grandha8@uwo.ca         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
len = length(sq);
    doublet = {'AA','AT','TA','AG','TT','TG','AC','TC','GA','CA','GT','GG','CT','GC','CG','CC'};
    alpha=0;
    numSeq = zeros(1,len,'double');
    kStrings=(2*alpha)+1;
    for K = 1:len
        if(alpha==0)
            if(K<len)
                t = sq(K:K+1);
            else
                t = strcat(sq(K:K),sq(1:1));
            end
            tp = find(strcmpi(doublet, t))-1;
            numSeq(K) = tp;
        else
            loc = 0;
            for index = K-alpha:K+alpha
                sPos = index;
                ePos = sPos+1;
                if(sPos<1)
                    sPos = len+sPos;
                elseif(sPos>len)
                    sPos = sPos-len;
                end
                if(ePos>len)
                    ePos = ePos-len;
                elseif(ePos<1)
                    ePos = ePos+len;
                end
                if(sPos==len && ePos==1)
                    t = strcat(sq(len:len),sq(1:1));
                else
                    t = sq(sPos:ePos);
                end
                loc = loc+(find(strcmpi(doublet, t))-1);
            end
            tp = loc/kStrings;
            numSeq(K) = tp;
        end
    end       
end

