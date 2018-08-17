function [ newCm ] = checkDimension( cm, aLabel, pLabel, n)
%function checks the dimensions of confusion matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gurjit Singh Randhawa  %
% Department of Computer Science,%
% Western University, Canada     %
% email: grandha8@uwo.ca         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB by default delete zero rows n colms from the confusion matrix
    [x,y] = size(cm);
    if(x==n && y==n)
        newCm = cm;
    else
    %newCm
    newCm = zeros(n,n);
    for i=1:n
        idx = find(aLabel==i);
        if(length(idx)==0)
            continue;            
        else
            temp = pLabel(idx);
            [uvals, ~, uidx] = unique(aLabel);
            output = [uvals, accumarray(uidx, 1)];
            [x,y]=size(output);
            for j=1:x
                 idx = output(j,1);
                 newCm(i,idx)=output(j,2);
            end
        end
    end
end

