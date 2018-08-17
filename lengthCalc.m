function [ maxLen, minLen, meanLen, medLen ] = lengthCalc( Seq )
    %function calculates length stats
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Gurjit Singh Randhawa  %
    % Department of Computer Science,%
    % Western University, Canada     %
    % email: grandha8@uwo.ca         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    len = cellfun('length',Seq);
    maxLen = max(len);
    minLen = min(len);
    meanLen = round(mean(len));
    medLen = round(median(len));    
end
