function txt = myupdatefcn(~,event_obj,Y, AcNmb)
%to show the information about selected point in 3D plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gurjit Singh Randhawa  %
% Department of Computer Science,%
% Western University, Canada     %
% email: grandha8@uwo.ca         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos = get(event_obj,'Position');
x = pos(1);
y = pos(2);
z = pos(3); 
for a=1:length(AcNmb)
   if(Y(a,1)==x && Y(a,2)==y && Y(a,3)==z)
       idx = a;
       break;
   end
end
txt = {['ANum:',AcNmb{idx}]}
end