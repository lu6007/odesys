function [ret ] = Is_OK( theta,rhs,texp,y0,options,group_index,gf_con,special_index)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for k=group_index
			theta(end-1)=gf_con(k);
            theta(end-2)=theta(end-2-k);
            [~,res]=forward_solver(rhs,texp{k},options,theta,y0);
% Estimate the error under this certain parameter.
            result{k}=res(:,8)/theta(end);
              if(length(result{k})~=length(texp{k}))
                ret=0;
                return;
            end
end
ret=0;
n=length(group_index);
time1=zeros(length(group_index));
time2=time1;
peak1=time1;
peak2=time1;
dp=time1;
for i=group_index
   t=texp{i};
   y=result{i};
   [num,crit]=Locate_crit(t,y,[t(1),t(end)]);
   if(num>=1)
       time1(i)=crit(1,1);
       peak1(i)=crit(1,2);
   end
   if(num>=3)
       time2(i)=crit(3,1);
       peak2(i)=crit(3,2);
       dp(i)=1;
   end
end



% yy=result{2};
% if(yy(end)>peak1(2))
% ans=2;
% return;
% end
if(dp(special_index)==1)
    ret=1;
end

if(length(find(dp==1))<=2)
    ret=2;
    return;
end


end

