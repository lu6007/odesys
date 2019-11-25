


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Name: Error_theta.m
% Author: Tongkai Li
% mail: ltk@pku.edu.cn
% Created Time: 2018ï¿?8ï¿?2ï¿?æ˜ŸæœŸï¿?10ï¿?4ï¿?4ï¿?%%%%%%%%%%%%%%%%%%%%%%
function[theta,ret,info]=Scaling(theta,texp,yexp,rhs,y0,options,weight, ...
    group_index,gf_con,turnon,special_index)
% Error_theta: calculate the error respect to theta only.
%   Input:
%       theta:  parameter.
%       remain: same as in other function.
%       tunron: scaling options.
%           case 0: just scaling, do not consider double peak constraint.
%           case 1: scaling while Latin_sampling. First match the first among the two peaks,
%           then calculate scaling factor.
%           case 2: scaling while optimizing. Calculate the scaling factor
%           while considering double peak constraint.
%           case 3: fix the scaling parameter while needed.
%   Output:
%       ret:    error.
ret=0;
info=1;
if(turnon==1)
   for k=special_index
    theta(end-1)=gf_con(k);
    theta(end-2)=theta(end-2-k);
            [t{k},result]=forward_solver(rhs,texp{k},options,theta,y0);
% Estimate the error under this certain parameter.
tt=texp{k};
       [~,locate_org]=Locate_crit(texp{k},yexp{k},[tt(1),tt(end)]);
                y{k}=result(:,8)/theta(end);
        [num_crit,locate]=Locate_crit(t{k},y{k},[tt(1),tt(end)]);
%         visualize(t,y,texp,yexp,k,showfigure);
% 		visualize(t,y,texp,yexp,k,2);
        if(num_crit==0)
            info=0;
            break;
        end
%         if(num_crit>=2)
%             fprintf('doublepeak');
%         end
            peak_time(1)=locate(1,1);
            peak_height(1)=locate(1,2);
            undimen_time=peak_time(1)/locate_org(1,1);
            undimen_conc=peak_height(1)/locate_org(1,2);
             theta=undimen(theta,undimen_time,1);
             theta(end)=theta(end)*undimen_conc;
        if(num_crit<3)
            info=0;
        else
            peak_time(2)=locate(3,1);
            peak_height(2)=locate(3,2);
            if(peak_height(2)<peak_height(1)&&peak_time(2)*peak_time(1)/locate_org(1,1)<tt(end))
                info=1;
            else
                info=0;
            end
        end
   end
end

    N=length(group_index);
    ret=0;
    ret3=0;
    ret1=0;
    ret2=0;
% Find the experiment time mesh that is in the interval.
     
for k=group_index
            rtemp=0;
            rtemp1=0;
            rtemp2=0;
            rtemp3=0;
			theta(end-1)=gf_con(k);
            theta(end-2)=theta(end-2-k);
            ttexp=texp{k};
            yyexp=yexp{k};
            [t{k},result]=forward_solver(rhs,ttexp,options,theta,y0);
% Estimate the error under this certain parameter.
            y=result(:,8)/theta(end);
            if(length(y)~=length(ttexp))
                ret=1e10;
                return;
            end
      for i=1:length(ttexp)-1
		rtemp3=rtemp3+(ttexp(i+1)-ttexp(i))*yyexp(i)^2;
        rtemp1=rtemp1+(ttexp(i+1)-ttexp(i))*y(i)^2;
        rtemp2=rtemp2+(ttexp(i+1)-ttexp(i))*y(i)*yyexp(i);
      end 
      if (turnon>=1 && k==special_index)
          tt=texp{k};
          [num_crit,~]=Locate_crit(t{k},y,[tt(1),tt(end)]);  
          if(num_crit<3)
            info=0;
          end
      end
			ret3=ret3+weight(k)/sum(weight(group_index))*N*rtemp3;
            ret1=ret1+weight(k)/sum(weight(group_index))*N*rtemp1;
            ret2=ret2+weight(k)/sum(weight(group_index))*N*rtemp2;
end
if (turnon==1||turnon==3)
    ret=ret3+ret1-2*ret2;
    return;
end
ret=(ret3*ret1-ret2^2)/ret1;
theta(end)=theta(end)*ret1/ret2;