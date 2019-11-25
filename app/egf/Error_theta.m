
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Name: Error_theta.m
% Author: Tongkai Li
% mail: ltk@pku.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[ret,info]=Error_theta(theta,texp,yexp,rhs,y0,options,weight,group_index,gf_con,turnon,special_index)
% Error_theta: calculate the error with respect to theta only.
%   Input:
%       theta:  parameter.
%       remain: same as in other function.
%   Output:
%       ret:    error.
ttt=texp{1};
tspan=[ttt(1),ttt(end)];
ret=0;
info=1;
num_group=length(group_index); % number of groups

%   Some weight selected when testing.
%     weight=[0,0,0,1,0];
%     weight=[1,1,3,1,1]/7;
%     weight=[1,3,1,1,1,1];
%     weight=[1,1,1,1,1];

y = cell(num_group, 1); 
sum_weight = sum(weight(group_index)); 
for k=group_index
    theta(end-1)=gf_con(k);
    theta(end-2)=theta(end-2-k);
    [~,result]=forward_solver(rhs,texp{k},options,theta,y0);
    
    y{k}=result(:,8); 
    % When the solver status error.            
    if(length(y{k})~=length(texp{k}))
        ret=1e10;
        return;
    end
    
    % Estimate the error under this certain parameter.
    ret=ret+weight(k)/sum_weight*num_group*Error_ODE2EXP(y{k},texp{k},yexp{k});
%     %%% Kathy: should stay with physiological units nM and sec
%     ret=ret+(weight(k)/sum_weight)*Error_ODE2EXP(y{k},texp{k},yexp{k}'*theta(end)); 
    if(k==special_index && turnon==1)
    [num_crit,~]=Locate_crit(texp{k},y{k},tspan);
    if(num_crit<3)
        % ret=ret+10;
        info=0;
    end
    end
end

if(max(theta(length(theta)-2-group_index))\min(theta(length(theta)-2-group_index))>1.05/0.95)
%     ret=ret+10;
    info=0;
end

return;