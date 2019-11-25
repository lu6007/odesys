%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Name: forward_solver.m
% Author: Tongkai Li
% mail: ltk@pku.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[t,result] =forward_solver(rhs,tspan,options,theta,y0)
%forward_solver: choose an ODE solver to solve the prolbem
%	input:
%		@rhs:	right hand site function .
%		tsapn:	time interval to solve in.
%		y0:		initial value of the variable.
%		theta:	parameter value.
%		option:	general options for the solver.
%	output: 
%		t:		time mesh 
%		result:	variable value with respect to time.

% tstart=tic;
[t,result]=ode15s(rhs,tspan,y0,options,theta);
% timeused=toc(tstart);
% if(timeused>20)
%     pause;
% end
end





















