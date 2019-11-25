%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Name: Error_ODE2EXP.m
% Author: Tongkai Li
% mail: ltk@pku.edu.cn
% Created Time: 2018年07月17日 星期二 20时51分46秒
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ret=Error_ODE2EXP(result,texp,yexp)
%Error_ODE2EXP: Error between the ODE simulated result and the experiment data.
%	Input:
%		t:		time mesh of the ODE.
%		result:	Variable value with respect to the mesh t.
%		texp:	time mesh of the experiment.
%		yexp:	Varialbe value with respect to the mesh texp.
%		tspan:	Interval for calculating the error.
%	Output:
%		ret:	the value of error function

ret=0;
% 	for i=1:length(texp)-1
% 		ret=ret+(texp(i+1)-texp(i))*(yode(i)-yexp(i))^2;
%  	end
% 
for i=1:length(texp)-1
    ret=ret+(texp(i+1)-texp(i))*(result(i)-yexp(i))^2;
end
%   ret=norm(result-yexp,2);

return;













