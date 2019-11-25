%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Name: rhs_2.m
% Author: Tongkai Li
% mail: ltk@pku.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f] = rhs_2(t, x, theta)
%  rhs:
%  rhs(t, x, theta) is the right hand side of the ODE dx/dt=f(t, x, theta),
%  including  function values respected to theta, x and time t.
%  where x represent the variables and theta are the set of unknown parameters.
%  The ODE system below is the nearest version.

	n = length(x);
	f = zeros(n, 1);
    if (t < 0)
		return ;
    end
    
% gf concentration are places at the 2nd place backward.
% volume parameter are picked at the 3nd place bacward.
    gf=theta(end-1);
    vol=theta(end-2);
    egfr_total=398.10/vol; % nM 
    fyn_total=174.98/vol; % nM
    ptp_total=296.58/vol; % nM
    sensor_total=500/vol; % nM
    
    gfgfr=x(1);
    dgfgfr=x(2);
    dgfgfrp=x(3);
    endo=x(4);
    deg=x(5);
    ptpa=x(6);
    fyna=x(7); % active fyn
    sensora=x(8);
    
    gfr=egfr_total-gfgfr-2*(dgfgfr+dgfgfrp+endo+deg);
    fyn=fyn_total-fyna;
    ptp=ptp_total-ptpa;
    sensor=sensor_total-sensora;
    
    kon1=theta(1);
    koff1=theta(2);
    kon2=theta(3);
    koff2=theta(4);
    kon3=theta(5);
    kcatoff3=theta(6);
    kdoff3=theta(7);
    kon4=theta(8);
    koff4=theta(9);
    kon5=theta(10);
    kcaton6=theta(11);
    kdon6=theta(12);
    Vmaxoff6=theta(13);
    kmoff6=theta(14);
    kcaton7=theta(15);
    kdon7=theta(16);
    kcatoff7=theta(17);
    kdoff7=theta(18);
    kcaton8=theta(19);
    kdon8=theta(20);
    kcatoff8=theta(21);
    kdoff8=theta(22);
    
%     kcaton9=kcaton7;
%     kdon9=kdon7;
%     Vmaxoff9=kcatoff7*ptpa;
%     Vmaxoff9=Vmaxoff9*ptpa;
%     kmoff9=kdoff7;
    
    v1=kon1*gf*gfr-koff1*gfgfr;
    v2=kon2*gfgfr^2-2*koff2*dgfgfr;
    v3=kon3*dgfgfr-kcatoff3*ptpa*dgfgfrp/(kdoff3+dgfgfrp);
    v4=kon4*dgfgfrp-koff4*endo;
    v5=kon5*endo;
    v6=kcaton6*dgfgfrp*ptp/(kdon6+ptp)-Vmaxoff6*ptpa/(kmoff6+ptpa);
    v7=kcaton7*dgfgfrp*fyn/(kdon7+fyn)-kcatoff7*fyna/(kdoff7+fyna);
    v8=kcaton8*(fyna)*sensor/(kdon8+sensor)-kcatoff8*sensora/(kdoff8+sensora);
    
    f(1)=v1-v2;
    f(2)=0.5*v2-v3;
    f(3)=v3-v4;
    f(4)=v4-v5;
    f(5)=v5;
    f(6)=v6;
    f(7)=v7;
    f(8)=v8;
    
end


