% function [sol, fval] = test_optimize()
% Example:
% >> test_optimize; 
% Ref: MATLAB documentation - Solve Constrained Nonlinear Optimization, Problem-Based
function [sol, fval] = test_optimize()
% change x, y to a single variable
x = optimvar('x', 2);
objective_fun = @objective_function; 
constraint_fun = @constraint_function; 
object = fcn2optimexpr(objective_fun, x);
problem = optimproblem('Objective',object);
% The constraint can be included in the function  
% because it is a rational function of the optimization variables.
tilt_ellipse = (constraint_fun(x)<= 0);
problem.Constraints.constr = tilt_ellipse; 
% Initial guess of the variables
x0.x = [-3; 3]; 
show(problem);

% Sovle the problem
[sol, fval] = solve(problem, x0);
disp('Solution 1:');
disp([sol.x' fval])

% Another initial value for a local solution
x0.x = [-1; 1];
[sol2,fval2] = solve(problem,x0);
disp('Solution 2:');
disp([sol2.x' fval2]);

% Visualize
f = @(x,y) objective_fun([x; y]);
g = @(x,y) constraint_fun([x;y]); % constraint = 0
range = [-5.5 -0.25 -0.25 7];
fimplicit(g,'k-')
axis(range);
hold on
fcontour(f,range,'LevelList',logspace(-1,1))
plot(sol.x(1),sol.x(2),'ro','LineWidth',2)
plot(sol2.x(1),sol2.x(2),'ko','LineWidth',2)
legend('Constraint','f Contours','Global Solution','Local Solution', ...
    'Location','northeast');
caxis([0,1]); colorbar; 
hold off
end


function f = objective_function(x)
f = exp(x(1)).*(4*x(1).^2 + 2*x(2).^2 + 4*x(1).*x(2) + 2*x(2) - 1);
end

function g = constraint_function(x)
g = (x(1).*x(2)/2 + (x(1)+2).^2 + (x(2)-2).^2/2 - 2); % <=0
end