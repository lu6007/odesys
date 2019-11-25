%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Name: Coor_Descent.m
% Author: Tongkai Li
% mail: ltk@pku.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Error,x_optimal,info]= Coor_Descent(x0,f,var_index,scaling)
% Coor_Descent: Apply coordinate descent method to find the optimal value of the @cost_function from x0.
%	Input:
%		x0:			initial value of the parameter.
%		f:	objective function to be optimized.	
%	Output:
%		Error:		The optimal value of the objective function.
%		x_optimal:	Parameter value to meet the optimal.
%		info:		Information of the results of the algorithm.
%	PS:
%		Here we apply the armijo line search.
%		The gradient are approximated by finite difference.
%	Inner parameter:
%		pho:		armijo line search parameter.
%		max_iter:	maximal iteration number.
%		max_armi:	maximal iteration in armijo line search.
		pho=0.5;
        tol=1e-1;
        %%% Kathy
		max_iter=2; % 50; % the max number of iterations needed
		max_armi=20;
        var_zindex=zeros(max(var_index),1);
		n=length(x0);
		delta=1e-6;
        x0=scaling(x0);
        xsave=x0;
		x=x0;
		Error=0;
		x_optimal=x0;
        info=0;

		for i=1:max_iter
			gradient=zeros(n,1);
%
			queue=randperm(length(var_index));
% Here we search over all directions, can be modified.
            xp=xsave;
%            tic
            for k=1:length(queue)
             ftemp=f(xp);
             x=xp;
             success=0;
            j=var_index(queue(k));
            [f_dr,isdouble]=f(xp+delta*[zeros(1,j-1),1,zeros(1,n-j)]);
            
 % If certain constraint is not satisfied, we set the gradient in the direction to zero
 % as 'bound' are reached.
            if(isdouble==0)
                gradient(j)=0;
                continue;
            else
                gradient(j)=(f_dr-ftemp)/delta;
            end
% When current position approaches zero and the search direction points out
% of the region, we reset the gradient and stops seraching in this
% direction as the gradient caculation is rough at this time.
            if(var_zindex(j)==1)
                if(gradient(j)<=0)
                  var_zindex(j)=0;
                else
                    gradient(j)=0;
                   continue;
                end
            end
            if(gradient(j)>0&& x(j)<tol*delta)
                gradient(j)=0;
                var_zindex(j)=1;
                continue;
            end
            
            if(isnan(gradient(j)))
                continue;
            end
                for l=0:max_armi
					x(j)=xp(j)-gradient(j)*pho^l;
% Bioparameters are all nonnegetive so we should hold these constraint when 
% we are doing line search.
					if(min(x)<0)
						continue;
					end
% line search principle, to be modified latter.
% isdouble here can be modified in the function to handle different needs
% or constraint in the resulted pattern.
                  [fvalue,isdouble]=f(x);
					if ((fvalue-ftemp)<0&&isdouble ==1)
					success=1;
					break;
					end
                end
               if(success==1)
				xp=x;
                end
            if(success==0)
                    gradient(j)=0;
                    x=xp;
                    continue;
            end

            end
% We do the analytic solution in the scaling parameter.
            [x,~]=scaling(x);
            % fprintf("norm of gradient is %10.4e\n",norm(gradient));
            if (norm(gradient)<1e-3)
                 % fprintf('Optmization complete as gradient < 1e-3.\n');
                 x_optimal=xsave;
		         Error=f(xsave);
                 info=1;
                return;
            end
            if(f(xsave)-f(x)<1e-3)
              % fprintf('Optmization finished as difference between two steps are small.\n');
                x_optimal=xsave;
		        Error=f(xsave);
                info=1;
                return;
            end
 %           toc
            xsave=x;
            % fprintf("norm of function is %f\n",f(x));
        end
        % fprintf("Optimization finished as Max_iter has been reached.\n");
		info=2;
		x_optimal=x;
		Error=f(x);


