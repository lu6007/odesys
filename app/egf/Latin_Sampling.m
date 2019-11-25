%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Name: Latin_Sampling.m
% Author: Tongkai Li
% mail: ltk@pku.edu.cn
% Created Time: 2018ï¿?7ï¿?7ï¿?æ˜ŸæœŸï¿?15ï¿?5ï¿?2ï¿?%%%%%%%%%%%%%%%%%%%%%%
% Latin hypercube sample

function[info,Error,Theta]= Latin_Sampling(M,scale,scaling1,scaling2,isok, ...
    var_index,vol_index, odesys_data)
root = odesys_data.root; 
exp_type = odesys_data.exp_type; 
% Latin_Sampling: Sample over the parameter space to find some 'good points' to start with.
%	Input:
%		M:		Devide number of the latin hypercube.
%		gf_con:	[gf] of different experiment group.
%		scale:	Specified value for each parameter.
%		texp:	time mesh for the experiment(cell type contain data of different group).
%		yexp:	varialbe value for the experiment(the same above).
%		@rhs et.al: Input for the forward solver.
%	Output:
%		Theta:	Sampled parameter value.
%		Error:	Error for each parameter value.
%       Info:   return information.
%	Inner parameter:
%		Max_sample:	Maximal sampling number.(some kind of early stop)
%		Max_step:	Maximal step number.
%		Step_sample:Number of latin queue in a certain latin hypercube.
%		Group_Index:Group index of the experiment that involved in this sampling.
%		Max_time:	Maximal runing time.(early abort)

    % Set up Max_time for simulation. By default, 12 hours.
    % For quick simulation, 1 h is enough to get a rough result.
	% For EGF simulation, we focous on the double peak pattern.
    
    % Max_time=3600*12; % seconds ; the optimization takes x60 time. 
    Max_time = odesys_data.max_total_time*60;  
%     Max_step=10000000;
%     tstart=tic;
    Max_sample = 10000; 
	sample_number=1;
        ng_number=0;
        bs_number=0;
        dp_number=0;
        bio_dp_number=0;
	% type=0;
	number_theta=length(var_index); 
    total_number=0;
    Min1=1e10;
    Min2=1e10;
   Theta_total(1,:)=zeros(length(scale),1);
   Error_total(1,:)=1e10;
   Theta_dp(1,:)=zeros(length(scale),1);
   Error_dp(1,:)=1e10;
   Theta_bs(1,:)=zeros(length(scale),1);
   Error_bs(1,:)=1e10;
   bio_Theta_dp(1,:)=zeros(length(scale),1);
   bio_Error_dp(1,:)=1e10;
% Each step we sample a latin hypercube. Then we may rescale the scale according to the error.
latin_number = zeros(number_theta, M); 
Theta = zeros(Max_sample, length(scale)); 
Error = zeros(Max_sample, 1); 
% for step=1:Max_step %%% Kathy, commented since this was not reached. 
step = 1; 
% Generate a latin hypercube.
    if(mod(step,1000)==0)
        fprintf("%d parameter combinations has been explored\n",step);
    end
    for i=1:number_theta
		latin_number(i,:)=randperm(M);
    end
    for i=1:M
        % For biochemical parameters, we may generate them in log space.
        % tic
        theta=odesys_data.theta_guess;
		theta(var_index)=exp(-M/2+latin_number(:,i)'-rand(1,number_theta)).*scale(var_index);
        % fix the volume factor: volfactor = 1
        theta(vol_index)=theta(intersect(var_index,vol_index)); 
        % theta=(latin_number(:,i)'+rand(1,number_theta)).*scale;
        [theta,error_temp,isdouble]=scaling1(theta);
        Theta(sample_number,:)=theta;
        Error(sample_number)=error_temp;
        switch exp_type
           case {'egf1', 'egf2'}
               % Theta(sample_number,:)=theta;
               % Sometimes we need to pick up parameters that lead to special patterns(double peak pattern for this case).
               % Then we choose another error estimator.
               if (isdouble==1)
               % theta,error_total]=scaling2(theta);
                   if(isok(theta)>=1)
                        dp_number=dp_number+1;
                        Theta_dp(dp_number,:)=theta;
                        Error_dp(dp_number)=error_temp;
                        if(theta(end)>=10)
                            bio_dp_number=bio_dp_number+1;
                            bio_Theta_dp(bio_dp_number,:)=theta;
                            bio_Error_dp(bio_dp_number)=error_temp;
                        end
                   end
                   if(isok(theta)==2)
                       dp_number=dp_number+1;
                       Theta_dp(dp_number,:)=theta;
                       Error_dp(dp_number)=error_temp;
                       bs_number=bs_number+1;
                       [theta,error_total]=scaling2(theta);
                       Theta_bs(bs_number,:)=theta;
                       Error_bs(bs_number)=error_total;
                   end
                    [theta,error_total]=scaling2(theta);
                    total_number=total_number+1;
                    Theta_total(total_number,:)=theta;
                    Error_total(total_number)=error_total;
                    if(Min2>error_total)
                        Min2=error_total;
                        fprintf(...
                            '   Min function value for total concentration renewed is %f\n',Min2);
                    end
                end % if(isdouble == 1)
           case {'pdgf1', 'pdgf2'}
           % Sometimes we need to pick up parameters that lead to 
           % special patterns(small total error here).Then we choose 
           % another error estimator.           
           if(error_temp<=10)
                [theta,error_total]=scaling2(theta);
                total_number=total_number+1;
                Theta_total(total_number,:)=theta;
                Error_total(total_number)=error_total;
                if(Min2>error_total)
                   Min2=error_total;
                    fprintf( ...
                        '   Min function value for total concentration renewed is %f\n',Min2);
                end
           end % if(error_temp <=10)

        end % switch case 'egf', 'pdgf'
        

         if (Min1>error_temp)
               Min1=error_temp;
               fprintf('      Min function value for single concentration renewed is %f.\n',Min1);
         end
         sample_number=sample_number+1;
         if sample_number > Max_sample
             fprintf('Latin_Sampling: Max_sample = %d has been reached. \n', Max_sample); 
             break; 
         end
    end % for i = 1:M
    
    % Save latin_sample.mat
    file_name = [root, 'latin_sample.mat'];
    fprintf('Latin_Sampling: saving %s\n', file_name); 
    %%% Kathy: clean up Theta and Error
    sample_number = sample_number - 1; % it is important to correct sample_number
    temp = Theta(1:sample_number, :); clear Theta; 
    Theta = temp; clear temp;
    temp = Error(1:sample_number); clear Error;
    Error = temp; clear temp;
    save(file_name, 'Theta', 'Error', 'Theta_total', 'Error_total', 'Theta_bs', ...
           'Error_bs', 'Theta_total', 'Error_total', 'Theta_dp', 'Error_dp', ...
           'bio_Theta_dp'); 
    
    % time_used=toc(tstart);
    % if(time_used>=Max_time)
    %	clear Theta Error
        fprintf( ...
            '\n ng_number = %d, bs_number = %d, total_number = %d, sample_number = %d\n\n', ...
            ng_number, bs_number, total_number, sample_number); 
        info=1;
return;



























