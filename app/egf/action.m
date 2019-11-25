% function data = action(exp_type)
% exp_type = 'egg2' (default) or 'pdgf2'
%
% Example:
% >> data = action();
% or, 
% >> exp_type = 'pdgf2';
% >> data = action(exp_type);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Name: Action.m
% Author: Tongkai Li, Shaoying Lu (shaoying.lu@gmail.com)
% mail: ltk@pku.edu.cn
% Created Time: 2018ï¿?7ï¿?0ï¿?æ˜ŸæœŸï¿?20ï¿?2ï¿?3ï¿?%%%%%%%%%%%%%%%%%%%%%%
% This script is used to plan to solve the parameter fiiting problem.
% 
%%% Kathy: 10/28/2019
% To-do list:
% * Fix scaling factor
% * Disable volume scale
% * Look for parameters within [1/10 10] * initial guess and select the
% ones closest to the intial guess. 
% * Optimize prediction from hypercubic Initial All and Double-peack value
% * Copy to a server computer to Run

% Done:
% * Change Prepare to a function : prepare()
% * Add '/' to the end of folders
% * Changed "results.mat" to "result.mat"
% * Change watch.m and Action.m to functions watch() and action() respectively. 
% * Changed the error to be in unit nM, independent of tspan and number of
% groups. Modified the functions Error_ODE2EXP and Error_theta()

% 
function data = action(exp_type)
if nargin == 0 || isempty(exp_type)
    exp_type = 'egf2';
end
   
%Set up some function handles and parameter initialization.
data = prepare(exp_type); 
scale = data.scale; 
f = data.f;
f_total = data.f_total;
Scaling4 = data.Scaling4;
Scaling_fixd = data.Scalingfixd; 
Scaling_fixt = data.Scalingfixt; 
Is_ok_db = data.Is_ok_db; 
theta = data.theta;
theta_guess = data.theta_guess;
total_index = data.total_index; 
special_index = data.special_index; 
error_threshold = data.error_threshold; 

M=6;
% Apply latin sampling to get nice initial values.
% [info,~,~]=Latin_Sampling(M,scale,Scaling4,Scaling23a, ...
%     Is_ok_db,[1:22,length(theta)-2-special_index], ...
%     (length(theta)-2-length(total_index):length(theta)-3), ...
%     data);
Latin_Sampling(M,scale,Scaling4,Scaling_fixt, ...
    Is_ok_db,[1:22,length(theta)-2-special_index], ...
    (length(theta)-2-length(total_index):length(theta)-3), data);

% load ../data/latin_sample.mat
file_name = [data.root, 'latin_sample.mat'];
fprintf('Action.m : loading data from %s\n', file_name);
temp = load(file_name); 
Error = temp.Error;
Theta = temp.Theta; clear temp;

% For egf simulations, we use error_dp and theta_dp as initials as we think
% the double peak patterns might be more important.
time_start = tic; 

% Error=Error_total; 
% Theta=Theta_total;
[~,sample_error_index]=sort(Error);
Error=Error(sample_error_index);
Theta=Theta(sample_error_index,:);

% Initialization for some parameters.
Time=fix(clock);
% reportaddress='../data/';
reportaddress = data.root; 
Opt_Theta=Theta;
Opt_Theta_t=Theta;
Init_Theta=Theta;
Opt_Thetas=Theta;
Error_Theta=zeros(size(Theta,1),1);
Error_Thetas=zeros(size(Theta,1),1);
Error_Theta_t=zeros(size(Theta,1),1);
% addressString1=sprintf('%s%d%d%d%d_data.mat',reportaddress,2018,Time(4),Time(5),Time(6));
% save(addressString1,'Theta_dp','Error_dp','Theta_bs','Error_bs');
addressString2=sprintf('%sresult.mat',reportaddress);

%% Main loop
i_temp=1;
fprintf('Index \t Single Value (nM) \t Initial All (nM) \t'); 
fprintf('All Value (nM)\t Ratio Initial/All \t Time (s) \t\n'); 
scale_func_double = Scaling_fixd; % fix scaling factor, % Scaling23 for not fixing scaling factor
scale_func_total = Scaling_fixt; % fix scaling factor, % Scaling23a for not fixing scaling factor
for i = 1:size(Theta, 1) %%% Kathy changed for running error 10/20
    % fprintf('Index %d starting.\n',i);
    theta=Theta(i,:);
    % fix scaling parameters
    theta(end)= theta_guess(end);
    init_total_error = sqrt(f_total(theta)); 
    % Mtheta is the theta value after fitting the double-peak pattern
    [error, Mtheta, ~]=Coor_Descent(theta,f,[1:22,length(theta)-2-special_index], ...
        scale_func_double);
    double_error = sqrt(f(Mtheta)); 
    % Function value of the double peak concentration 
    fprintf('%5d \t %17.4f %17.4e\t', i, double_error, init_total_error);  
    Opt_Thetas(i_temp,:) = Mtheta;
    Error_Thetas(i_temp) = error;
    if double_error > error_threshold % nM
        fprintf('\n');
        continue;
    end
    %
    [error, Mtheta, ~]=Coor_Descent(Mtheta,f_total,(23:length(Mtheta)-3),scale_func_total);
    Opt_Theta(i_temp,:)=Mtheta;
    Error_Theta(i_temp)=error;
    %
    [error, Mtheta, ~]=Coor_Descent(Mtheta,f_total,(1:length(Mtheta)-3),scale_func_total);
    time_used = toc(time_start); 
    time_start = tic; 
    total_error = sqrt(f_total(Mtheta)); 
    % Function value after volumne scaling, and time_used
    fprintf('%17.4f \t %9.4e \t %8.2f \n', total_error, ...
        init_total_error/total_error, time_used);
    %
    Opt_Theta_t(i_temp,:)=Mtheta;
    Error_Theta_t(i_temp)=error;
    save(addressString2,'Opt_Theta','Error_Theta','Init_Theta','Opt_Thetas', ...
        'Error_Thetas','Opt_Theta_t','Error_Theta_t');
    % fprintf("The %d initial parameter has been optimized.\n",i_temp);
    %%% Kathy: what's the difference between i_temp and i? 
    %%% Ans: Not much difference for now. 
    i_temp=i_temp+1;
    
end
return

