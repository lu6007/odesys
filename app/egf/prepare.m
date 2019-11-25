% function data = prepare(exp_tye)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Name: prepare.m
% Author: Tongkai Li
% mail: ltk@pku.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to load data and initiate parameter.
function data = prepare(exp_type)
root = '/Users/kathylu/Documents/sof/odesys/data/1019_2019/';
data.exp_type = exp_type; 
% 2 sec of hypercube and 120 sec of optimization
data.max_total_time = 20; % 60;  % sec % need 1200 sec to get double-peak
mt = data.max_total_time/60; 
fprintf('prepare(): Experiment type = %s\n', exp_type); 
fprintf('Estimated latin Sampling time = %d sec, %10.3f hr\n', ...
    mt, mt/3600.0);
fprintf('    Estimated total optimization time = %d sec, %10.3f hr\n', ...
    mt*60, mt/3600.0 * 60);

switch data.exp_type
    case 'egf1'
        data.root = [root, 'egf/1/']; 
        data.var_index = (1:27);  % rhs % this index defines the parameters that will be optimized. 
        theta = ones(1,30);
        gf_con = 80.65/500*[500,100,25,10,5,1];
        data.conc_str = {'500 ng/ml', '100 ng/ml',  '25 ng/ml', '10 ng/ml', ...
            '5 ng/ml', '1 ng/ml'};
        data.para_nametag={'k_{on1}','k_{off1}','k_{on2}','k_{off2}','k_{on3}','k_{catoff3}', ...
            'k_{doff3}','k_{on4}','k_{off4}','k_{on5}','k_{caton6}','k_{don6}', ...
            'V_{maxoff6}','k_{moff6}','k_{caton7}','k_{don7}','k_{catoff7}',...
            'k_{doff7}','k_{caton8}','k_{don8}','k_{catoff8}','k_{doff8}'};
        theta_guess=theta;
        theta_guess(end)=0.08659;
        total_index=[1,2,3,4,5]; 
    case 'egf2'
        data.root = [root, 'egf/2/']; 
        data.var_index = (1:27);  % rhs % this index defines the parameters that will be optimized. 
        theta = ones(1,30);
        gf_con = 80.65/500*[500,100,25,10,5,1];
        data.conc_str = {'500 ng/ml', '100 ng/ml',  '25 ng/ml', '10 ng/ml', ...
            '5 ng/ml', '1 ng/ml'};
        data.para_nametag={'k_{on1}','k_{off1}','k_{on2}','k_{off2}','k_{on3}','k_{catoff3}', ...
            'k_{doff3}','k_{on4}','k_{off4}','k_{on5}','k_{caton6}','k_{don6}', ...
            'V_{maxoff6}','k_{moff6}','k_{caton7}','k_{don7}','k_{catoff7}',...
            'k_{doff7}','k_{caton8}','k_{don8}','k_{catoff8}','k_{doff8}'};
        theta_guess=theta;
        theta_guess(end)=0.08659;
        total_index=[1,2,3,4,5]; 
        error_threshold = 6; 

    case 'egf2_copy'
        data.root = [root, 'egf/2_copy_1029_2019/']; 
        data.var_index = (1:27);  % rhs % this index defines the parameters that will be optimized. 
        theta = ones(1,30);
        gf_con = 80.65/500*[500,100,25,10,5,1];
        data.conc_str = {'500 ng/ml', '100 ng/ml',  '25 ng/ml', '10 ng/ml', ...
            '5 ng/ml', '1 ng/ml'};
        data.para_nametag={'k_{on1}','k_{off1}','k_{on2}','k_{off2}','k_{on3}','k_{catoff3}', ...
            'k_{doff3}','k_{on4}','k_{off4}','k_{on5}','k_{caton6}','k_{don6}', ...
            'V_{maxoff6}','k_{moff6}','k_{caton7}','k_{don7}','k_{catoff7}',...
            'k_{doff7}','k_{caton8}','k_{don8}','k_{catoff8}','k_{doff8}'};
        theta_guess=theta;
        theta_guess(end)=0.08659;
        total_index=[1,2,3,4,5]; 
        % data.error_threshold = 600; 

    case 'pdgf1'
        data.root = [root, 'pdgf/1/'];
        data.var_index = (1:28); % rhs
        theta = ones(1,31); 
        gf_con = 80.65/500*[1,5,10,25,50,100]; % something is wrong
        data.conc_str = {'1 ng/ml', '5 ng/ml', '10 ng/ml', '25 ng/ml', ...
            '50 ng/ml', '100 ng/ml', '50+ ng/ml', '100+ ng/ml'}; 
        data.para_nametag={'k_{on1}','k_{off1}','k_{on2}','k_{off2}','k_{on3}', ...
            'k_{catoff3}','k_{doff3}','k_{on4}','k_{off4}','k_{on5}', ...
            'k_{caton6}','k_{don6}','V_{maxoff6}','k_{moff6}', ...
            'k_{caton7}','k_{don7}','k_{catoff7}','k_{doff7}', ...
            'k_{caton8}','k_{don8}','k_{catoff8}','k_{doff8}'};
        theta_guess=theta;
        theta_guess(end)=0.08659;
        total_index=[1,2,3,4,5,6];
    case 'pdgf2'
        data.root = [root, 'pdgf/2/'];
        data.var_index = (1:28); % rhs
        theta = ones(1,31); 
        gf_con = 80.65/500*[1,5,10,25,50,100]; % something is wrong
        data.conc_str = {'1 ng/ml', '5 ng/ml', '10 ng/ml', '25 ng/ml', ...
            '50 ng/ml', '100 ng/ml', '50+ ng/ml', '100+ ng/ml'}; 
        data.para_nametag={'k_{on1}','k_{off1}','k_{on2}','k_{off2}','k_{on3}', ...
            'k_{catoff3}','k_{doff3}','k_{on4}','k_{off4}','k_{on5}', ...
            'k_{caton6}','k_{don6}','V_{maxoff6}','k_{moff6}', ...
            'k_{caton7}','k_{don7}','k_{catoff7}','k_{doff7}', ...
            'k_{caton8}','k_{don8}','k_{catoff8}','k_{doff8}'};
        theta_guess=theta;
        theta_guess(end)=0.08659;
        total_index=[1,2,3,4,5,6];
        error_threshold = 5; 
end % switch data.exp_type
%%% Kathy
% theta_guess(end) = 865.9; % nM/AU converting from FRET o concentration
fprintf('Function prepare(): set theta_guess(end) to %f. \n', theta_guess(end));

% Load experiment data. 
file_name = [data.root, '../input/exp.mat'];
fprintf('prepare(): loading %s\n', file_name); 
temp = load(file_name); 
Tq = temp.Tq; Yq = temp.Yq; YqL = temp.YqL; YqU = temp.YqU; clear temp; 
rand('seed',sum(100*clock));

% Set up interval and select data point according to it.
tspan=[0;2350];
for i=1:length(Tq)
    temp_index=find(Tq{i}>=tspan(1)&Tq{i}<=tspan(2));
    Tq{i}=Tq{i}(temp_index);
    Yq{i}=Yq{i}(temp_index);
    YqL{i}=YqL{i}(temp_index);
    YqU{i}=YqU{i}(temp_index);    
end
texp=Tq; 
yexp=Yq;

%set up r.h.s function, initial value, variable, index, nametags, options
%and so on.
rhs=@rhs_2;
y0=[0,0,0,0,0,0,0,0];
% var_index=[1:27];
% theta=zeros(1,30);
special_index=3;

% gf_con = [80.65, 16.13, 8.065, 4.033, 1.515, 0.76, 0.1515];
% gf_con = 80.65/500*[500,100,50,25,10,5,1];
% gf_con = 80.65/500*[500,100,25,10,5,1];
% gf_con = 80.65/500*[1,5,10,25,50,100,200,500];

options=odeset('AbsTol',1e-16,'RelTol',1e-8,'MaxStep',1,'Stats','off');
data.line_type = {'r-', 'g-', 'b-', 'k-', 'y-', 'm-', 'c-'};

% conc_str = {'500 ng/ml', '100 ng/ml', '50 ng/ml', '25 ng/ml', '10 ng/ml', '5 ng/ml', '1 ng/ml'};
% conc_str = {'500 ng/ml', '100 ng/ml',  '25 ng/ml', '10 ng/ml', '5 ng/ml', '1 ng/ml'};
% conc_str = {'1 ng/ml', '5 ng/ml', '10 ng/ml', '25 ng/ml', '50 ng/ml', '100 ng/ml', '50+ ng/ml', '100+ ng/ml'}; 
%set up initial scale guess for each variable.

% scale=[3*1e-3,1e-3,1e-3,3*1e-3,1e-3,0.5,3,0.5,3,1e-3,0.5,3,0.5,3,1e-3,...
% 0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5];
% scale=[0.003,0.006,0.0004,0.07,0.002,0.001,0.001,0.001,0.001,30,0.45,...
% 100,1.8,100,0.1,0.1,0.1,0.1,0.005,9.3,0.005,240,1];
data.scale=0.36*[ones(1,length(theta)-1),100];
% para_nametag={'k_{on1}','k_{off1}','k_{on2}','k_{off2}','k_{on3}','k_{catoff3}','k_{doff3}',...
%      'k_{on4}','k_{off4}','k_{on5}','k_{caton6}','k_{don6}','V_{maxoff6}','k_{moff6}',...
%      'k_{caton7}','k_{don7}','k_{catoff7}','k_{doff7}','k_{caton8}','k_{don8}',...
%      'k_{catoff8}','k_{doff8}'};
% 22 kinetic parameter; 6 volume related parameter; 1*[gf]; 1* scaling factor. 
% PDGF prepare.m should have different number of parameters. 
 
% total_index=[1,2,3,4,5];
weight=[1,1,1,1,1,1,1];
%function handle for simplicity.
% f:Fucntion handle for the loss function of specific pattern(group_index).
% f_show: show pictures during calculation(only when testing algorithm).
% f_show_total: show pictures of all concentrations.
% f_total: Function hanlde for loss function summation of all concentrations.
% scaling4: Scaling for the specific pattern(double peak group_index).
% Scaling23: Scaling for the specific pattern (no need to be doublepeak).
% Scaling23a: Scaling for all concentrations.(solve for the
% readout/concentration scaling factor).
% Is_ok_db: Testing double peak and bell-shape pattern. (Adding constraints
% for specific pattern, pattern selection).
data.Scalingfixd=@(theta)(Scaling(theta,texp,yexp,rhs,y0,options,weight,special_index,gf_con,3,special_index));
data.Scalingfixt=@(theta)(Scaling(theta,texp,yexp,rhs,y0,options,weight,total_index,gf_con,3,special_index));
switch data.exp_type
    case {'egf1', 'egf2'}
        turn_on = 0;
    case {'pdgf1', 'pdgf2'}
        turn_on = 1;
end
data.f_show=@(theta)(Error_theta(theta,texp,yexp,rhs,y0,options,weight,special_index,gf_con,turn_on, special_index));
data.f_show_total=@(theta)(Error_theta(theta,texp,yexp,rhs,y0,options,weight,total_index,gf_con,turn_on, special_index));
data.f=@(theta)(Error_theta(theta,texp,yexp,rhs,y0,options,weight,special_index,gf_con,turn_on, special_index));
data.f_total=@(theta)(Error_theta(theta,texp,yexp,rhs,y0,options,weight,total_index,gf_con,turn_on, special_index));
data.Scaling4=@(theta)(Scaling(theta,texp,yexp,rhs,y0,options,weight, ...
    special_index,gf_con,turn_on,special_index));
turn_on = 0; 

data.Scaling23=@(theta)(Scaling(theta,texp,yexp,rhs,y0,options,weight,special_index,gf_con,turn_on,special_index));
data.Scaling23a=@(theta)(Scaling(theta,texp,yexp,rhs,y0,options,weight,total_index,gf_con,turn_on,special_index));
data.Is_ok_db=@(theta) (Is_OK(theta,rhs,texp,y0,options,total_index,gf_con,special_index));

data.gf_con = gf_con;
data.texp = texp; 
data.option = options; 
data.y0 = y0; 
data.weight = weight;
data.yexp = yexp; 
data.tspan = tspan; 
data.special_index = special_index; 
data.Tq = Tq; 
data.Yq = Yq; 
data.YqU = YqU; 
data.YqL = YqL; 
data.theta_guess=theta_guess;
data.theta = theta; 
data.total_index = total_index; 
data.error_threshold = error_threshold; 

return







