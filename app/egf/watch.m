% function data = watch(exp_type)
% exp_type = 'egg1' (default) or 'pdgf1'
%
% Example:
% >> data = watch();
% or, 
% >> exp_type = 'pdgf1';
% >> data = watch(exp_type);

% Author: Tongkai Li (ltk@pku.edu.cn), Shaoying Lu (shaoying.lu@gmail.com)
function data = watch(exp_type)
if nargin == 0 || isempty(exp_type)
    exp_type = 'egf1';
end

% This script is used to show the curve of the ODE or the experiment data.
data = prepare(exp_type); 
% root = data.root; 
gf_con = data.gf_con; 
texp = data.texp; 
option = data.option; 
y0 = data.y0;
weight = data.weight; 
yexp = data.yexp; 
tspan = data.tspan; 
special_index = data.special_index; 
Tq = data.Tq; 
Yq = data.Yq; 
YqU = data.YqU; 
YqL = data.YqL; 
conc_str = data.conc_str;
para_nametag = data.para_nametag; 
line_type = data.line_type;

% load the results of our algorithm: result.mat 
file_name = [data.root, 'result.mat'];
fprintf('watch: loading %s\n', file_name); 
temp = load(file_name); 
field_name = fieldnames(temp);
for k=1:numel(field_name)
    cmd = [field_name{k}, '=temp.', field_name{k}, ';'];
    eval(cmd); clear cmd; 
end
clear temp;

%%% Kathy, we don't seem to need to load the latin_sample results. 
% %
% file_name = [data.root, 'latin_sample.mat'];
% fprintf('watch: loading %s\n', file_name); 
% temp = load(file_name); 
% field_name = fieldnames(temp);
% for k=1:numel(field_name)
%     cmd = [field_name{k}, '=temp.', field_name{k}, ';'];
%     eval(cmd); clear cmd; 
% end
% clear temp; 

% temp_index is the index we want to watch. By default, it is [1:5] in EGF
% simulation and [1:6] in PDGF simulation.
switch data.exp_type
    case {'egf1', 'egf2'}
        temp_index=(1:5);
    case {'pdgf1', 'pdgf2'}
        temp_index = (1:6); 
end
num_group=length(temp_index);

% make a dir to save the results.
path=[data.root, 'output/'];
if(exist(path,'dir')==0)
    mkdir(path);
end

% set up plot index.By default, we choose the Opt_Theta & Error_Theta, which is the
% optimized results under different initial values.
% Other Options listed as follows:
% Init_Theta: Parameter picked initially before the optimizing algorithm.
% Theta_dp and Error_dp: parameter of double peak concentration in simulation.
% Theta_total and Error_total: parameters of total error in simulation.
% Theta_bs and Error_bs: parameters of bell-shape pattern in simulation¡£
% Opt_Thetas and Error_Thetas: parameters optimized in special pattern(group_index).
% Opt_Theta and Error_Theta: parameters optmized after special patter, the
% variables here are volume factors.
% Opt_Theta_t and Error_Theta_t: parameters optimized with respect to all
% concentrations and variables.


% watch_theta=Theta_dp;
% Error_Temp=Error_dp;
watch_theta=Opt_Theta_t;
Error_Temp=Error_Theta_t;
[~,watch_index]=sort(Error_Temp(Error_Temp>0));
% [~,watch_index]=sort(Error_dp(Error_dp>0));
% watch_index=[1:length(watch_theta)];
% tspan=[0;3000];

%show the results one by one, you can also go to one index directly by
%changing the number in the loop.
for i=1:length(watch_index)
    
path=sprintf('%s%d',path,i);
if(exist(path,'dir')==0)
    mkdir(path);
end   
theta=watch_theta(watch_index(i),:);
%If you want to adjust the time or concentration scale manually, change the
%line following.
theta=undimen(theta,1,1);

    ret =0;
    positive=1;
    Error=zeros(1,length(temp_index));
    t = cell(num_group, 1);
    y = cell(num_group, 1);
    yy = cell(num_group, 1);
    index = 7; % watch actin fyn kinase
    sum_weight = sum(weight(temp_index)); 
    for k=temp_index
        theta(end-1)=gf_con(k);
        theta(end-2)=theta(end-2-k);
        [t{k},result]=forward_solver(@rhs_2,texp{k},option,theta,y0);
        % Estimate the error under this certain parameter.
        % y{k}=result(:,8)/theta(end); % active biosensor
        y{k} = result(:,8); % active biosensor
        yy{k} = result(:,index); 
        % As the system is nonlinear and may be stiffness, we should mark those
        % solutions that are not strictly non-negative.
        if(min(min(result))<0)
           fprintf('negative concentration in %s\n',conc_str{k});
           positive=0;
        end
        %%% independent of number of groups
        % Error(k)=weight(k)/sum(weight(temp_index))*num_group*Error_ODE2EXP(y{k},texp{k},yexp{k});
        Error(k)=(weight(k)/sum_weight)*Error_ODE2EXP(y{k},texp{k},yexp{k}'*theta(end));
        ret=ret+Error(k);
%        [num_crit,~]=Locate_crit(t{k},y{k},tspan);
%         if(k==special_index&&num_crit<3)
%             info=0;
%         end
    end
    % if failed to get positive solution, continue to show the next one.
    if(positive==0)       
        continue;
    end
    fprintf('Now we are showing the results of the %d-th parameter combination:\n',watch_index(i));
    fprintf('  L2 error in double peak pattern is %e.\n',Error(special_index));
    fprintf('  Total L2 error is %e.\n',ret);

    savename=sprintf('%spara.txt',path);
    savename_txt=fopen(savename,'w');
    fprintf(savename_txt,'  L2 error in double peak pattern is %e.\n',Error(special_index));
    fprintf(savename_txt,'  Total L2 error is %e.\n',ret);
    fprintf(savename_txt, '  The value of parameters are showed as follows:\n');
    index_nametag=(1:22);
    for si=index_nametag
       fprintf(savename_txt,'%s(Theta%d): %e.\n',para_nametag{si},si,theta(si));
    end

    %plot of the fit-experiment curve, uncomment the line if you do not need lower and
    %upper bounds.    
    for k=temp_index       
        figure(k);
        ki=fix((k-0.01)/3);
        kj=mod(k-1,3);
          set(k,'unit','normalized','Position',[0.3,0.3,0.4,0.4]);
        plot(t{k},y{k},'r');
        hold on
        plot(Tq{k},Yq{k},'b');
        hold on
        plot(Tq{k},3*YqU{k}-2*Yq{k},'b--');
        hold on
        plot(Tq{k},3*YqL{k}-2*Yq{k},'b--');
        hold off
        
        title(conc_str{k});
		xlabel('Time (sec)');
		ylabel('Biosensor Signal');
        legend('model','expermiment','mean+3\sigma','mean-3\sigma');
%         % save corresponding curves by index name.
%         savename=sprintf('%s%d.pdf',path,k);
%         saveas(gcf,savename);
        legend('off');
        set(k,'unit','normalized','Position',[0.01+kj/2.5,1.25-0.6*ki,0.8/2,0.4]/2);
    end
        
    fprintf('The scaling factor concentration/readout is %f.\n',theta(end)); 
    fprintf(savename_txt,'The scaling factor concentration/readout is %f.\n',theta(end));
    for k=temp_index
        fprintf('The volume factor vol/base for species %s is %f.\n',conc_str{k},1/theta(end-2-k)); 
    end
    
    figure(8);
    ki=2;
    kj=0;
    set(8,'unit','normalized','Position',[0.3,0.3,0.4,0.4]);
	for k = temp_index
		plot(t{k},y{k},line_type{k});
        hold on;
    end
    legend(conc_str{temp_index});
    title('showtotal');
    xlabel('Time (sec)');
    ylabel('Biosensor Signal');
    savename=sprintf('%stotal.pdf',path);
    saveas(gcf,savename);
    set(8,'unit','normalized','Position',[0.01+kj/2.5,1.25-0.6*ki,0.8/2,0.4]/2);
    hold off    
    fprintf('Press enter to watch next parameter combination.\n\n');
    fclose(savename_txt);
    
    % draw other variables
    figure(9);
    set(8,'unit','normalized','Position',[0.3,0.3,0.4,0.4]);
	for k = temp_index
		plot(t{k},y{k},line_type{k});
        hold on;
    end
    legend(conc_str{temp_index});
    title('Other Variables');
    xlabel('Time (sec)');
    ylabel('yy [Fyn\_act] (nM)');
    hold off    
    fprintf('Press enter to watch next parameter combination.\n\n');
    
    pause();
end

return