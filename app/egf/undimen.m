%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Name: undimen.m
% Author: Tongkai Li
% mail: ltk@pku.edu.cn
% Created Time: 2018�?8�?1�?星期�?19�?8�?6�?%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[result]=undimen(theta,time,conc)
%	undimen: used to apply undimensionalization.
%	Input:
%		theta:	original theta;
%		time:	time multiple number.
%		conc:	conc multiple number.
%	Output:
%		result:	results after undimensionalization.


%     	result=theta;
% 	for i=[1:6,8,10,11,13,15]
% 		result(i)=result(i)*time;
%     end                                                                                                                                                                                                                                                                                                                                                                                                              
% 	for i=[6:9,11,12,13,14,17:18]
% 		result(i)=result(i)/conc;
% 	end
%     
%     
%     result(23)=result(23)/time;
    	result=theta;
    for i=[1,2,3,4,5,6,8,9,10,11,13,15,17,19,21]
		result(i)=result(i)*time;
    end
    
%     for i=[15,19:22]
% 		result(i)=result(i)/conc;
%     end
%     for i=[2]
%         result(i)=result(i)*conc; 
%     end
    
    
