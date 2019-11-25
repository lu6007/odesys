%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Name: Locate_peaks.m
% Author: Tongkai Li
% mail: ltk@pku.edu.cn
% Created Time: 2018年07月18日 星期三 09时41分49秒
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[num_crit,locate]= Locate_crit(t,y,tspan)
%Locate_peaks: Locate peaks and wells in the current data.
%	Input:
%		t:		time mesh.
%		y:		variable value.
%		tspan:	time interval to locate in.
%	Output:
%		num_crit:	number of critical points in current data.
%		location:	time points where the peaks appear.
%	PS:
%		This alogrithm locates peaks slightly different from the real peaks.
%		Mainly because that small peaks can be neglect in our case.
%		So we have two creterions:
%			Time between a peak and a well must not be small.
%			Falling distance between a peak and a well must not be small.

	tstart=tspan(1);
	tend=tspan(2);

%Find the time mesh that is in the interval.
	Index=find(t<=tend & t>=tstart);

	t=t(Index);
	y=y(Index);
	Max=max(abs(y));
	
	n=length(t);
	num_crit=0;
    locate=[];
%Collected time that the data remains monotonous.
	time_collect=0;	
%Collected distance that the data remains monotonous.
	dis_collect=0;
%Time threshhold to decide wheather a monotonous interval is legal.
	time_thresh=6;
%Distance threshhold to decide wheather monotonous interval is legal.
	dis_thresh=1e-2*Max;
%Start and end of a monotonous interval.	
	mono_end=1;
	mono_start=1;
%Sign of current interval(1 represent increase, -1 vise versa).
	msign=sign(y(2)-y(1));
%Sign of the origin interval(ignore small pertubations)
	osign=0;

	for i=1:n-1
		change=y(i+1)-y(i);
		if(change*msign>=0&&i<n-1)
			time_collect=time_collect+t(i+1)-t(i);
			dis_collect=dis_collect+change;
			mono_end=i+1;
		else
			if (time_collect>=time_thresh && abs(dis_collect)>=dis_thresh)
				if(-msign==osign)
					num_crit=num_crit+1;
					locate(num_crit,1)=t(mono_start);
					locate(num_crit,2)=y(mono_start);
				end
				osign=msign;
			end
				time_collect=0;
				dis_collect=0;
				mono_start=i;
				mono_end=i;
				msign=-msign;
		end
	 end
