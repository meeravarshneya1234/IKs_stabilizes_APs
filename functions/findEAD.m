% function EAD = findEAD(t,V,t_cutoff,flag)
% EAD = 0;
% if ~exist('t_cutoff','var')
%     t_cutoff = 5;
% end 
% 
% if ~exist('flag','var')
%     flag =0;
% end 
% 
% 
% [~,t_start] = findpeaks(V);
% if flag ==1
%     t_start = t(t_start(2));
% else
%     t_start = t(t_start(1));
% end
% t_end = find(t>t_start & V<-75,1);
% 
% time = t((find(t==t_start(1))):t_end(1));
% Voltage = V((find(t==t_start(1))):t_end(1));
% dVdt = diff(Voltage)./diff(time);
% slope_threshold = 0.009; % note changed this from 0.1 to get rid of EADs in the Grandi model 
% %t_cutoff = 5 ;
% above_threshold = dVdt > slope_threshold;
% time_between_depols = diff(time(above_threshold));
% if any(time_between_depols > t_cutoff)
%     EAD = 1;
% end
% 
% end

function EAD = findEAD(t,V,t_cutoff,flag)
EAD = 0;
if ~exist('t_cutoff','var')
    t_cutoff = 5;
end 

if ~exist('flag','var')
    flag =0;
end 


[~,t_start] = findpeaks(V);
if flag ==1
    t_start = t(t_start(2));
else
    t_start = t(t_start(1));
end
t_end = find(t>t_start & V<-75,1);

time = t((find(t==t_start(1))):t_end(1));
Voltage = V((find(t==t_start(1))):t_end(1));
dVdt = diff(Voltage)./diff(time);
slope_threshold = 0.009; % note changed this from 0.1 to get rid of EADs in the Grandi model 
%t_cutoff = 5 ;
above_threshold = dVdt > slope_threshold;
% figure
% plot(t,V)
time_between_depols = diff(time(above_threshold));
% nonEADs = find(time_between_depols > 3);
% if ~isempty(nonEADs)
%     time_between_depols = time_between_depols(nonEADs(end)+1:end);
% end

if any(time_between_depols > t_cutoff)
    %if (time_between_depols(1) > time_between_depols(end))
    figure
    plot(t,V)
    hold on
    plot(time,Voltage)
    plot(time(above_threshold),Voltage(above_threshold),'go')
    inds = find(time_between_depols > t_cutoff);
    temp = find(above_threshold == 1);
    plot(time(temp(inds)),Voltage(temp(inds)),'mo')
    EAD = 1;
    %end
end

close(gcf)
end

% Old 
% dVdt = diff(V)./diff(t);
% slope_threshold = 0.1; %mV/ms, rate of change considered to be a significant depolarization or repolarization
% t_cutoff = 10 ; %ms, the time between depolarization periods above which it is considered two separate depolarizations
% V_slope = diff(V)./diff(t);
% dices_above = V_slope > slope_threshold;
% time_between_depols = diff(t(dices_above));
% 
% if sum(time_between_depols > t_cutoff) > 1
%     EAD = 1;
% end