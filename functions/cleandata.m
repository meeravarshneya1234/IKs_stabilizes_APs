function [APfails,EADs] = cleandata(APDs,times,volts,t_cutoff,flag)
                            %% -- cleandata.m -- %%
% Description: Finds index of APs with EADs or fail to repolarize. 

% Inputs:
% --> APDs -- population's APDs 
% --> times -- population's time matricies  
% --> volts -- population's voltage matricies  
% --> t_cutoff -- time between EAD cutoff 
% --> flag -- set to 1 if Heijman model 


% Outputs:
% --> APfails - index of APs that failed to repolarize. 
% --> EADs - index of APs that formed EADs. 
%--------------------------------------------------------------------------
%% Run Function 
if ~exist('t_cutoff','var')
    t_cutoff = 5;
end

if ~exist('flag','var')
    flag =0;
end 

APfails = zeros(1,length(times));
EADs = zeros(1,length(times));

for i = 1:length(times)
    
    t = times{i};
    V = volts{i};
    
    % Cell fails if APD is NaN, AP does not begin at rest, AP does not reach full potential 
    if isnan(APDs(i)) || V(1) > -70 || max(V) < 0
        APfails(i) = 1;
    else
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
        slope_threshold = 0.009;
        above_threshold = dVdt > slope_threshold;
        time_between_depols = diff(time(above_threshold));
        if any(time_between_depols > t_cutoff)
            EADs(i) = 1;
        end
    end

end


