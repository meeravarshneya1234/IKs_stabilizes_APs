function [ic,Vind] = ICs(model_name,steady_state,pcl)
if steady_state
    load('ICs_SS.mat')
    if sum(pcl == [500,1000,5000]) > 0
        str =['PCL_' num2str(pcl)];
        ic = ICs.(model_name).(str);
    else 
        disp(['Initial conditions for steady state PCL = ' num2str(pcl) ' do not exist. Using the non-steady state inital conditions.'])
        ic = initial_conditions(model_name);
    end 
end 

% Determine index of voltage in the state variables. 
if strcmp(model_name,'Shannon')
    Vind = 32;
elseif strcmp(model_name,'Grandi') || strcmp(model_name,'Morotti')
    Vind = 39;
else 
    Vind = 1;
end 

