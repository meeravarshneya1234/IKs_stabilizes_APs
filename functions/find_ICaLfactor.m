function final_factor = find_ICaLfactor(ICaL_Factor,settings,flags)

total = [];
n = 0.1;

while isempty(total)
    settings.ICaLB = ICaL_Factor;
    [currents,State,Ti,APD]=mainHRdBA(settings,flags);
    [time,volts,APDs] = split_data(Ti,State,settings);
    
    [APfails,nEADs] = cleandata(APDs,time,volts);
    [~,ind_failed] = find(APfails ==1); % number of failed to repolarize
    [~,ind_EADs] = find(nEADs==1); % number of EADs
    total = [ind_EADs ind_failed];
    
    figure 
    plot(Ti,State(:,1),'linewidth',2)
    
    if isempty(total)
        disp(['Ca_scale: ' num2str(ICaL_Factor) ' no EADs.'])
        ICaL_Factor = ICaL_Factor + n;
    else
        beats = (settings.freq - (settings.storeLast - 1)):1:settings.freq;
        disp(['EADs occur at a ICaL Factor of ' num2str(ICaL_Factor) ' on beat ' num2str(beats(total(1))) '.'])
        final_factor = ICaL_Factor;
    end
    close(gcf)
end

% total = [];
% n = 0.1;
% 
% while isempty(total)
%     settings.Ca_scale = ICaL_Factor;
%     [currents,State,Ti,APD]=mainHRdBA(settings,flags);
%     [time,volts,APDs] = split_data(Ti,State,settings);
%     
%     [APfails,nEADs] = cleandata(APDs,time,volts);
%     [~,ind_failed] = find(APfails ==1); % number of failed to repolarize
%     [~,ind_EADs] = find(nEADs==1); % number of EADs
%     total = [ind_EADs ind_failed];
%     
%     if isempty(total)
%         disp(['Ca_scale: ' num2str(ICaL_Factor) ' no EADs.'])
%         ICaL_Factor = ICaL_Factor + n;
%     else
%         disp(['Ca_scale: ' num2str(ICaL_Factor) ' EADs.'])
%         final_factor = ICaL_Factor;
%     end
% end



