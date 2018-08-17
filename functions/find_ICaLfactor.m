function final_factor = find_ICaLfactor(ICaL_Factor,settings,flags,sims)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--- "Slow delayed rectifier current protects ventricular myocytes from
% arrhythmic dynamics across multiple species: a computational study" ---%

% By: Varshneya,Devenyi,Sobie 
% For questions, please contact Dr.Eric A Sobie -> eric.sobie@mssm.edu 
% or put in a pull request or open an issue on the github repository:
% https://github.com/meeravarshneya1234/IKs_stabilizes_APs.git. 

%--- Note:
% Results displayed in manuscript were run using MATLAB 2016a on a 64bit
% Intel Processor. For exact replication of figures it is best to use these
% settings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
                            %% -- find_ICaLfactor.m -- %%
% Description: Determines the factor required to induce an EAD (increase in
% ICaL current) in different variations of the Heijman model (i.e. no ISO,
% ISO, block of PKA phosphorylation of IKs)

% Inputs:
% --> ICaL_factor - [double] value to begin increase of ICaL to induce EADs
% --> settings - [struct array] AP stimulation protocol (PCL,nBeats,...)
% --> flags - [struct array] block or no block of PKA targets, 0 - block, 1 - no block  

% Outputs:
% --> final_factor - [double] first ICaL factor that induces EADs 

%---: Functions required to run this function :---%
% ** mainHRdBA.m - run Heijman model 
% ** split_data.m - splits data based on number of beats saved
% ** cleandata.m - determines index of APs with EADs 
%--------------------------------------------------------------------------
%%

total = [];
n = 0.1; 
while isempty(total)
    settings.ICaLB = ICaL_Factor;
    [currents,State,Ti,APD]=mainHRdBA(settings,flags); % run simulation 
    [time,volts,APDs] = split_data(Ti,State,settings); % split data based on number of beats kept 
    
    [APfails,nEADs] = cleandata(APDs,time,volts);
    [~,ind_failed] = find(APfails ==1); % number of failed to repolarize
    [~,ind_EADs] = find(nEADs==1); % number of EADs
    total = [ind_EADs ind_failed];
    
    figure 
    plot(Ti,State(:,1),'linewidth',2)
    
    if isempty(total) % did not find factor - keep going 
        disp([sims ' - Ca_scale: ' num2str(ICaL_Factor) ' no EADs.'])
        ICaL_Factor = ICaL_Factor + n;    
        close(gcf)

    else % found factor, leave function now 
        beats = (settings.freq - (settings.storeLast - 1)):1:settings.freq;
        disp([sims ' - EADs occur at a ICaL Factor of ' num2str(ICaL_Factor) ' on beat ' num2str(beats(total(1))) '.'])
        final_factor = ICaL_Factor;
        title([sims ' : ICaL Factor - ' num2str(final_factor)])    
    end
    
end

