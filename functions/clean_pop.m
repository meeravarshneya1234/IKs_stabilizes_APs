function Xnew = clean_pop(settings,X)

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

                            %% -- clean_pop.m --%%
% Description: removes APs with EADs or that fail to repolarize from a 
% population of APs. If there are no EADs, it returns the original population. 

% Inputs:
% --> settings - [struct array] AP stimulation protocol (PCL,nBeats,...)
% --> X - [struct array] Original population of APs that needs to be cleaned 

% Outputs:
% --> Xnew - [struct array] APs with EADs or failed to repolarize are removed (clean
% population) 

%---: Functions required to run this function :---%
% ** cleandata.m - Determines index of APs with EADs 
% ** rerunAPs.m - If APs with EADs are present, this function reruns new
% APs to replace the ones with EADs. (example: original population is 300,
% 5 with EADs,function reruns 5 APs so output returns population of 300) 

%% 
%--- Determine if the population has any EADs ---%
[APfails,nEADs] = cleandata(cell2mat(X.APDs(:,1)),X.times(:,1),X.V(:,1),settings.t_cutoff,settings.flag); 
[~,ind_failed] = find(APfails ==1); % number of failed to repolarize
[~,ind_EADs] = find(nEADs==1); % number of EADs
indexs = [ind_EADs ind_failed]; %total 

%--- Remove APs with EADs and rerun the number of APs that were removed  ---%
if isempty(indexs) % no EADs leave function 
%%% uncomment the following section to graph the population    
%     x = X.times(:,1);
%     y = X.V(:,1);
    
%     figure
%     fig = gcf;
%     
%     figure(fig) % plot original population
%     hold on
%     cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
%     title([settings.model_name ' Original'])
    Xnew = X;
    disp(['No EADs in ' settings.model_name ' population.'])
else % EADs present 
    
    % plot original population
    x = X.times(:,1); y = X.V(:,1);
    figure
    fig = gcf;
    
    figure(fig) 
    subplot(1,3,1)
    hold on
    cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
    title([settings.model_name ' Original'])
    
    % plot APs with EADs in population 
    x = X.times(indexs,1);
    y = X.V(indexs,1);

    figure(fig)
    subplot(1,3,2)
    hold on
    cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
    title([settings.model_name ' APs with EADs'])
    
    % new matrix without the APs with EADs 
    clean_datatable = [];
    clean_datatable.times = X.times(~(nEADs' + APfails'));
    clean_datatable.V= X.V(~(nEADs' + APfails'));
    clean_datatable.states = X.states(~(nEADs' + APfails'));
    clean_datatable.APDs = X.APDs(~(nEADs' + APfails'));
    clean_datatable.scaling = X.scalings(~(nEADs' + APfails'),:);
    
    x = clean_datatable.times(:,1);
    y = clean_datatable.V(:,1);
    
    % plot population without EADs "cleaned"
    figure(fig) 
    subplot(1,3,3)
    hold on
    cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
    title([settings.model_name ' EADs removed'])
    
    set(fig,'Position',[20,20,600,300])
    
    EADs = settings.variations-length(clean_datatable.APDs); % number of APs with EADs 

    % rerun population so no EADs 
    settings.variations = EADs;
    Xnew = rerunAPs(settings,clean_datatable);
    x = Xnew.times(:,1);
    y = Xnew.V(:,1);
    
    figure
    hold on
    cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
    title([settings.model_name ' final'])
    disp(['Number of EADs in ' settings.model_name ' population: ' num2str(EADs)])

end

