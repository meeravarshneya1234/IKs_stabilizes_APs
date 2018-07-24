function Xnew = clean_pop(settings,X)


%% Determine if the population has any EADs.
[APfails,nEADs] = cleandata(X.APDs(:,1),X.times(:,1),X.V(:,1),settings.t_cutoff,settings.flag);
[~,ind_failed] = find(APfails ==1); % number of failed to repolarize
[~,ind_EADs] = find(nEADs==1); % number of EADs
indexs = [ind_EADs ind_failed];

if isempty(indexs) % no EADs leave function 
    x = X.times(:,1);
    y = X.V(:,1);
    
    figure
    fig = gcf;
    
    figure(fig) % plot original population
    hold on
    cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
    title([settings.model_name ' Original'])

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

