function X = rerunAPs(settings,clean_datatable)
EADs = 1;

%% Run the Population 
while EADs > 0
    X = pop_program(settings);
    X = reformat_data(X, settings.variations);
    
    n = settings.variations+1;
    X.times(n:settings.totalvars,1)=clean_datatable.times;
    X.V(n:settings.totalvars,1)=clean_datatable.V;
    X.states(n:settings.totalvars,1)=clean_datatable.states;
    X.APDs(n:settings.totalvars,1)=clean_datatable.APDs;
    X.scalings(n:settings.totalvars,:) = clean_datatable.scaling;
    
    [APfails,nEADs] = cleandata(X.APDs(:,1),X.times(:,1),X.V(:,1),settings.t_cutoff,settings.flag);
    [number_of_failed,~] = find(APfails ==1); % number of failed to repolarize
    [number_of_EADs,~] = find(nEADs==1); % number of EADs
    
    clean_datatable.times = X.times(~(nEADs' + APfails'));
    clean_datatable.V = X.V(~(nEADs' + APfails'));
    clean_datatable.states = X.states(~(nEADs' + APfails'));
    clean_datatable.APDs = X.APDs(~(nEADs' + APfails'));
    clean_datatable.scaling = X.scalings(~(nEADs' + APfails'),:);
    
    EADs = length(number_of_failed) + length(number_of_EADs);
    settings.variations = EADs;
end

%% check to make sure no duplicates in scalings after rerun 
x = cell2mat(X.scalings);
[u,I,J] = unique(x, 'rows', 'first');
hasDuplicates = size(u,1) < size(x,1);
ixDupRows = setdiff(1:size(x,1), I);
dupRowValues = x(ixDupRows,:);
if ~isempty(dupRowValues)
    disp('Population contains duplicate scalings')
end 