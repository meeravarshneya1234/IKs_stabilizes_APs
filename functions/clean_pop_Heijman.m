function Xnew = clean_pop_Heijman(settings,flags,X)

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
                        %% -- clean_pop_Heijman.m --%%
% Description: Removes APs with EADs or that fail to repolarize from the
% Heijman population. If there are no EADs, it returns the original
% population. If APs with EADs are present, this function reruns new APs to
% replace the ones with EADs. (example: original population is 300, 5 with
% EADs,function reruns 5 APs so output returns population of 300)

% Inputs:
% --> settings - [struct array] AP stimulation protocol (PCL,nBeats,...)
% --> flags - [struct array] indicates PKA targets to block: 0 - block; 1 -
% no block 
% --> X - [struct array] original population of APs that needs to be
% cleaned

% Outputs:
% --> Xnew - [struct array] APs with EADs or failed to repolarize are
% removed (clean population)

%---: Functions used in this script :---%
% ** cleandata.m - Determines index of APs with EADs
%% 
%--- Determine if the population has any EADs ---%
[APfails,nEADs] = cleandata(cell2mat(X.APDs),X.times(:,1),X.V(:,1),settings.t_cutoff,settings.flag);
[~,ind_failed] = find(APfails ==1); % number of failed to repolarize
[~,ind_EADs] = find(nEADs==1); % number of EADs
indexs = [ind_EADs ind_failed];

%--- Remove APs with EADs and rerun the number of APs that were removed  ---%
if isempty(indexs) % no EADs leave function
%     x = X.times(:,1);
%     y = X.V(:,1);
%     
%     figure
%     fig = gcf;
%     
%     figure(fig) % plot original population
%     hold on
%     cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
%     title('Heijman Original')
    
    Xnew = X;
    disp('No EADs in Heijman population.')
else % EADs present
    
    % plot original population
    x = X.times(:,1); y = X.V(:,1);
    figure
    fig = gcf;
    figure(fig)
    subplot(1,3,1)
    hold on
    cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
    title('Heijman Original')
    
    % plot APs with EADs in population
    x = X.times(indexs,1); y = X.V(indexs,1);
    figure(fig)
    subplot(1,3,2)
    hold on
    cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
    title('Heijman APs with EADs')
    
    % new matrix without the APs with EADs
    clean_datatable = [];
    clean_datatable.times = X.times(~(nEADs' + APfails'));
    clean_datatable.V= X.V(~(nEADs' + APfails'));
    clean_datatable.states = X.states(~(nEADs' + APfails'));
    clean_datatable.APDs = X.APDs(~(nEADs' + APfails'));
    clean_datatable.scaling = X.scalings(~(nEADs' + APfails'),:);
    clean_datatable.currents = X.currents(~(nEADs' + APfails'),:);
    
    x = clean_datatable.times(:,1);
    y = clean_datatable.V(:,1);
    
    % plot population without EADs "cleaned"
    figure(fig)
    subplot(1,3,3)
    hold on
    cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
    title('Heijman EADs removed')
    set(fig,'Position',[20,20,600,300])
    
    loop = length(indexs);% number of APs to rerun
    
    while loop > 0 % loop through keep only APs that do not form EADs
        scalings = exp(settings.sigma*randn(18,length(indexs)))' ;
        variations = loop;
        parfor i = 1:variations
            scaling = scalings(i,:);
            [currents,State,Ti,APD]=mainHRdBA(settings,flags,scaling);
            X1(i).times = Ti;
            X1(i).V =  State(:,1);
            X1(i).states =  State;
            X1(i).currents= currents;
            X1(i).APDs = APD;
            X1(i).scalings = scaling;
            
            disp(['Model Variant # ', num2str(i)]);
        end
        Xnew = reformat_data(X1, variations);
        
        n = loop + 1;
        Xnew.times(n:settings.totalvars,1)=clean_datatable.times;
        Xnew.V(n:settings.totalvars,1)=clean_datatable.V;
        Xnew.states(n:settings.totalvars,1)=clean_datatable.states;
        Xnew.APDs(n:settings.totalvars,1)=clean_datatable.APDs;
        Xnew.scalings(n:settings.totalvars,:) = clean_datatable.scaling;
        Xnew.currents(n:settings.totalvars,:) = clean_datatable.currents;

        [APfails,nEADs] = cleandata(cell2mat(Xnew.APDs(:,1)),Xnew.times(:,1),Xnew.V(:,1),settings.t_cutoff,settings.flag);
        [number_of_failed,~] = find(APfails ==1); % number of failed to repolarize
        [number_of_EADs,~] = find(nEADs==1); % number of EADs
        
        clean_datatable.times = Xnew.times(~(nEADs' + APfails'));
        clean_datatable.V = Xnew.V(~(nEADs' + APfails'));
        clean_datatable.states = Xnew.states(~(nEADs' + APfails'));
        clean_datatable.APDs = Xnew.APDs(~(nEADs' + APfails'));
        clean_datatable.scaling = Xnew.scalings(~(nEADs' + APfails'),:);
        clean_datatable.currents = Xnew.currents(~(nEADs' + APfails'),:);

        loop = length(number_of_failed) + length(number_of_EADs);
        
        %check to make sure no duplicates in scalings after rerun
        xnew = cell2mat(Xnew.scalings);
        [u,I,J] = unique(xnew, 'rows', 'first');
        hasDuplicates = size(u,1) < size(xnew,1);
        ixDupRows = setdiff(1:size(xnew,1), I);
        dupRowValues = xnew(ixDupRows,:);
        if ~isempty(dupRowValues)
            disp('Population contains duplicate scalings')
        end
    end
    Xnew = clean_datatable;

    figure
    hold on
    cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
    title('Heijman final')
    disp(['Number of EADs in Heijman population: ' num2str(length(indexs))])
    
end

