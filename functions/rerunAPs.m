function X = rerunAPs(settings,clean_datatable)

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
                            %% -- rerunAPs.m --%%
% Description: Run new APs to replace ones removed with EADs.

% Inputs:
% --> settings - [struct array] AP stimulation protocol (PCL,nBeats,...)
% --> clean_datatable - [struct array] cleaned population of APs 

% Outputs:
% --> X - [struct array] final population 

%---: Functions required to run this function :---%
% ** pop_program.m - runs population for each model using parloop
% ** reformat_data.m - restructures output from the pop_program.m 
% ** cleandata.m - Determines index of APs with EADs 
%% 
%--- Run Population ---%
EADs = 1;
while EADs > 0 % keep running population until all APs have no EADs
    X = pop_program(settings); % run population simulation 
    X = reformat_data(X, settings.variations);% reformat data that is in the parallel format
    
    n = settings.variations+1;
    X.times(n:settings.totalvars,1)=clean_datatable.times;
    X.V(n:settings.totalvars,1)=clean_datatable.V;
    X.states(n:settings.totalvars,1)=clean_datatable.states;
    X.APDs(n:settings.totalvars,1)=clean_datatable.APDs;
    X.scalings(n:settings.totalvars,:) = clean_datatable.scaling;
    
    [APfails,nEADs] = cleandata(cell2mat(X.APDs(:,1)),X.times(:,1),X.V(:,1),settings.t_cutoff,settings.flag);
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

%--- Check to make sure no duplicate model variants after rerun ---%
x = cell2mat(X.scalings);
[u,I,J] = unique(x, 'rows', 'first');
hasDuplicates = size(u,1) < size(x,1);
ixDupRows = setdiff(1:size(x,1), I);
dupRowValues = x(ixDupRows,:);
if ~isempty(dupRowValues)
    disp('Population contains duplicate scalings')
end 