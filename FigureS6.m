%% Figure S6: Population Variability across multiple species.

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
%--- Description of Figure: 
% Build heterogenous populations of APs in multiple species models by 
% imposeing random variation in parameters controlling ionic current levels. 

%---: Functions required to run this part :---%
% pop_program.m - runs population simulation using parfor loop
% reformat_data.m - reformats data collected during parallel loop into
% easier format
% clean_pop.m - inspects for AP population for EADs and returns population with no EADs 
%--------------------------------------------------------------------------
%% 
%---- Set Up Simulation Protocol ----%
modelnames = {'Fox','Hund','Shannon','Livshitz','Devenyi','TT06'};

celltypes = {'','','','','','endo','endo','endo'}; % size should be same as modelnames, enter one for each model
% options only available for human models as follows:
% TT04, TT06, Ohara -> 'epi', 'endo', or 'mid' 
% Grandi -> 'epi' or 'endo' 
% remaining models -> ''

settings.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
amps = [36.4,32.2,35,35,30.9,22.6]; % Stimulus amplitude for each model, same order as model names (See Table S1)
settings.nBeats = 100 ; % Number of beats to simulate 
settings.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 1;% Run models using the steady state initial conditions 

% Ks, Kr, Ca scaling must have the same length vector. Mainly change
% Ks_scale and Kr_scale when you want to test what different ratios of the
% two look like 
settings.Ks_scale = 1; % perturb IKs, set to 1 for baseline
settings.Kr_scale = 1; % perturb IKr, set to 1 for baseline
settings.Ca_scale = 1; % perturb ICaL, set to 1 for baseline

settings.variations = 300; % number of variants in the population 
settings.sigma = 0.2; % standard deviation assigned for population 

% settings for clean_pop.m function 
settings.totalvars = settings.variations;
t_cutoffs = [5 7 5 5 5 5]; % time between EAD cutoff, order based on order for modelnames
settings.flag = 0; % set to 1 only if the model has notch like Heijman 

%% Run Simulation 

for i = 1:length(modelnames)
    settings.model_name = modelnames{i};
    settings.celltype = celltypes{i};
    settings.stim_amp = amps(i);
    X = pop_program(settings);
    X = reformat_data(X, settings.variations);
    
    settings.t_cutoff = t_cutoffs(i);
    Xnew = clean_pop(settings,X);
       
    % histogram calculations
    APDs = cell2mat(Xnew.APDs);
    normAPDs = APDs - median(APDs);
    temp = min(normAPDs):25:max(normAPDs);
    bins = linspace(min(normAPDs),max(normAPDs),length(temp));
    medians(i) = median(APDs);
    
    pert1 = prctile(APDs,90);
    pert2 = prctile(APDs,10);
    Xnew.APDSpread =(pert1 - pert2)/ median(APDs);
    
    figure
    histoutline(normAPDs,bins,'linewidth',4);
    hold on
    title(modelnames(i))
    set(gca,'FontSize',12,'FontWeight','bold')
    xlabel('APD (ms)','FontSize',12,'FontWeight','bold')
    ylabel('Count','FontSize',12,'FontWeight','bold')
    xlim([-300 300])
    ylim([0 180])
    
    X1.(modelnames{i}) = Xnew;
    Xnew = []; X =[];% reset for next model 
end

X1.Heijman  = FigureS6Heijman();

all_names = {'Fox','Hund','Heijman','Livshitz','Devenyi','Shannon','TT06'};
figure 
summary_barplot = gcf;
ax_summary = axes('parent', summary_barplot);
Spread = cell2mat(cellfun(@(x) X1.(x).APDSpread,all_names,'UniformOutput',0));
bar(Spread,0.5)
set(ax_summary, 'xticklabel',all_names)
ylabel('APD Spread')

% Note in mauscript the APD Spread bar plot includes TT04 and Grandi APD 
% Spreads which were taken from data collected in Figure3.m file and Ohara
% APD Spread taken from Figure1.m 
