%% Set Up Simulation Protocol 
modelnames = {'Fox','Hund','Shannon','Livshitz','Devenyi','TT06','Ohara'};
% options - 'Fox', 'Hund', 'Livshitz',
% 'Devenyi','Shannon','TT04','TT06','Grandi','Ohara'

celltypes = {'','','','','','endo','endo','endo','endo'}; % size should be same as model_name, enter one for each model
% options only available for human models - 'epi', 'endo', 'mid'

settings.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
amps = [36.4,32.2,35,35,30.9,22.6,32.2]; % Stimulus amplitude for each model, same order as model names
settings.nBeats = 100 ; % Number of beats to simulate 
settings.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 1;% Run models using the steady state initial conditions 

% Ks, Kr, Ca scaling must have the same length vector. Mainly change
% Ks_scale and Kr_scale when you want to test what different ratios of the
% two look like 
settings.Ks_scale = 1; % perturb IKs, set to 1 for baseline
settings.Kr_scale = 1; % perturb IKr, set to 1 for baseline
settings.Ca_scale = 1; % perturb ICaL, set to 1 for baseline

settings.variations = 300;
settings.totalvars = settings.variations;
settings.sigma = 0.2;
t_cutoffs = [5 7 5 5 5 5 5];
settings.flag = 0;

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
    normAPDs = Xnew.APDs - median(Xnew.APDs);
    temp = min(normAPDs):25:max(normAPDs);
    bins = linspace(min(normAPDs),max(normAPDs),length(temp));
    medians(i) = median(Xnew.APDs);
    
    pert1 = prctile(Xnew.APDs,90);
    pert2 = prctile(Xnew.APDs,10);
    Xnew.APDSpread =(pert1 - pert2)/ median(Xnew.APDs);
    
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

run FigureS2Heijman.m 
X1.Heijman = X_Heijman; 

all_names = {'Fox','Hund','Heijman','Livshitz','Devenyi','Shannon','TT06','Ohara'};
figure 
summary_barplot = gcf;
ax_summary = axes('parent', summary_barplot);
Spread = cell2mat(cellfun(@(x) X1.(x).APDSpread,all_names,'UniformOutput',0));
bar(Spread,0.5)
set(ax_summary, 'xticklabel',all_names)
ylabel('APD Spread')

% Note in mauscript the APD Spread bar plot includes TT04 and Grandi APD 
% Spreads which were taken from data collected in Figure3.m file. 
