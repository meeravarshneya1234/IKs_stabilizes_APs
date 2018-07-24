%% Figure 4:Susceptibility to Early Afterdepolarizations (EADs) in Grandi and TT04 models. 
%--- Note: 
% Results displayed in manuscript were run using MATLAB 2016a on a 64bit
% Intel Processor. For exact replication of figures it is best to use these settings.

%--- Description of Figures: 
% Calcium perturbations and calculated ICaL for Grandi and TT04

%--- Functions used in this script:
% main_program.m - runs simulation 
% * find_APD.m - determines APD when AP returns to -75 mV
% * parameters.m - lists the model parameters 
% * ICs - returns initial conditions - either steady state or ones from
% paper
%   * inital_conditions.m - lists the model initial conditions from paper
%   if user does not want to use steady state 
% * currents.m - calculates ICaL, IKs, IKr 
% dydt functions for each model 
% plotting.m - graphs the figures
% * mtit - create main title (from MATLAB file exchange) 

%% Set Up Simulation Protocol 
modelnames = {'TT04','Grandi'};
% options - 'Fox', 'Hund', 'Livshitz',
% 'Devenyi','Shannon','TT04','TT06','Grandi','Ohara'

celltypes = {'endo','endo'}; % size should be same as model_name, enter one for each model
% options only available for human models as follows
% TT04, TT06, Ohara - 'epi', 'endo', 'mid'
% Grandi - 'epi', 'endo'

settings.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
amps = [25,20.6]; %Stimulus amplitude for each model, same order as model names
nBeats = [91 91]; %number of beats to stimulate first EAD 

settings.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 1;
% Ks, Kr, Ca scaling must have the same length vector. Mainly change
% Ks_scale and Kr_scale when you want to test what different ratios of the
% two look like 
settings.Ks_scale = 1; % perturb IKs, set to 1 for baseline
settings.Kr_scale = 1; % perturb IKr, set to 1 for baseline
%settings.Ca_scale = ; % perturb ICaL, set to 1 for baseline
ICa_scales = { [1,14.4,27.7,27.8]
    [1,1.3,1.5,1.6]}; %calcium perturbation for each model based on order of modelnames

%% Run Simulation
for i = 1:length(modelnames) 
    settings.model_name = modelnames{i};
    settings.celltype = celltypes{i};
    settings.stim_amp = amps(i);
    settings.Ca_scale = ICa_scales{i};
    settings.nBeats = nBeats(i);
    X = main_program(settings);     
    
    figure
    colors = hsv(length(ICa_scales{i}));
    ca_scale = ICa_scales{i};
    
    for ii = 1:length(ca_scale)
        subplot(1,2,1)
        plot(X.times{ii},X.V{ii},'color',colors(ii,:),'linewidth',2);
        hold on
        xlabel('time (ms)')
        ylabel('V (mv)')
        set(gcf,'Position',[20,20,600,300])
        ylim([-100 100])
        
        if ii == 1 || ii == 3
            subplot(1,2,2)
            plot(X.times{ii},X.ICaL{ii},'color',colors(ii,:),'linewidth',2);
            hold on
            xlabel('time (ms)')
            ylabel('ICaL (A/F)')
            set(gcf,'Position',[20,20,600,300])
        end
        
    end
    %legend(figurestrings);
    mtit(modelnames{i} ,...
        'fontsize',14);
    
    X1.(modelnames{i}) = X;  
end 

ICaL= [X1.Grandi.Area_Ca(1) X1.TT04.Area_Ca(1);
X1.Grandi.Area_Ca(3) X1.TT04.Area_Ca(3)]';

figure
handle = gcf;
ax = axes('parent', handle);
bar(abs(ICaL),'stacked')
set(ax, 'xticklabel',({'Grandi','TT04'}))
legend('baseline ICaL','ICaL below threshold')
