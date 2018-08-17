%% Figure 4:Susceptibility to Early Afterdepolarizations (EADs) in Grandi 
%% TT04, and Ohara models. 

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
%% Figure 4
%--- Description of Figures: 
% Calcium perturbations and calculated ICaL for Grandi, Ohara, TT04

%---: Functions required to run this script :---%
% main_program.m - runs simulation 
% mtit - create main title (from MATLAB file exchange) 
%--------------------------------------------------------------------------
%%
%---- Set Up Simulation ----%
modelnames = {'TT04','Grandi','Ohara'};
% options -> 'Fox','Hund','Livshitz','Devenyi','Shannon','TT04','TT06','Grandi','Ohara'

celltypes = {'endo','endo','endo'}; % size should be same as modelnames, enter one for each model
% options only available for human models as follows:
% TT04, TT06, Ohara -> 'epi', 'endo', or 'mid' 
% Grandi -> 'epi' or 'endo' 
% remaining models -> ''

settings.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
amps = [25,20.6,32.2]; % Stimulus amplitude (see Supplement Table 1 for details) 
nBeats = [91 91 91]; %number of beats to stimulate first EAD 
settings.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 1;% 1 - run steady state conditions 0 - do not run steady state conditions 

% Ks, Kr, Ca scaling must have the same length vector. 
settings.Ks_scale = 1; % perturb IKs, set to 1 for baseline
settings.Kr_scale = 1; % perturb IKr, set to 1 for baseline
%settings.Ca_scale = ; % perturb ICaL, set to 1 for baseline
ICa_scales = { [1,14.4,27.7,27.8]
    [1,1.3,1.5,1.6]
    [1 8 15.1 15.2]}; %calcium perturbation for each model based on order of modelnames

%---- Run Simulation ----%
%% Plot Figure 4A & 4B
for i = 1:length(modelnames) 
    settings.model_name = modelnames{i};
    settings.celltype = celltypes{i};
    settings.stim_amp = amps(i);
    settings.Ca_scale = ICa_scales{i};
    settings.nBeats = nBeats(i);
    X = main_program(settings);     
    
    figure % Figure 4A - APs 
    colors = hsv(length(ICa_scales{i}));
    ca_scale = ICa_scales{i};
    
    for ii = 1:length(ca_scale) % Figure 4B - ICaL levels at baseline & before EAD
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

ICaL= [X1.Grandi.Area_Ca(1) X1.TT04.Area_Ca(1) X1.Ohara.Area_Ca(1);
X1.Grandi.Area_Ca(3) X1.TT04.Area_Ca(3) X1.Ohara.Area_Ca(3)]';

%% Plot Figure 4C
figure
handle = gcf;
ax = axes('parent', handle);
bar(abs(ICaL),'stacked')
set(ax, 'xticklabel',({'Grandi','TT04','Ohara'}))
legend('baseline ICaL','ICaL below threshold')
