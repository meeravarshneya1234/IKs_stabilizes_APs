%% Figure S5: Baseline analysis of models 

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
% APs + IKr IKs Waveforms for 8 additional models + IKs Fractions 

%---: Functions required to run this part :---%
% main_program.m - runs single AP simulation
%--------------------------------------------------------------------------
%% 
%---- Set Up Simulation Protocol ----%
modelnames = {'Fox','Hund','Shannon','Livshitz','Devenyi','TT04','TT06','Ohara','Grandi'};

celltypes = {'','','','','','endo','endo','endo','endo'}; % size should be same as modelnames, enter one for each model
% options only available for human models as follows:
% TT04, TT06, Ohara -> 'epi', 'endo', or 'mid' 
% Grandi -> 'epi' or 'endo' 
% remaining models -> ''

settings.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
amps = [36.4,32.2,35,35,30.9,25,22.6,32.2,20.6]; % Stimulus amplitude for each model, same order as model names
settings.nBeats = 1 ; % Number of beats to simulate 
settings.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 1;% Run models using the steady state initial conditions 

% Ks and Kr scaling must have the same length vector. 
settings.Ks_scale = 1; % perturb IKs, set to 1 for baseline
settings.Kr_scale = 1; % perturb IKr, set to 1 for baseline
settings.Ca_scale = 1; % perturb ICaL, set to 1 for baseline

%---- Run Simulation ----%
models_to_graph = {'Fox','Hund','Shannon','Livshitz','Devenyi','TT06','Ohara'}; % only plot these 

% for loop going through each model 
for i = 1:length(modelnames)    
    settings.model_name = modelnames{i};
    settings.celltype = celltypes{i};
    settings.stim_amp = amps(i);
    X = main_program(settings);
    
    if sum(strcmp(modelnames{i},models_to_graph)) == 1
        figure
        subplot(1,2,1) % AP plot
        plot(X.times{:},X.V{:},'linewidth',2)
        hold on
        xlabel('time (ms)')
        ylabel('V (mv)')
        set(gcf,'Position',[20,20,600,300])
        
        subplot(1,2,2) % IKs and IKr plots
        plot(X.times{:},X.IKs{:},'color','b','linewidth',2);
        hold on
        plot(X.times{:},X.IKr{:},'color','r','linewidth',2);
        xlabel('time (ms)')
        ylabel('current (A/F)')
        legend('IKs','IKr')
        set(gcf,'Position',[20,20,600,300])
        
        mtit(modelnames{i},...
            'fontsize',14);
    end
    
    X1.(modelnames{i}) = X;
end 

X_Heijman = FigureS5Heijman();
X1.Heijman_ISO_0 = X_Heijman.ISO_0;

models = {'Fox','Hund','Heijman_ISO_0','Livshitz','Devenyi','Shannon','TT04','TT06','Ohara','Grandi'};
figure % IKs Fractions for all models (Supplement 1)
summary_barplot = gcf;
ax_summary = axes('parent', summary_barplot);
IKs_Fractions = cell2mat(cellfun(@(x) X1.(x).IKs_Fraction,models,'UniformOutput',0));
bar(IKs_Fractions,0.5)
set(ax_summary, 'xticklabel',models)
ylabel('IKs Fraction')
xtickangle(90)
set(gcf,'Position',[20,20,600,300])
