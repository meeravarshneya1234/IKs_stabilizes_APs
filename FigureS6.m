%% Figure S6:Susceptibility to EADs across multiple species.

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
%% Figure S6
%--- Description of Figures: 
% Calcium perturbations and calculated ICaL for all models. 

%--- Functions used in this script:
% main_program.m - runs simulation 
% mtit - create main title (from MATLAB file exchange) 
%--------------------------------------------------------------------------
%%
%---- Set Up Simulation ----%
modelnames = {'Fox','Hund','Shannon','Livshitz','Devenyi','TT04','TT06','Ohara','Grandi'};
% options - 'Fox', 'Hund', 'Livshitz','Devenyi','Shannon','TT04','TT06','Grandi','Ohara'

celltypes = {'','','','','','endo','endo','endo','endo'}; % size should be same as model_name, enter one for each model
%celltypes = {''}; % size should be same as model_name, enter one for each model
% options only available for human models as follows
% TT04, TT06, Ohara - 'epi', 'endo', 'mid'
% Grandi - 'epi', 'endo'

settings.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
amps = [36.4,32.2,35,35,30.9,25,22.6,32.2,20.6]; %Stimulus amplitude for each model, same order as model names
nBeats = [100 99 91 93 92 91 91 91 91]; %number of beats to stimulate first EAD , see Figure S5 for details
settings.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 1;% 1 - run steady state conditions 0 - do not run steady state conditions 

% Ks and Kr scaling must have the same length vector.
settings.Ks_scale = 1; % perturb IKs, set to 1 for baseline
settings.Kr_scale = 1; % perturb IKr, set to 1 for baseline
ICa_scales = {[1 20 38.9 39] %calcium perturbation for each model based on order of modelnames
    [1,34.4,67.7,67.8]
    [1,2,2.9,3]
    [1,29.6,58,58.1]
    [1,2,2.9,3]
    [1,14.4,27.7,27.8]
    [1,11.7,22.3,22.4]
    [1,8.1,15.1,15.2]
    [1,1.3,1.5,1.6]};

%---- Run Simulation ----%
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
    
    mtit(modelnames{i} ,...
        'fontsize',14);
    
    X1.(modelnames{i}) = X;  

end 

X1.Heijman = FigureS6Heijman();

models = {'Fox','Hund','Heijman','Livshitz','Devenyi','Shannon','TT04','TT06','Ohara','Grandi'};
figure % plot Baseline AUC of ICaL for all models 
summaryBL_barplot = gcf;
axBL_summary = axes('parent', summaryBL_barplot);
ICaL_Baseline = abs(cell2mat(cellfun(@(x) X1.(x).Area_Ca(1),models,'UniformOutput',0)));
bar(ICaL_Baseline,0.5)
set(axBL_summary, 'xticklabel',models)
ylabel('Baseline ICaL')
xtickangle(90)
set(gca,'FontSize',12,'FontWeight','bold')

figure % plot Before First EAD AUC of ICaL for all models 
summary_barplot = gcf;
ax_summary = axes('parent', summary_barplot);
ICaL_Before_EAD = abs(cell2mat(cellfun(@(x) X1.(x).Area_Ca(3),models,'UniformOutput',0)));
bar(ICaL_Before_EAD,0.5)
set(ax_summary, 'xticklabel',models)
ylabel('ICaL Before EAD')
xtickangle(90)
set(gca,'FontSize',12,'FontWeight','bold')


