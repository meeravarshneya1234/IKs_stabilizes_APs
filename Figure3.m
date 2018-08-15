%% Figure 3: Dramatic differences seen between three human myocyte models in 
%% K+ currents and population variability.

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
%% Figure 3A & 3B
%--- Description of Figure: 
% APs + IKr and IKs Waveforms for human models - Grandi, Ohara, and TT04 
% Models + IKs Fraction Barplot. All panels with Ohara data taken from the 
% code in Figure 1. 

%---: Functions required to run this script :---%
% main_program.m - runs single AP simulation 
%--------------------------------------------------------------------------
%%
%---- Set Up High and Low IKs Simulation ----%
modelnames = {'TT04','Grandi'};
% options -> 'Fox','Hund','Livshitz','Devenyi','Shannon','TT04','TT06','Grandi','Ohara'

celltypes = {'endo','endo'}; % size should be same as modelnames, enter one for each model
% options only available for human models as follows:
% TT04, TT06, Ohara -> 'epi', 'endo', or 'mid' 
% Grandi -> 'epi' or 'endo' 
% remaining models -> ''

settings1.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings1.stim_delay = 100 ; % Time the first stimulus, [ms]
settings1.stim_dur = 2 ; % Stimulus duration
amps = [25,20.6]; % Stimulus amplitude for each model, same order as model names
settings1.nBeats = 1 ; % Number of beats to simulate 
settings1.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings1.steady_state = 1;% Run models using the steady state initial conditions 

% Ks, Kr, Ca scaling must have the same length vector. 
settings1.Ks_scale = 1; % perturb IKs, set to 1 for baseline
settings1.Kr_scale = 1; % perturb IKr, set to 1 for baseline
settings1.Ca_scale = 1; % perturb ICaL, set to 1 for baseline

%---- Run Simulation ----%
%% Plot Figure 3A - Baseline AP plots, IKs and IKr waveforms of each model 
for i = 1:length(modelnames)% for loop going through each model 
    settings1.model_name = modelnames{i};
    settings1.celltype = celltypes{i};
    settings1.stim_amp = amps(i);
    X = main_program(settings1);
    
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
    
    title_var = strjoin({'Figure 3A', modelnames{i}});
    
    mtit(title_var,...
        'fontsize',14);
    
    str = modelnames{i};
    X1.(str) = X;
end

%% Figure 3B - IKs Fractions 
figure
barplot = gcf;
ax_summary = axes('parent', barplot);
IKs_Fractions = [X1.Grandi.IKs_Fraction X1.TT04.IKs_Fraction];
bar(IKs_Fractions,0.5)
set(ax_summary, 'xticklabel',{'Grandi','TT04'})
ylabel('IKs Fraction')
title('Figure 3B')
xtickangle(90)

%--------------------------------------------------------------------------
%% Figures 3C & 3D 
%--- Description of Figures: 
% Population variability and calculated APD Spreads for Grandi, Ohara, TT04.
% All panels with Ohara data taken from the code in Figure 1. 

%--- Functions used in this script:
% pop_program.m - runs simulation 
% histoutline.m - create histograms of distribution (from MATLAB file exchange) 
%--------------------------------------------------------------------------
%%
%---- Set Up Population Simulations ----%
modelnames = {'Grandi','TT04'};
% options -> 'Fox','Hund','Livshitz','Devenyi','Shannon','TT04','TT06','Grandi','Ohara'

celltypes = {'endo','endo'}; % size should be same as modelnames, enter one for each model
% options only available for human models as follows:
% TT04, TT06, Ohara -> 'epi', 'endo', or 'mid' 
% Grandi -> 'epi' or 'endo' 
% remaining models -> ''

settings2.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings2.stim_delay = 100 ; % Time the first stimulus, [ms]
settings2.stim_dur = 2 ; % Stimulus duration
amps = [20.6,25]; % Stimulus amplitude 
settings2.nBeats = 100 ; % Number of beats to simulate 
settings2.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings2.steady_state = 1;% Run models using the steady state initial conditions 

% Ks, Kr, Ca scaling must have the same length vector. 
settings2.Ks_scale = 1; % perturb IKs, set to 1 for baseline
settings2.Kr_scale = 1; % perturb IKr, set to 1 for baseline
settings2.Ca_scale = 1; % perturb ICaL, set to 1 for baseline

settings2.variations = 300;
settings2.totalvars = settings2.variations;
settings2.sigma = 0.2;
settings2.t_cutoff = 3;
settings2.flag = 0;

%---- Run Simulation ----%
%% Plot Figure 3c & 3D - Populations and Histograms of populations
figure % figure for histograms of each population 
histos = gcf;
colors = 'rb';

for i = 1:length(modelnames)
    settings2.model_name = modelnames{i};
    settings2.celltype = celltypes{i};
    settings2.stim_amp = amps(i);
    Xtemp = pop_program(settings2);
    X = reformat_data(Xtemp,settings2.variations);
    %X.APDs = cell2mat(X.APDs);
    
    Xnew = clean_pop(settings2,X);
    
    figure % plot population: only first 20 of 300 
    for ii = 1:20
        plot(Xnew.times{ii},Xnew.V{ii},'linewidth',2)
        hold on
    end
    title_var = strjoin({'Figure 3C', modelnames{i}});
    title(title_var)
    xlabel('time (ms)')
    ylabel('voltage (mV)')
    ylim([-100 60])
    
    % histogram calculations
    APDs = cell2mat(Xnew.APDs);
    temp = min(APDs):25:max(APDs);
    bins = linspace(min(APDs),max(APDs),length(temp));
    
    figure(histos)
    histoutline(APDs,bins,colors(i),'linewidth',4);
    hold on
    
    X2.(modelnames{i}) = Xnew;
    Xnew = []; X =[];% reset for next model 
end

figure(histos) 
set(gca,'FontSize',12,'FontWeight','bold')
xlabel('APD (ms)','FontSize',12,'FontWeight','bold')
ylabel('Count','FontSize',12,'FontWeight','bold')
hold off
box('off')
legend(modelnames)
title('Figure 3D')
