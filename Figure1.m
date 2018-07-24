%% Figure 1: Altering the contribution of IKs within the same model influences population variability. 
%--- Note: 
% Results displayed in manuscript were run using MATLAB 2016a on a 64bit
% Intel Processor. For exact replication of figures it is best to use these settings.

%% Figure 1A 
%--- Description of Figure: 
% Simulations in high(10x) and low(0.1x) IKs models performed in Ohara Model  

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

%% Set Up Simulation 
settings1.model_name = 'Ohara';
% options - 'Fox', 'Hund', 'Livshitz',
% 'Devenyi','Shannon','TT04','TT06','Grandi','Ohara'

settings1.celltype ='endo'; % size should be same as modelnames, enter one for each model
% options only available for human models as follows
% TT04, TT06, Ohara - 'epi', 'endo', 'mid'
% Grandi - 'epi', 'endo'

settings1.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings1.stim_delay = 100 ; % Time of first stimulus, [ms]
settings1.stim_dur = 2 ; % Stimulus duration
settings1.stim_amp = 32.2; % Stimulus amplitude (see Supplement Table 1 for details) 
settings1.nBeats = 100 ; % Number of beats to simulate 
settings1.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings1.steady_state = 1; % 1 - run steady state conditions 0 - do not run steady state conditions 

% Ks, Kr, Ca scaling must have the same length vector. Mainly change
% Ks_scale and Kr_scale when you want to test what different ratios of the
% two look like 
settings1.Ks_scale = [0.1 10]; % perturb IKs, set to 1 for baseline
settings1.Kr_scale = [1.08 0.39]; % perturb IKr, set to 1 for baseline
settings1.Ca_scale = 1; % perturb ICaL, set to 1 for baseline

%% Run Simulation
X1 = main_program(settings1);
figure
plot(X1.times{1,1},X1.V{1,1},'linewidth',2,'color','b') %-low ohara iks model
hold on 
plot(X1.times{2,1},X1.V{2,1},'linewidth',2,'color','r') %-high phara iks model
xlabel('time(ms)')
ylabel('Voltage(mV)')
legend('low iks model','high iks model')
%note the code to run normal Ohara AP is located in FigureS1.m 

%% Figures 1C & 1D 
%--- Description of Figures: 
% Population variability and calculated APD Spreads for high iks and low iks
% versions of the Ohara model. 

%--- Functions used in this script:
% pop_program.m - runs simulation 
% * find_APD.m - determines APD when AP returns to -75 mV
% * parameters.m - lists the model parameters 
% * ICs - returns initial conditions - either steady state or ones from paper
%   * inital_conditions.m - lists the model initial conditions from paper
%   if user does not want to use steady state 
% * currents.m - calculates ICaL, IKs, IKr 
% histoutline.m - create histograms of distribution (from MATLAB file exchange) 

%% Set Up Simulation 
settings2.model_name = 'Ohara';
% options - 'Fox', 'Hund', 'Livshitz',
% 'Devenyi','Shannon','TT04','TT06','Grandi','Ohara'

settings2.celltype = 'endo'; % size should be same as model_name, enter one for each model
% options only available for human models - 'epi', 'endo', 'mid',
% TT04, TT06, Ohara - 'epi', 'endo', 'mid'
% Grandi - 'epi', 'endo'

settings2.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings2.stim_delay = 100 ; % Time the first stimulus, [ms]
settings2.stim_dur = 2 ; % Stimulus duration
settings2.stim_amp = 32.2; % Stimulus amplitude (see Supplement Table 1 for details) 
settings2.nBeats = 100 ; % Number of beats to simulate 
settings2.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings2.steady_state = 1;% 1 - run steady state conditions 0 - do not run steady state conditions 

% Ks, Kr, and Ca scaling must have the same length vector. Mainly change
% Ks_scale and Kr_scale when you want to test what different ratios of the
% two look like 
settings2.Ks_scale = [0.1 10]; % perturb IKs, set to 1 for baseline
settings2.Kr_scale = [1.08 0.39]; % perturb IKr, set to 1 for baseline
settings2.Ca_scale =  1; % perturb ICaL, set to 1 for baseline
settings2.variations = 300;
settings2.sigma = 0.2;

%% Run Simulation 
X2 = pop_program(settings2);

figure
pop1 = gcf;
colors = 'br';
for i = 1:2    
    % histogram calculations 
    APDs = X2.APDs(:,i);
    temp = min(APDs):25:max(APDs);
    bins = linspace(min(APDs),max(APDs),length(temp));
    
    figure(pop1)
    histoutline(newAPDs,bins,'linewidth',4,'color',colors(i));
    hold on
    set(gca,'FontSize',12,'FontWeight','bold')
    xlabel('APD (ms)','FontSize',12,'FontWeight','bold')
    ylabel('Count','FontSize',12,'FontWeight','bold')
    
    pert1 = prctile(APDs,90);
    pert2 = prctile(APDs,10);
    spreads(i) =(pert1 - pert2 )/ median(APDs);
end

figure 
summary_barplot = gcf;
ax_summary = axes('parent', summary_barplot);
bar(spreads,0.5)
ylabel('APD Spread')
xticklabels({'Low','High'})
%note the code to run normal Ohara AP is located in FigureS2.m 
