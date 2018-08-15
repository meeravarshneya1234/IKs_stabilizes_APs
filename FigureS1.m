%% Figure S1: Altering the contribution of IKs within the same model 
%% influences population variability - TT04 and Grandi 

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
%% Figure S1A, S1B, S1C
%--- Description of Figure:
% Simulations in low(0.05x), baseline IKs models performed in TT04 Model

%---: Functions required to run this part :---%
% main_program.m - runs single AP simulation
%--------------------------------------------------------------------------

%---- Set Up Low IKs Simulation ----%
settings1.model_name = 'TT04';
settings1.celltype = 'endo';
settings1.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings1.stim_delay = 100 ; % Time the first stimulus, [ms]
settings1.stim_dur = 2 ; % Stimulus duration
settings1.stim_amp = 25; % Stimulus amplitude
settings1.nBeats = 100 ; % Number of beats to simulate 
settings1.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings1.steady_state = 1;

% Ks, Kr, Ca scaling must have the same length vector. 
settings1.Ks_scale = 0.05; % perturb IKs, set to 1 for baseline
settings1.Kr_scale = 3.13; % perturb IKr, set to 1 for baseline
settings1.Ca_scale = 1; % perturb ICaL, set to 1 for baseline

%---- Run Simulation ----%
X_TT04_highlow = main_program(settings1);

%---- Set Up and Run Baseline IKs Simulation ----%
settings1.Ks_scale = 1; % perturb IKs, set to 1 for baseline
settings1.Kr_scale = 1; % perturb IKr, set to 1 for baseline
settings1.nBeats = 1; % Number of beats to simulate 
X_TT04_base = main_program(settings1);

%---- Plot Figures ----%
figure
plot(X_TT04_base.times{1,1},X_TT04_base.V{1,1},'linewidth',2,'color','k') %baseline tt04 model
hold on 
plot(X_TT04_highlow.times{1,1},X_TT04_highlow.V{1,1},'linewidth',2,'color','b') %-low tt04 iks model
xlabel('time(ms)')
ylabel('Voltage(mV)')
title('Figure S1A - TT04')
legend('baseline model', 'low iks model')
xlim([0 600])

figure
subplot(1,2,1) % plot the baseline tt04 iks and ikr waveforms 
plot(X_TT04_base.times{1,1},X_TT04_base.IKs{1,1},'linewidth',2,'color','k')
hold on 
plot(X_TT04_base.times{1,1},X_TT04_base.IKr{1,1},'--k','linewidth',2)
ylim([-0.1 1.6])
xlim([0 600])
title('TT04 BL Model')
xlabel('time(ms)')
ylabel('Current (A/F)')
legend('IKs','IKr')

subplot(1,2,2) % plot the low tt04 iks and ikr waveforms 
plot(X_TT04_highlow.times{1,1},X_TT04_highlow.IKs{1,1},'linewidth',2,'color','b')
hold on 
plot(X_TT04_highlow.times{1,1},X_TT04_highlow.IKr{1,1},'--b','linewidth',2)
ylim([-0.1 1.6])
xlim([0 600])
title('TT04 Low Model')
xlabel('time(ms)')
ylabel('Current (A/F)')
legend('IKs','IKr')

mtit('Figure S1B' ,...
	     'fontsize',14);

% plot Iks fractions 
figure
IKs_Fractions = [X_TT04_highlow.IKs_Fraction X_TT04_base.IKs_Fraction];
bar(IKs_Fractions,0.5)
xticklabels({'Low','Baseline'})
ylabel('IKs Fraction')
title('Figure S1C - TT04')

%% Figure S1D & S1E
%--- Description of Figure:
% Simulations in low(0.05x), baseline IKs models performed in TT04 Model

%--- Description of Figure: 
% Build heterogenous populations of APs in TT04 low and baseline models

%---: Functions required to run this part :---%
% pop_program.m - runs population simulation using parfor loop
% reformat_data.m - reformats data collected during parallel loop into
% easier format
% clean_pop.m - inspects for AP population for EADs and returns population with no EADs 
%--------------------------------------------------------------------------

%---- Set Up Low and Base IKs TT04 Population Simulation ----%
settings2.model_name = 'TT04';
settings2.celltype = 'endo'; 
settings2.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings2.stim_delay = 100 ; % Time the first stimulus, [ms]
settings2.stim_dur = 2 ; % Stimulus duration
settings2.stim_amp = 25;
settings2.nBeats = 100 ; % Number of beats to simulate 
settings2.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings2.steady_state = 1;

% Ks and Kr scaling must have the same length vector. 
Ks_scale = [1 0.05]; % perturb IKs, set to 1 for baseline
Kr_scale = [1 3.13]; % perturb IKr, set to 1 for baseline
settings2.Ca_scale = 1; % perturb ICaL, set to 1 for baseline
settings2.variations = 300;
settings2.sigma = 0.2;
settings2.flag = 0;
settings2.totalvars = settings2.variations;

%---- Run Simulation ----%
strs = {'base','low'};
for i = 1:2    
    settings2.Ks_scale = Ks_scale(i);
    settings2.Kr_scale = Kr_scale(i);
    X = pop_program(settings2);
    X = reformat_data(X, settings2.variations);
    
    settings2.t_cutoff = 5;
    TT04_Pop.(strs{i}) = clean_pop(settings2,X);
   
end

figure  % plot the baseline tt04 histogram 
pop1 = gcf;
APDs = cell2mat(TT04_Pop.base.APDs);
temp = min(APDs):25:max(APDs);
bins = linspace(min(APDs),max(APDs),length(temp));
histoutline(APDs,bins,'linewidth',4,'color','k');
hold on

pert1 = prctile(APDs,90);
pert2 = prctile(APDs,10);
spreads(1) =(pert1 - pert2 )/ median(APDs); % baseline TT04 APD spread  

% histogram calculations
APDs = cell2mat(TT04_Pop.low.APDs);
temp = min(APDs):25:max(APDs);
bins = linspace(min(APDs),max(APDs),length(temp));

figure(pop1)
histoutline(APDs,bins,'linewidth',4,'color','b');
hold on
set(gca,'FontSize',12,'FontWeight','bold')
xlabel('APD (ms)','FontSize',12,'FontWeight','bold')
ylabel('Count','FontSize',12,'FontWeight','bold')
legend('Baseline', 'Low')

pert1 = prctile(APDs,90);
pert2 = prctile(APDs,10);
spreads(2) =(pert1 - pert2 )/ median(APDs); % low tt04 APD spread
title('Figure S1D  - TT04')

figure % plot all APD Spreads 
bar(spreads,0.5)
ylabel('APD Spread')
xticklabels({'Baseline','Low'})
title('Figure S1E - TT04')

%% Figure S1F, S1G, S1H
%--- Description of Figure:
% Simulations in high (20x), baseline IKs models performed in Grandi Model

%---: Functions required to run this part :---%
% main_program.m - runs single AP simulation
%--------------------------------------------------------------------------

%---- Set Up High IKs Simulation ----%
settings3.model_name = 'Grandi';
settings3.celltype = 'endo'; 
settings3.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings3.stim_delay = 100 ; % Time the first stimulus, [ms]
settings3.stim_dur = 2 ; % Stimulus duration
settings3.stim_amp  = 20.6; % Stimulus amplitude
settings3.nBeats = 100 ; % Number of beats to simulate 
settings3.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings3.steady_state = 1;

% Ks, Kr, Ca scaling must have the same length vector. 
settings3.Ks_scale = 20; % perturb IKs, set to 1 for baseline
settings3.Kr_scale = 0.04; % perturb IKr, set to 1 for baseline
settings3.Ca_scale = 1; % perturb ICaL, set to 1 for baseline

%---- Run Simulation ----%
X_Grandi_highlow = main_program(settings3);

%---- Set Up and Run Baseline IKs Simulation ----%
settings3.Ks_scale = 1; % perturb IKs, set to 1 for baseline
settings3.Kr_scale = 1; % perturb IKr, set to 1 for baseline
settings3.nBeats = 1; % Number of beats to simulate 

X_Grandi_base = main_program(settings3);

figure
plot(X_Grandi_base.times{1,1},X_Grandi_base.V{1,1},'linewidth',2,'color','k') %baseline grandi model
hold on 
plot(X_Grandi_highlow.times{1,1},X_Grandi_highlow.V{1,1},'linewidth',2,'color','r') %-high grandi iks model
xlabel('time(ms)')
ylabel('Voltage(mV)')
title('Grandi - Figure S1F')
legend('baseline model', 'high iks model')
xlim([0 600])

figure
subplot(1,2,1) % plot the baseline grandi iks and ikr waveforms 
plot(X_Grandi_base.times{1,1},X_Grandi_base.IKs{1,1},'linewidth',2,'color','k')
hold on 
plot(X_Grandi_base.times{1,1},X_Grandi_base.IKr{1,1},'--k','linewidth',2)
ylim([-0.1 0.3])
xlim([0 600])
title('Grandi BL Model')
xlabel('time(ms)')
ylabel('Current (A/F)')
legend('IKs','IKr')

subplot(1,2,2) % plot the high grandi iks and ikr waveforms 
plot(X_Grandi_highlow.times{1,1},X_Grandi_highlow.IKs{1,1},'linewidth',2,'color','r')
hold on 
plot(X_Grandi_highlow.times{1,1},X_Grandi_highlow.IKr{1,1},'--r','linewidth',2)
ylim([-0.1 0.3])
xlim([0 600])
title('Grandi High Model')
xlabel('time(ms)')
ylabel('Current (A/F)')
legend('IKs','IKr')

mtit('Figure S1G' ,...
	     'fontsize',14);

% plot Iks fractions 
figure
IKs_Fractions = [X_Grandi_highlow.IKs_Fraction X_Grandi_base.IKs_Fraction];
bar(IKs_Fractions,0.5)
xticklabels({'High','Baseline'})
ylabel('IKs Fraction')
title('Figure S1H - Grandi')

%% Figure S1I, S1J
%--- Description of Figure:
% Simulations in high (20x), baseline IKs models performed in Grandi Model

%--- Description of Figure: 
% Build heterogenous populations of APs in Grandi high and baseline models

%---: Functions required to run this part :---%
% pop_program.m - runs population simulation using parfor loop
% reformat_data.m - reformats data collected during parallel loop into
% easier format
% clean_pop.m - inspects for AP population for EADs and returns population with no EADs 
%--------------------------------------------------------------------------

%---- Set Up High and Base IKs Grandi Population Simulation ----%
settings4.model_name = 'Grandi';
settings4.celltype = 'endo'; 
settings4.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings4.stim_delay = 100 ; % Time the first stimulus, [ms]
settings4.stim_dur = 2 ; % Stimulus duration
settings4.stim_amp = 20.6;
settings4.nBeats = 100 ; % Number of beats to simulate 
settings4.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings4.steady_state = 1;

% Ks and Kr scaling must have the same length vector. 
Ks_scale = [1 20]; % perturb IKs, set to 1 for baseline
Kr_scale = [1 0.03]; % perturb IKr, set to 1 for baseline
settings4.Ca_scale = 1; % perturb ICaL, set to 1 for baseline
settings4.variations = 300;
settings4.sigma = 0.2;
settings4.flag = 0;
settings4.t_cutoff = 5;
settings4.totalvars = settings4.variations;

%---- Run Simulation ----%
strs = {'base','high'};
for i = 1:2    
    settings4.Ks_scale = Ks_scale(i);
    settings4.Kr_scale = Kr_scale(i);
    X = pop_program(settings4);
    X = reformat_data(X, settings4.variations);
    
    Grandi_Pop.(strs{i}) = clean_pop(settings4,X);
end

figure  % plot the baseline grandi histogram 
pop1 = gcf;
APDs = cell2mat(Grandi_Pop.base.APDs);
temp = min(APDs):25:max(APDs);
bins = linspace(min(APDs),max(APDs),length(temp));
histoutline(APDs,bins,'linewidth',4,'color','k');
hold on

pert1 = prctile(APDs,90);
pert2 = prctile(APDs,10);
spreads(1) =(pert1 - pert2 )/ median(APDs); % baseline grandi APD spread  

% histogram calculations
APDs = cell2mat(Grandi_Pop.high.APDs);
temp = min(APDs):25:max(APDs);
bins = linspace(min(APDs),max(APDs),length(temp));

figure(pop1)
histoutline(APDs,bins,'linewidth',4,'color','r');
hold on
set(gca,'FontSize',12,'FontWeight','bold')
xlabel('APD (ms)','FontSize',12,'FontWeight','bold')
ylabel('Count','FontSize',12,'FontWeight','bold')
legend('Baseline', 'High')

pert1 = prctile(APDs,90);
pert2 = prctile(APDs,10);
spreads(2) =(pert1 - pert2 )/ median(APDs); % high grandi APD spread
title('Figure S1I - Grandi')


figure % plot all APD Spreads 
summary_barplot = gcf;
ax_summary = axes('parent', summary_barplot);
bar(spreads,0.5)
ylabel('APD Spread')
xticklabels({'Baseline','High'})
title('Figure S1J - Grandi')
