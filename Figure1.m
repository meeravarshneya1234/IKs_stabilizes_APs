%% Figure 1: Altering the contribution of IKs within the same model 
%% influences population variability. 

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
%% Figure 1A & 1B
%--- Description of Figure:
% Simulations in high(10x),low(0.1x), baseline IKs models performed in Ohara
% Model

%---: Functions required to run this part :---%
% main_program.m - runs single AP simulation
%--------------------------------------------------------------------------
%%
%---- Set Up High and Low IKs Simulation ----%
settings1.model_name = 'Ohara';
% options -> 'Fox','Hund','Livshitz','Devenyi','Shannon','TT04','TT06','Grandi','Ohara'

settings1.celltype ='endo'; % size should be same as modelnames, enter one for each model
% options only available for human models as follows:
% TT04, TT06, Ohara -> 'epi', 'endo', or 'mid' 
% Grandi -> 'epi' or 'endo' 
% remaining models -> ''

settings1.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings1.stim_delay = 100 ; % Time of first stimulus, [ms]
settings1.stim_dur = 2 ; % Stimulus duration
settings1.stim_amp = 32.2; % Stimulus amplitude (see Supplement Table 1 for details) 
settings1.nBeats = 100; % Number of beats to simulate 
settings1.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings1.steady_state = 1; % 1 - run steady state conditions 0 - do not run steady state conditions 

% Ks, Kr, Ca scaling must have the same length vector. 
settings1.Ks_scale = [0.1 10]; % perturb IKs, set to 1 for baseline
settings1.Kr_scale = [1.08 0.39]; % perturb IKr, set to 1 for baseline
settings1.Ca_scale = 1; % perturb ICaL, set to 1 for baseline

%---- Run Simulation ----%
X1 = main_program(settings1);

%---- Set Up and Run Baseline IKs Simulation ----%
settings1.Ks_scale = 1; % perturb IKs, set to 1 for baseline
settings1.Kr_scale = 1; % perturb IKr, set to 1 for baseline
settings1.nBeats = 1; % Number of beats to simulate 

X1_base = main_program(settings1);

%% Plot Figure 1A - Ohara High, Low APs
figure
plot(X1.times{1,1},X1.V{1,1},'linewidth',2,'color','b') %-low ohara iks model
hold on 
plot(X1_base.times{1,1},X1_base.V{1,1},'linewidth',2,'color','k') %-baseline ohara iks model
plot(X1.times{2,1},X1.V{2,1},'linewidth',2,'color','r') %-high ohara iks model
xlabel('time(ms)')
ylabel('Voltage(mV)')
legend('low iks model','baseline iks model','high iks model')

%% Plot Figure 1B - Ohara High, Low IKs & IKr Waveforms + IKs Fraction
figure
subplot(1,3,1) % plot the baseline ohara iks and ikr waveforms 
plot(X1_base.times{1,1},X1_base.IKs{1,1},'linewidth',2,'color','k')
hold on 
plot(X1_base.times{1,1},X1_base.IKr{1,1},'--k','linewidth',2)
ylim([-0.1 1])
xlim([0 600])
title('Ohara BL Model')
xlabel('time(ms)')
ylabel('Current (A/F)')

subplot(1,3,2) % plot the low ohara iks and ikr waveforms 
plot(X1.times{1,1},X1.IKs{1,1},'linewidth',2,'color','b')
hold on 
plot(X1.times{1,1},X1.IKr{1,1},'--b','linewidth',2)
ylim([-0.1 1])
xlim([0 600])
title('Ohara Low Model')
xlabel('time(ms)')
ylabel('Current (A/F)')

subplot(1,3,3) % plot the low ohara iks and ikr waveforms 
plot(X1.times{2,1},X1.IKs{2,1},'linewidth',2,'color','r')
hold on 
plot(X1.times{2,1},X1.IKr{2,1},'--r','linewidth',2)
ylim([-0.1 1])
xlim([0 600])
title('Ohara High Model')
xlabel('time(ms)')
ylabel('Current (A/F)')

% plot Iks fractions
figure
IKs_Fractions = [X1.IKs_Fraction(1) X1_base.IKs_Fraction X1.IKs_Fraction(2)];
bar(IKs_Fractions,0.5)
xticklabels({'Low','Baseline','High'})
ylabel('IKs Fraction')

% zoomed in on low iks model iks
figure
bar(X1.IKs_Fraction(1))
ylim([0 0.02])
ylabel('low IKs Fraction zoomed in')

%--------------------------------------------------------------------------
%% Figures 1C & 1D 
%--- Description of Figures:
% Population variability and calculated APD Spreads for high iks and low
% iks versions of the Ohara model.

%--- Functions used in this part:
% pop_program.m - runs population simulation 
% histoutline.m - create histograms of distribution (from MATLAB file exchange)
%--------------------------------------------------------------------------
%%
%---- Set Up Population Simulations for High and Low IKs Simulation ----%
settings2.model_name = 'Ohara';
% options -> 'Fox','Hund','Livshitz','Devenyi','Shannon','TT04','TT06','Grandi','Ohara'

settings2.celltype = 'endo'; % size should be same as modelnames, enter one for each model
% options only available for human models as follows:
% TT04, TT06, Ohara -> 'epi', 'endo', or 'mid' 
% Grandi -> 'epi' or 'endo' 
% remaining models -> ''

settings2.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings2.stim_delay = 100 ; % Time the first stimulus, [ms]
settings2.stim_dur = 2 ; % Stimulus duration
settings2.stim_amp = 32.2; % Stimulus amplitude (see Supplement Table 1 for details) 
settings2.nBeats = 100 ; % Number of beats to simulate 
settings2.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings2.steady_state = 1;% 1 - run steady state conditions 0 - do not run steady state conditions 

% Ks, Kr, and Ca scaling must have the same length vector. 
Ks_scale = [0.1 1 10]; % perturb IKs, set to 1 for baseline
Kr_scale = [1.08 1 0.39]; % perturb IKr, set to 1 for baseline
settings2.variations = 300;
settings2.sigma = 0.2;

%---- Run Simulation ----%
for i = 1:length(Ks_scale) % for each ratio of Iks:Ikr
    settings2.Ks_scale = Ks_scale(i);
    settings2.Kr_scale = Kr_scale(i);
    
    disp(['GKs = GKs*' num2str(Ks_scale(i)) '; GKr = GKr*' num2str(Kr_scale(i))])
    
    Xtemp = pop_program(settings2);    
    F = fieldnames(Xtemp);
    for ind = 1:settings2.variations
        for iF = 1:length(F)
            aF = F{iF};
            X2.(aF){ind,i} = Xtemp(ind).(aF);
        end
    end
 
end 

%% Plot Figure 1C - Ohara High, Baseline, Low IKs histogram of populations
figure
pop1 = gcf;
colors = 'bkr';
orig_APDs = cell2mat(X2.APDs);
for i = 1:3    
    % histogram calculations
    APDs = orig_APDs(:,i);
    temp = min(APDs):25:max(APDs);
    bins = linspace(min(APDs),max(APDs),length(temp));
    
    figure(pop1)
    histoutline(APDs,bins,'linewidth',4,'color',colors(i));
    hold on
    set(gca,'FontSize',12,'FontWeight','bold')
    xlabel('APD (ms)','FontSize',12,'FontWeight','bold')
    ylabel('Count','FontSize',12,'FontWeight','bold')
    
    pert1 = prctile(APDs,90);
    pert2 = prctile(APDs,10);
    spreads(i) =(pert1 - pert2)/ median(APDs);
end

%% Plot Figure 1D - Ohara High, Baseline, Low IKs APD Spreads
figure 
summary_barplot = gcf;
ax_summary = axes('parent', summary_barplot);
bar(spreads,0.5)
ylabel('APD Spread')
xticklabels({'Low','Baseline','High'})

for i = 1:3
    figure % plot population: only first 20 of 300
    for ii = 1:300
        plot(X2.times{ii,i},X2.V{ii,i},'linewidth',2)
        hold on
    end
end