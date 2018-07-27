%% Figure 2:Altering the contribution of IKs within the same model influences 
%% arrhythmia susceptibility. 
%--- Note: 
% Results displayed in manuscript were run using MATLAB 2016a on a 64bit
% Intel Processor. For exact replication of figures it is best to use these settings.
%% Figures 2A,2B,2C
%--- Description of Figures: 
% Calcium perturbations and calculated ICaL for high iks and low iks
% versions of Ohara model. 

%---: Functions required to run this script :---%
% main_program.m - runs single AP simulation 
%--------------------------------------------------------------------------
%% Set Up Simulation 
settings1.model_name = 'Ohara';
% options - 'Fox', 'Hund', 'Livshitz',
% 'Devenyi','Shannon','TT04','TT06','Grandi','Ohara'

settings1.celltype = 'endo'; % size should be same as model_name, enter one for each model
% options only available for human models as follows
% TT04, TT06, Ohara - 'epi', 'endo', 'mid'
% Grandi - 'epi', 'endo'

settings1.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings1.stim_delay = 100 ; % Time the first stimulus, [ms]
settings1.stim_dur = 2 ; % Stimulus duration
settings1.stim_amp = 32.2; % Stimulus amplitude (see Supplement Table 1 for details) 
nBeats = [92 91 91]; % Number of beats to simulate (see Methods) 
settings1.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings1.steady_state = 1; % 1 - run steady state conditions 0 - do not run steady state conditions 

% Ks, Kr, Ca scaling must have the same length vector. Mainly change
% Ks_scale and Kr_scale when you want to test what different ratios of the
% two look like 
Ks_scale = [0.1 1 10]; % perturb IKs, set to 1 for baseline
Kr_scale = [1.08 1 0.39]; % perturb IKr, set to 1 for baseline
settings1.Ca_scale = [1 7.5 15.5]; % perturb ICaL, set to 1 for baseline

%% Run Simulation
strs = {'low','base','high'};
for i = 1:length(strs)
    settings1.Ks_scale = Ks_scale(i);
    settings1.Kr_scale = Kr_scale(i);
    settings1.nBeats = nBeats(i);
    X = main_program(settings1);
    
    figure
    plot(X.times{1,1},X.V{1,1},'linewidth',2)
    hold on
    plot(X.times{1,2},X.V{1,2},'linewidth',2)
    plot(X.times{1,3},X.V{1,3},'linewidth',2)
    ylim([-100 60])
    title([strs{i} ' IKs Model'])
    legend('ICaL*1','ICaL*7.5','ICaL*15.5')
    xlabel('time (ms)')
    ylabel('Voltage (mV)')
    
    highlow_EAD.(strs{i}) = X;
end

%% Figure 2D
%--- Description of Figures: 
% Calcium perturbations (wide range) and calculated ICaL for high iks and low iks
% models using the Ohara model. 

%% Set Up Simulation 
settings2.model_name = 'Ohara';
% options - 'Fox', 'Hund', 'Livshitz',
% 'Devenyi','Shannon','TT04','TT06','Grandi','Ohara'

settings2.celltype = {'endo'}; % size should be same as model_name, enter one for each model
% options only available for human models as follows
% TT04, TT06, Ohara - 'epi', 'endo', 'mid'
% Grandi - 'epi', 'endo'

settings2.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings2.stim_delay = 100 ; % Time the first stimulus, [ms]
settings2.stim_dur = 2 ; % Stimulus duration
settings2.stim_amp = 32.2;  % Stimulus amplitude (see Supplement Table 1 for details) 
nBeats = [92 91 91] ; % Number of beats to simulate 
settings2.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings2.steady_state = 1;% 1 - run steady state conditions 0 - do not run steady state conditions 

% Ks, Kr, Ca scaling must have the same length vector. Mainly change
% Ks_scale and Kr_scale when you want to test what different ratios of the
% two look like 
Ks_scale = [0.1 1 10]; % perturb IKs, set to 1 for baseline
Kr_scale = [1.08 1 0.39]; % perturb IKr, set to 1 for baseline
settings2.Ca_scale = 1:0.1:22; % perturb ICaL, set to 1 for baseline

%% Run Simulation
strs = {'low','base','high'};
for i = 1:length(Ks_scale)
    settings2.Ks_scale = Ks_scale(i);
    settings2.Kr_scale = Kr_scale(i);
    settings2.nBeats = nBeats(i);
    X = main_program(settings2);
    highlow_multipleEAD.(strs{i}) = X;
end

% clean the data - includes APDs with EADs etc 
for i = 1:3
    indxs = find(highlow_multipleEAD.(strs{i}).APDs>500 );
    highlow_multipleEAD.(strs{i}).APDs(1,indxs) = NaN;
    
    indxs = find(highlow_multipleEAD.(strs{i}).APDs<250);
    highlow_multipleEAD.(strs{i}).APDs(1,indxs) = NaN;
end 

figure
plot(settings2.Ca_scale,highlow_multipleEAD.low.APDs,'linewidth',2,'color','b')
hold on
plot(settings2.Ca_scale,highlow_multipleEAD.base.APDs,'linewidth',2,'color','k')
plot(settings2.Ca_scale,highlow_multipleEAD.high.APDs,'linewidth',2,'color','r')
xlabel('ICaL Factor')
ylabel('APDs (ms)')