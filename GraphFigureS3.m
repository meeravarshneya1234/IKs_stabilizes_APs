settings.model_name = 'Ohara';
% options - 'Fox', 'Hund', 'Livshitz',
% 'Devenyi','Shannon','TT04','TT06','Grandi','Ohara'
% 'TT06_opt','Grandi_opt','Ohara_opt'

settings.celltype = 'endo'; % size should be same as model_name, enter one for each model
% options only available for human models as follows
% TT04, TT06, Ohara - 'epi', 'endo', 'mid'
% Grandi - 'epi', 'endo'

settings.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
settings.stim_amp = 32.2;
settings.nBeats = 100;
settings.numbertokeep = 10;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 1;

% Ks and Kr scaling must have the same length vector. Mainly change
% Ks_scale and Kr_scale when you want to test what different ratios of the
% two look like 
settings.Ks_scale = 1; % perturb IKs, set to 1 for baseline
settings.Kr_scale = 1; % perturb IKr, set to 1 for baseline
settings.Ca_scale = 15.2;% perturb ICaL, set to 1 for baseline

%% Run Simulation
datatable= main_program(settings);

%% Graph Figure
load('ohara_ead_10APs.mat')
figure 
plot(datatable.Ohara.times{1},datatable.Ohara.V{1},'linewidth',2,'color','b')
