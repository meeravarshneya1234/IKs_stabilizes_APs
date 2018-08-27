%% Figure S1:Procedure to determine EAD Threshold

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
%% Figure S1
%--- Description of Figures: 
% AP waveforms of beats 91-100 of the O’Hara Human Model under calcium
% perturbation at 1 Hz pacing rate. EAD appearance was not consistent at
% every beat. Thus, in order to address this variability, beats 91-100 were
% analyzed and the first ICaL Factor that caused an EAD between beats 
% 91-100 was considered the threshold for triggering EADs in the model.
% This procedure was repeated for each model. The O’Hara Model displayed an
% EAD first on the 91st beat at an ICaL Factor of 15.2. This analysis was
% replicated in each model.

%---: Functions required to run this script :---%
% main_program.m - runs simulation 
%--------------------------------------------------------------------------
%%
%---- Set Up Simulation ----%
settings.model_name = 'Ohara'; %model name 
settings.celltype = 'endo'; % cell type 'endo', 'epi', or 'mid'

settings.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
settings.stim_amp = 32.2;% Stimulus amplitude for each model
settings.nBeats = 100; % Number of beats to simulate 
settings.numbertokeep = 10;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 1;% Run models using the steady state initial conditions 

% Ks, Kr, Ca scaling must have the same length vector. 
settings.Ks_scale = 1; % perturb IKs, set to 1 for baseline
settings.Kr_scale = 1; % perturb IKr, set to 1 for baseline
settings.Ca_scale = 15.2;% perturb ICaL, set to 1 for baseline

%---- Run Simulation ----%
datatable= main_program(settings);

%---- Graph Figure ----%
figure 
plot(datatable.times{1},datatable.V{1},'linewidth',2,'color','b')
xlabel('Time (ms)')
ylabel('Voltage (mV)')
set(gca,'FontSize',12,'FontWeight','bold')
title('Figure S5')

