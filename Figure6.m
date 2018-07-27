%% Figure 6:AP clamp simulations illustrate changes in K+ currents with 
%% alterations in AP shape. 
%--- Note: 
% Results displayed in manuscript were run using MATLAB 2016a on a 64bit
% Intel Processor. For exact replication of figures it is best to use these settings.

%--- Description of Figure: 
% AP Clamp simulations performed in the O'Hara model

%---: Functions required to run this script :---%
% main_program.m - runs single AP simulation 
% main_APClamp.m - runs the AP clamp simulation using parfor loop 
% plotting_APClamp.m - plots the waveforms from the main_APClamp.m function 

%% Run under baseline conditions to get input voltage and time for AP Clamp
settings.model_name = 'Ohara';
settings.celltype = 'endo';

settings.PCL = 1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
settings.stim_amp = 32.2; % Stimulus amplitude
settings.nBeats = 1 ; % Number of beats to simulate
settings.numbertokeep = 1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 1;

% Ks and Kr scaling must have the same length vector. Mainly change
% Ks_scale and Kr_scale when you want to test what different ratios of the
% two look like
settings.Ks_scale = 1; % perturb IKs, set to 1 for baseline
settings.Kr_scale = 1; % perturb IKr, set to 1 for baseline
settings.Ca_scale =1; % perturb ICaL, set to 1 for baseline

base_datatable= main_program(settings);
statevars_input = base_datatable.states{1,1}; % used as input for the AP clamp 
t_input = base_datatable.times{1,1}; % used as input for the AP clamp 

%% Figure 6A & 6B - Run AP Clamp for multiple APDs - Linear Scale
apclamp_settings.model_name = 'Ohara';
apclamp_settings.celltype = 'endo';
[~,Vind] = ICs(apclamp_settings.model_name,1,1000);
statevar_i_original = statevars_input(1:end ~= Vind); % all state variables besides V
apclamp_settings.statevar_i_original = cellfun(@(v) v(end), statevar_i_original);
apclamp_settings.volts = statevars_input{1,Vind};
apclamp_settings.times = t_input;
apclamp_settings.repol_change =  [0.5,1,1.5,2,2.5]; % scaling
apclamp_settings.numberofAPs = 100;
apclamp_settings.numberkeep = 1;
apclamp_settings.PCL = 1000;
apclamp_settings.vals_graph = [0.5,1,1.5,2,2.5];

APClamp_data = main_APClamp(apclamp_settings); % run simulation 
indxs = find(ismember(roundx(apclamp_settings.repol_change,2),roundx(apclamp_settings.vals_graph,2)));
data_to_graph = APClamp_data(:,indxs);
plotting_APClamp(apclamp_settings,data_to_graph) % plot results 

%% Figure 6C - Run AP Clamp for multiple APDs - Log Scale
apclamp_settings.model_name = 'Ohara';
apclamp_settings.celltype = 'endo';
[~,Vind] = ICs(apclamp_settings.model_name,1,1000);
statevar_i_original = statevars_input(1:end ~= Vind); % all state variables besides V
apclamp_settings.statevar_i_original = cellfun(@(v) v(end), statevar_i_original);
apclamp_settings.volts = statevars_input{1,Vind};
apclamp_settings.times = t_input;
logscalefactors = linspace(log(1/3),log(3),11) ;
apclamp_settings.repol_change =  exp(logscalefactors) ; % scaling
apclamp_settings.numberofAPs = 100;
apclamp_settings.numberkeep = 1;
apclamp_settings.PCL = 1000;

log_APClamp_data = main_APClamp(apclamp_settings);

BL = find(roundx(apclamp_settings.repol_change,2) ==1.00);
AKr = cell2mat(arrayfun(@(x) x.Area_Kr(1), log_APClamp_data,'UniformOutput', 0));
AKs = cell2mat(arrayfun(@(x) x.Area_Ks(1), log_APClamp_data,'UniformOutput', 0));
Volts = cell2mat(arrayfun(@(x) x.Voltage(:,1), log_APClamp_data,'UniformOutput', 0));
times = cell2mat(arrayfun(@(x) x.time(:,1), log_APClamp_data,'UniformOutput', 0));
IKr = cell2mat(arrayfun(@(x) x.IKr(:,1), log_APClamp_data,'UniformOutput', 0));
IKs = cell2mat(arrayfun(@(x) x.IKs(:,1), log_APClamp_data,'UniformOutput', 0));

figure
loglog(apclamp_settings.repol_change,AKr/AKr(BL),'k-x', 'linewidth',2)
hold on
loglog(apclamp_settings.repol_change,AKs/AKs(BL),'k-o','linewidth',2)
xlabel('APD Factor')
ylabel('norm AUC')
legend('IKr','IKs')
set(gca,'FontSize',12,'FontWeight','bold')

figure
bar([AKr(end)/AKr(BL),AKs(end)/AKs(BL)],0.5)
ylabel('norm current (A/F)')
set(gca,'XTickLabels',{'IKr','IKs'})

