%% Figure S10: AP Clamp Simulations across Models

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
%% Figure S10
%--- Description of Figure:
% Runs AP clamps simulations in all models as done in Figure6.m 

%---: Functions required to run this part :---%
% main_program.m - runs baseline simulation using parfor loop
% main_APClamp.m - runs ap clamp simulation using parfor loop
% reformat_data.m - reformats the data collected from inject_current_program
% FigureS10Heijman.m - runs AP Clamp simulation for Heijman model 
%--------------------------------------------------------------------------
%%
%---- Set Up Simulation ----%
modelnames = {'Fox','Hund','Shannon','Livshitz','Devenyi','TT04','TT06','Grandi'};
celltypes = {'','','','','','endo','endo','endo','endo'}; % size should be same as model_name, enter one for each model
amps = [36.4,32.2,35,35,30.9,25,22.6,32.2,20.6 ]; % Stimulus amplitude 
n = [11 11 11 11 11 10 10 10]; % number of APD factors to plot because TT04 and TT06 cannot be increased by 3

% graph the log scaled IKs vs IKr 
figure 
fig1 = gcf;
% graph the levels of IKs vs IKr at the last scaling 
figure 
fig2 = gcf;

for index = 1:length(modelnames)
    %% Step 1 - Run under baseline conditions to get input voltage and time for AP Clamp    
    settings.model_name = modelnames{index};
    settings.celltype = celltypes{index};
    
    settings.PCL = 1000 ;  % Interval bewteen stimuli,[ms]
    settings.stim_delay = 100 ; % Time the first stimulus, [ms]
    settings.stim_dur = 2 ; % Stimulus duration
    settings.stim_amp = amps(index); % Stimulus amplitude
    settings.nBeats = 1 ; % Number of beats to simulate
    settings.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
    settings.steady_state = 1;
    
    % Ks, Kr scaling must have the same length vector. Mainly change
    % Ks_scale and Kr_scale when you want to test what different ratios of the
    % two look like
    settings.Ks_scale = 1; % perturb IKs, set to 1 for baseline
    settings.Kr_scale = 1; % perturb IKr, set to 1 for baseline
    settings.Ca_scale =1; % perturb ICaL, set to 1 for baseline
    
    datatable= main_program(settings);
    statevars_input = datatable.states{1,1};
    t_input = datatable.times{1,1};

    %% Run AP Clamp for multiple APDs - Log Scale
    apclamp_settings.model_name = modelnames{index};
    apclamp_settings.celltype = celltypes{index};     
    [~,Vind] = ICs(apclamp_settings.model_name,1,1000);
    statevar_i_original = statevars_input(1:end ~= Vind); % all state variables besides V
    apclamp_settings.statevar_i_original = cellfun(@(v) v(end), statevar_i_original);
    apclamp_settings.volts = statevars_input{1,Vind};
    apclamp_settings.times = t_input;
    logscalefactors = linspace(log(1/3),log(3),11) ;
    apclamp_settings.repol_change =  exp(logscalefactors(1:n(index))); % scaling
    apclamp_settings.numberofAPs = 100;
    apclamp_settings.numberkeep = 1;
    apclamp_settings.PCL = 1000;
    
    log_APClamp_data = main_APClamp(apclamp_settings);
    
    BL = find(round(apclamp_settings.repol_change,2) ==1.00);
    AKr = cell2mat(arrayfun(@(x) x.Area_Kr(1), log_APClamp_data,'UniformOutput', 0));
    AKs = cell2mat(arrayfun(@(x) x.Area_Ks(1), log_APClamp_data,'UniformOutput', 0));
    
    figure(fig1)
    subplot(3,3,index+1)
    loglog(apclamp_settings.repol_change(1:n(index)),AKr(1:n(index))/AKr(BL),'k-x', 'linewidth',2)
    hold on
    loglog(apclamp_settings.repol_change(1:n(index)),AKs(1:n(index))/AKs(BL),'k-o','linewidth',2)
    set(gca,'FontSize',12,'FontWeight','bold')
    title(modelnames{index})
    xlim([0,3])
    ylim([0.001 100])
    
    
    figure(fig2)
    subplot(3,3,index+1)
    bar([AKr(n(index))/AKr(BL),AKs(n(index))/AKs(BL)],0.5)
    set(gca,'XTickLabels',{'IKr','IKs'})
    title(modelnames{index})      
end

%---- Run Simulation in Heijman Model ----%
Heijman_APClamp = FigureS10Heijman();
Heijman_APClamp = reformat_data(Heijman_APClamp,11);

logscalefactors = linspace(log(1/3),log(3),11) ;
apclamp_settings.repol_change =  exp(logscalefactors(1:11)); % scaling

BL = find(round(apclamp_settings.repol_change,2) == 1.00);
AKr = cell2mat(Heijman_APClamp.Area_Kr);
AKs = cell2mat(Heijman_APClamp.Area_Ks);

figure(fig1)
subplot(3,3,1)
loglog(apclamp_settings.repol_change,AKr/AKr(BL),'k-x', 'linewidth',2)
hold on
loglog(apclamp_settings.repol_change,AKs/AKs(BL),'k-o','linewidth',2)
title('Heijman')

figure(fig2)
subplot(3,3,1)
bar([AKr(end)/AKr(BL),AKs(end)/AKs(BL)],0.5)
set(gca,'XTickLabels',{'IKr','IKs'})
title('Heijman')