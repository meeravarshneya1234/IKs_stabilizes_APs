%% Figure S4:Altering the contribution of IKs within the same model  
%% influences arrhythmia susceptibility. 

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
%% Figures S4A,S4B
%--- Description of Figures: 
% Calcium perturbations in low(0.05x),baseline IKs models performed in TT04 Model

%---: Functions required to run this part :---%
% main_program.m - runs single AP simulation 
%--------------------------------------------------------------------------
%%
%---- Set Up Simulation ----%
settings.model_name = 'TT04';
settings.celltype = 'endo'; 
settings.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
settings.stim_amp = 25; % Stimulus amplitude
settings.nBeats = 91; % Number of beats to simulate 
settings.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 1;

% Ks, Kr, Ca scaling must have the same length vector. 
Ks_scale = [1 0.05]; % perturb IKs, set to 1 for baseline
Kr_scale = [1 3.13]; % perturb IKr, set to 1 for baseline
settings.Ca_scale = [1 8.1 16.1] ; % perturb ICaL, set to 1 for baseline

%---- Run and Plot Simulation ----%
strs = {'base', 'low'}; 
titlestrings = {'Figure S4A - TT04 Baseline IKs','Figure S4B - TT04 Low IKs'};
for i = 1:length(Ks_scale)
    settings.Ks_scale = Ks_scale(i);
    settings.Kr_scale = Kr_scale(i);
    X = main_program(settings);
    TT04.(strs{i}) = X;    
    
    figure % plot figures 
    plot(TT04.(strs{i}).times{1,1},TT04.(strs{i}).V{1,1},'linewidth',2)
    hold on
    plot(TT04.(strs{i}).times{1,2},TT04.(strs{i}).V{1,2},'linewidth',2)
    plot(TT04.(strs{i}).times{1,3},TT04.(strs{i}).V{1,3},'linewidth',2)
    ylim([-100 100])
    title(titlestrings(i))
    xlabel('time (ms)')
    ylabel('voltage (mV)')
    set(gca,'FontSize',12,'FontWeight','bold')
    legend('ICaL*1','ICaL*8.1','ICaL*16.1')

end

%% Figure S4C
%--- Description of Figures: 
% Calcium perturbations (wide range) in low(0.05x),baseline IKs models performed
% in TT04 Model

%---: Functions required to run this part :---%
% main_program.m - runs single AP simulation 
%--------------------------------------------------------------------------

%---- Set Up Simulation ----%
settings.model_name= 'TT04';
settings.celltype = 'endo'; 
settings.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
settings.stim_amp = 25; % Stimulus amplitude
settings.nBeats = 91 ; % Number of beats to simulate 
settings.numbertokeep = 1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 1;
% Ks,Kr,Ca scaling must have the same length vector.
Ks_scale = [1 0.05]; % perturb IKs, set to 1 for baseline
Kr_scale = [1 3.13]; % perturb IKr, set to 1 for baseline
settings.Ca_scale = 1:0.1:27.8; % perturb ICaL, set to 1 for baseline

%---- Run Simulation ----%
strs = {'base','low'}; 
for i = 1:length(strs)
    settings.Ks_scale = Ks_scale(i);
    settings.Kr_scale = Kr_scale(i);
    X = main_program(settings);   
    TT04_highlow.(strs{i}) = X;
end

% check for presence of EADs and find the index right before the first EAD
% (baseline iks model)
t_cutoff = 5;
flag = 0;
[APfails1,EADs1] = cleandata(TT04_highlow.base.APDs,TT04_highlow.base.times,TT04_highlow.base.V,t_cutoff,flag);
[~,ind_failed] = find(APfails1 ==1); % number of failed to repolarize
[~,ind_EADs] = find(EADs1 ==1); % number of EADs
total = [ind_EADs ind_failed];
base_APDs = TT04_highlow.base.APDs;
base_APDs(min(total):end) = NaN;

% check for presence of EADs and find the index right before the first EAD
% (low iks model)
[APfails2,EADs2] = cleandata(TT04_highlow.low.APDs,TT04_highlow.low.times,TT04_highlow.low.V,t_cutoff,flag);
[~,ind_failed] = find(APfails2 ==1); % number of failed to repolarize
[~,ind_EADs] = find(EADs2==1); % number of EADs
total = [ind_EADs ind_failed];
low_APDs = TT04_highlow.low.APDs;
low_APDs(min(total):end) = NaN;

%---- Plot Simulation ----%
figure
plot(settings.Ca_scale,base_APDs,'linewidth',2,'color','k')
hold on 
plot(settings.Ca_scale,low_APDs,'linewidth',2,'color','b')
xlabel('ICaL Factor')
ylabel('APDs (ms)')
title('Figure S4C')
set(gca,'FontSize',12,'FontWeight','bold')

%% Figures S2D,S2E
%--- Description of Figures: 
% Calcium perturbations in high(20x),baseline IKs models performed
% in Grandi Model

%---: Functions required to run this part :---%
% main_program.m - runs single AP simulation 
%--------------------------------------------------------------------------

%---- Set Up Simulation ----%
settings.model_name = 'Grandi';
settings.celltype = 'endo'; 
settings.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
settings.stim_amp = 20.6; % Stimulus amplitude
settings.nBeats = 91; % Number of beats to simulate 
settings.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 1;
% Ks,Kr,Ca scaling must have the same length vector. 
Ks_scale = [1 20]; % perturb IKs, set to 1 for baseline
Kr_scale = [1 0.04]; % perturb IKr, set to 1 for baseline
settings.Ca_scale = [1 1.6 2]; % perturb ICaL, set to 1 for baseline

%---- Run and Plot Simulation ----%
strs = {'base','high'}; 
titlestrings = {'Figure S4D - Grandi Baseline IKs','Figure S4E - Grandi High IKs'};
for i = 1:length(Ks_scale)
    settings.Ks_scale = Ks_scale(i);
    settings.Kr_scale = Kr_scale(i);
    X = main_program(settings);
    Grandi.(strs{i}) = X;    
    
    figure % plot figures 
    plot(Grandi.(strs{i}).times{1,1},Grandi.(strs{i}).V{1,1},'linewidth',2)
    hold on
    plot(Grandi.(strs{i}).times{1,2},Grandi.(strs{i}).V{1,2},'linewidth',2)
    plot(Grandi.(strs{i}).times{1,3},Grandi.(strs{i}).V{1,3},'linewidth',2)
    ylim([-100 100])
    title(titlestrings(i))
    xlabel('time (ms)')
    ylabel('voltage (mV)')
    set(gca,'FontSize',12,'FontWeight','bold')
    legend('ICaL*1','ICaL*1.6','ICaL*2')

end
%% Figures S4F
%--- Description of Figures: 
% Calcium perturbations (wide range) in high(20x) and baseline IKs models performed
% in Grandi Model

%---: Functions required to run this part :---%
% main_program.m - runs single AP simulation 
%--------------------------------------------------------------------------

%---- Set Up Simulation ----%
settings.model_name= 'Grandi';
settings.celltype = 'endo';
settings.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
settings.stim_amp = 20.6; % Stimulus amplitude
settings.nBeats = 91 ; % Number of beats to simulate 
settings.numbertokeep = 1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 1;

% Ks, Kr, Ca scaling must have the same length vector. 
Ks_scale = [1 20]; % perturb IKs, set to 1 for baseline
Kr_scale = [1 0.04]; % perturb IKr, set to 1 for baseline
settings.Ca_scale = 1:0.1:3.2; % perturb ICaL, set to 1 for baseline

%---- Run Simulation ----%
strs = {'base','high'}; 
for i = 1:length(strs)
    settings.Ks_scale = Ks_scale(i);
    settings.Kr_scale = Kr_scale(i);
    X = main_program(settings);   
    Grandi_highlow.(strs{i}) = X;
end

% check for presence of EADs and find the index right before the first EAD
% (baseline iks model)
t_cutoff = 5;
flag = 0;
[APfails1,EADs1] = cleandata(Grandi_highlow.base.APDs,Grandi_highlow.base.times,Grandi_highlow.base.V,t_cutoff,flag);
[~,ind_failed] = find(APfails1 ==1); % number of failed to repolarize
[~,ind_EADs] = find(EADs1 ==1); % number of EADs
total = [ind_EADs ind_failed];
base_APDs = Grandi_highlow.base.APDs;
base_APDs(min(total):end) = NaN;

% check for presence of EADs and find the index right before the first EAD
% (high iks model)
[APfails2,EADs2] = cleandata(Grandi_highlow.high.APDs,Grandi_highlow.high.times,Grandi_highlow.high.V,t_cutoff,flag);
[number_of_failed,ind_failed] = find(APfails2 ==1); % number of failed to repolarize
[number_of_EADs,ind_EADs] = find(EADs2==1); % number of EADs
total = [ind_EADs ind_failed];
high_APDs = Grandi_highlow.high.APDs;
high_APDs(min(total):end) = NaN;

%---- Run Simulation ----%
figure
plot(settings.Ca_scale,base_APDs,'linewidth',2,'color','k')
hold on 
plot(settings.Ca_scale,high_APDs,'linewidth',2,'color','r')
xlabel('ICaL Factor')
ylabel('APDs (ms)')
set(gca,'FontSize',12,'FontWeight','bold')
title('Figure 24F')

