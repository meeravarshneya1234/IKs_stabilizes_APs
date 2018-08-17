%% Figure S8: Additional correlation plots 

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
%% Figure S8A
%--- Description of Figure:
% The scaling factor (ICaL factor) by which the calcium current was
% increased does not correlate with the levels of Iks. Since the baseline
% ICaL magnitude differed between models, the ICaL factor did not provide a
% fair basis of comparison. To compare susceptibility to EADs, we
% integrated the ICaL waveform at the level just under the EAD threshold as
% depicted in Figures 4 & 5.

%---: Functions required to run this part :---%
% inject_current_program.m - runs injection simulation using parfor loop
% reformat_data.m - reformats the data collected from inject_current_program
%--------------------------------------------------------------------------
%%
modelnames = {'Fox','Hund','Heijman','Shannon','Livshitz','Devenyi','TT04','TT06','Ohara','Grandi'};
load('data.mat')
IKs_Fraction = cell2mat(cellfun(@(x) datatable.(x).IKs_Fraction,modelnames,'UniformOutput',0));
ICaL_Factor = [39 67.8 2 3 58.1 3 27.8 22.4 15.2 1.6]; % values obtained from Figure4.m and FigureS6.m.

figure 
plot(IKs_Fraction,ICaL_Factor,'o','MarkerSize',10,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
text(IKs_Fraction,ICaL_Factor, modelnames);
h=lsline;
xlabel('IKs Fraction')
ylabel('ICaL Factor')
title('Figure S8A')
set(gca,'FontSize',12,'FontWeight','bold')

%% Figure S8B
%--- Description of Figure:
% The percent change in APD between the APD when 0.1 A/F inward current was
% injected into the models and the baseline APD was calculated and compared
% with the levels of IKs. This shows that models with low IKs are much more
% susceptible to APD prolongation when an arbitrary inward current is
% injected as compared to models with high IKs.

%---: Functions required to run this part :---%
% inject_current_program.m - runs injection simulation using parfor loop
% reformat_data.m - reformats the data collected from inject_current_program
%--------------------------------------------------------------------------
%%
%---- Set Up Simulation ----%
modelnames = {'Fox','Hund','Livshitz','Devenyi','Shannon','TT04','TT06','Ohara','Grandi'};
% options - 'Fox', 'Hund', 'Livshitz','Devenyi','Shannon','TT04','TT06','Grandi','Ohara'

celltypes = {'','','','','','endo','endo','endo','endo'}; 
% size should be same as model_name, enter one for each model
% options only available for human models as follows
% TT04, TT06, Ohara - 'epi', 'endo', 'mid'
% Grandi - 'epi', 'endo'

settings.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
amps = [36.4,32.2,35,30.9,35,25,22.6,32.2,20.6]; % Stimulus amplitude (see Table S1 for details) 
settings.nBeats = 100 ; % Number of beats to simulate 
settings.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 1;% 1 - run steady state conditions 0 - do not run steady state conditions 
settings.Inject = [0 0.1]; %amount of current to inject 

%---- Run Simulation ----%

for i = 1:length(modelnames) 
    settings.model_name = modelnames{i};
    settings.celltype = celltypes{i};
    settings.stim_amp = amps(i);
    X = inject_current_program(settings);    
    Xnew = reformat_data(X,length(settings.Inject));  
    
    APDs = cell2mat(Xnew.APDs);
    change(i) =(APDs(2) - APDs(1))/APDs(1);
end 

% Heijman model sims 
X = FigureS8Heijman();
Xnew = reformat_data(X,1);
APDs = cell2mat(Xnew.APDs);
change_Heijman = (APDs(2) - APDs(1))/APDs(1);

%---- Plot Correlation ----%
change = [change(1:2) change_Heijman change(3:end)]; %percent change in APD between the APD when 0.1 A/F inward current and baseline APD
models = {'Fox','Hund','Heijman','Livshitz','Devenyi','Shannon','TT04','TT06','Ohara','Grandi'};
figure 
plot(IKs_Fraction,change,'o','MarkerSize',10,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
text(IKs_Fraction,change, models);
h=lsline;
xlabel('IKs Fraction')
title('Figure S8B')
ylabel('%Change in APD at 0.1 A/F Inject')
set(gca,'FontSize',12,'FontWeight','bold')