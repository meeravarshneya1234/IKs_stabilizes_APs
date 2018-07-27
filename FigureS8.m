%% Figure S8: Altering the contribution of IKs within the same model influences population variability. 
%--- Note: 
% Results displayed in manuscript were run using MATLAB 2016a on a 64bit
% Intel Processor. For exact replication of figures it is best to use these settings.

%% Figure 1A 
%--- Description of Figure: 
% Simulations in high(10x) and low(0.1x) IKs models performed in Ohara Model  

%---: Functions required to run this script :---%
% main_program.m - runs single AP simulation 
%--------------------------------------------------------------------------
%% Run Current Injection Simulation 
modelnames = {'Fox','Hund','Livshitz','Devenyi','Shannon','TT04','TT06','Ohara','Grandi'};
% options - 'Fox', 'Hund', 'Livshitz',
% 'Devenyi','Shannon','TT04','TT06','Grandi','Ohara'

celltypes = {'','','','','','endo','endo','endo','endo'}; 
% size should be same as model_name, enter one for each model
% options only available for human models as follows
% TT04, TT06, Ohara - 'epi', 'endo', 'mid'
% Grandi - 'epi', 'endo'

settings.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
amps = [36.4,32.2,35,30.9,35,25,22.6,32.2,20.6]; % Stimulus amplitude 
settings.nBeats = 100 ; % Number of beats to simulate 
settings.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 1;
inject_scales = {[0 0.1 0.4]
[0 0.4 0.9]
[0 1.3 2.7]
[0 0.4 0.9]
[0 0.4 0.8]
[0 0.6 1.3]
[0 1 1.9]
[0 0.3 0.7]
[0 0.1 0.3]};
%% Run Simulation 

for i = 1:length(modelnames) 
    settings.model_name = modelnames{i};
    settings.celltype = celltypes{i};
    settings.stim_amp = amps(i);
    settings.Inject = inject_scales{i};
    X = inject_current_program(settings);    
    Xnew = reformat_data(X,length(settings.Inject));  
    colors = hsv(length(inject_scales{i}));
    
    figure
    for ii = 1:length(settings.Inject)
        plot(Xnew.times{ii},Xnew.V{ii},'color',colors(ii,:),'linewidth',2);
        hold on
        xlabel('time (ms)')
        ylabel('V (mv)')
        ylim([-100 80])
        set(gcf,'Position',[20,20,300,300])           
    end
    title(modelnames{i})
    
    Xnew.inject_factor = settings.Inject(end);
    str = modelnames{i};
    inject_datatable.(str) = Xnew;     
end 

inject_datatable.Heijman = FigureS8Heijman;
