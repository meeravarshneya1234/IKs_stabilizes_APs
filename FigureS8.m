%% Figure S8: Inducing proarrhythmic behavior through constant inward current 
%% inject across multiple species.

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
%% Figure S8
%--- Description of Figure:
% Injecting a constant inward current while voltage is greater than -60 mV
% to induce proarrhythmic behavior in each of the models.  

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
inject_scales = {[0 0.1 0.4] %amount of current to inject 
[0 0.4 0.9]
[0 1.3 2.7]
[0 0.4 0.9]
[0 0.4 0.8]
[0 0.6 1.3]
[0 1 1.9]
[0 0.3 0.7]
[0 0.1 0.3]};

%---- Run Simulation ----%

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
    set(gca,'FontSize',12,'FontWeight','bold')
    title(modelnames{i})
    
    Xnew.inject_factor = settings.Inject(end);
    str = modelnames{i};
    inject_datatable.(str) = Xnew;     
end 

inject_datatable.Heijman = FigureS8Heijman();

%---- Plot Summary Plot ----%
models = {'Fox','Hund','Heijman','Livshitz','Devenyi','Shannon','TT04','TT06','Ohara','Grandi'};
figure 
summary_barplot = gcf;
ax_summary = axes('parent', summary_barplot);
Injects = abs(cell2mat(cellfun(@(x) inject_datatable.(x).inject_factor,models,'UniformOutput',0)));
bar(Injects,0.5)
set(ax_summary, 'xticklabel',models)
xtickangle(90)
ylabel('Current Injected (A/F)')
set(gca,'FontSize',12,'FontWeight','bold')

