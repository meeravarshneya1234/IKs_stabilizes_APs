modelnames = {'Fox','Hund','Shannon','Livshitz','Devenyi','TT04','TT06','Ohara','Grandi'};
% options - 'Fox', 'Hund', 'Livshitz',
% 'Devenyi','Shannon','TT04','TT06','Grandi','Ohara'

celltypes = {'','','','','','endo','endo','endo','endo'}; % size should be same as model_name, enter one for each model
%celltypes = {'endo'}; % size should be same as model_name, enter one for each model
% options only available for human models as follows
% TT04, TT06, Ohara - 'epi', 'endo', 'mid'
% Grandi - 'epi', 'endo'

settings.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
amps = [36.4,32.2,35,35,30.9,25,22.6,32.2,20.6]; % Stimulus amplitude 
settings.nBeats = 100 ; % Number of beats to simulate 
settings.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 1;
settings.Inject = 0:0.1:3; %current to inject into the models 
%% Run Simulation 

for i = 1:length(modelnames) 
    settings.model_name = modelnames{i};
    settings.celltype = celltypes{i};
    settings.stim_amp = amps(i);
    X = inject_current_program(settings);    
    Xnew = reformat_data(X,length(settings.Inject));  
    APDs = cell2mat(Xnew.APDs);
    
    % find which 
    for ii = 1:length(APDs)
        if APDs(ii) > 500
            APDs(ii) = NaN;
        end
    end
    index = find(isnan(APDs),1);
    if isempty(index)
        index = length(APDs);
    end
    
    Xnew.change_APD = ((APDs-APDs(1))/APDs(1))*100;
    Xnew.inject_factor = settings.Inject(index);
    str = modelnames{i};
    inject_datatable.(str) = Xnew;     
end 

load('baseline_data.mat')
IKs_Fraction = cell2mat(cellfun(@(x) base_datatable.(x).IKs_Fraction,modelnames,'UniformOutput',0));

inject_factor = cell2mat(cellfun(@(x) inject_datatable.(x).inject_factor,modelnames,'UniformOutput',0));
change = cell2mat(cellfun(@(x) inject_datatable.(x).change_APD(2),modelnames,'UniformOutput',0));  

figure
plot(IKs_Fraction,inject_factor,'o')
text(IKs_Fraction,inject_factor, modelnames);
h=lsline;
xlabel('IKs Fraction')
ylabel('Current Injection (A/F)')


figure 
plot(IKs_Fraction,change,'o')
text(IKs_Fraction,change, modelnames);
h=lsline;
xlabel('IKs Fraction')
ylabel('Change in APD')

