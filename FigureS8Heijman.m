function datatable = FigureS8Heijman
%% Figure S8 Heijman: Inducing proarrhythmic behavior through constant inward current 
%% inject in Heijman model. 

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
                            %% -- FigureS8Heijman.m -- %%
% Description: Injecting a constant inward current while voltage is greater than -60 mV
% to induce proarrhythmic behavior in each of the models. 

% Outputs:
% --> X_Heijman - struct that outputs the APDs, time, voltage, and state variables 

%---: Functions required to run this script :---%
% mainHRdBA_Inject.m - runs injection simulation for Heijman only 
%--------------------------------------------------------------------------
%%
% settings
settings.PCL = 1000;% Interval bewteen stimuli,[ms]
settings.freq = 100; %number of beats to stimulate first EAD 
settings.storeLast = 2; % Determine how many beats to keep. 1 = last beat, 2 = last two beats 
settings.stimdur = 2;% Stimulus duration
settings.Istim = -36.7;%Stimulus amplitude for each model, see Table S1 for details
settings.ISO = 0; %concentration of ISO to add; 0 = no ISO 
flags.SS = 1; % 1 - run steady state conditions 0 - do not run steady state conditions 
settings.showProgress = 0;

% option to block PKA targets: no block = 1; 100% block = 0 
flags.ICaL = 1; flags.IKs = 1; flags.PLB = 1; flags.TnI = 1; flags.INa = 1;
flags.INaK = 1; flags.RyR = 1; flags.IKur = 1;

Injects = [0 0.1 0.2]; %amount of current to inject into cell [A/F]
colors = hsv(length(Injects));

figure
fig = gcf;
disp('Heijman Model')
for ii = 1:length(Injects)
    
    settings.Inject = Injects(ii);
    [currents,State,Ti,APDs,settings]=mainHRdBA_Inject(settings,flags);
       
    datatable.APDs(:,ii) = APDs;
    datatable.times{:,ii} = Ti;
    datatable.V{:,ii} =  State(:,1);
    datatable.statevars{:,ii} =  State;
    datatable.currents{:,ii} = currents;       
   
    V = State(:,1);
   
    figure(fig)
    plot(Ti,V,'linewidth',2,'color',colors(ii,:))
    xlim([900,2000]) % adding the "time delay" in each plot
    xlabel('time (ms)')
    ylabel('V (mv)')
    ylim([-100 80])
    set(gcf,'Position',[20,20,300,300]) 
    hold on
    
    disp(['Current Inject = ' num2str(Injects(ii))])

end
datatable.inject_factor = Injects(end); %added here to collect the value of the current required to form EAD

