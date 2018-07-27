function datatable = FigureS7Heijman
%% Figure S7 Heijman: Inducing proarrhythmic behavior through constant inward current 
%% inject in Heijman model. 
%--- Note: 
% Results displayed in manuscript were run using MATLAB 2016a on a 64bit
% Intel Processor. For exact replication of figures it is best to use these settings.

%--- Description of Figure: 
% Injecting a constant inward current while voltage is greater than -60 mV
% to induce proarrhythmic behavior in each of the models.  

%---: Functions required to run this script :---%
% mainHRdBA_Inject.m - runs injection simulation for Heijman only 
%--------------------------------------------------------------------------

% settings
settings.PCL = 1000;
settings.freq = 100;
settings.storeLast = 2; % save both beats 99 and 100
settings.stimdur = 2;
settings.Istim = -36.7;
settings.showProgress = 1;
settings.Inject = 0;
Injects = [0 0.1 0.2];

settings.ISO = 0;
flags.SS = 1;

flags.ICaL = 1;
flags.IKs = 1;
flags.PLB = 1;
flags.TnI = 1;
flags.INa = 1;
flags.INaK = 1;
flags.RyR = 1;
flags.IKur = 1;

figure
fig = gcf;

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
    plot(Ti,V,'linewidth',2)
    xlim([900,2000]) % adding the "time delay" in each plot
    xlabel('time (ms)')
    ylabel('V (mv)')
    ylim([-100 80])
    set(gcf,'Position',[20,20,300,300]) 
    hold on

end
