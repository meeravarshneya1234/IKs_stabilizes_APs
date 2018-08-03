function X_Heijman = FigureS3Heijman()

settings.PCL = 1000;
settings.freq =1;
settings.storeLast = 2; % save both beats 99 and 100
settings.stimdur = 2;
settings.Istim = -36.7;
settings.showProgress = 1;

%colors = repmat('krgbmc',1,1000);
iso_concs = 0;
settings.SS = 1;

%
flags.ICaL = 1;
flags.IKs = 1;
flags.PLB = 1;
flags.TnI = 1;
flags.INa = 1;
flags.INaK = 1;
flags.RyR = 1;
flags.IKur = 1;

for ii = 1:length(iso_concs)
    settings.ISO = iso_concs(ii);
    
    [currents,State,Ti,APDs,settings]=mainHRdBA(settings,flags);
    
    str = ['ISO_' num2str(iso_concs(ii))];
    
    X_Heijman.(str).APDs = APDs;
    X_Heijman.(str).times = Ti;
    X_Heijman.(str).V =  State(:,1);
    X_Heijman.(str).statevars =  State;
    X_Heijman.(str).currents = currents;       
   
    % plot figure 
    V = State(:,1);
    IKs = currents.iks; 
    IKr = currents.ikr; 
    
    figure
    subplot(1,2,1)
    plot(Ti,V,'linewidth',2)
    xlim([900,2000]) % adding the "time delay" in each plot
    xlabel('Time (ms)','FontSize',12,'FontWeight','bold')
    ylabel('Voltage (mV)','FontSize',12,'FontWeight','bold')
    set(gcf,'Position',[20,20,600,300])
    
    
    subplot(1,2,2)
    plot(Ti,IKs,'linewidth',2)
    hold on
    plot(Ti,IKr,'linewidth',2)
    xlabel('Time (ms)','FontSize',12,'FontWeight','bold')
    xlim([900,2000])% adding the "time delay" in each plot
    set(gcf,'Position',[20,20,600,300])
    mtit('Heijman','fontsize',14);

    
    % get information about the last AP (beat 100)
    stop = length(Ti); %beginning of AP
    start = find(Ti(end)-settings.bcl==Ti,1,'last'); % end of AP
    
    Tlast = Ti(start:stop); 
    Vlast = State(start:stop,1);
    ICaL_last = currents.ical(start:stop,1);
    IKs_last = currents.iks(start:stop,1);
    IKr_last = currents.ikr(start:stop,1);
    [~,indV] = max(Vlast);
    x2 = find(floor(Vlast)==floor(Vlast(1)) & Tlast > Tlast(indV),1); % and ends
    
    X_Heijman.(str).Area_Ca = trapz(Tlast(1:x2),ICaL_last(1:x2));
    X_Heijman.(str).Area_Ks = trapz(Tlast(1:x2),IKs_last(1:x2));
    X_Heijman.(str).Area_Kr = trapz(Tlast(1:x2),IKr_last(1:x2));
    X_Heijman.(str).IKs_Fraction = X_Heijman.(str).Area_Ks/(X_Heijman.(str).Area_Kr+X_Heijman.(str).Area_Ks);   

end
