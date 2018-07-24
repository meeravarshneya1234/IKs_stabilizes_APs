%%

% settings 
settings.freq =92;
settings.storeLast =3;
settings.stimdur = 2;
settings.Istim = -36.7;
settings.showProgress = 1;
settings.bcl =1000;

colors = hsv(4);iso_concs = 0;
flags.SS = 1;

%
flags.ICaL = 1;
flags.IKs = 1;
flags.PLB = 1;
flags.TnI = 1;
flags.INa = 1;
flags.INaK = 1;
flags.RyR = 1;
flags.IKur = 1;
Ca_scale = [1,1.5,1.9,2];

figure 
fig1 = gcf;
xlabel('Time (ms)','FontSize',12,'FontWeight','bold')
ylabel('Voltage (mV)','FontSize',12,'FontWeight','bold')

for ii = 1:length(iso_concs)
    settings.ISO = iso_concs(ii);   
    
    for i = 1:length(Ca_scale)
        settings.ICaLB = Ca_scale(i);

        [currents,State,Ti,APDs,settings]=mainHRdBA(settings,flags);
           
        datatable.APDs{ii,i} = APDs;
        datatable.times{ii,i} = Ti;
        datatable.V{ii,i} =  State(:,1);
        datatable.statevars{ii,i} =  State;
        datatable.currents{ii,i} = currents;
        %ICs.(str) = State(end,:);
        [~,indV] = max(State(:,1));
        V = State(:,1);
        x2 = find(floor(V)==floor(V(1)) & Ti > Ti(indV),1); % and ends
        
        datatable.Area_Ca{ii,i} = trapz(Ti(1:x2),currents.ical(1:x2));
        %datatable.Area_Ks(ii,i) = trapz(Ti(1:x2),currents.iks(1:x2));
        %datatable.Area_Kr(ii,i) = trapz(Ti(1:x2),currents.ikr(1:x2));
        %datatable.IKs_Fraction = datatable.Area_Ks(ii,i)/(datatable.Area_Kr(ii,i)+datatable.Area_Ks(ii,i));
        
%         stop = length(Ti);
%         start = find(Ti(end)-settings.bcl==Ti,1,'last');
%         T = Ti(start:stop); V = State(start:stop,1);
%         
%         T = T - min(T);
%         T = T +100;
        figure(fig1)
        plot(Ti,V,colors(i),'linewidth',2)
        
        hold on       
        %test_SS(Ti,State,settings,APDs,num_beats)
    end 
        
%         disp(['Calculated ISO Concentration = ', num2str(iso_concs(ii))]);
%         figurelegend{ii} = ['Iso = ' num2str(iso_concs(ii))];
              
%         filename = ['Heijman_SS_2000beats' num2str(iso_concs(ii)) '.mat'];
%         save(filename,'datatable')
    
end
