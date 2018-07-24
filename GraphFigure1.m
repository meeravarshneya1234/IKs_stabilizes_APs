%% Graph Figure 1 based on saved data.  
% Altering the contribution of IKs within the same model influences population variability. 

%--- Notes: 
% Results displayed in manuscript were run using MATLAB 2016a on a 64bit
% Intel Processor. For exact replication of figures it is best to use these settings.

%% Plot Figure 1A - Ohara High, Low APs
load('baseline_data.mat','base_datatable')%-original ohara model data 
load('ohara_highlow_BL.mat','highlow_BL')%-high and low ohara model data 

figure
plot(base_datatable.Ohara.times{1,1},base_datatable.Ohara.V{1,1},'linewidth',2,'color','k') %baseline ohara model
hold on 
plot(highlow_BL.times{1,1},highlow_BL.V{1,1},'linewidth',2,'color','b') %-low ohara iks model
plot(highlow_BL.times{2,1},highlow_BL.V{2,1},'linewidth',2,'color','r') %-high phara iks model
xlabel('time(ms)')
ylabel('Voltage(mV)')

%% Plot Figure 1B - Ohara High, Low IKs & IKr Waveforms + IKs Fraction
figure
subplot(1,3,1) % plot the baseline ohara iks and ikr waveforms 
plot(base_datatable.Ohara.times{1,1},base_datatable.Ohara.IKs{1,1},'linewidth',2,'color','k')
hold on 
plot(base_datatable.Ohara.times{1,1},base_datatable.Ohara.IKr{1,1},'--k','linewidth',2)
ylim([-0.1 1])
xlim([0 600])
title('Ohara BL Model')
xlabel('time(ms)')
ylabel('Current (A/F)')

subplot(1,3,2) % plot the low ohara iks and ikr waveforms 
plot(highlow_BL.times{1,1},highlow_BL.IKs{1,1},'linewidth',2,'color','b')
hold on 
plot(highlow_BL.times{1,1},highlow_BL.IKr{1,1},'--b','linewidth',2)
ylim([-0.1 1])
xlim([0 600])
title('Ohara Low Model')
xlabel('time(ms)')
ylabel('Current (A/F)')

subplot(1,3,3) % plot the low ohara iks and ikr waveforms 
plot(highlow_BL.times{2,1},highlow_BL.IKs{2,1},'linewidth',2,'color','r')
hold on 
plot(highlow_BL.times{2,1},highlow_BL.IKr{2,1},'--r','linewidth',2)
ylim([-0.1 1])
xlim([0 600])
title('Ohara High Model')
xlabel('time(ms)')
ylabel('Current (A/F)')

% plot Iks fractions 
figure
IKs_Fractions = [highlow_BL.IKs_Fraction(1) base_datatable.Ohara.IKs_Fraction highlow_BL.IKs_Fraction(2)];
bar(IKs_Fractions,0.5)
xticklabels({'Low','Baseline','High'})
ylabel('IKs Fraction')

% zoomed in on low iks model iks 
figure
bar(highlow_BL.IKs_Fraction(1),0.5)
ylim([0 0.01])
ylabel('low IKs Fraction zoomed in')

%% Plot Figures 1C & 1D
load('population_data.mat','pop_datatable') %-original ohara model data 
load('ohara_highlow_pop.mat','highlow_pop') %-high and low ohara model data 

figure  % plot the baseline ohara histogram 
pop1 = gcf;
APDs = pop_datatable.Ohara.APDs;
temp = min(APDs):25:max(APDs);
bins = linspace(min(APDs),max(APDs),length(temp));
histoutline(APDs,bins,'linewidth',4,'color','k');
hold on

pert1 = prctile(APDs,90);
pert2 = prctile(APDs,10);
spreads(1) =(pert1 - pert2 )/ median(APDs); % baseline ohara APD spread  

colors = 'br';
for i = 1:2    % plot the high low ohara histogram 
    % histogram calculations 
    APDs = highlow_pop.APDs(:,i);
    temp = min(APDs):25:max(APDs);
    bins = linspace(min(APDs),max(APDs),length(temp));
    
    figure(pop1)
    histoutline(APDs,bins,'linewidth',4,'color',colors(i));
    hold on
    set(gca,'FontSize',12,'FontWeight','bold')
    xlabel('APD (ms)','FontSize',12,'FontWeight','bold')
    ylabel('Count','FontSize',12,'FontWeight','bold')

    pert1 = prctile(APDs,90);
    pert2 = prctile(APDs,10);
    spreads(i+1) =(pert1 - pert2 )/ median(APDs); % high low ohara APD spread  
end

figure % plot all APD Spreads 
summary_barplot = gcf;
ax_summary = axes('parent', summary_barplot);
bar(spreads,0.5)
ylabel('APD Spread')
xticklabels({'Baseline','Low','High'})
