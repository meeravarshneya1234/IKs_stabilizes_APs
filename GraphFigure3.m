%% Graph Figure 2 based on saved data.  
% Dramatic differences seen between two human myocyte models in K+ currents and population variability.

%--- Notes: 
% Results displayed in manuscript were run using MATLAB 2016a on a 64bit
% Intel Processor. For exact replication of figures it is best to use these settings.

%% Plot Figure 3A - Grandi and TT04 APs + IKs/IKr Waveforms 

load('baseline_data.mat','base_datatable') %saved data from the simulation 
models = {'Grandi','Ohara','TT04'};

for ii = 1:length(models)   
    figure
    subplot(1,2,1) % AP plot
    plot(base_datatable.(models{ii}).times{:},base_datatable.(models{ii}).V{:},'linewidth',2)
    hold on
    xlabel('time (ms)')
    ylabel('V (mv)')
    set(gcf,'Position',[20,20,600,300])  
    xlim([0 700])    
    ylim([-100 50])
    
    subplot(1,2,2) % IKs and IKr plots
    plot(base_datatable.(models{ii}).times{:},base_datatable.(models{ii}).IKs{:},'color','b','linewidth',2);
    hold on
    plot(base_datatable.(models{ii}).times{:},base_datatable.(models{ii}).IKr{:},'color','r','linewidth',2);
    xlabel('time (ms)')
    ylabel('current (A/F)')
    legend('IKs','IKr')
    set(gcf,'Position',[20,20,600,300])
    xlim([0 700])


    mtit(models{ii},...
        'fontsize',14);
end

%% Plot Figure 3B - IKs Fraction Barplot Grandi and TT04 
figure 
barplot = gcf;
ax_summary = axes('parent', barplot);
IKs_Fractions = [base_datatable.Grandi.IKs_Fraction base_datatable.Ohara.IKs_Fraction base_datatable.TT04.IKs_Fraction];
bar(IKs_Fractions,0.5)
set(ax_summary, 'xticklabel',{'Grandi','Ohara','TT04'})
ylabel('IKs Fraction')
xtickangle(90)

%% Plot Figure 3C - Population variability for Grandi and TT04 
load('population_data.mat')

figure
histos = gcf;
colors = 'rbg';
model_names = {'Grandi','Ohara','TT04'};
for i = 1:length((model_names))
    figure
    for ii = 1:20 % only plot the first 20 for the figure
        plot(pop_datatable.(model_names{i}).times{ii},pop_datatable.(model_names{i}).V{ii},'linewidth',2)
        hold on
    end
    title(model_names{i})
    xlabel('time (ms)')
    ylabel('voltage (mV)')
    ylim([-100 60])
    xlim([0 700])
      
    % histogram calculations 
    APDs = pop_datatable.(model_names{i}).APDs;
    temp = min(APDs):25:max(APDs);
    bins = linspace(min(APDs),max(APDs),length(temp));
    
    figure(histos)
    histoutline(APDs,bins,colors(i),'linewidth',4);
    hold on  
end

figure(histos) 
set(gca,'FontSize',12,'FontWeight','bold')
xlabel('APD (ms)','FontSize',12,'FontWeight','bold')
ylabel('Count','FontSize',12,'FontWeight','bold')
hold off
box('off')
legend(model_names)