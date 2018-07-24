%% Recreate Figure 1C - Population variability for Grandi and TT04 
load('population_data.mat')

figure
histos = gcf;
colors = 'rgb';
model_names = {'Grandi','Ohara','TT04'};
for i = 1:length((model_names))
    figure
    for ii = 1:20 % only plot the first 20 for the figure
        plot(new_datatable.(model_names{i}).times{ii},new_datatable.(model_names{i}).V{ii},'linewidth',2)
        hold on
    end
    title(model_names{i})
    xlabel('time (ms)')
    ylabel('voltage (mV)')
    ylim([-100 60])
      
    % histogram calculations 
    APDs = new_datatable.(model_names{i}).APDs;
%     newAPDs = APDs - median(APDs);
%     medians(i) = median(APDs);
    newAPDs = APDs; 
    temp = min(newAPDs):25:max(newAPDs);
    bins = linspace(min(newAPDs),max(newAPDs),length(temp));
    
    figure(histos)
    histoutline(newAPDs,bins,colors(i),'linewidth',4);
    hold on  
end

figure(histos) 
set(gca,'FontSize',12,'FontWeight','bold')
xlabel('APD (ms)','FontSize',12,'FontWeight','bold')
ylabel('Count','FontSize',12,'FontWeight','bold')
hold off
box('off')
legend(model_names)


%% Create Supplement 2 - Histograms of all populations except Grandi, TT04
model_names = {'Fox','Hund','Livshitz','Devenyi','Shannon','TT06','Ohara','Heijman_ISO_0'};
for i = 1:length((model_names))     
    % histogram calculations 
    APDs = new_datatable.(model_names{i}).APDs;
    newAPDs = APDs - median(APDs);
    temp = min(newAPDs):25:max(newAPDs);
    bins = linspace(min(newAPDs),max(newAPDs),length(temp));
    medians(end+1) = median(APDs);

    figure
    histoutline(newAPDs,bins,'linewidth',4);
    hold on
    title(model_names(i))
    set(gca,'FontSize',12,'FontWeight','bold')
    xlabel('APD (ms)','FontSize',12,'FontWeight','bold')
    ylabel('Count','FontSize',12,'FontWeight','bold')
    xlim([-300 300])
    ylim([0 180])

end

%% Create Supplement 2 - Bar Graphs of the APD Spread 
all_names = {'Fox','Hund','Heijman_ISO_0','Livshitz','Devenyi','Shannon','TT04','TT06','Ohara','Grandi',};

figure 
summary_barplot = gcf;
ax_summary = axes('parent', summary_barplot);
Spread = cell2mat(cellfun(@(x) new_datatable.(x).APDSpread,all_names,'UniformOutput',0));
bar(Spread,0.5)
set(ax_summary, 'xticklabel',all_names)
ylabel('APD Spread')


