%% ---- Simulations for Calcium Perturbation Analysis of all 9 models ---- %%
% Heijman model was run in: 
%% Plot Figures 4 & Supplement 4

%% Plot Figure 4A & 4B  - Grandi and TT04 ICaL Perturbations 
load('ead_data.mat','datatable')
modelnames = {'TT04','Grandi','Ohara'};
ICa_scales = {[1 14.4 27.7 27.8] %calcium perturbation for each model based on order of modelnames
[1 1.3 1.5 1.6]
[1 8 14.9 15]};

for ii = 1:length(modelnames)
    figure 
    colors = hsv(length(ICa_scales{ii}));
    name = modelnames{ii};
    ca_scale = ICa_scales{ii};

    for i = 1:length(ca_scale)
        subplot(1,2,1)
        plot(datatable.(name).times{i},datatable.(name).V{i},'color',colors(i,:),'linewidth',2);
        hold on
        xlabel('time (ms)')
        ylabel('V (mv)')
        set(gcf,'Position',[20,20,600,300])
        ylim([-100 100])

        if i == 1 || i == 3
            subplot(1,2,2)
            plot(datatable.(name).times{i},datatable.(name).ICaL{i},'color',colors(i,:),'linewidth',2);
            hold on
            xlabel('time (ms)')
            ylabel('ICaL (A/F)')
            set(gcf,'Position',[20,20,600,300])
        end
       
    end
    %legend(figurestrings);
    mtit(name ,...
        'fontsize',14);
end

%% Plot Figure 4C - TT04 and Grandi ICaL Before First EAD 
ICaL= [datatable.Grandi.Area_Ca(1) datatable.Ohara.Area_Ca(1) datatable.TT04.Area_Ca(1);
datatable.Grandi.Area_Ca(3) datatable.Ohara.Area_Ca(3) datatable.TT04.Area_Ca(3)]';

figure
handle = gcf;
ax = axes('parent', handle);
bar(abs(ICaL),'stacked')
set(ax, 'xticklabel',({'Grandi','Ohara','TT04'}))

%% Plot Supplement 4 - APs + ICaLs Waveforms 

modelnames = {'Fox','Hund','Shannon','Livshitz','Devenyi','TT06','Ohara',};
ICa_scales = {[1 21 40.8 40.9] %calcium perturbation for each model based on order of modelnames
[1 34.5 67.9 68]
[1 2.15 3.2 3.3]
[1 21.3 41.5 41.6]
[1 1.95 2.8 2.9]
[1 11.65 22.2 22.3]
[1 8 14.9 15]};

for ii = 1:length(modelnames)
    figure 
    colors = hsv(length(ICa_scales{ii}));
    name = modelnames{ii};
    ca_scale = ICa_scales{ii};

    for i = 1:length(ca_scale)
        subplot(1,2,1)
        plot(datatable.(name).times{i},datatable.(name).V{i},'color',colors(i,:),'linewidth',2);
        hold on
        xlabel('time (ms)')
        ylabel('V (mv)')
        set(gcf,'Position',[20,20,600,300])
        ylim([-100 100])

        if i == 1 || i == 3
            subplot(1,2,2)
            plot(datatable.(name).times{i},datatable.(name).ICaL{i},'color',colors(i,:),'linewidth',2);
            hold on
            xlabel('time (ms)')
            ylabel('ICaL (A/F)')
            set(gcf,'Position',[20,20,600,300])
        end
       
    end
    %legend(figurestrings);
    mtit(name ,...
        'fontsize',14);
end

%% Plot Supplement 4 - APs + ICaLs Waveforms for Heijman Model 
figure
colors = hsv(4);
ca_scale = [1,1.5,1.9,2];
model_name = 'Heijman_ISO_0';
for i = 1:length(ca_scale)
    subplot(1,2,1)
    plot(datatable.(model_name).times{i},datatable.(model_name).V{i},'color',colors(i,:),'linewidth',2);
    hold on
    xlabel('time (ms)')
    ylabel('V (mv)')
    set(gcf,'Position',[20,20,600,300])
    xlim([1900 3000])
    
    if i == 1 || i ==3
        subplot(1,2,2)
        plot(datatable.(model_name).times{i},datatable.(model_name).currents{i}.ical,'color',colors(i,:),'linewidth',2);
        hold on
        xlabel('time (ms)')
        ylabel('ICaL (A/F)')
        set(gcf,'Position',[20,20,600,300])
        xlim([1900 3000])
    end
end

%% Create Supplememt 4 - ICaL BL and Before EAD Bar Plot 
all_names = {'Fox','Hund','Heijman_ISO_0','Livshitz','Devenyi','Shannon','TT04','TT06','Ohara','Grandi'};
datatable.Heijman_ISO_0.Area_Ca = cell2mat(datatable.Heijman_ISO_0.Area_Ca); 

figure % plot Baseline AUC of ICaL for all models 
summaryBL_barplot = gcf;
axBL_summary = axes('parent', summaryBL_barplot);
ICaL_Baseline = abs(cell2mat(cellfun(@(x) datatable.(x).Area_Ca(1),all_names,'UniformOutput',0)));
bar(ICaL_Baseline,0.5)
set(axBL_summary, 'xticklabel',all_names)
ylabel('Baseline ICaL')
savefig('BL_Ca.fig')


figure % plot Before First EAD AUC of ICaL for all models 
summary_barplot = gcf;
ax_summary = axes('parent', summary_barplot);
ICaL_Before_EAD = abs(cell2mat(cellfun(@(x) datatable.(x).Area_Ca(3),all_names,'UniformOutput',0)));
bar(ICaL_Before_EAD,0.5)
set(ax_summary, 'xticklabel',all_names)
ylabel('ICaL Before EAD')
savefig('EAD_Ca.fig')
