%% ---- Simulations for Baseline Analysis of all 10 models ---- %%
%% Plot Figures 3A, 3B, Supplement 1 

%% Plot Figure 3A & Supplement Figure 1 
load('baseline_data.mat','datatable') %saved data from the simulation 
models = {'Fox','Hund','Livshitz','Devenyi','Shannon','TT04','TT06','Ohara','Grandi'};

for ii = 1:length(models)   
    figure
    subplot(1,2,1) % AP plot
    plot(datatable.(models{ii}).times{:},datatable.(models{ii}).V{:},'linewidth',2)
    hold on
    xlabel('time (ms)')
    ylabel('V (mv)')
    set(gcf,'Position',[20,20,600,300])    
    
    subplot(1,2,2) % IKs and IKr plots
    plot(datatable.(models{ii}).times{:},datatable.(models{ii}).IKs{:},'color','b','linewidth',2);
    hold on
    plot(datatable.(models{ii}).times{:},datatable.(models{ii}).IKr{:},'color','r','linewidth',2);
    xlabel('time (ms)')
    ylabel('current (A/F)')
    legend('IKs','IKr')
    set(gcf,'Position',[20,20,600,300])

    mtit(models{ii},...
        'fontsize',14);
end

%% Plot Supplement Figure 1 - Heijman model 
figure
subplot(1,2,1) % AP plot
plot(datatable.Heijman_ISO_0.times,datatable.Heijman_ISO_0.V,'linewidth',2)
hold on
xlabel('time (ms)')
ylabel('V (mv)')
xlim([900,2000])% adding the "time delay" in each plot

subplot(1,2,2) % IKs and IKr plots
plot(datatable.Heijman_ISO_0.times,datatable.Heijman_ISO_0.currents.iks,'color','b','linewidth',2);
hold on
plot(datatable.Heijman_ISO_0.times,datatable.Heijman_ISO_0.currents.ikr,'color','r','linewidth',2);
xlabel('time (ms)')
ylabel('current (A/F)')
legend('IKs','IKr')
xlim([900,2000])% adding the "time delay" in each plot

mtit('Heijman no ISO',...
    'fontsize',14);
set(gcf,'Position',[20,20,600,300])

%% Plot Figure 3B - IKs Fraction Barplot Grandi and TT04 
figure % IKs Fractions for two models (Figure 3B)
barplot = gcf;
ax_summary = axes('parent', barplot);
IKs_Fractions = [datatable.Grandi.IKs_Fraction datatable.Ohara.IKs_Fraction datatable.TT04.IKs_Fraction];
bar(IKs_Fractions,0.5)
set(ax_summary, 'xticklabel',{'Grandi','Ohara','TT04'})
ylabel('IKs Fraction')
xtickangle(90)

%% Plot Supplemental Figure 1 - IKs Fraction Barplot all models 
models = {'Fox','Hund','Heijman_ISO_0','Livshitz','Devenyi','Shannon','TT04','TT06','Ohara','Grandi'};
figure % IKs Fractions for all models (Supplement 1)
summary_barplot = gcf;
ax_summary = axes('parent', summary_barplot);
IKs_Fractions = cell2mat(cellfun(@(x) datatable.(x).IKs_Fraction,models,'UniformOutput',0));
bar(IKs_Fractions,0.5)
set(ax_summary, 'xticklabel',models)
ylabel('IKs Fraction')
xtickangle(90)
set(gcf,'Position',[20,20,600,300])
