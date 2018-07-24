%% Graph Figure 2 based on saved data.  
% Altering the contribution of IKs within the same model influences arrhythmia susceptibility. 

%--- Notes: 
% Results displayed in manuscript were run using MATLAB 2016a on a 64bit
% Intel Processor. For exact replication of figures it is best to use these settings.

%% Plot Figure 2A,B,C - Ohara High, Base, Low APs with Calcium Perturbation 

load('ohara_highlow_EAD.mat')
strings = {'low','base','high'};
for i = 1:length(strings)
    figure
    plot(highlow_EAD.(strings{i}).times{1,1},highlow_EAD.(strings{i}).V{1,1},'linewidth',2)
    hold on
    plot(highlow_EAD.(strings{i}).times{1,2},highlow_EAD.(strings{i}).V{1,2},'linewidth',2)
    plot(highlow_EAD.(strings{i}).times{1,3},highlow_EAD.(strings{i}).V{1,3},'linewidth',2)
    ylim([-100 100])
    title([strings{i} ' IKs Model'])
    legend('ICaL*1','ICaL*7.5','ICaL*15.5')
    xlabel('time (ms)')
    ylabel('Voltage (mV)')
end

%% Figure 2D - Ohara High, Base, Low APD vs Calcium Perturbations 
Ca_scale = 1:0.1:22;
load('ohara_highlow_multipleEAD.mat')
figure
plot(Ca_scale,highlow_multipleEAD.low.APDs,'linewidth',2,'color','b')
hold on
plot(Ca_scale,highlow_multipleEAD.base.APDs,'linewidth',2,'color','k')
plot(Ca_scale,highlow_multipleEAD.high.APDs,'linewidth',2,'color','r')
xlabel('ICaL Factor')
ylabel('APDs (ms)')
