function plotting_APClamp(settings,datatable)

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
                         %% -- plotting_APClamp.m -- %%
% Description: Plot figures for Figure6.m O'Hara model AP clamp
% simualtions.

% Inputs:
% --> settings - [struct array] AP stimulation protocol (PCL,nBeats,...)
% --> datatable - [struct array] AP Clamp data

% Outputs:
% --> Xnew - [struct array] APs with EADs or failed to repolarize are
% removed (clean population)

%---: Functions required to run this function :---%
% ** mtit.m - Creates a title for subplots. (from MATLAB File Exchange)
%--------------------------------------------------------------------------

%--- Set Up Figures ---%
figure('name','Voltage')
Volt = gcf;

figure('name','IKs')
Ks_Figure = gcf;

figure('name','IKr')
Kr_Figure = gcf;

colors = jet(length(settings.vals_graph));

%-- AP Clamp Variables
time_clamp = cell2mat(arrayfun(@(x) x.time, datatable,'UniformOutput', 0));
volt_clamp = cell2mat(arrayfun(@(x) x.Voltage, datatable,'UniformOutput', 0));

IKs_clamp = cell2mat(arrayfun(@(x) x.IKs, datatable,'UniformOutput', 0));
xs1_clamp = cell2mat(arrayfun(@(x) x.statevars(:,35), datatable,'UniformOutput', 0));
xs2_clamp = cell2mat(arrayfun(@(x) x.statevars(:,36), datatable,'UniformOutput', 0));

IKr_clamp = cell2mat(arrayfun(@(x) x.IKr, datatable,'UniformOutput', 0));
xrf_clamp = cell2mat(arrayfun(@(x) x.statevars(:,33), datatable,'UniformOutput', 0));
xrs_clamp = cell2mat(arrayfun(@(x) x.statevars(:,34), datatable,'UniformOutput', 0));
Axrf=1.0./(1.0+exp((volt_clamp+54.81)/38.21));
Axrs=1.0-Axrf;
xr_clamp=Axrf.*xrf_clamp+Axrs.*xrs_clamp;
rkr_clamp=1.0./(1.0+exp((volt_clamp+55.0)/75.0))*1.0./(1.0+exp((volt_clamp-10.0)/30.0));

%% Figure 6A
set(Volt,'defaultaxescolororder',colors) %black and gray
figure(Volt)
plot(time_clamp,volt_clamp,'linewidth',2)
hold on
ylabel('voltage (mV)')
xlabel('time (ms)')
set(gca,'FontSize',12,'FontWeight','bold')
title('Figure 6A')

%% Figure 6B IKs 
set(Ks_Figure,'defaultaxescolororder',colors) %black and gray
figure(Ks_Figure)
subplot(2,2,1)
plot(time_clamp,IKs_clamp,'linewidth',2)
hold on
title('IKs')
set(gca,'FontSize',12,'FontWeight','bold')


subplot(2,2,3)
plot(time_clamp,xs1_clamp,'linewidth',2)
hold on
title('xs1')
set(gca,'FontSize',12,'FontWeight','bold')


subplot(2,2,4)
plot(time_clamp,xs2_clamp,'linewidth',2)
hold on
title('xs2')
set(gca,'FontSize',12,'FontWeight','bold')
mtit('Figure 6B - IKs')

%% Figure 6B - IKr 
set(Kr_Figure,'defaultaxescolororder',colors) %black and gray
figure(Kr_Figure)
subplot(2,2,1)
plot(time_clamp,IKr_clamp,'linewidth',2)
hold on
title('IKr')
set(gca,'FontSize',12,'FontWeight','bold')


subplot(2,2,3)
plot(time_clamp,xr_clamp,'linewidth',2)
hold on
title('xr')
set(gca,'FontSize',12,'FontWeight','bold')


subplot(2,2,4)
plot(time_clamp,rkr_clamp,'linewidth',2)
hold on
title('rkr')
set(gca,'FontSize',12,'FontWeight','bold')
mtit('Figure 6B - IKr')






