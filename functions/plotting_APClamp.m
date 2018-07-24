function plotting_APClamp(settings,datatable)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('name','Voltage')
Volt = gcf;

figure('name','IKs')
Ks_Figure = gcf;

figure('name','IKr')
Kr_Figure = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

set(Volt,'defaultaxescolororder',colors) %black and gray
figure(Volt)
plot(time_clamp,volt_clamp,'linewidth',2)
hold on
ylabel('voltage (mV)')
xlabel('time (ms)')

set(Ks_Figure,'defaultaxescolororder',colors) %black and gray
figure(Ks_Figure)
subplot(2,2,1)
plot(time_clamp,IKs_clamp,'linewidth',2)
hold on
title('IKs')

subplot(2,2,3)
plot(time_clamp,xs1_clamp,'linewidth',2)
hold on
title('xs1')


subplot(2,2,4)
plot(time_clamp,xs2_clamp,'linewidth',2)
hold on
title('xs2')


set(Kr_Figure,'defaultaxescolororder',colors) %black and gray
figure(Kr_Figure)
subplot(2,2,1)
plot(time_clamp,IKr_clamp,'linewidth',2)
hold on
title('IKr')

subplot(2,2,3)
plot(time_clamp,xr_clamp,'linewidth',2)
hold on
title('xr')


subplot(2,2,4)
plot(time_clamp,rkr_clamp,'linewidth',2)
hold on
title('rkr')





