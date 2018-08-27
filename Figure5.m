%% Figure 5:Comparison of population variability & susceptibility to EADs
%% across models from multiple species.

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
%% Figure 5
%--- Description of Figure:
% To test whether the protective effect of high IKs levels can be considered
% a general phenomenon, we examined population variability and the susceptibility
% to EADs across models.

%---: Functions required to run this script :---%
% main_program.m - runs simulation
% mtit - create main title (from MATLAB file exchange)
%--------------------------------------------------------------------------
%%
%---- Data saved from earlier figures, plotted here: ----%
% (1) IKs Fraction calculated in FigureS3.m
% (2) APD Spread calculated in Figure3.m and FigureS4.m
% (3) EAD calculated in Figure4.m and FigureS6.m
% (4) Inject calculated in FigureS7.m

modelnames = {'Fox','Hund','Heijman','Shannon','Livshitz','Devenyi','TT04','TT06','Ohara','Grandi'};
load('data.mat')
IKs_Fraction = cell2mat(cellfun(@(x) datatable.(x).IKs_Fraction,modelnames,'UniformOutput',0));
APDSpread = cell2mat(cellfun(@(x) datatable.(x).APDSpread,modelnames,'UniformOutput',0));
EAD = cell2mat(cellfun(@(x) datatable.(x).EAD,modelnames,'UniformOutput',0));
Inject = cell2mat(cellfun(@(x) datatable.(x).Injection,modelnames,'UniformOutput',0));

%% Plot Figure 5A - IKs Fraction vs APD Spread
figure
plot(IKs_Fraction,APDSpread,'o','MarkerSize',10,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
text(IKs_Fraction,APDSpread, modelnames);
h=lsline;
xlabel('IKs Fraction')
ylabel('APD Spread')
title('Figure 5A')
set(gca,'FontSize',12,'FontWeight','bold')


%% Plot Figure 5B - IKs Fraction vs ICaL before EAD
figure
plot(IKs_Fraction,abs(EAD),'o','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor',[.6 .6 1])
text(IKs_Fraction,abs(EAD), modelnames);
h=lsline;
xlabel('IKs Fraction')
ylabel('AUC ICaL Before EAD')
title('Figure 5B')
set(gca,'FontSize',12,'FontWeight','bold')


%% Plot Figure 5C - IKs Fraction vs Current Injected
figure
plot(IKs_Fraction,Inject,'o','MarkerSize',10,...
    'MarkerEdgeColor','green',...
    'MarkerFaceColor',[.6 1 .6])
text(IKs_Fraction,Inject, modelnames);
h=lsline;
xlabel('IKs Fraction')
ylabel('Current Injected')
title('Figure 5C')
set(gca,'FontSize',12,'FontWeight','bold')


