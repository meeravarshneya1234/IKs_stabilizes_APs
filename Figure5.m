%% Figure 5:Comparison of population variability & susceptibility to EADs 
%% across models from multiple species. 
%--- Note: 
% Results displayed in manuscript were run using MATLAB 2016a on a 64bit
% Intel Processor. For exact replication of figures it is best to use these 
% settings.

%--- Description of Figures: 
% To test whether the protective effect of high IKs levels can be considered 
% a general phenomenon, we examined population variability and the susceptibility 
% to EADs across models. 

%--- Functions used in this script:
% main_program.m - runs simulation 
% mtit - create main title (from MATLAB file exchange) 
%--------------------------------------------------------------------------
%% Figure 5A & 5B
% Data saved from earlier figures, plotted here:
% (1) IKs Fraction calculated in FigureS3.m
% (2) APD Spread calculated in Figure3.m and FigureS4.m 
% (3) EAD calculated in Figure4.m and FigureS5.m 
% (4) Inject calculated in FigureS8.m 
modelnames = {'Fox','Hund','Heijman','Shannon','Livshitz','Devenyi','TT04','TT06','Ohara','Grandi'};
load('data.mat')
IKs_Fraction = cell2mat(cellfun(@(x) datatable.(x).IKs_Fraction,modelnames,'UniformOutput',0));
APDSpread = cell2mat(cellfun(@(x) datatable.(x).APDSpread,modelnames,'UniformOutput',0));
EAD = cell2mat(cellfun(@(x) datatable.(x).EAD,modelnames,'UniformOutput',0));
Inject = cell2mat(cellfun(@(x) datatable.(x).Injection,modelnames,'UniformOutput',0));

% Figure 5A
figure 
plot(IKs_Fraction,APDSpread,'o','MarkerSize',10,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
text(IKs_Fraction,APDSpread, modelnames);
h=lsline;
xlabel('IKs Fraction')
ylabel('APD Spread')

% Figure 5B
figure 
plot(IKs_Fraction,abs(EAD),'o','MarkerSize',10,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor',[.6 .6 1])
text(IKs_Fraction,abs(EAD), modelnames);
h=lsline;
xlabel('IKs Fraction')
ylabel('AUC ICaL Before EAD')

% Figure 5C
figure 
plot(IKs_Fraction,Inject,'o','MarkerSize',10,...
    'MarkerEdgeColor','green',...
    'MarkerFaceColor',[.6 1 .6])
text(IKs_Fraction,Inject, modelnames);
h=lsline;
xlabel('IKs Fraction')
ylabel('Current Injected')

