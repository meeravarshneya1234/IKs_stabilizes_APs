function X_Heijman = FigureS6Heijman()
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
                            %% -- FigureS6Heijman.m -- %%
% Description: Runs the Heijman model population simulation. 

% Outputs:
% --> X_Heijman - outputs the APDs, time, voltage, and state variables 

%---: Functions required to run this script :---%
% mainHRdBA.m - runs Heijman model simulation (downloaded code online)
% clean_pop_Heijman.m - 
% reformat_data.m - 
%--------------------------------------------------------------------------
%%
% settings
settings.bcl = 1000;
settings.freq =100;
settings.storeLast = 1;
settings.stimdur = 2;
settings.Istim = -36.7;
settings.showProgress = 1;

% parameters to vary to create population, set to 1 for default
c.ICaLB = 1; c.IKsB = 1; c.IKrB = 1; c.INaKB =1; c.INaCaB = 1; c.IKpB = 1;
c.IK1B = 1; c.INabB = 1; c.ITo1B = 1; c.ITo2B = 1; c.INaB = 1; c.INaLB = 1;
c.IClB = 1; c.IpCaB = 1; c.ICabB = 1; c.IrelB = 1; c.IupB = 1; c.IleakB = 1;

variations = 300;
settings.sigma = 0.2;
scalings = exp(settings.sigma*randn(length(fieldnames(c)),variations))' ;

settings.totalvars = variations;
settings.t_cutoff = 3;
settings.flag = 1;

settings.ISO = 0;
settings.SS = 1;

% option to block PKA targets: no block = 1; 100% block = 0 
flags.ICaL = 1; flags.IKs = 1; flags.PLB = 1; flags.TnI = 1; flags.INa = 1;
flags.INaK = 1; flags.RyR = 1; flags.IKur = 1;

X = [];

parfor i = 1:variations
    scaling = scalings(i,:);
    [currents,State,Ti,APD]=mainHRdBA(settings,flags,scaling);
    X(i).times = Ti;
    X(i).V =  State(:,1);
    X(i).states =  State;
    X(i).currents= currents;
    X(i).APDs = APD;
    X(i).scalings = scaling;
    
    disp(['Model Variant # ', num2str(i)]);
end

X = reformat_data(X,variations);
X_Heijman = clean_pop_Heijman(settings,flags,X);
normAPDs = X_Heijman.APDs - median(X_Heijman.APDs);
temp = min(normAPDs):25:max(normAPDs);
bins = linspace(min(normAPDs),max(normAPDs),length(temp));
medians_Heijman = median(APDs);

pert1 = prctile(APDs,90);
pert2 = prctile(APDs,10);
X_Heijman.APDSpread =(pert1 - pert2)/ median(APDs);

figure
histoutline(normAPDs,bins,'linewidth',4);
hold on
title('Heijman')
set(gca,'FontSize',12,'FontWeight','bold')
xlabel('APD (ms)','FontSize',12,'FontWeight','bold')
ylabel('Count','FontSize',12,'FontWeight','bold')
xlim([-300 300])
ylim([0 180])
