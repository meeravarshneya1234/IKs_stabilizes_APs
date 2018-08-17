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
% Description: Runs the Heijman model calcium perturbation simulation. 

% Outputs:
% --> X_Heijman - struct that outputs the APDs, time, voltage, and state variables 

%---: Functions required to run this script :---%
% mainHRdBA.m - runs Heijman model simulation (downloaded code online)
% mtit.m - create main title (from MATLAB file exchange) 
%--------------------------------------------------------------------------

% settings 
settings.freq =100; %number of beats to stimulate first EAD , see Figure S5 for details
settings.storeLast =3; % Determine how many beats to keep. 1 = last beat, 2 = last two beats 
settings.stimdur = 2;% Stimulus duration
settings.Istim = -36.7; %Stimulus amplitude for each model, see Table S1 for details
settings.showProgress = 0;
settings.bcl =1000; % Interval bewteen stimuli,[ms]

colors = hsv(4);
settings.ISO = 0; % concentration of ISO; 0 = no ISO 
settings.SS = 1; % 1 - run steady state conditions 0 - do not run steady state conditions 

% option to block PKA targets: no block = 1; 100% block = 0 
flags.ICaL = 1; flags.IKs = 1; flags.PLB = 1; flags.TnI = 1; flags.INa = 1;
flags.INaK = 1; flags.RyR = 1; flags.IKur = 1;

Ca_scale = [1,1.5,1.9,2];

disp('Heijman Model')
figure
fig1 = gcf;

for i = 1:length(Ca_scale)
    settings.ICaLB = Ca_scale(i);
    
    [currents,State,Ti,APDs,settings]=mainHRdBA(settings,flags);
    
    X_Heijman.APDs{i} = APDs;
    X_Heijman.times{i} = Ti;
    X_Heijman.V{i} =  State(:,1);
    X_Heijman.statevars{i} =  State;
    X_Heijman.currents{i} = currents;
    [~,indV] = max(State(:,1));
    V = State(:,1);
    x2 = find(floor(V)==floor(V(1)) & Ti > Ti(indV),1); % and ends
    
    X_Heijman.Area_Ca(i) = trapz(Ti(1:x2),currents.ical(1:x2));
%     X_Heijman.Area_Ks(i) = trapz(Ti(1:x2),currents.iks(1:x2));
%     X_Heijman.Area_Kr(i) = trapz(Ti(1:x2),currents.ikr(1:x2));
%     X_Heijman.IKs_Fraction = datatable.Area_Ks(ii,i)/(datatable.Area_Kr(ii,i)+datatable.Area_Ks(ii,i));
    
    figure(fig1)
    subplot(1,2,1)
    plot(X_Heijman.times{i},X_Heijman.V{i},'color',colors(i,:),'linewidth',2);
    hold on
    xlabel('time (ms)')
    ylabel('V (mv)')
    set(gcf,'Position',[20,20,600,300])
    xlim([1900 3000])
    
    if i == 1 || i ==3
        subplot(1,2,2)
        plot(X_Heijman.times{i},X_Heijman.currents{i}.ical,'color',colors(i,:),'linewidth',2);
        hold on
        xlabel('time (ms)')
        ylabel('ICaL (A/F)')
        set(gcf,'Position',[20,20,600,300])
        xlim([1900 3000])
    end
    disp(['Completed GKs*1 & PCa*', num2str(Ca_scale(i))])
    
end
mtit('Heijman' ,...
    'fontsize',14);
