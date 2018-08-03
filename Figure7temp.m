% %% Figure 7: IKs stabilizes the AP during beta-adrenergic stimulation.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %--- "Slow delayed rectifier current protects ventricular myocytes from
% % arrhythmic dynamics across multiple species: a computational study" ---%
% 
% % By: Varshneya,Devenyi,Sobie 
% % For questions, please contact Dr.Eric A Sobie -> eric.sobie@mssm.edu 
% % or put in a pull request or open an issue on the github repository:
% % https://github.com/meeravarshneya1234/IKs_stabilizes_APs.git. 
% 
% %--- Note:
% % Results displayed in manuscript were run using MATLAB 2016a on a 64bit
% % Intel Processor. For exact replication of figures it is best to use these
% % settings.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %--------------------------------------------------------------------------
% %% Figure 7A 
% %--- Description of Figure: 
% % Population simulations of the following conditions:
% % (1) no ISO 
% % (2) ISO 
% % (3) ISO + IKs Phosp Block 
% % (4) ISO + ICaL Phosp Block (Results presented in SFigure11)
% % (5) ISO + PLB Phosp Block (Results presented in SFigure11)
% % (6) ISO + TnI Phosp Block (Results presented in SFigure11)
% % (7) ISO + INa Phosp Block (Results presented in SFigure11)
% % (8) ISO + INaK Phosp Block (Results presented in SFigure11)
% % (9) ISO + RyR Phosp Block (Results presented in SFigure11)
% % (10) ISO + IKur Phosp Block (Results presented in SFigure11)
% 
% %---: Functions required to run this script :---%
% % mainHRdBA.m - runs Heijman model simulation 
% % cleandata.m - remove APs with EADs
% %--------------------------------------------------------------------------
% %% 
% %---- Set Up Simulation ----%
% % parameters to vary in the model
% c.ICaLB = 1; c.IKsB = 1; c.IKrB = 1; c.INaKB =1; c.INaCaB = 1; c.IKpB = 1;
% c.IK1B = 1;  c.INabB = 1; c.ITo1B = 1; c.ITo2B = 1; c.INaB = 1; c.INaLB = 1;
% c.IClB = 1; c.IpCaB = 1; c.ICabB = 1; c.IrelB = 1; c.IupB = 1; c.IleakB = 1;
% 
% % settings
% settings.bcl = 1000;
% settings.freq = 100;
% settings.storeLast = 1;
% settings.stimdur = 2;
% settings.Istim = -36.7;
% settings.showProgress = 0;
% settings.SS = 1; % run steady state initial conditions 
% 
% variations = 20;
% settings.sigma = 0.2;
% scalings = exp(settings.sigma*randn(length(fieldnames(c)),variations))' ;
% 
% settings.totalvars = variations;
% settings.t_cutoff = 3;
% settings.flag = 1;
% 
% sims = {'noISO','ISO','IKs','ICaL','PLB','TnI','INa','INaK','RyR','IKur'};
% % following variables should be the same length as sims matrix 
% ISO = [0 1 1 1 1 1 1 1 1 1]; 
% % change to 0 if blocking target 
% IKs = [1 1 0 1 1 1 1 1 1 1]; % block IKs Phosphorylation in third iteration 
% ICaL = [1 1 1 0 1 1 1 1 1 1];
% PLB = [1 1 1 1 0 1 1 1 1 1];
% TnI = [1 1 1 1 1 0 1 1 1 1];
% INa = [1 1 1 1 1 1 0 1 1 1];
% INaK = [1 1 1 1 1 1 1 0 1 1];
% RyR = [1 1 1 1 1 1 1 1 0 1];
% IKur = [1 1 1 1 1 1 1 1 1 0];
% 
% %---- Run Simulation ----%
% for ii = 1:length(sims)
%     flags.ISO = ISO(ii);
%     flags.IKs = IKs(ii); % block IKs Phosphorylation in third iteration
%     flags.ICaL = ICaL(ii);
%     flags.PLB = PLB(ii);
%     flags.TnI = TnI(ii);
%     flags.INa = INa(ii);
%     flags.INaK = INaK(ii);
%     flags.RyR = RyR(ii);
%     flags.IKur = IKur(ii);
% 
%     parfor i = 1:variations
%         scaling = scalings(i,:);
%         [currents,State,Ti,APD]=mainHRdBA(settings,flags,scaling);
%         datatable(i).times = Ti;
%         datatable(i).V =  State(:,1);
%         datatable(i).states =  State;
%         datatable(i).currents= currents;
%         datatable(i).APDs = APD;
%         
%         disp(['Model Variant # ', num2str(i), ' ', sims{ii}]);
%     end
%     
%     X.(sims{ii}) = reformat_data(datatable,variations);    
% end
% 
% %% Plot Figure 7A,7B,7C
% 
% %---- Check for the presence of EADs in no ISO population. ----%
% [APfails,nEADs] = cleandata(cell2mat(X.noISO.APDs),X.noISO.times(:,1),X.noISO.V(:,1),settings.t_cutoff,settings.flag);
% [~,ind_failed] = find(APfails ==1); % number of failed to repolarize
% [~,ind_EADs] = find(nEADs==1); % number of EADs
% indexs = [ind_EADs ind_failed];
% 
% if isempty(indexs)
%     disp('No EADs in no ISO population, so no APs removed.')
%     for ii = 1:length(sims)       
%         APDs = cell2mat(X.(sims{ii}).APDs);
%         pert1 = prctile(APDs,90);
%         pert2 = prctile(APDs,10);
%         X.(sims{ii}).APDSpread =(pert1 - pert2 )/ median(APDs);
%     end
%     
% else
%     clean_datatable = [];
%     for ii = 1:length(sims)
%         clean_datatable.(sims{ii}).times =  X.(sims{ii}).times(~(nEADs' + APfails'));
%         clean_datatable.(sims{ii}).V =  X.(sims{ii}).V(~(nEADs' + APfails'));
%         clean_datatable.(sims{ii}).states =  X.(sims{ii}).states(~(nEADs' + APfails'));
%         clean_datatable.(sims{ii}).APDs =  X.(sims{ii}).APDs(~(nEADs' + APfails'));
%         clean_datatable.(sims{ii}).scaling =  scalings(~(nEADs' + APfails'),:);
%         
%         APDs = cell2mat(clean_datatable.(sims{ii}).APDs);
%         pert1 = prctile(APDs,90);
%         pert2 = prctile(APDs,10);
%         clean_datatable.(sims{ii}).APDSpread =(pert1 - pert2 )/ median(APDs);
%         
%         figure
%         for i = 1:20
%             plot(X.(sims{ii}).times{i,1},X.(sims{ii}).V{i,1},'linewidth',2)
%             hold on
%         end
%         xlabel('time(ms)')
%         ylabel('voltage(mV)')
%         title(sims{ii})
%     end
%     X = clean_datatable;
%     disp([num2str(length(indexs)) ' APs removed from each population.'])
% end
% 
% %% Plot Figure 7D 
% spreads = cell2mat(cellfun(@(x) X.(x).APDSpread,sims,'UniformOutput',0));
% figure 
% summary_barplot = gcf;
% ax_summary = axes('parent', summary_barplot);
% bar(spreads,0.5)
% ylabel('APD Spread')
% xticklabels(sims)
% 
% disp(['Final Population has: ' num2str(variations - length(indexs)) ' AP variants.'])

%% Figure 7E 
%--- Description of Figure: 
% Simulations in high(10x) and low(0.1x) IKs models performed in Ohara Model  

%---: Functions required to run this script :---%
% find_ICaLfactor.m - finds the ICaL factor that induces EAD in AP 
%--------------------------------------------------------------------------
%% 
%---- Set up Simulation ----%
% settings 
settings.freq =100;
settings.storeLast = 10;
settings.stimdur = 2;
settings.Istim = -36.7;
settings.showProgress = 0;
settings.bcl =1000;

%ISO = [0 1 1];
settings.SS = 1;

flags.ICaL = 1; flags.IKs = 1; flags.PLB = 1; flags.TnI = 1; flags.INa = 1;
flags.INaK = 1; flags.RyR = 1; flags.IKur = 1;
starts = [2 36.9 12.9 5.6 42.5 1 1 1 1 1]; % if you want to change where to begin calcium perturbation

sims = {'noISO','ISO','IKs','ICaL','PLB','TnI','INa','INaK','IKur'};
% following variables should be the same length as sims matrix 
ISO = [0 1 1 1 1 1 1 1 1 1]; 
% change to 0 if blocking target 
IKs = [1 1 0 1 1 1 1 1 1 1]; % block IKs Phosphorylation in third iteration 
ICaL = [1 1 1 0 1 1 1 1 1 1];
PLB = [1 1 1 1 0 1 1 1 1 1];
TnI = [1 1 1 1 1 0 1 1 1 1];
INa = [1 1 1 1 1 1 0 1 1 1];
INaK = [1 1 1 1 1 1 1 0 1 1];
%RyR = [1 1 1 1 1 1 1 1 1 1];
IKur = [1 1 1 1 1 1 1 1 1 0];

% starts = 50;
% sims = {'PLB'};
% PLB = 0;
% ISO = 1;
parfor index = 1:length(sims)
    %     flags.IKs = IKs(index);
    %     flags.ICaL = ICaL(index);
    %     flags.PLB = PLB(index);
    %     flags.TnI = TnI(index);
    %     flags.INa = INa(index);
    %     flags.INaK = INaK(index);
    %     flags.RyR = RyR(index);
    %     flags.IKur = IKur(index);
    %
    %     settings.ISO = ISO(index);
    
    inputs = [ISO(index),IKs(index),ICaL(index),PLB(index),TnI(index)...
        INa(index),INaK(index),1,IKur(index)];    
    
    flags = setparameters(inputs);
    start = starts(index);
    %disp(sims{index})
    ICaL_Factors(index) = find_ICaLfactor(start,settings,flags,sims{index});
%     figure(gcf)
%     title([sims{index} ' : ICaL Factor - ' num2str(ICaL_Factors(index))])
end

figure 
summary_barplot = gcf;
ax_summary = axes('parent', summary_barplot);
bar(ICaL_Factors,0.5)
ylabel('ICaL Factor')
xticklabels(sims)


% EADs occur at a ICaL Factor of 2 on beat 91. -- no ISO 
% EADs occur at a ICaL Factor of 36.9 on beat 98. -- ISO 
% EADs occur at a ICaL Factor of 12.9 on beat 91. -- ISO block IKs
% EADs occur at a ICaL Factor of 5.6 on beat 91. -- ISO block ICaL
