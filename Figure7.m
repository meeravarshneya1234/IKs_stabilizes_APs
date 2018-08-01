% %% Figure 7: IKs stabilizes the AP during ?-adrenergic stimulation.
% %--- Note: 
% % Results displayed in manuscript were run using MATLAB 2016a on a 64bit
% % Intel Processor. For exact replication of figures it is best to use these settings.
% 
% %% Figure 7A 
% %--- Description of Figure: 
% % Simulations in high(10x) and low(0.1x) IKs models performed in Ohara Model  
% 
% %---: Functions required to run this script :---%
% % main_program.m - runs single AP simulation 
% % cleandata.m 
% %--------------------------------------------------------------------------
% 
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
% variations = 300;
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
% for ii = 1:length(sims)
%     settings.ISO = ISO(ii);
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
%     figure
%     for i = 1:20
%         plot(X.(sims{ii}).times{i,1},X.(sims{ii}).V{i,1},'linewidth',2)
%         hold on
%     end 
%     xlabel('time(ms)')
%     ylabel('voltage(mV)')
%     title(sims{ii})   
% end
% 
% %% Check for the presence of EADs in no ISO population.
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
%     end
%     X = clean_datatable;
%     disp([num2str(length(indexs)) ' APs removed from each population.'])
% end
% 
% 
% spreads = cell2mat(cellfun(@(x) X.(x).APDSpread,sims,'UniformOutput',0));
% figure 
% summary_barplot = gcf;
% ax_summary = axes('parent', summary_barplot);
% bar(spreads,0.5)
% ylabel('APD Spread')
% xticklabels(sims)
% 
% disp(['Final Population has: ' num2str(variations - length(indexs)) ' AP variants.'])

%% Figure 7B 
%--- Description of Figure: 
% Simulations in high(10x) and low(0.1x) IKs models performed in Ohara Model  

%---: Functions required to run this script :---%
% main_program.m - runs single AP simulation 
% cleandata.m 
%--------------------------------------------------------------------------

% settings 
settings.freq =100;
settings.storeLast = 10;
settings.stimdur = 2;
settings.Istim = -36.7;
settings.showProgress = 0;
settings.bcl =1000;

colors = hsv(4);
ISO = [0 1 1];
settings.SS = 1;

flags.ICaL = 1; flags.IKs = 1; flags.PLB = 1; flags.TnI = 1; flags.INa = 1;
flags.INaK = 1; flags.RyR = 1; flags.IKur = 1;
starts = [2 36.9 12.9];

sims = {'noISO','ISO','IKsBP'};
IKsBP = [1 1 0];

for index = 1:length(sims)
    flags.IKs = IKsBP(index);
    settings.ISO = ISO(index);
    start = starts(index);
    ICaL_Factors(index) = find_ICaLfactor(start,settings,flags);
    disp(sims{index})
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
