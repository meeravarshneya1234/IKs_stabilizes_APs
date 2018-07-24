% parameters to vary in the model
c.ICaLB = 1;
c.IKsB = 1;
c.IKrB = 1;
c.INaKB =1;
c.INaCaB = 1;
c.IKpB = 1;
c.IK1B = 1;
c.INabB = 1;
c.ITo1B = 1;
c.ITo2B = 1;
c.INaB = 1;
c.INaLB = 1;
c.IClB = 1;
c.IpCaB = 1;
c.ICabB = 1;
c.IrelB = 1;
c.IupB = 1;
c.IleakB = 1;

% settings
settings.bcl = 1000;
settings.freq = 100;
settings.storeLast = 10;
settings.stimdur = 2;
settings.Istim = -36.7;
settings.showProgress = 0;

variations = 300;
settings.sigma = 0.2;
%scalings = exp(settings.sigma*randn(length(fieldnames(c)),variations))' ;
load('scalings.mat') %use same scaling matrix for every population

settings.totalvars = variations;
settings.t_cutoff = 3;
settings.flag = 1;

colors = repmat('krgbmc',1,1000);
%settings.ISO = 0;
ISO = [0 1 1]; 
flags.SS = 1;

flags.ICaL = 1;
% flags.IKs = 1;
IKsBP = [1 1 0]; % block IKs Phosphorylation in third iteration 
flags.PLB = 1;
flags.TnI = 1;
flags.INa = 1;
flags.INaK = 1;
flags.RyR = 1;
flags.IKur = 1;

sims = {'noISO','ISO','IKsBP'};
ICa_scales = { [1:0.1:4] [35:0.1:40] [10:0.1:15]}; %calcium perturbation for each model based on order of modelnames

for index = 1:length(sims)    
    flags.IKs = IKsBP(index);
    settings.ISO = ISO(index);   
    scale = ICa_scales{index};
        
    for ii = 1:length(scale) % for each calcium factor 
        settings.Ca_scale = scale(ii);
        parfor i = 1:variations % for each model variant 
            scaling = scalings(i,:);
            [currents,State,Ti,APD]=mainHRdBA_EAD(settings,flags,scaling);
            
            X(i).times = Ti;
            X(i).V =  State(:,1);
            X(i).states =  State;
            X(i).currents= currents;
            X(i).APDs = APD;
        end
        X_all{ii} = X;
        X =[];
        disp([sims{index} ' : ' num2str(ii) ' of ' num2str(length(scale))])
    end
    X_all_EAD.(sims{index}) = X_all; % save each population in a single struct
    X_all =[];
end

%     [APfails,nEADs] = cleandata(APD,times,V,settings.t_cutoff,settings.flag);
%     [~,ind_failed] = find(APfails ==1); % number of failed to repolarize
%     [~,ind_EADs] = find(nEADs==1); % number of EADs
%     total = [ind_EADs ind_failed];
%     ICaL_factors(i) = Ca_scale;

