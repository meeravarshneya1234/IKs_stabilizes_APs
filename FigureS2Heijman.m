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
settings.freq =100;
settings.storeLast = 1;
settings.stimdur = 2;
settings.Istim = -36.7;
settings.showProgress = 1;

variations = 300;
settings.sigma = 0.2;
scalings = exp(settings.sigma*randn(length(fieldnames(c)),variations))' ;

settings.totalvars = variations;
settings.t_cutoff = 3;
settings.flag = 1;

colors = repmat('krgbmc',1,1000);
settings.ISO = 0;
flags.SS = 1;

flags.ICaL = 1;
flags.IKs = 1;
flags.PLB = 1;
flags.TnI = 1;
flags.INa = 1;
flags.INaK = 1;
flags.RyR = 1;
flags.IKur = 1;

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
