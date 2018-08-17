function datatable = FigureS10Heijman()

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

                            %% -- FigureS10Heijman.m -- %%
% Description: Run AP Clamp in Heijman model 

% Outputs:
% --> datatable - struct that outputs the APDs, time, voltage, and state variables 

%---: Functions required to run this script :---%
% mainHRdBA.m - runs baseline simulation for Heijman  
% funHRdBA_APclamp - runs AP Clamp for Heijman model 
%--------------------------------------------------------------------------
%% Run under baseline conditions to get input voltage and time for AP Clamp
settings.freq =1;%number of beats to stimulate 
settings.storeLast =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats 
settings.stimdur = 2;% Stimulus duration
settings.Istim = -36.7;%Stimulus amplitude for each model, see Table S1 for details
settings.bcl =1000;% Interval bewteen stimuli,[ms]
settings.ISO = 0; %concentration of ISO to add; 0 = no ISO 
settings.SS = 1;% 1 - run steady state conditions 0 - do not run steady state conditions 
settings.showProgress = 0;

% option to block PKA targets: no block = 1; 100% block = 0 
flags.ICaL = 1; flags.IKs = 1; flags.PLB = 1; flags.TnI = 1; flags.INa = 1;
flags.INaK = 1; flags.RyR = 1; flags.IKur = 1;

[currents,State,Ti,APDs,settings]=mainHRdBA(settings,flags);
t_input = Ti; 
state_vars_input = State;

%% Run AP Clamp for multiple APDs - Log Scale
statevar_i_original_orig = state_vars_input(end,2:end);
logscalefactors = linspace(log(1/3),log(3),11) ;
repol_change = exp(logscalefactors) ; % scaling
numberofAPs = 100;
numberkeep = 1;
PCL = 1000;
volts = state_vars_input(:,1);

start_point = find(volts>-75,1);
[~,ind_V] = max(volts);
end_point = find( t_input>t_input(ind_V) & volts<-75,1);

[APD_base,repol_index,tinit] = find_APD(t_input,volts);


% code to include a parfor loop 
tic 
parfor ind = 1:length(repol_change)
    settings = [];
    settings.ISO = 0;
    settings.dataElectrophysiol = constants_Electrophysiol;
    settings.dataSignaling = constants_Signaling;
    settings = setDefaultSettings(settings);

    statevar_i_original = statevar_i_original_orig;

    State = zeros(0,length(statevar_i_original_orig));
    Ti = zeros(0,1);
    Voltages = zeros(0,1);
    IKr = zeros(0,1);
    IKs = zeros(0,1);
    Area_Kr = zeros(0,1);
    Area_Ks = zeros(0,1);
    Ratio = zeros(0,1);
    sv_matrix = zeros(0,145);
    
    new_APD = (APD_base*repol_change(ind))+tinit;
    new_end = find(new_APD>t_input,1,'last');
    
    % new V & t values between -75 mV back to -75 mV
    %l1 = (repol_index-start_point)+1;
    l1 = 350;
    t_interp1 = linspace(t_input(start_point),t_input(end_point),l1); % original data
    t_repol = linspace(t_input(start_point),t_input(new_end),l1);
    v_repol = interp1(t_input(start_point:end_point),volts(start_point:end_point),t_interp1);
    
    % new V & t during relaxation
    %l2 = 120;%(length(t_input)- (length(v_repol)+start_point))+1;
    l2 = 350;
    t_interp2 = linspace(t_input(end_point+1),t_input(numel(t_input)),l2);
    t_relax = linspace((t_input(new_end)+1),t_input(numel(t_input)),l2);
    v_relax = interp1(t_input(end_point+1:numel(t_input)),volts(end_point+1:numel(volts)),t_interp2);
    
    t = [t_input(1:start_point-1); t_repol';t_relax'];
    V = [volts(1:start_point-1); v_repol';v_relax'];
    
    %set up segements for ode solver
    simints = length(t)-1 ;
    intervals = cell(simints,1) ;
    voltage = zeros(simints,1) ;
    for i=1:length(intervals)
        intervals{i} = [t(i),t(i+1)] ;
        voltage(i) = V(i) ;
    end
    
    index = 0;
    for ii = 1:numberofAPs
        time = 0 ;
        Voltage =voltage(1) ;
        statevar_i = statevar_i_original;
        statevars = statevar_i ;
        tic
        for i=1:simints            
            [post,posstatevars] = ode15s(@funHRdBA_APclamp,intervals{i},statevar_i,[],flags,voltage(i),settings) ;
            time = [time;post(2:end)] ;
            statevars = [statevars;posstatevars(2:end,:)] ;
            statevar_i = posstatevars(end,:) ;
            Voltage = [Voltage;voltage(i)*ones(length(post)-1,1)] ;            
        end
        toc
        true_vals=ismember(time,t);
        true_statevars = statevars(true_vals,:);
        outputcell = num2cell(true_statevars,1) ;
        statevar_i_original = cellfun(@(v) v(end), outputcell); % inital conditions for next AP
        
        figure 
        plot(time(true_vals),Voltage(true_vals),'linewidth',2)
        
        disp(['AP # ' num2str(ii)])
        
        if ii > (numberofAPs - numberkeep) 
            Vind = 1;
            index = index + 1;
            [nRows,~] = size(true_statevars);
            State(end+1:end+nRows,:)  = true_statevars;
            Ti(end+1:end+nRows,1) = time(true_vals)+PCL*(index-1);
            Voltages(end+1:end+nRows,1) = Voltage(true_vals);
            sv_matrix = [true_statevars(:,1:(Vind-1)) Voltage(true_vals) true_statevars(:,(Vind:end))];
            % calculate IKr and IKs 
            sigdata = settings.dataSignaling;
            K_i = sv_matrix(:, 10);
            Na_i = sv_matrix(:,8);
            data = constants_Electrophysiol;
            EK=log(settings.K_o./K_i)/data.frt;
            EKs = log((settings.K_o+data.prnak*settings.Na_o)./(K_i+data.prnak*Na_i))/data.frt;
            Ca_i = sv_matrix(:,2);
            MM_IKsNP_O1 = sv_matrix(:,47); MM_IKsNP_O2 = sv_matrix(:,48);
            MM_IKsP_O1 = sv_matrix(:,72); MM_IKsP_O2 = sv_matrix(:,73);
            fIKsP = min(max(((sv_matrix(:,129) + sigdata.IKs_AKAP_PKA) / sigdata.IKs_tot - (0.0306 + sigdata.IKs_AKAP_PKA / sigdata.IKs_tot)) / (0.7850 - (0.0306 + sigdata.IKs_AKAP_PKA / sigdata.IKs_tot)), 0), 1);
            IKsNP = settings.IKsNPParams(end) .* (1+0.6./(1+(((3.8E-05)./Ca_i).^1.4))) .* (MM_IKsNP_O1 + MM_IKsNP_O2) .* (Voltages-EKs);
            IKsP = settings.IKsPParams(end) .* (1+0.6./(1+(((3.8E-05)./Ca_i).^1.4))) .* (MM_IKsP_O1 + MM_IKsP_O2) .* (Voltages-EKs);
            iKs = (settings.IKsB) .* (fIKsP .* IKsP + (1 - fIKsP) .* IKsNP);
            xr = sv_matrix(:,22);
            iKr = (settings.IKrB) * (data.GKrmax*sqrt(settings.K_o/5.4)) .* xr .* (1./(1+exp((V+10)/15.4))) .* (V-EK);
            
            IKr(end+1:end+nRows,1)= iKr;
            IKs(end+1:end+nRows,1)= iKs;
            Area_Kr(end+1:end+index,1) = trapz(time(true_vals),iKr);
            Area_Ks(end+1:end+index,1) = trapz(time(true_vals),iKs);
            Ratio(end+1:end+index,1) =  trapz(time(true_vals),iKs)/(trapz(time(true_vals),iKr) + trapz(time(true_vals),iKs));
        end       
    end       
    datatable(ind).statevars = sv_matrix;
    datatable(ind).Voltage = Voltages;
    datatable(ind).time = Ti;
    datatable(ind).IKs = IKs;
    datatable(ind).IKr = IKr;
    datatable(ind).Area_Kr = Area_Kr;
    datatable(ind).Area_Ks = Area_Ks;
    datatable(ind).Ratio =  Ratio;
    
    disp(['APD ' num2str(repol_change(ind))])
end
toc


