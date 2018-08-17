function [datatable] = inject_current_program(settings)

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
                    %% -- inject_current_program.m -- %%
% Description: Runs the "inject current" simulation in the specified model.

% Inputs:
% --> settings - [struct array] AP stimulation protocol (PCL,nBeats,...)

% Outputs:
% --> datatable - [struct array] AP info including: APD, time matrix
% and state variables, AUC of ICaL,IKs,IKr waveforms, ICaL,IKs,IKr currents, 
% IKs Fraction. 

%---: Functions used in this script :---%
% ** parameters.m - Extracts the model parameters for the model called. 
% ** ICs.m - Extracts the initial conditions for the model called. 
% ** dydt functions - one for each model, calculates state variables.
%--------------------------------------------------------------------------

%% 1--- Load Model Parameters 
% First check if running the optimized models, the celltype selected is endo.
[p,c] = parameters(settings.model_name,settings.celltype);

%% 2--- Load Model Initial Conditions of State Variables 
[y0,V_ind] = ICs(settings.model_name,settings.steady_state,settings.PCL);
%y0 = cell2mat(struct2cell(ic))';

%% 3--- Define Simulation Protocol 
odefcn = str2func(strcat('dydt_',char(settings.model_name),'_inject'));

stim_starts = settings.stim_delay + settings.PCL*(0:settings.nBeats-1)  ;
stim_ends = stim_starts + settings.stim_dur ;

% Create intervals for each beat 
simints = 3*settings.nBeats ;
for i=1:settings.nBeats
    intervals(3*i-2,:) = [settings.PCL*(i-1),stim_starts(i)] ; %beginning 
    intervals(3*i-1,:) = [stim_starts(i),stim_ends(i)] ; %stimulus 
    intervals(3*i,:) = [stim_ends(i),settings.PCL*i] ; % stimulus ends 
end
tend = settings.nBeats*settings.PCL ;              % end of simulation, ms
intervals(end,:) = [stim_ends(end),tend] ;

% Determine when to apply stim_amp or 0 amp  
Istim = zeros(simints,1) ;
stimindices = 3*(1:settings.nBeats) - 1 ; % apply stimulus on second part of intervals
Istim(stimindices) = -settings.stim_amp ; 

%% 4--- Run Simulation
disp(strjoin({char(settings.model_name),'Model Test'}))
parfor I=1:length(settings.Inject) % for every CaL perturbation
    Inject = settings.Inject(I);
    
    statevar_i = y0;
    %%% stimulate cell
    if (settings.nBeats > settings.numbertokeep)
        for i=1:simints-3*settings.numbertokeep
            options = odeset('RelTol',1e-3,'AbsTol',1e-6);
            [post,posstatevars] = ode15s(odefcn,intervals(i,:),statevar_i,options,Istim(i),p,c,Inject) ;
            statevar_i = posstatevars(end,:) ;
            t = post(end) ;
        end % for
        statevars = statevar_i ;
        for i=simints-3*settings.numbertokeep+1:simints
            options = odeset('RelTol',1e-3,'AbsTol',1e-6);
            [post,posstatevars] = ode15s(odefcn,intervals(i,:),statevar_i,options,Istim(i),p,c,Inject) ;
            t = [t;post(2:end)] ;
            statevars = [statevars;posstatevars(2:end,:)] ;
            statevar_i = posstatevars(end,:) ;
        end % for
    else
        t = 0 ;
        statevars = statevar_i ;
        for i=1:simints
            options = odeset('RelTol',1e-3,'AbsTol',1e-6);
            [post,posstatevars] = ode15s(odefcn,intervals(i,:),statevar_i,options,Istim(i),p,c,Inject) ;
            t = [t;post(2:end)] ;
            statevars = [statevars;posstatevars(2:end,:)] ;
            statevar_i = posstatevars(end,:) ;
        end % for
    end % if
    
    t = t - min(t) ;
    V = statevars(:,V_ind);
    outputcell = num2cell(statevars,1) ;
    state_variables = outputcell;
    time_variables = t;
    
    % Determine APD
    APD = find_APD(t,V);
    
    % Calculate Currents
    [ICaL,IKr,IKs] = currents(settings.model_name,outputcell,p,c);
    
    % determine Integrated AUCs of currents
    x1 = find(t==settings.stim_delay); % find interval where the AP begins
    x2 = find(floor(V)==floor(V(x1))+3 & t > t(x1)+10,1); % and ends
    Area_Ks = trapz(t(x1:x2),IKs(x1:x2));
    Area_Kr = trapz(t(x1:x2),IKr(x1:x2));
    Area_Ca = trapz(t(x1:x2),ICaL(x1:x2));
    
    % determine fraction of IKs based on ratio of IKs and IKr
    IKs_Fraction = Area_Ks/(Area_Kr+Area_Ks);   
    disp([int2str(I),' of ',int2str((length(settings.Inject))),' simulations completed.'])
    
    datatable(I).APDs = APD;
    datatable(I).Area_Ca = Area_Ca;
    datatable(I).Area_Ks = Area_Ks;
    datatable(I).Area_Kr = Area_Kr;
    datatable(I).ICaL= ICaL;
    datatable(I).IKs = IKs;
    datatable(I).IKr = IKr;
    datatable(I).IKs_Fraction = IKs_Fraction;
    datatable(I).times = time_variables;
    datatable(I).states= state_variables;
    datatable(I).V = V;
    
end
