function datatable = pop_program(settings)
%% -- pop_program.m -- %%
% Description: Runs a heterogenous population.

% Inputs:
% --> settings - AP stimulation protocol (PCL,nBeats,...)

% Outputs:
% --> datatable - struct array with AP info including: APD, time matrix
% and state variables.

%---: Functions used in this script :---%
% --* parameters.m - Extracts the model parameters for the model called.
% --* ICs.m - Extracts the initial conditions for the model called.
% --* scaling_factors.m - Varies the parameters to create the population.
% --* dydt functions - one for each model, calculates state variables.
%--------------------------------------------------------------------------
%% 1--- Load Model Parameters
[p,c] = parameters(settings.model_name,settings.celltype);
% p --> main model parameters
% c --> model parameters that can be varied to build population

%% 2--- Load Model Initial Conditions of State Variables
[y0,V_ind] = ICs(settings.model_name,settings.steady_state,settings.PCL);
% y0 --> initial conditions
% Vind --> index of voltage in the initial conditions matrix

%% 3--- Define Simulation Protocol
odefcn = str2func(strcat('dydt_',char(settings.model_name)));
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

%% 4--- Define Perturbation Protocol
GKs_baseline = c.GKs_; %save baseline value
GKr_baseline = c.GKr_; %save baseline value
baselineparameters = c;
%% 5--- Run Simulation
disp(strjoin({char(settings.model_name),'Model Test'}))

F = fieldnames(c);

baselineparameters.GKr_ = GKr_baseline*settings.Kr_scale;
baselineparameters.GKs_ = GKs_baseline*settings.Ks_scale;

if ~isfield(settings,'scalings') %if you want to specify a certain scaling matrix ahead of time
    scalings = exp(settings.sigma*randn(length(F),settings.variations))';
else
    scalings = settings.scalings;
end

parfor ii=1:settings.variations
    
    scaling = scalings(ii,:);
    c = scaling_factors(scaling,baselineparameters);
    statevar_i = y0;
    
    %%% stimulate cell
    if (settings.nBeats > settings.numbertokeep)
        for i=1:simints-3*settings.numbertokeep
            options = odeset('RelTol',1e-3,'AbsTol',1e-6);
            [post,posstatevars] = ode15s(odefcn,intervals(i,:),statevar_i,options,Istim(i),p,c) ;
            statevar_i = posstatevars(end,:) ;
            t = post(end) ;
        end % for
        statevars = statevar_i ;
        for i=simints-3*settings.numbertokeep+1:simints
            options = odeset('RelTol',1e-3,'AbsTol',1e-6);
            [post,posstatevars] = ode15s(odefcn,intervals(i,:),statevar_i,options,Istim(i),p,c) ;
            t = [t;post(2:end)] ;
            statevars = [statevars;posstatevars(2:end,:)] ;
            statevar_i = posstatevars(end,:) ;
        end % for
    else
        t = 0 ;
        statevars = statevar_i ;
        for i=1:simints
            options = odeset('RelTol',1e-3,'AbsTol',1e-6);
            [post,posstatevars] = ode15s(odefcn,intervals(i,:),statevar_i,options,Istim(i),p,c) ;
            t = [t;post(2:end)] ;
            statevars = [statevars;posstatevars(2:end,:)] ;
            statevar_i = posstatevars(end,:) ;
        end % for
    end % if
    
    time = t - min(t) ;
    V = statevars(:,V_ind);
    
    APD = find_APD(time,V);
    disp([int2str(ii),' of ',int2str(settings.variations),' simulations completed.'])
    
    datatable(ii).APDs = APD;
    datatable(ii).times = time;
    datatable(ii).states = statevars;
    datatable(ii).V = V;
    datatable(ii).scalings = scaling;
    
end


