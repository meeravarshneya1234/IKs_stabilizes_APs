function [datatable] = main_program(settings)

%% 1--- Load Model Parameters 
% First check if running the optimized models, the celltype selected is endo.
[p,c] = parameters(settings.model_name,settings.celltype);

%% 2--- Load Model Initial Conditions of State Variables 
[y0,V_ind] = ICs(settings.model_name,settings.steady_state,settings.PCL);
%y0 = cell2mat(struct2cell(ic))';

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
PCa_baseline = c.PCa_; %save baseline value

%% 5--- Run Simulation
n = 0; % number of simulations completed 
disp(strjoin({char(settings.model_name),'Model Test'}))

for Ks=1:length(settings.Ks_scale) % for every Ks perturbation 
    c.GKs_=GKs_baseline*settings.Ks_scale(Ks);
    c.GKr_=GKr_baseline*settings.Kr_scale(Ks);    

    for Ca=1:length(settings.Ca_scale) % for every CaL perturbation 
        c.PCa_=PCa_baseline*settings.Ca_scale(Ca);
        
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
        
%         figurestrings{Ca} = ['P_C_a* ' num2str(settings.Ca_scale(Ca))]; % for legend
        
        n = n+1;               
        disp(['Completed GKs*', num2str(settings.Ks_scale(Ks)),' & PCa*', num2str(settings.Ca_scale(Ca))])
        disp([int2str(n),' of ',int2str((length(settings.Ca_scale)*length(settings.Ks_scale))),' simulations completed.'])
        
        datatable.APDs(Ks,Ca) = APD;
        datatable.Area_Ca(Ks,Ca) = Area_Ca;
        datatable.Area_Ks(Ks,Ca) = Area_Ks;
        datatable.Area_Kr(Ks,Ca) = Area_Kr;
        datatable.ICaL{Ks,Ca} = ICaL;
        datatable.IKs{Ks,Ca} = IKs;
        datatable.IKr{Ks,Ca} = IKr;
        datatable.IKs_Fraction(Ks,Ca) = IKs_Fraction;
        datatable.times{Ks,Ca} = time_variables;
        datatable.states{Ks,Ca} = state_variables;
        datatable.V{Ks,Ca} = V;
        
    end
end
