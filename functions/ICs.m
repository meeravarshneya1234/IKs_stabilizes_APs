function [ic,Vind] = ICs(model_name,steady_state,pcl)

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
                            %% -- ICs.m -- %%
% Description: Extracts either the steady state or original paper initial
% conditions for the model called. 

% Note: Steady state conditions were collected for 500 ms, 1000 ms, 5000 ms
% for each model. If a different pacing length is desired please set 
% 'steady_state' variable to 0 and run for extended number of beats to
% reach steady state. For more information on steady state analysis, please see
% Methods and Table S1. 

% Inputs:
% --> model_name - [string] name of model for which the parameters are needed.
% options: 'Fox'(dog),'Hund'(dog),'Livshitz'(guinea pig),'Devenyi'(guinea pig),
% 'Shannon'(rabbit),'TT04'(human),'TT06'(human),'Grandi'(human),'Ohara'(human)
% --> steady_state - [binary] 1 (steady state conditions); 0 (non-steady state
% conditions) 
% --> pcl - [binary] the pacing length for the simulation to be run. 

% Outputs:
% --> ic - [struct array] with initial conditions
% --> Vind - [double] index of Voltage in the matrix of initial conditions

%---: Functions used in this script :---%
% ** initial_conditions.m - Extracts the initial conditions from the paper if 
% steady_state = 0 or if PCL is not 500,1000,5000.
% -------------------------------------------------------------------------
%% 

if steady_state == 1 % use steady state initial conditions 
    load('ICs_SS.mat')
    if sum(pcl == [500,1000,5000]) > 0
        str =['PCL_' num2str(pcl)];
        ic = ICs.(model_name).(str);
    else % steady state initial conditions do not exist for that PCL 
        disp(['Initial conditions for steady state PCL = ' num2str(pcl) ' do not exist. Using the non-steady state inital conditions.'])
        ic = initial_conditions(model_name);
    end
else % use initial conditions from original manuscript
    ic = initial_conditions(model_name);
    
end

% Determine index of voltage in the state variables. 
if strcmp(model_name,'Shannon')
    Vind = 32;
elseif strcmp(model_name,'Grandi') || strcmp(model_name,'Morotti')
    Vind = 39;
else 
    Vind = 1;
end 

