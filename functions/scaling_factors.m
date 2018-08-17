function c = scaling_factors(scaling,baselineparameters)

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
                            %% -- scaling_factors.m -- %%
% Description: Randomly vary the model parameters that control ion transfer rates. 

% Inputs:
% --> scaling - [array] random multiplicative factor to scale each parameter
% --> baselineparameters - [struct array] save the original value of the parameter 

% Outputs:
% --> c - [struct array] altered parameters to create population. 
%--------------------------------------------------------------------------
%% 
F = fieldnames(baselineparameters);
for iF = 1:length(F)
    aF = F{iF};
    c.(aF) = baselineparameters.(aF) * scaling(iF);
end