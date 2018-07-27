function c = scaling_factors(scaling,baselineparameters)
                            %% -- scaling_factors.m -- %%
% Description: Randomly vary the model parameters that control ion transfer rates. 

% Inputs:
% --> scaling - random multiplicative factor to scale each parameter
% --> baselineparameters - save the original value of the parameter 

% Outputs:
% --> c - altered parameters to create population. 
%--------------------------------------------------------------------------
%% 
F = fieldnames(baselineparameters);
for iF = 1:length(F)
    aF = F{iF};
    c.(aF) = baselineparameters.(aF) * scaling(iF);
end