function X2 = reformat_data(X1, variations) 

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
                        %% -- reformat_data.m -- %%
% Description: Reformats data that was run using the parfor loop. 

% Inputs:
% --> X1 - [struct array] with AP info including: APD, time matrix
% and state variables.
% --> variations - [double] number of variants in the population 

% Outputs:
% --> X2 - [struct array] reformatted 
%--------------------------------------------------------------------------
%%
F = fieldnames(X1);
for ind = 1:variations
    for iF = 1:length(F)
        aF = F{iF};
        X2.(aF){ind,1} = X1(ind).(aF);
    end
end
