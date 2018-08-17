function [APD,repoldex,tinit] = find_APD(time, Voltage)

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
                            %% -- find_APD.m -- %%
% Description: Calculates a single AP's APD. Determined as the interval
% between the stimulus and when the membrane potential return to below -75
% mV.

% Inputs:
% --> time - [array] single AP time vector in ms 
% --> Voltage - [array] single AP voltage vector in mV

% Outputs:
% --> APD - [double] action potential duration 
% --> repoldex  - [double] index when voltage returns to -75 mV (used in AP
% clamp code)
% --> tinit -  [double] time at beginning of action potential (used in AP
% clamp code)
%--------------------------------------------------------------------------
%%

dVdt = diff(Voltage)./diff(time) ;
[~,dexmax] = max(dVdt) ;

%t of maximum dV/dt, consider this beginning of action potential
tinit = time(dexmax) ;

% Then determine peak V of action potential, for two reasons,
% 1) Because repolarization must, by definition, occur after this
% 2) To compute 50%, 90%, etc., must have this value
[~,peakdex] = max(Voltage) ;
tpeak = time(peakdex) ;
repoldex = find(time > tpeak & Voltage < -75,1) ; % APD based on repolarization back to -75 mV instead of -60 mV
if isempty(repoldex)
    APD = NaN;
else
    repoltime = time(repoldex(1)) ;
    APD = round((repoltime - tinit)) ;
end