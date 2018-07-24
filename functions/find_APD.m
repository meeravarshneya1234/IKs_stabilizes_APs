function [APD,repoldex,tinit] = find_APD(time, Voltage)
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

end