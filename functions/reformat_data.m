% function X2 = reformat_data(X1, variations) 
% for ind = 1:variations
%     X2.APDs(ind,1) = X1(ind).APDs;
%     X2.times{ind,1} = X1(ind ).times;
%     X2.V{ind,1} = X1(ind ).V;
%     X2.states{ind,1} = X1(ind ).states;
%     X2.scalings{ind,1} = X1(ind).scalings;
%     
%     if isfield(X1,'currents')
%         X2.currents{ind,1} = X1(ind).currents;
%     end
% end

function X2 = reformat_data(X1, variations) 
F = fieldnames(X1);
for ind = 1:variations
    for iF = 1:length(F)
        aF = F{iF};
        X2.(aF){ind,1} = X1(ind).(aF);
    end
end
