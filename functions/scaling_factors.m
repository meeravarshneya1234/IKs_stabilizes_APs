function c = scaling_factors(scaling,baselineparameters)
F = fieldnames(baselineparameters);
for iF = 1:length(F)
    aF = F{iF};
    c.(aF) = baselineparameters.(aF) * scaling(iF);
end